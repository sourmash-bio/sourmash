import weakref

from ._lowlevel import ffi, lib
from .exceptions import exceptions_by_code, SourmashError

attached_refs = weakref.WeakKeyDictionary()


class RustObject:
    __dealloc_func__ = None
    _objptr = None
    _shared = False

    def __init__(self):
        raise TypeError("Cannot instanciate %r objects" % self.__class__.__name__)

    @classmethod
    def _from_objptr(cls, ptr, shared=False):
        rv = object.__new__(cls)
        rv._objptr = ptr
        rv._shared = shared
        return rv

    def _methodcall(self, func, *args):
        return rustcall(func, self._get_objptr(), *args)

    def _get_objptr(self):
        if not self._objptr:
            raise RuntimeError("Object is closed")
        return self._objptr

    def __del__(self):
        if self._objptr is None or self._shared:
            return
        f = self.__class__.__dealloc_func__
        if f is not None:
            rustcall(f, self._objptr)
            self._objptr = None


def decode_str(s):
    """Decodes a SourmashStr"""
    try:
        if s.len == 0:
            return u""
        return ffi.unpack(s.data, s.len).decode("utf-8", "replace")
    finally:
        if s.owned:
            lib.sourmash_str_free(ffi.addressof(s))


def encode_str(s):
    """Encodes a SourmashStr"""
    rv = ffi.new("SourmashStr *")
    if isinstance(s, str):
        s = s.encode("utf-8")
    rv.data = ffi.from_buffer(s)
    rv.len = len(s)
    # we have to hold a weak reference here to ensure our string does not
    # get collected before the string is used.
    attached_refs[rv] = s
    return rv


def rustcall(func, *args):
    """Calls rust method and does some error handling."""
    lib.sourmash_err_clear()
    rv = func(*args)
    err = lib.sourmash_err_get_last_code()
    if not err:
        return rv
    msg = lib.sourmash_err_get_last_message()
    cls = exceptions_by_code.get(err, SourmashError)
    exc = cls(decode_str(msg))
    backtrace = decode_str(lib.sourmash_err_get_backtrace())
    if backtrace:
        exc.rust_info = backtrace
    raise exc

###
### data class utils - from https://gist.github.com/mikeholler/4be180627d3f8fceb55704b729464adb
###


from dataclasses import is_dataclass
from typing import TypeVar, Type, Callable, List, Dict, Any

_T = TypeVar("_T")
_Self = TypeVar("_Self")
_VarArgs = List[Any]
_KWArgs = Dict[str, Any]


def _kwarg_only_init_wrapper(
        self: _Self,
        init: Callable[..., None],
        *args: _VarArgs,
        **kwargs: _KWArgs
) -> None:
    if len(args) > 0:
        raise TypeError(
            f"{type(self).__name__}.__init__(self, ...) only allows keyword arguments. Found the "
            f"following positional arguments: {args}"
        )
    init(self, **kwargs)


def _positional_arg_only_init_wrapper(
        self: _Self,
        init: Callable[..., None],
        *args: _VarArgs,
        **kwargs: _KWArgs
) -> None:
    if len(kwargs) > 0:
        raise TypeError(
            f"{type(self).__name__}.__init__(self, ...) only allows positional arguments. Found "
            f"the following keyword arguments: {kwargs}"
        )
    init(self, *args)


def require_kwargs_on_init(cls: Type[_T]) -> Type[_T]:
    """
    Force a dataclass's init function to only work if called with keyword arguments.
    If parameters are not positional-only, a TypeError is thrown with a helpful message.
    This function may only be used on dataclasses.
    This works by wrapping the __init__ function and dynamically replacing it. Therefore,
    stacktraces for calls to the new __init__ might look a bit strange. Fear not though,
    all is well.
    Note: although this may be used as a decorator, this is not advised as IDEs will no longer
    suggest parameters in the constructor. Instead, this is the recommended usage::
        from dataclasses import dataclass
        @dataclass
        class Foo:
            bar: str
        require_kwargs_on_init(Foo)
    """
    if cls is None:
        raise TypeError("Cannot call with cls=None")
    if not is_dataclass(cls):
        raise TypeError(
            f"This decorator only works on dataclasses. {cls.__name__} is not a dataclass."
        )

    original_init = cls.__init__

    def new_init(self: _Self, *args: _VarArgs, **kwargs: _KWArgs) -> None:
        _kwarg_only_init_wrapper(self, original_init, *args, **kwargs)

    # noinspection PyTypeHints
    cls.__init__ = new_init  # type: ignore

    return cls

