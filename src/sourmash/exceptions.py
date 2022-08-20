from ._lowlevel import lib


__all__ = ['SourmashError']
exceptions_by_code = {}


class SourmashError(Exception):
    code = None

    def __init__(self, msg):
        Exception.__init__(self)
        self.message = msg
        self.rust_info = None

    def __str__(self):
        rv = self.message
        if self.rust_info is not None:
            return u'%s\n\n%s' % (rv, self.rust_info)
        return rv


class IndexNotSupported(SourmashError):
    def __init__(self):
        SourmashError.__init__(self, "This index format is not supported in this version of sourmash")


def _make_error(error_name, base=SourmashError, code=None):
    class Exc(base):
        pass

    Exc.__name__ = Exc.__qualname__ = error_name
    if code is not None:
        Exc.code = code
    globals()[Exc.__name__] = Exc
    __all__.append(Exc.__name__)
    return Exc


def _get_error_base(error_name):
    pieces = error_name.split("Error", 1)
    if len(pieces) == 2 and pieces[0] and pieces[1]:
        base_error_name = pieces[0] + "Error"
        base_class = globals().get(base_error_name)
        if base_class is None:
            base_class = _make_error(base_error_name)
        return base_class
    return SourmashError


def _make_exceptions():
    for attr in dir(lib):
        if not attr.startswith('SOURMASH_ERROR_CODE_'):
            continue

        code = getattr(lib, attr)
        if code == 1104:
            exceptions_by_code[code] = ValueError
        elif code < 100 or code > 10000:
            error_name = attr[20:].title().replace("_", "")
            base = _get_error_base(error_name)
            exc = _make_error(error_name, base=base, code=getattr(lib, attr))
            exceptions_by_code[exc.code] = exc
        else:
            exceptions_by_code[code] = ValueError

_make_exceptions()
