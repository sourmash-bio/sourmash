from ._compat import implements_to_string
from ._lowlevel import lib


__all__ = ['SourmashError']
exceptions_by_code = {}


@implements_to_string
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


def _make_exceptions():
    for attr in dir(lib):
        if not attr.startswith('SOURMASH_ERROR_CODE_'):
            continue

        class Exc(SourmashError):
            pass

        code = getattr(lib, attr)
        if code < 100 or code > 10000:
            Exc.__name__ = attr[20:].title().replace('_', '')
            Exc.code = getattr(lib, attr)
            globals()[Exc.__name__] = Exc
            Exc.code = code
            exceptions_by_code[code] = Exc
            __all__.append(Exc.__name__)
        else:
            exceptions_by_code[code] = ValueError

_make_exceptions()
