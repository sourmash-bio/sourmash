import re
import pathlib

import cffi

_directive_re = re.compile(r'^\s*#.*?$(?m)')

header = (pathlib.Path(__file__).parent.parent.parent / pathlib.Path("include/sourmash.h")).read_text()
header = _directive_re.sub('', header)

ffibuilder = cffi.FFI()
ffibuilder.set_source("sourmash._lowlevel__ffi", None)
ffibuilder.cdef(header)

if __name__ == "__main__":
    ffibuilder.compile()
