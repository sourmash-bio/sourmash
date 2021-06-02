__all__ = ['lib', 'ffi']

import os
from sourmash._lowlevel__ffi import ffi

lib = ffi.dlopen(os.path.join(os.path.dirname(__file__), '_lowlevel__lib.abi3.so'), 2)
del os
