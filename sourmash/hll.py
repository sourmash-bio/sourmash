# -*- coding: UTF-8 -*-

import sys
from tempfile import NamedTemporaryFile

from ._lowlevel import ffi, lib
from .utils import RustObject, rustcall, decode_str
from .exceptions import SourmashError


class HLL(RustObject):
    __dealloc_func__ = lib.hll_free

    def __init__(self, error_rate, ksize):
        self._objptr = lib.hll_with_error_rate(error_rate, ksize)

    def ksize(self):
        return self._methodcall(lib.hll_ksize)
