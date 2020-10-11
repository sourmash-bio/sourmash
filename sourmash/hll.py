# -*- coding: UTF-8 -*-

import sys
from tempfile import NamedTemporaryFile

from ._lowlevel import ffi, lib
from .utils import RustObject, rustcall, decode_str
from .exceptions import SourmashError
from .minhash import to_bytes


class HLL(RustObject):
    __dealloc_func__ = lib.hll_free

    def __init__(self, error_rate, ksize):
        self._objptr = lib.hll_with_error_rate(error_rate, ksize)

    def __len__(self):
        return self.cardinality()

    def cardinality(self):
        return self._methodcall(lib.hll_cardinality)

    @property
    def ksize(self):
        return self._methodcall(lib.hll_ksize)

    def add_sequence(self, sequence, force=False):
        "Add a sequence into the sketch."
        self._methodcall(lib.hll_add_sequence, to_bytes(sequence), len(sequence), force)

    def add_kmer(self, kmer):
        "Add a kmer into the sketch."
        if len(kmer) != self.ksize:
            raise ValueError("kmer to add is not {} in length".format(self.ksize))
        self.add_sequence(kmer)

    def add(self, h):
        if isinstance(h, str):
            return self.add_kmer(h)
        return self._methodcall(lib.hll_add_hash, h)
