# -*- coding: UTF-8 -*-

import sys
from tempfile import NamedTemporaryFile

from ._lowlevel import ffi, lib
from .utils import RustObject, rustcall, decode_str
from .exceptions import SourmashError
from .minhash import to_bytes, MinHash


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

    def update(self, other):
        if isinstance(other, HLL):
            return self._methodcall(lib.hll_merge, other._objptr)
        elif isinstance(other, MinHash):
            return self._methodcall(lib.hll_update_mh, other._objptr)
        else:
            # FIXME: we could take sets here too (or anything that can be
            # converted to a list of ints...)
            raise TypeError("Must be a HyperLogLog or MinHash")

    def similarity(self, other):
        if isinstance(other, HLL):
            return self._methodcall(lib.hll_similarity, other._objptr)
        else:
            # FIXME: we could take sets here too (or anything that can be
            # converted to a list of ints...)
            raise TypeError("other must be a HyperLogLog")

    def containment(self, other):
        if isinstance(other, HLL):
            return self._methodcall(lib.hll_containment, other._objptr)
        else:
            # FIXME: we could take sets here too (or anything that can be
            # converted to a list of ints...)
            raise TypeError("other must be a HyperLogLog")

    def intersection(self, other):
        if isinstance(other, HLL):
            return self._methodcall(lib.hll_intersection_size, other._objptr)
        else:
            # FIXME: we could take sets here too (or anything that can be
            # converted to a list of ints...)
            raise TypeError("other must be a HyperLogLog")

    @staticmethod
    def load(filename):
        hll_ptr = rustcall(lib.hll_from_path, to_bytes(filename))
        return HLL._from_objptr(hll_ptr)

    @staticmethod
    def from_buffer(buf):
        hll_ptr = rustcall(lib.hll_from_buffer, buf, len(buf))
        return HLL._from_objptr(hll_ptr)

    def save(self, filename):
        self._methodcall(lib.hll_save, to_bytes(filename))

    def to_bytes(self, compression=1):
        size = ffi.new("uintptr_t *")
        rawbuf = self._methodcall(lib.hll_to_buffer, size)
        size = size[0]

        rawbuf = ffi.gc(rawbuf, lambda o: lib.nodegraph_buffer_free(o, size), size)
        buf = ffi.buffer(rawbuf, size)

        return buf

    def count(self, h):
        self.add(h)

    def get(self, h):
        raise NotImplementedError("HLL doesn't support membership query")

    def matches(self, mh):
        if not isinstance(mh, MinHash):
            # FIXME: we could take sets here too (or anything that can be
            # converted to a list of ints...)
            raise ValueError("mh must be a MinHash")

        return self._methodcall(lib.hll_matches, mh._objptr)
