# -*- coding: UTF-8 -*-
# cython: language_level=3, c_string_type=str, c_string_encoding=ascii

from __future__ import unicode_literals

from cython.operator cimport dereference as deref, address

from libcpp cimport bool
from libc.stdint cimport uint32_t

from ._minhash cimport KmerMinHash, KmerMinAbundance, _hash_murmur


cdef uint32_t MINHASH_DEFAULT_SEED = 42


cdef bytes to_bytes(s):
    if not isinstance(s, (basestring, bytes)):
        raise TypeError("Requires a string-like sequence")

    if isinstance(s, unicode):
        s = s.encode('utf-8')
    return s


def hash_murmur(kmer, uint32_t seed=MINHASH_DEFAULT_SEED):
    "hash_murmur(string, [,seed])\n\n"
    "Compute a hash for a string, optionally using a seed (an integer). "
    "The current default seed is returned by hash_seed()."

    return _hash_murmur(to_bytes(kmer), seed)


cdef class MinHash(object):

    def __init__(self, unsigned int n, unsigned int ksize,
                       bool is_protein=False,
                       bool track_abundance=False,
                       uint32_t seed=MINHASH_DEFAULT_SEED,
                       HashIntoType max_hash=0):
        self.track_abundance = track_abundance

        cdef KmerMinHash *mh = NULL
        if track_abundance:
            mh = new KmerMinAbundance(n, ksize, is_protein, seed, max_hash)
        else:
            mh = new KmerMinHash(n, ksize, is_protein, seed, max_hash)

        self._this.reset(mh)

    def __copy__(self):
        a = MinHash(deref(self._this).num, deref(self._this).ksize,
                    deref(self._this).is_protein, self.track_abundance,
                    deref(self._this).seed, deref(self._this).max_hash)
        a.merge(self)
        return a

    def add_sequence(self, sequence, bool force=False):
        deref(self._this).add_sequence(to_bytes(sequence), force)

    def __len__(self):
        return deref(self._this).num

    cpdef get_mins(self, bool with_abundance=False):
        cdef KmerMinAbundance *mh = <KmerMinAbundance*>address(deref(self._this))
        if with_abundance and self.track_abundance:
            return mh.mins
        elif self.track_abundance:
            return [it.first for it in mh.mins]
        else:
            return [it for it in deref(self._this).mins]

    @property
    def seed(self):
        return deref(self._this).seed

    def is_protein(self):
        return deref(self._this).is_protein

    def add_hash(self, uint64_t h):
        deref(self._this).add_hash(h)

    def count_common(self, MinHash other):
        cdef KmerMinAbundance *mh = NULL
        cdef KmerMinAbundance *other_mh = NULL
        cdef uint64_t n = 0

        if self.track_abundance:
            mh = <KmerMinAbundance*>address(deref(self._this))
            if other.track_abundance:
                other_mh = <KmerMinAbundance*>address(deref(other._this))
                n = mh.count_common(deref(other_mh))
            else:
                n = mh.count_common(deref(other._this))
        else:
            if other.track_abundance:
                other_mh = <KmerMinAbundance*>address(deref(other._this))
                n = other_mh.count_common(deref(self._this))
            else:
                n = deref(self._this).count_common(deref(other._this))

        return n

    def compare(self, MinHash other):
        n = self.count_common(other)

        combined = MinHash(deref(self._this).num, deref(self._this).ksize,
                           deref(self._this).is_protein, self.track_abundance,
                           deref(self._this).seed, deref(self._this).max_hash)
        combined += self
        combined += other

        size = max(deref(combined._this).size(), 1)
#        size = max(deref(self._this).size(), 1)
        return n / size

    def __iadd__(self, MinHash other):
        cdef KmerMinAbundance *mh = <KmerMinAbundance*>address(deref(self._this))
        cdef KmerMinAbundance *other_mh = <KmerMinAbundance*>address(deref(other._this))
        if self.track_abundance:
             mh.merge(deref(other_mh))
        else:
            deref(self._this).merge(deref(other._this))

        return self
    merge = __iadd__

    cpdef set_abundances(self, dict values):
        if self.track_abundance:
            for k, v in values.items():
                (<KmerMinAbundance*>address(deref(self._this))).mins[k] = v
        else:
            raise RuntimeError("Use track_abundance=True when constructing "
                               "the MinHash to use set_abundances.")

    def add_protein(self, sequence):
        cdef uint32_t ksize = deref(self._this).ksize // 3
        if len(sequence) < ksize:
            return

        if not deref(self._this).is_protein:
            raise ValueError("cannot add amino acid sequence to DNA MinHash!")

        for i in range(0, len(sequence) - ksize + 1):
            deref(self._this).add_word(to_bytes(sequence[i:i + ksize]))
