# -*- coding: UTF-8 -*-
# cython: c_string_type=str, c_string_encoding=ascii

from __future__ import unicode_literals

from cython.operator cimport dereference as deref

from libcpp cimport bool
from libc.stdint cimport uint32_t

from _minhash cimport KmerMinHash, KmerMinAbundance, _hash_murmur

cdef uint32_t MINHASH_DEFAULT_SEED = 42

def hash_murmur(str kmer, uint32_t seed=MINHASH_DEFAULT_SEED):
    "hash_murmur(string, [,seed])\n\n"
    "Compute a hash for a string, optionally using a seed (an integer). "
    "The current default seed is returned by hash_seed()."

    return _hash_murmur(kmer, seed)


cdef class MinHash(object):

    def __cinit__(self, unsigned int n, unsigned int ksize,
                        bool is_protein=False,
                        bool track_abundance=False,
                        uint32_t seed=MINHASH_DEFAULT_SEED):
        self.track_abundance = track_abundance
        if track_abundance:
            self._this = new KmerMinAbundance(n, ksize, is_protein, seed)
        else:
            self._this = new KmerMinHash(n, ksize, is_protein, seed)


    cpdef add_sequence(self, str sequence, bool force=True):
        # TODO: add exception handling
        self._this.add_sequence(sequence, force)

    def __len__(self):
        return self._this.num

    cpdef get_mins(self, bool with_abundance=False):
        if with_abundance and self.track_abundance:
            return dict((<KmerMinAbundance*>self._this).mins)
        elif self.track_abundance:
            return list(sorted((<KmerMinAbundance*>self._this).mins.keys()))
        else:
            return list(sorted(self._this.mins))

    def is_protein(self):
        return self._this.is_protein

    cpdef add_hash(self, uint64_t h):
        self._this.add_hash(h)

    cpdef uint64_t count_common(self, MinHash other):
        # TODO: add exception handling
        cdef KmerMinAbundance *mh = NULL
        cdef KmerMinAbundance *other_mh = NULL
        cdef uint64_t n = 0

        if self.track_abundance:
            mh = <KmerMinAbundance*>(self._this)
            if other.track_abundance:
                other_mh = <KmerMinAbundance*>(other._this)
                n = mh.count_common(deref(other_mh))
            else:
                n = mh.count_common(deref(other._this))
        else:
            if other.track_abundance:
                other_mh = <KmerMinAbundance*>(other._this)
                n = other_mh.count_common(deref(self._this))
            else:
                n = self._this.count_common(deref(other._this))

        return n

    def compare(self, MinHash other):
        # TODO: add exception handling
        n = self.count_common(other)
        if self.track_abundance:
            size = (<KmerMinAbundance*>self._this).mins.size()
        else:
            size = self._this.mins.size()
        return n / size

    cpdef set_abundances(self, dict values):
        # TODO: need to check if track_abundance is True
        for k, v in values.items():
            (<KmerMinAbundance*>self._this).mins[k] = v

    cpdef add_protein(self, str sequence):
        cdef uint32_t ksize = self._this.ksize / 3
        if len(sequence) < ksize:
            return

        if not self._this.is_protein:
            raise ValueError("cannot add amino acid sequence to DNA MinHash!")

        for i in range(0, len(sequence) - ksize + 1):
            self._this.add_word(sequence[i:i + ksize])

    #cdef add_protein()
    #cdef __copy__()
    #cdef merge()
