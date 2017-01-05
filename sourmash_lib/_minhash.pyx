from cython.operator cimport dereference as deref

from libcpp cimport bool
from libc.stdint cimport uint32_t

from _minhash cimport KmerMinHash, KmerMinAbundance, MINHASH_DEFAULT_SEED, _hash_murmur


#cpdef uint64_t hash_murmur(string kmer, uint32_t seed=MINHASH_DEFAULT_SEED):
#    "hash_murmur(string, [,seed])\n\n"
#    "Compute a hash for a string, optionally using a seed (an integer). "
#    "The current default seed is returned by hash_seed()."
#
#    return _hash_murmur(kmer, seed)


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
        self._this.add_sequence(sequence.encode('utf-8'), force)

    def __len__(self):
        return self._this.num

    cpdef get_mins(self, bool with_abundance=False):
        if with_abundance and self.track_abundance:
            return (<KmerMinAbundance*>self._this).mins
        elif self.track_abundance:
            return (<KmerMinAbundance*>self._this).mins.keys()
        else:
            return self._this.mins

    def is_protein(self):
        return self._this.is_protein

    cpdef add_hash(self, uint64_t h):
        self._this.add_hash(h)

    def count_common(self, MinHash other):
        # TODO: add exception handling
        # TODO: check abundance
        return self._this.count_common(deref(other._this))

    def compare(self, MinHash other):
        # TODO: add exception handling
        n = self.count_common(other)
        #size = self._this.mins.size()
        size = self._this.num
        return n / size

    cpdef set_abundances(self, dict values):
        self._this.mins = values

    #cdef add_protein()
    #cdef get_mins()
    #cdef __copy__()
    #cdef merge()
    #cdef is_protein()
