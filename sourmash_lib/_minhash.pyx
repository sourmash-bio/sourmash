# -*- coding: UTF-8 -*-
# cython: language_level=3, c_string_type=str, c_string_encoding=ascii

from __future__ import unicode_literals

from cython.operator cimport dereference as deref, address

from libcpp cimport bool
from libc.stdint cimport uint32_t

from ._minhash cimport KmerMinHash, KmerMinAbundance, _hash_murmur
import math
import copy


# default MurmurHash seed
cdef uint32_t MINHASH_DEFAULT_SEED = 42


def get_minhash_default_seed():
    return MINHASH_DEFAULT_SEED


# we use the 64-bit hash space of MurmurHash only
cdef uint64_t MINHASH_MAX_HASH = 2**64 - 1


def get_minhash_max_hash():
    return MINHASH_MAX_HASH


def get_max_hash_for_scaled(scaled):
    if scaled == 0:
        return 0
    elif scaled == 1:
        return get_minhash_max_hash()

    return int(round(get_minhash_max_hash() / scaled, 0))


def get_scaled_for_max_hash(max_hash):
    if max_hash == 0:
        return 0
    return int(round(get_minhash_max_hash() / max_hash, 0))


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


def dotproduct(a, b, normalize=True):
    """
    Compute the dot product of two dictionaries {k: v} where v is
    abundance.
    """

    if normalize:
        norm_a = math.sqrt(sum([ x*x for x in a.values() ]))
        norm_b = math.sqrt(sum([ x*x for x in b.values() ]))

        if norm_a == 0.0 or norm_b == 0.0:
            return 0.0
    else:
        norm_a = 1.0
        norm_b = 1.0

    prod = 0.
    for k, abundance in a.items():
        prod += (float(abundance) / norm_a) * (b.get(k, 0) / norm_b)

    return prod


cdef class MinHash(object):

    def __init__(self, unsigned int n, unsigned int ksize,
                       bool is_protein=False,
                       bool track_abundance=False,
                       uint32_t seed=MINHASH_DEFAULT_SEED,
                       HashIntoType max_hash=0,
                       mins=None, HashIntoType scaled=0):
        self.track_abundance = track_abundance

        if max_hash and scaled:
            raise ValueError('cannot set both max_hash and scaled')
        elif scaled:
            max_hash = get_max_hash_for_scaled(scaled)

        if max_hash and n:
            raise ValueError('cannot set both n and max_hash')

        cdef KmerMinHash *mh = NULL
        if track_abundance:
            mh = new KmerMinAbundance(n, ksize, is_protein, seed, max_hash)
        else:
            mh = new KmerMinHash(n, ksize, is_protein, seed, max_hash)

        self._this.reset(mh)

        if mins:
            if track_abundance:
                self.set_abundances(mins)
            else:
                self.add_many(mins)


    def __copy__(self):
        a = MinHash(deref(self._this).num, deref(self._this).ksize,
                    deref(self._this).is_protein, self.track_abundance,
                    deref(self._this).seed, deref(self._this).max_hash)
        a.merge(self)
        return a

    def __getstate__(self):             # enable pickling
        with_abundance = False
        if self.track_abundance:
            with_abundance = True

        return (deref(self._this).num,
                deref(self._this).ksize,
                deref(self._this).is_protein,
                self.get_mins(with_abundance=with_abundance),
                None, self.track_abundance, deref(self._this).max_hash,
                deref(self._this).seed)

    def __setstate__(self, tup):
        (n, ksize, is_protein, mins, _, track_abundance, max_hash, seed) =\
          tup

        self.track_abundance = track_abundance

        cdef KmerMinHash *mh = NULL
        if track_abundance:
            mh = new KmerMinAbundance(n, ksize, is_protein, seed, max_hash)
            self._this.reset(mh)
            self.set_abundances(mins)
        else:
            mh = new KmerMinHash(n, ksize, is_protein, seed, max_hash)
            self._this.reset(mh)
            self.add_many(mins)

    def __reduce__(self):
        return (MinHash,
               (deref(self._this).num,
                deref(self._this).ksize,
                deref(self._this).is_protein,
                self.track_abundance,
                deref(self._this).seed,
                deref(self._this).max_hash,
                self.get_mins(with_abundance=self.track_abundance),
                0))

    def __richcmp__(self, other, op):
        if op == 2:
            return self.__getstate__() == other.__getstate__()
        raise Exception("undefined comparison")

    def copy_and_clear(self):
        a = MinHash(deref(self._this).num, deref(self._this).ksize,
                    deref(self._this).is_protein, self.track_abundance,
                    deref(self._this).seed, deref(self._this).max_hash)
        return a

    def add_sequence(self, sequence, bool force=False):
        deref(self._this).add_sequence(to_bytes(sequence), force)

    def add(self, kmer):
        "Add kmer into sketch."
        self.add_sequence(kmer)

    def add_many(self, hashes):
        "Add many hashes in at once."
        for hash in hashes:
            self.add_hash(hash)

    def update(self, other):
        "Update this estimator from all the hashes from the other."
        self.add_many(other.get_mins())

    def __len__(self):
        return deref(self._this).mins.size()

    cpdef get_mins(self, bool with_abundance=False):
        cdef KmerMinAbundance *mh = <KmerMinAbundance*>address(deref(self._this))
        if with_abundance and self.track_abundance:
            return dict(zip(mh.mins, mh.abunds))
        else:
            return [it for it in sorted(deref(self._this).mins)]

    def get_hashes(self):
        return self.get_mins()

    def subtract_mins(self, other):
        a = set(self.get_mins())
        b = set(other.get_mins())
        return a - b

    @property
    def seed(self):
        return deref(self._this).seed

    @property
    def num(self):
        return deref(self._this).num

    @property
    def scaled(self):
        if self.max_hash:
            return get_scaled_for_max_hash(self.max_hash)
        return 0

    @property
    def is_protein(self):
        return deref(self._this).is_protein

    @property
    def ksize(self):
        return deref(self._this).ksize

    @property
    def max_hash(self):
        mm = deref(self._this).max_hash

        return mm

    def add_hash(self, uint64_t h):
        deref(self._this).add_hash(h)

    def count_common(self, MinHash other):
        return deref(self._this).count_common(deref(other._this))

    def downsample_n(self, new_num):
        if self.num and self.num < new_num:
            raise ValueError('new sample n is higher than current sample n')

        a = MinHash(new_num, deref(self._this).ksize,
                    deref(self._this).is_protein, self.track_abundance,
                    deref(self._this).seed, 0)
        if self.track_abundance:
            a.set_abundances(self.get_mins(with_abundance=True))
        else:
            a.add_many(self.get_mins())

        return a

    def downsample_max_hash(self, *others):
        max_hashes = [ x.max_hash for x in others ]
        new_max_hash = min(self.max_hash, *max_hashes)
        new_scaled = get_scaled_for_max_hash(new_max_hash)

        return self.downsample_scaled(new_scaled)

    def downsample_scaled(self, new_num):
        if self.num:
            raise ValueError('num != 0 - cannot downsample a standard MinHash')

        max_hash = self.max_hash
        if max_hash is None:
            raise ValueError('no max_hash available - cannot downsample')

        old_scaled = get_scaled_for_max_hash(self.max_hash)
        if old_scaled > new_num:
            raise ValueError('new scaled is lower than current sample scaled')

        new_max_hash = get_max_hash_for_scaled(new_num)

        a = MinHash(0, deref(self._this).ksize,
                    deref(self._this).is_protein, self.track_abundance,
                    deref(self._this).seed, new_max_hash)
        if self.track_abundance:
            a.set_abundances(self.get_mins(with_abundance=True))
        else:
            a.add_many(self.get_mins())

        return a

    def intersection(self, MinHash other):
        if self.num != other.num:
            err = 'must have same num: {} != {}'.format(self.num,
                                                            other.num)
            raise TypeError(err)
        else:
            num = self.num

        if self.track_abundance and other.track_abundance:
            combined_mh = new KmerMinAbundance(num,
                                          deref(self._this).ksize,
                                          deref(self._this).is_protein,
                                          deref(self._this).seed,
                                          deref(self._this).max_hash)

        else:
            combined_mh = new KmerMinHash(num,
                                          deref(self._this).ksize,
                                          deref(self._this).is_protein,
                                          deref(self._this).seed,
                                          deref(self._this).max_hash)

        combined_mh.merge(deref(self._this))
        combined_mh.merge(deref(other._this))

        common = set(self.get_mins())
        common.intersection_update(other.get_mins())
        common.intersection_update(combined_mh.mins)

        return common, max(combined_mh.size(), 1)

    def compare(self, MinHash other):
        common, size = self.intersection(other)
        n = len(common)
        return n / size

    def jaccard(self, MinHash other):
        return self.compare(other)

    def similarity(self, other, ignore_abundance=False):
        """\
        Calculate similarity of two sketches.

        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate a distance metric
        based on the cosine similarity.

        Note, because the term frequencies (tf-idf weights) cannot be negative,
        the angle will never be < 0deg or > 90deg.

        See https://en.wikipedia.org/wiki/Cosine_similarity
        """

        # if either signature is flat, calculate Jaccard only.
        if not (self.track_abundance and other.track_abundance) or \
          ignore_abundance:
            return self.jaccard(other)
        else:
            # can we merge? if not, raise exception.
            aa = copy.copy(self)
            aa.merge(other)

            a = self.get_mins(with_abundance=True)
            b = other.get_mins(with_abundance=True)

            prod = dotproduct(a, b)
            prod = min(1.0, prod)

            distance = 2*math.acos(prod) / math.pi
            return 1.0 - distance

    def contained_by(self, other):
        """\
        Calculate how much of self is contained by other.
        """
        return self.count_common(other) / len(self.get_mins())

    def similarity_ignore_maxhash(self, MinHash other):
        a = set(self.get_mins())
        if not a:
            return 0.0

        b = set(other.get_mins())

        overlap = a.intersection(b)
        return float(len(overlap)) / float(len(a))

    def __iadd__(self, MinHash other):
        cdef KmerMinAbundance *mh = <KmerMinAbundance*>address(deref(self._this))
        cdef KmerMinAbundance *other_mh = <KmerMinAbundance*>address(deref(other._this))

        if self.track_abundance and other.track_abundance:
            deref(mh).merge(deref(other_mh))
        else:
            deref(self._this).merge(deref(other._this))

        return self
    merge = __iadd__

    cpdef set_abundances(self, dict values):
        if self.track_abundance:
            added = 0

            for k, v in sorted(values.items()):
                if not self.max_hash or k <= self.max_hash:
                    deref(self._this).mins.push_back(k)
                    (<KmerMinAbundance*>address(deref(self._this))).abunds.push_back(v)
                    added += 1
                    if self.num > 0 and added >= self.num:
                        break
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

    def is_molecule_type(self, molecule):
        if molecule.upper() == 'DNA' and not self.is_protein:
            return True
        if molecule == 'protein' and self.is_protein:
            return True
        return False
