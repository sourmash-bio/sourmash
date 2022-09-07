# -*- coding: utf-8 -*-
"""
sourmash submodule that provides MinHash class and utility functions.

class MinHash - core MinHash class.
class FrozenMinHash - read-only MinHash class.
"""
from __future__ import unicode_literals, division
from .distance_utils import jaccard_to_distance, containment_to_distance, set_size_exact_prob
from .logging import notify

import numpy as np


__all__ = ['get_minhash_default_seed',
           'get_minhash_max_hash',
           'hash_murmur',
           'MinHash',
           'FrozenMinHash']

from collections.abc import Mapping

from . import VERSION
from ._lowlevel import ffi, lib
from .utils import RustObject, rustcall
from .exceptions import SourmashError
from deprecation import deprecated

# default MurmurHash seed
MINHASH_DEFAULT_SEED = 42


def get_minhash_default_seed():
    "Return the default seed value used for the MurmurHash hashing function."
    return MINHASH_DEFAULT_SEED


# we use the 64-bit hash space of MurmurHash only
# this is 2 ** 64 - 1 in hexadecimal
MINHASH_MAX_HASH = 0xFFFFFFFFFFFFFFFF


def get_minhash_max_hash():
    "Return the maximum hash value."
    return MINHASH_MAX_HASH


def _get_max_hash_for_scaled(scaled):
    "Convert a 'scaled' value into a 'max_hash' value."
    if scaled == 0:
        return 0
    elif scaled == 1:
        return get_minhash_max_hash()

    return min(
        int(round(get_minhash_max_hash() / scaled, 0)),
        MINHASH_MAX_HASH
    )


def _get_scaled_for_max_hash(max_hash):
    "Convert a 'max_hash' value into a 'scaled' value."
    if max_hash == 0:
        return 0
    return min(
        int(round(get_minhash_max_hash() / max_hash, 0)),
        MINHASH_MAX_HASH
    )


def to_bytes(s):
    # Allow for strings, bytes or int
    # Single item of byte string = int

    if isinstance(s, bytes):
        return s

    if not isinstance(s, (str, bytes, int)):
        raise TypeError("Requires a string-like sequence")

    if isinstance(s, str):
        s = s.encode("utf-8")
    elif isinstance(s, int):
        s = bytes([s])

    return s


def hash_murmur(kmer, seed=MINHASH_DEFAULT_SEED):
    "hash_murmur(string, [,seed])\n\n"
    "Compute a hash for a string, optionally using a seed (an integer). "
    "The current default seed is returned by hash_seed()."

    return lib.hash_murmur(to_bytes(kmer), seed)


def translate_codon(codon):
    "Translate a codon into an amino acid."
    try:
        return rustcall(lib.sourmash_translate_codon,
                        to_bytes(codon)).decode('utf-8')
    except SourmashError as e:
        raise ValueError(e.message)


def flatten_and_downsample_scaled(mh, *scaled_vals):
    "Flatten MinHash object and downsample to max of scaled values."
    assert mh.scaled
    assert all( (x > 0 for x in scaled_vals) )

    mh = mh.flatten()
    scaled = max(scaled_vals)
    if scaled > mh.scaled:
        return mh.downsample(scaled=scaled)
    return mh


def flatten_and_downsample_num(mh, *num_vals):
    "Flatten MinHash object and downsample to min of num values."
    assert mh.num
    assert all( (x > 0 for x in num_vals) )

    mh = mh.flatten()
    num = min(num_vals)
    if num < mh.num:
        return mh.downsample(num=num)
    return mh


def flatten_and_intersect_scaled(mh1, mh2):
    "Flatten and downsample two scaled MinHash objs, then return intersection."
    scaled = max(mh1.scaled, mh2.scaled)
    mh1 = mh1.flatten().downsample(scaled=scaled)
    mh2 = mh2.flatten().downsample(scaled=scaled)

    return mh1 & mh2


class _HashesWrapper(Mapping):
    "A read-only view of the hashes contained by a MinHash object."
    def __init__(self, h):
        self._data = h

    def __getitem__(self, key):
        return self._data[key]

    def __repr__(self):
        return repr(self._data)

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __eq__(self, other):
        return list(self.items()) == list(other.items())

    def __setitem__(self, k, v):
        raise RuntimeError("cannot modify hashes directly; use 'add' methods")


class MinHash(RustObject):
    """\
    The core sketch object for sourmash.

    MinHash objects store and provide functionality for subsampled hash values
    from DNA, RNA, and amino acid sequences. MinHash also supports both the
    standard MinHash behavior (bounded size or ``num``) and a non-standard
    MinHash, called "modulo hash" behavior, or ``scaled``. Please see
    the API examples at

        https://sourmash.readthedocs.io/en/latest/api-example.html#sourmash-minhash-objects-and-manipulations

    for more information.

    Basic usage:

    >>> from sourmash import MinHash
    >>> mh1 = MinHash(n=20, ksize=3)
    >>> mh1.add_sequence('ATGAGAGACGATAGACAGATGAC')

    >>> mh2 = MinHash(n=20, ksize=3)
    >>> mh2.add_sequence('ATGAGActCGATAGaCAGATGAC')

    >>> round(mh1.similarity(mh2), 2)
    0.85
    """
    __dealloc_func__ = lib.kmerminhash_free

    def __init__(
        self,
        n,
        ksize,
        *,
        is_protein=False,
        dayhoff=False,
        hp=False,
        track_abundance=False,
        seed=MINHASH_DEFAULT_SEED,
        max_hash=0,
        mins=None,
        scaled=0,
    ):
        """\
        Create a sourmash.MinHash object.

        To create a standard (``num``) MinHash, use:
           ``MinHash(<num>, <ksize>, ...)``

        To create a ``scaled`` MinHash, use
            ``MinHash(0, <ksize>, scaled=<int>, ...)``

        Optional arguments:
           * is_protein (default False) - aa k-mers
           * dayhoff (default False) - dayhoff encoding
           * hp (default False) - hydrophilic/hydrophobic aa
           * track_abundance (default False) - track hash multiplicity
           * mins (default None) - list of hashvals, or (hashval, abund) pairs
           * seed (default 42) - murmurhash seed
        """
        # support max_hash in constructor, for now.
        if max_hash:
            if scaled:
                raise ValueError("cannot set both max_hash and scaled")
            scaled = _get_scaled_for_max_hash(max_hash)

        if scaled and n:
            raise ValueError("cannot set both n and max_hash")

        if not n and not scaled:
            raise ValueError("cannot omit both n and scaled")

        if dayhoff or hp:
            is_protein = False

        if dayhoff:
            hash_function = lib.HASH_FUNCTIONS_MURMUR64_DAYHOFF
            ksize = ksize*3
        elif hp:
            hash_function = lib.HASH_FUNCTIONS_MURMUR64_HP
            ksize = ksize*3
        elif is_protein:
            hash_function = lib.HASH_FUNCTIONS_MURMUR64_PROTEIN
            ksize = ksize*3
        else:
            hash_function = lib.HASH_FUNCTIONS_MURMUR64_DNA

        self._objptr = lib.kmerminhash_new(
            scaled, ksize, hash_function, seed, track_abundance, n
        )

        if mins:
            if track_abundance:
                self.set_abundances(mins)
            else:
                self.add_many(mins)

    def __copy__(self):
        "Create a new copy of this MinHash."
        a = MinHash(
            self.num,
            self.ksize,
            is_protein=self.is_protein,
            dayhoff=self.dayhoff,
            hp=self.hp,
            track_abundance=self.track_abundance,
            seed=self.seed,
            max_hash=self._max_hash,
        )
        a.merge(self)
        return a

    copy = __copy__

    def __getstate__(self):
        "support pickling via __getstate__/__setstate__"

        # note: we multiple ksize by 3 here so that
        # pickle protocols that bypass __setstate__ <coff numpy coff>
        # get a ksize that makes sense to the Rust layer. See #2262.
        return (
            self.num,
            self.ksize if self.is_dna else self.ksize*3,
            self.is_protein,
            self.dayhoff,
            self.hp,
            self.hashes,
            None,
            self.track_abundance,
            self._max_hash,
            self.seed,
        )

    def __setstate__(self, tup):
        "support pickling via __getstate__/__setstate__"
        (n, ksize, is_protein, dayhoff, hp, mins, _, track_abundance,
         max_hash, seed) = tup

        self.__del__()

        hash_function = (
            lib.HASH_FUNCTIONS_MURMUR64_DAYHOFF if dayhoff else
            lib.HASH_FUNCTIONS_MURMUR64_HP if hp else
            lib.HASH_FUNCTIONS_MURMUR64_PROTEIN if is_protein else
            lib.HASH_FUNCTIONS_MURMUR64_DNA
        )

        scaled = _get_scaled_for_max_hash(max_hash)
        self._objptr = lib.kmerminhash_new(
            scaled, ksize, hash_function, seed, track_abundance, n
        )
        if track_abundance:
            self.set_abundances(mins)
        else:
            self.add_many(mins)

    def __eq__(self, other):
        "equality testing via =="
        return self.__getstate__() == other.__getstate__()

    def copy_and_clear(self):
        "Create an empty copy of this MinHash."
        a = MinHash(
            self.num,
            self.ksize,
            is_protein=self.is_protein,
            dayhoff=self.dayhoff,
            hp=self.hp,
            track_abundance=self.track_abundance,
            seed=self.seed,
            max_hash=self._max_hash,
        )
        return a

    def add_sequence(self, sequence, force=False):
        "Add a sequence into the sketch."
        self._methodcall(lib.kmerminhash_add_sequence, to_bytes(sequence),
                         force)

    def seq_to_hashes(self, sequence, *, force=False, bad_kmers_as_zeroes=False, is_protein=False):
        """Convert sequence to hashes without adding to the sketch.

        If input sequence is DNA and this is a protein, dayhoff, or hp
        MinHash, translate the DNA appropriately before hashing.

        If input sequence is protein, set is_protein=True.

        If `force = True` and `bad_kmers_as_zeroes = True`,
        invalid kmers hashes will be represented as `0`.
        """

        if is_protein and self.moltype not in ("protein", "dayhoff", "hp"):
            raise ValueError("cannot add protein sequence to DNA MinHash")

        if bad_kmers_as_zeroes and not force:
            raise ValueError("cannot represent invalid kmers as 0 while force is not set to True")

        size = ffi.new("uintptr_t *")
        hashes_ptr = self._methodcall(lib.kmerminhash_seq_to_hashes, to_bytes(sequence), len(sequence), force, bad_kmers_as_zeroes, is_protein, size)
        size = size[0]

        try:
            return ffi.unpack(hashes_ptr, size)

        finally:
            lib.kmerminhash_slice_free(hashes_ptr, size)

    def kmers_and_hashes(self, sequence, *, force=False, is_protein=False):
        """Convert sequence into (k-mer, hashval) tuples without adding
        it to the sketch.

        If input sequence is DNA and this is a protein, dayhoff, or hp
        MinHash, translate the DNA appropriately before hashing.

        If input sequence is protein, set is_protein=True.

        If 'force' is True, invalid k-mers will be represented with 'None'.
        """
        import screed

        bad_kmers_as_zeroes = False
        if force:
            bad_kmers_as_zeroes = True

        sequence = sequence.upper()
        hashvals = self.seq_to_hashes(sequence,
                                      force=force, is_protein=is_protein,
                                      bad_kmers_as_zeroes=bad_kmers_as_zeroes)

        if bad_kmers_as_zeroes:
            hashvals = [ None if h == 0 else h for h in hashvals ]

        ksize = self.ksize
        translate = False
        if self.moltype == 'DNA':
            pass
        elif is_protein:
            pass
        else:                   # translate input DNA sequence => aa
            assert self.moltype in ('protein', 'dayhoff', 'hp')
            translate = True
            ksize = self.ksize * 3

        # special code for translation -
        if translate:
            # forward AND reverse complement => twice the k-mers
            n_kmers = (len(sequence) - ksize + 1) * 2
            assert n_kmers == len(hashvals)

            # generate reverse complement of sequence
            seqrc = screed.rc(sequence)

            hash_i = 0
            for frame in (0, 1, 2):
                # get forward k-mers
                for start in range(0, len(sequence) - ksize + 1 - frame, 3):
                    kmer = sequence[start + frame:start + frame + ksize]
                    yield kmer, hashvals[hash_i]
                    hash_i += 1

                # get rc k-mers
                for start in range(0, len(seqrc) - ksize + 1 - frame, 3):
                    kmer = seqrc[start + frame:start + frame + ksize]
                    yield kmer, hashvals[hash_i]
                    hash_i += 1
        else:
            # otherwise, all very straightforward :)
            n_kmers = len(sequence) - ksize + 1
            assert n_kmers == len(hashvals)
            for i, hashval in zip(range(0, n_kmers), hashvals):
                kmer = sequence[i:i+ksize]
                yield kmer, hashval

    def add_kmer(self, kmer):
        "Add a kmer into the sketch."
        if self.is_dna:
            if len(kmer) != self.ksize:
                raise ValueError("kmer to add is not {} in length".format(self.ksize))
        else:
            if len(kmer) != self.ksize*3:
                raise ValueError("kmer to add is not {} in length".format(self.ksize*3))
        self.add_sequence(kmer)

    def add_many(self, hashes):
        """Add many hashes to the sketch at once.

        ``hashes`` can be either an iterable (list, set, etc.), or another
        ``MinHash`` object.
        """
        if isinstance(hashes, MinHash):
            self._methodcall(lib.kmerminhash_add_from, hashes._objptr)
        else:
            self._methodcall(lib.kmerminhash_add_many, list(hashes), len(hashes))

    def remove_many(self, hashes):
        """Remove many hashes from a sketch at once.

        ``hashes`` can be either an iterable (list, set, etc.), or another
        ``MinHash`` object.
        """
        if isinstance(hashes, MinHash):
            self._methodcall(lib.kmerminhash_remove_from, hashes._objptr)
        else:
            self._methodcall(lib.kmerminhash_remove_many, list(hashes), len(hashes))

    def __len__(self):
        "Number of hashes."
        return self._methodcall(lib.kmerminhash_get_mins_size)

    @deprecated(deprecated_in="3.5", removed_in="5.0",
                current_version=VERSION,
                details='Use .hashes property instead.')
    def get_mins(self, with_abundance=False):
        """Return list of hashes or if ``with_abundance`` a list
        of (hash, abund).
        """
        mins = self.hashes
        if not with_abundance:
            return mins.keys()
        return mins


    @deprecated(deprecated_in="3.5", removed_in="5.0",
                current_version=VERSION,
                details='Use .hashes property instead.')
    def get_hashes(self):
        "Return the list of hashes."
        return self.hashes.keys()

    @property
    def hashes(self):
        size = ffi.new("uintptr_t *")
        mins_ptr = self._methodcall(lib.kmerminhash_get_mins, size)
        size = size[0]

        try:
            if self.track_abundance:
                size_abunds = ffi.new("uintptr_t *")
                abunds_ptr = self._methodcall(lib.kmerminhash_get_abunds, size_abunds)
                size_abunds = size_abunds[0]
                assert size == size_abunds
                result = dict(zip(ffi.unpack(mins_ptr, size), ffi.unpack(abunds_ptr, size)))
                lib.kmerminhash_slice_free(abunds_ptr, size)
                return _HashesWrapper(result)
            else:
                d = ffi.unpack(mins_ptr, size)
                return _HashesWrapper({ k : 1 for k in d })

        finally:
            lib.kmerminhash_slice_free(mins_ptr, size)


    @property
    def seed(self):
        return self._methodcall(lib.kmerminhash_seed)

    @property
    def num(self):
        return self._methodcall(lib.kmerminhash_num)

    @property
    def scaled(self):
        mx = self._methodcall(lib.kmerminhash_max_hash)
        if mx:
            return _get_scaled_for_max_hash(mx)
        return 0

    @property
    def is_dna(self):
        return not (self.is_protein or self.dayhoff or self.hp)

    @property
    def is_protein(self):
        return self._methodcall(lib.kmerminhash_is_protein)

    @property
    def dayhoff(self):
        return self._methodcall(lib.kmerminhash_dayhoff)

    @property
    def hp(self):
        return self._methodcall(lib.kmerminhash_hp)

    @property
    def ksize(self):
        k = self._methodcall(lib.kmerminhash_ksize)
        if not self.is_dna:
            assert k % 3 == 0
            k = int(k / 3)
        return k

    @property
    @deprecated(deprecated_in="3.5", removed_in="5.0",
                current_version=VERSION,
                details='Use scaled instead.')
    def max_hash(self):
        return self._methodcall(lib.kmerminhash_max_hash)

    # a non-deprecated `max_hash` property for internal testing purposes only
    @property
    def _max_hash(self):
        return self._methodcall(lib.kmerminhash_max_hash)

    @property
    def track_abundance(self):
        return self._methodcall(lib.kmerminhash_track_abundance)

    @track_abundance.setter
    def track_abundance(self, b):
        if self.track_abundance == b:
            return

        if b is False:
            self._methodcall(lib.kmerminhash_disable_abundance)
        elif len(self) > 0:
            raise RuntimeError("Can only set track_abundance=True if the MinHash is empty")
        else:
            self._methodcall(lib.kmerminhash_enable_abundance)

    def add_hash(self, h):
        "Add a single hash value."
        return self._methodcall(lib.kmerminhash_add_hash, h)

    def add_hash_with_abundance(self, h, a):
        "Add a single hash value with an abundance."
        if self.track_abundance:
            return self._methodcall(lib.kmerminhash_add_hash_with_abundance, h, a)
        else:
            raise RuntimeError(
                "Use track_abundance=True when constructing "
                "the MinHash to use add_hash_with_abundance."
            )

    def clear(self):
        "Clears all hashes and abundances."
        return self._methodcall(lib.kmerminhash_clear)

    def count_common(self, other, downsample=False):
        """\
        Return the number of hashes in common between ``self`` and ``other``.

        Optionally downsample ``scaled`` objects to highest ``scaled`` value.
        """
        if not isinstance(other, MinHash):
            raise TypeError("Must be a MinHash!")
        return self._methodcall(lib.kmerminhash_count_common, other._get_objptr(), downsample)

    def intersection_and_union_size(self, other):
        "Calculate intersection and union sizes between `self` and `other`."
        if not isinstance(other, MinHash):
            raise TypeError("Must be a MinHash!")
        if not self.is_compatible(other):
            raise TypeError("incompatible MinHash objects")

        usize = ffi.new("uint64_t *")
        common = self._methodcall(lib.kmerminhash_intersection_union_size,
                                  other._get_objptr(), usize)

        usize = ffi.unpack(usize, 1)[0]
        return common, usize

    def downsample(self, *, num=None, scaled=None):
        """Copy this object and downsample new object to either `num` or
        `scaled`.
        """
        # first, evaluate provided parameters --

        # at least one must be specified!
        if num is None and scaled is None:
            raise ValueError('must specify either num or scaled to downsample')

        # both cannot be specified
        if num is not None and scaled is not None:
            raise ValueError('cannot specify both num and scaled')

        if num is not None:
            # cannot downsample a scaled MinHash with num:
            if self.scaled:
                raise ValueError("cannot downsample a scaled MinHash using num")
            # cannot upsample
            if self.num < num:
                raise ValueError("new sample num is higher than current sample num")

            # acceptable num value? make sure to set max_hash to 0.
            max_hash = 0
            
        elif scaled is not None:
            # cannot downsample a num MinHash with scaled
            if self.num:
                raise ValueError("cannot downsample a num MinHash using scaled")
            if self.scaled > scaled:
                raise ValueError(f"new scaled {scaled} is lower than current sample scaled {self.scaled}")

            # acceptable scaled value? reconfigure max_hash, keep num 0.
            max_hash = _get_max_hash_for_scaled(scaled)
            num = 0

        # end checks! create new object:
        a = MinHash(
            num, self.ksize,
            is_protein=self.is_protein, dayhoff=self.dayhoff, hp=self.hp,
            track_abundance=self.track_abundance, seed=self.seed,
            max_hash=max_hash
        )
        # copy over hashes:
        if self.track_abundance:
            a.set_abundances(self.hashes)
        else:
            a.add_many(self)

        return a

    def flatten(self):
        """If track_abundance=True, return a new flattened MinHash."""
        if self.track_abundance:
            # create new object:
            a = MinHash(
                self.num, self.ksize,
                is_protein=self.is_protein, dayhoff=self.dayhoff, hp=self.hp,
                track_abundance=False, seed=self.seed, max_hash=self._max_hash
            )
            a.add_many(self)

            return a
        return self

    def jaccard(self, other, downsample=False):
        "Calculate Jaccard similarity of two MinHash objects."
        if self.num != other.num:
            err = "must have same num: {} != {}".format(self.num, other.num)
            raise TypeError(err)
        return self._methodcall(lib.kmerminhash_similarity, other._get_objptr(), True, downsample)

    def jaccard_ani(self, other, *, downsample=False, jaccard=None, prob_threshold=1e-3,  err_threshold=1e-4):
        "Use jaccard to estimate ANI between two MinHash objects."
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate ANI for scaled MinHashes")
        self_mh = self
        other_mh = other
        scaled = self.scaled
        if downsample:
            scaled = max(self_mh.scaled, other_mh.scaled)
            self_mh = self.downsample(scaled=scaled)
            other_mh = other.downsample(scaled=scaled)
        if jaccard is None:
            jaccard = self_mh.similarity(other_mh, ignore_abundance=True)
        avg_sketch_kmers = (len(self_mh) + len(other_mh))/2
        avg_n_kmers = round(avg_sketch_kmers * scaled)  # would be better if hll estimate - see #1798
        j_aniresult = jaccard_to_distance(jaccard, self_mh.ksize, scaled,
                                          n_unique_kmers=avg_n_kmers,
                                          prob_threshold = prob_threshold,
                                          err_threshold = err_threshold)
        # null out ANI if either mh size estimation is inaccurate
        if not self.size_is_accurate() or not other.size_is_accurate():
            j_aniresult.size_is_inaccurate = True
        return j_aniresult

    def similarity(self, other, ignore_abundance=False, downsample=False):
        """Calculate similarity of two sketches.

        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate the angular
        similarity, a distance metric based on the cosine similarity.

        Note, because the term frequencies (tf-idf weights) cannot be negative,
        the angle will never be < 0deg or > 90deg.

        See https://en.wikipedia.org/wiki/Cosine_similarity
        """
        return self._methodcall(lib.kmerminhash_similarity,
                                other._get_objptr(),
                                ignore_abundance, downsample)

    def angular_similarity(self, other):
        "Calculate the angular similarity."
        if not (self.track_abundance and other.track_abundance):
            raise TypeError("Error: Angular (cosine) similarity requires both sketches to track hash abundance.")
        return self._methodcall(lib.kmerminhash_angular_similarity,
                                other._get_objptr())

    def is_compatible(self, other):
        return self._methodcall(lib.kmerminhash_is_compatible, other._get_objptr())

    def contained_by(self, other, downsample=False):
        """
        Calculate how much of self is contained by other.
        """
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate containment for scaled MinHashes")
        if not len(self):
            return 0.0
        return self.count_common(other, downsample) / len(self)
        # with bias factor
        #return self.count_common(other, downsample) / (len(self) * (1- (1-1/self.scaled)^(len(self)*self.scaled)))


    def contained_by_debiased(self, other, downsample=False):
        """
        Calculate how much of self is contained by other.
        """
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate containment for scaled MinHashes")
        denom = len(self)
        if not denom:
            return 0.0
        total_denom = float(denom * self.scaled) # would be better if hll estimate - see #1798
        bias_factor = 1.0 - (1.0 - 1.0/self.scaled) ** total_denom
        containment = self.count_common(other, downsample) / (denom * bias_factor)
        # debiasing containment can lead to vals outside of 0-1 range!?
        if containment >= 1:
            return 1
        elif containment <= 0:
            return 0
        else:
            return containment


    def containment_ani(self, other, *, downsample=False, containment=None, confidence=0.95, estimate_ci = False, prob_threshold=1e-3):
        "Use self contained by other to estimate ANI between two MinHash objects."
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate ANI for scaled MinHashes")
        self_mh = self
        other_mh = other
        scaled = self.scaled
        if downsample:
            scaled = max(self_mh.scaled, other_mh.scaled)
            self_mh = self.downsample(scaled=scaled)
            other_mh = other.downsample(scaled=scaled)
        containment = self_mh.contained_by_debiased(other_mh) # recalc debiased containment
        #containment = self_mh.contained_by(other_mh) # recalc debiased containment
        n_kmers = len(self_mh) * scaled # would be better if hll estimate - see #1798

        c_aniresult = containment_to_distance(containment, self_mh.ksize, self_mh.scaled,
                                                        n_unique_kmers=n_kmers, confidence=confidence,
                                                        estimate_ci = estimate_ci, prob_threshold=prob_threshold)
        # null out ANI if either mh size estimation is inaccurate
        if not self.size_is_accurate() or not other.size_is_accurate():
            c_aniresult.size_is_inaccurate = True
        return c_aniresult


    def max_containment(self, other, downsample=False):
        """
        Calculate maximum containment.
        """
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate containment for scaled MinHashes")
        min_denom = min((len(self), len(other)))
        if not min_denom:
            return 0.0

        return self.count_common(other, downsample) / min_denom


    def max_containment_debiased(self, other, downsample=False):
        """
        Calculate maximum containment.
        """
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate containment for scaled MinHashes")
        min_denom = min((len(self), len(other)))
        if not min_denom:
            return 0.0
        total_denom =  float(min_denom * self.scaled) # would be better if hll estimate - see #1798
        bias_factor = 1.0 - (1.0 - 1.0/self.scaled) ** total_denom
        max_containment = self.count_common(other, downsample) / (min_denom * bias_factor)
        # debiasing containment can lead to vals outside of 0-1 range!?
        if max_containment >= 1:
            return 1
        elif max_containment <= 0:
            return 0
        else:
            return max_containment

    def max_containment_ani(self, other, *, downsample=False, max_containment=None, confidence=0.95, estimate_ci=False, prob_threshold=1e-3):  
        "Use max_containment to estimate ANI between two MinHash objects."
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate ANI for scaled MinHashes")
        self_mh = self
        other_mh = other
        scaled = self.scaled
        if downsample:
            scaled = max(self_mh.scaled, other_mh.scaled)
            self_mh = self.downsample(scaled=scaled)
            other_mh = other.downsample(scaled=scaled)
        max_containment = self_mh.max_containment_debiased(other_mh) #recalc debiased max containment
        min_n_kmers = min(len(self_mh), len(other_mh))
        n_kmers = min_n_kmers * scaled  # would be better if hll estimate - see #1798

        c_aniresult = containment_to_distance(max_containment, self_mh.ksize, scaled,
                                           n_unique_kmers=n_kmers,confidence=confidence,
                                           estimate_ci = estimate_ci, prob_threshold=prob_threshold)
        # null out ANI if either mh size estimation is inaccurate
        if not self.size_is_accurate() or not other.size_is_accurate():
            c_aniresult.size_is_inaccurate = True
        return c_aniresult

    def avg_containment(self, other, *, downsample=False):
        """
        Calculate average containment.
        Note: this is average of the containments, *not* count_common/ avg_denom
        """
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate containment for scaled MinHashes")

        c1 = self.contained_by(other, downsample)
        c2 = other.contained_by(self, downsample)

        return (c1 + c2)/2

    def avg_containment_debiased(self, other, downsample=False):
        """
        Calculate average containment, debiased.
        Note: this is average of the containments, *not* count_common/ avg_denom
        """
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate containment for scaled MinHashes")

        c1 = self.contained_by_debiased(other, downsample)
        c2 = other.contained_by_debiased(self, downsample)

        return (c1 + c2)/2

    def avg_containment_ani(self, other, *, downsample=False, prob_threshold=1e-3):
        """
        Calculate average containment ANI.
        Note: this is average of the containment ANI's, *not* ANI using count_common/ avg_denom
        """
        if not (self.scaled and other.scaled):
            raise TypeError("Error: can only calculate ANI for scaled MinHashes")
        a1 = self.containment_ani(other, downsample=downsample, prob_threshold=prob_threshold).ani
        a2 = other.containment_ani(self, downsample=downsample, prob_threshold=prob_threshold).ani
        if any([a1 is None, a2 is None]):
            return None
        return (a1 + a2)/2

    def __add__(self, other):
        if not isinstance(other, MinHash):
            raise TypeError("can only add MinHash objects to MinHash objects!")

        if self.num and other.num:
            if self.num != other.num:
                raise TypeError(f"incompatible num values: self={self.num} other={other.num}")

        new_obj = self.to_mutable()
        new_obj += other
        return new_obj
    __or__ = __add__

    def __iadd__(self, other):
        if not isinstance(other, MinHash):
            raise TypeError("can only add MinHash objects to MinHash objects!")
        self._methodcall(lib.kmerminhash_merge, other._get_objptr())
        return self

    def merge(self, other):
        if not isinstance(other, MinHash):
            raise TypeError("can only add MinHash objects to MinHash objects!")
        self._methodcall(lib.kmerminhash_merge, other._get_objptr())

    def intersection(self, other):
        if not isinstance(other, MinHash):
            raise TypeError("can only intersect MinHash objects")
        if self.track_abundance or other.track_abundance:
            raise TypeError("can only intersect flat MinHash objects")

        ptr = self._methodcall(lib.kmerminhash_intersection, other._get_objptr())
        return MinHash._from_objptr(ptr)
    __and__ = intersection

    def set_abundances(self, values, clear=True):
        """Set abundances for hashes from ``values``, where
        ``values[hash] = abund``

        If ``abund`` value is set to zero, the ``hash`` will be removed from the sketch.
        ``abund`` cannot be set to a negative value.
        """
        if self.track_abundance:
            hashes = []
            abunds = []

            for h, v in values.items():
                hashes.append(h)                
                if v < 0:
                    raise ValueError("Abundance cannot be set to a negative value.")
                abunds.append(v)

            self._methodcall(lib.kmerminhash_set_abundances, hashes, abunds, len(hashes), clear)
        else:
            raise RuntimeError(
                "Use track_abundance=True when constructing "
                "the MinHash to use set_abundances."
            )

    def add_protein(self, sequence):
        "Add a protein sequence."
        self._methodcall(lib.kmerminhash_add_protein, to_bytes(sequence))

    @property
    def moltype(self):                    # TODO: test in minhash tests
        if self.is_protein:
            return 'protein'
        elif self.dayhoff:
            return 'dayhoff'
        elif self.hp:
            return 'hp'
        else:
            return 'DNA'

    def to_mutable(self):
        "Return a copy of this MinHash that can be changed."
        return self.__copy__()

    def to_frozen(self):
        "Return a frozen copy of this MinHash that cannot be changed."
        new_mh = self.__copy__()
        new_mh.into_frozen()
        return new_mh

    def into_frozen(self):
        "Freeze this MinHash, preventing any changes."
        self.__class__ = FrozenMinHash

    def inflate(self, from_mh):
        """return a new MinHash object with abundances taken from 'from_mh'

        note that this implicitly does an intersection: hashes that have
        no abundance in 'from_mh' are set to abundance 0 and removed from
        'self'.
        """
        if not self.track_abundance and from_mh.track_abundance:
            orig_abunds = from_mh.hashes
            abunds = { h: orig_abunds.get(h, 0) for h in self.hashes }

            abund_mh = from_mh.copy_and_clear()

            abund_mh.downsample(scaled=self.scaled)
            abund_mh.set_abundances(abunds)

            return abund_mh
        else:
            raise ValueError("inflate operates on a flat MinHash and takes a MinHash object with track_abundance=True")

    @property
    def sum_abundances(self):
        if self.track_abundance:
            return sum(v for v in self.hashes.values())
        return None

    @property
    def mean_abundance(self):
        if self.track_abundance:
            return np.mean(list(self.hashes.values()))
        return None

    @property
    def median_abundance(self):
        if self.track_abundance:
            return np.median(list(self.hashes.values()))
        return None

    @property
    def std_abundance(self):
        if self.track_abundance:
            return np.std(list(self.hashes.values()))
        return None

    @property
    def unique_dataset_hashes(self):
        """
        Approximate total number of hashes (num_hashes *scaled).
        """
        if not self.scaled:
            raise TypeError("can only approximate unique_dataset_hashes for scaled MinHashes")
        # TODO: replace set_size with HLL estimate when that gets implemented
        return len(self) * self.scaled # + (self.ksize - 1) for bp estimation

    def size_is_accurate(self, relative_error=0.20, confidence=0.95):
        """
        Computes the probability that the estimate: sketch_size * scaled deviates from the true
        set_size by more than relative_error. This relies on the fact that the sketch_size
        is binomially distributed with parameters sketch_size and 1/scaled. The two-sided Chernoff
        bounds are used.
        Returns True if probability is greater than or equal to the desired confidence.
        """
        if not self.scaled:
            raise TypeError("Error: can only estimate dataset size for scaled MinHashes")
        if any([not (0 <= relative_error <= 1), not (0 <= confidence <= 1)]):
            raise ValueError("Error: relative error and confidence values must be between 0 and 1.")
        # to do: replace unique_dataset_hashes with HLL estimation when it gets implemented 
        probability = set_size_exact_prob(self.unique_dataset_hashes, self.scaled, relative_error=relative_error)
        return probability >= confidence


class FrozenMinHash(MinHash):
    def add_sequence(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def add_kmer(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def add_many(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def remove_many(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def add_hash(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def add_hash_with_abundance(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def clear(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def set_abundances(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def add_protein(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def downsample(self, *, num=None, scaled=None):
        if scaled and self.scaled == scaled:
            return self
        if num and self.num == num:
            return self

        down_mh = MinHash.downsample(self, num=num, scaled=scaled)
        down_mh.into_frozen()
        return down_mh

    def flatten(self):
        if not self.track_abundance:
            return self
        flat_mh = MinHash.flatten(self)
        flat_mh.into_frozen()
        return flat_mh

    def __iadd__(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def merge(self, *args, **kwargs):
        raise TypeError('FrozenMinHash does not support modification')

    def to_mutable(self):
        "Return a copy of this MinHash that can be changed."
        mut = MinHash.__new__(MinHash)
        state_tup = self.__getstate__()

        mut.__setstate__(state_tup)
        return mut

    def to_frozen(self):
        "Return a frozen copy of this MinHash that cannot be changed."
        return self

    def into_frozen(self):
        "Freeze this MinHash, preventing any changes."
        pass

    def __setstate__(self, tup):
        "support pickling via __getstate__/__setstate__"
        (n, ksize, is_protein, dayhoff, hp, mins, _, track_abundance,
         max_hash, seed) = tup

        self.__del__()

        hash_function = (
            lib.HASH_FUNCTIONS_MURMUR64_DAYHOFF if dayhoff else
            lib.HASH_FUNCTIONS_MURMUR64_HP if hp else
            lib.HASH_FUNCTIONS_MURMUR64_PROTEIN if is_protein else
            lib.HASH_FUNCTIONS_MURMUR64_DNA
        )

        scaled = _get_scaled_for_max_hash(max_hash)
        self._objptr = lib.kmerminhash_new(
            scaled, ksize, hash_function, seed, track_abundance, n
        )
        if track_abundance:
            MinHash.set_abundances(self, mins)
        else:
            MinHash.add_many(self, mins)

    def __copy__(self):
        return self
    copy = __copy__
