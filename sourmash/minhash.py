# -*- coding: UTF-8 -*-
from __future__ import unicode_literals, division

import math
import copy
import collections

from . import VERSION
from ._lowlevel import ffi, lib
from .utils import RustObject, rustcall, decode_str
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

    return int(round(get_minhash_max_hash() / scaled, 0))


def _get_scaled_for_max_hash(max_hash):
    "Convert a 'max_hash' value into a 'scaled' value."
    if max_hash == 0:
        return 0
    return int(round(get_minhash_max_hash() / max_hash, 0))


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


class _HashesWrapper(collections.Mapping):
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

        # ok, for Rust API, go from scaled back to max_hash
        max_hash = _get_max_hash_for_scaled(scaled)
        self._objptr = lib.kmerminhash_new(
            n, ksize, is_protein, dayhoff, hp, seed, int(max_hash), track_abundance
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
            max_hash=self.max_hash,
        )
        a.merge(self)
        return a

    def __getstate__(self):
        "support pickling via __getstate__/__setstate__"
        return (
            self.num,
            self.ksize,
            self.is_protein,
            self.dayhoff,
            self.hp,
            self.get_mins(with_abundance=self.track_abundance),
            None,
            self.track_abundance,
            self.max_hash,
            self.seed,
        )

    def __setstate__(self, tup):
        "support pickling via __getstate__/__setstate__"
        (n, ksize, is_protein, dayhoff, hp, mins, _, track_abundance,
         max_hash, seed) = tup

        self.__del__()
        self._objptr = lib.kmerminhash_new(
            n, ksize, is_protein, dayhoff, hp, seed, max_hash, track_abundance
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
            self.is_protein,
            self.dayhoff,
            self.hp,
            self.track_abundance,
            self.seed,
            self.max_hash,
        )
        return a

    def add_sequence(self, sequence, force=False):
        "Add a sequence into the sketch."
        self._methodcall(lib.kmerminhash_add_sequence, to_bytes(sequence),
                         force)

    @deprecated(deprecated_in="3.5", removed_in="4.0",
                current_version=VERSION,
                details='Use add_kmer instead.')
    def add(self, kmer):
        "Add a kmer into the sketch."
        self.add_sequence(kmer)

    def add_kmer(self, kmer):
        "Add a kmer into the sketch."
        if len(kmer) != self.ksize:
            raise ValueError("kmer to add is not {} in length".format(self.ksize))
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
        "Remove many hashes at once; ``hashes`` must be an iterable."
        self._methodcall(lib.kmerminhash_remove_many, list(hashes), len(hashes))

    @deprecated(deprecated_in="3.5", removed_in="4.0",
                current_version=VERSION,
                details='Use add_many instead.')
    def update(self, other):
        "Update this sketch from all the hashes in the other."
        self.add_many(other)

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
        size = ffi.new("uintptr_t *")
        mins_ptr = self._methodcall(lib.kmerminhash_get_mins, size)
        size = size[0]

        try:
            if with_abundance and self.track_abundance:
                size_abunds = ffi.new("uintptr_t *")
                abunds_ptr = self._methodcall(lib.kmerminhash_get_abunds, size_abunds)
                size_abunds = size_abunds[0]
                assert size == size_abunds
                result = dict(zip(ffi.unpack(mins_ptr, size), ffi.unpack(abunds_ptr, size)))
                lib.kmerminhash_slice_free(abunds_ptr, size)
            else:
                result = ffi.unpack(mins_ptr, size)
        finally:
            lib.kmerminhash_slice_free(mins_ptr, size)

        return result

    @deprecated(deprecated_in="3.5", removed_in="5.0",
                current_version=VERSION,
                details='Use .hashes property instead.')
    def get_hashes(self):
        "Return the list of hashes."
        return self.get_mins()

    @property
    def hashes(self):
        if self.track_abundance:
            return _HashesWrapper(self.get_mins(with_abundance=True))
        else:
            d = self.get_mins()
            return _HashesWrapper({ k : 1 for k in d })

    @deprecated(deprecated_in="3.5", removed_in="4.0",
                current_version=VERSION)
    def subtract_mins(self, other):
        """Get the list of mins in this MinHash, after removing the ones in
        ``other``.
        """
        a = set(self.get_mins())
        b = set(other.get_mins())
        return a - b

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
        return self._methodcall(lib.kmerminhash_ksize)

    @property
    @deprecated(deprecated_in="3.5", removed_in="5.0",
                current_version=VERSION,
                details='Use scaled instead.')
    def max_hash(self):
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

    @deprecated(deprecated_in="3.5", removed_in="4.0",
                current_version=VERSION,
                details='Use translate_codon function at module level instead.')
    def translate_codon(self, codon):
        "Translate a codon into an amino acid."
        try:
            return rustcall(lib.sourmash_translate_codon,
                            to_bytes(codon)).decode('utf-8')
        except SourmashError as e:
            raise ValueError(e.message)

    def count_common(self, other, downsample=False):
        """\
        Return the number of hashes in common between ``self`` and ``other``.

        Optionally downsample ``scaled`` objects to highest ``scaled`` value.
        """
        if not isinstance(other, MinHash):
            raise TypeError("Must be a MinHash!")
        return self._methodcall(lib.kmerminhash_count_common, other._get_objptr(), downsample)

    def downsample(self, num=None, scaled=None):
        """Copy this object and downsample new object to either `num` or
        `scaled`.
        """
        if num is None and scaled is None:
            raise ValueError('must specify either num or scaled to downsample')
        elif num is not None:
            if self.num and self.num < num:
                raise ValueError("new sample num is higher than current sample num")
            max_hash=0
        elif scaled is not None:
            if self.num:
                raise ValueError("num != 0 - cannot downsample a standard MinHash")
            max_hash = self.max_hash
            if max_hash is None:
                raise ValueError("no max_hash available - cannot downsample")

            old_scaled = _get_scaled_for_max_hash(self.max_hash)
            if old_scaled > scaled:
                raise ValueError(
                    "new scaled {} is lower than current sample scaled {}".format(
                        scaled, old_scaled
                    )
                )

            max_hash = _get_max_hash_for_scaled(scaled)
            num = 0
        ###

        # create new object:
        a = MinHash(
            num, self.ksize, self.is_protein, self.dayhoff, self.hp,
            self.track_abundance, self.seed, max_hash
        )
        # copy over hashes:
        if self.track_abundance:
            a.set_abundances(self.get_mins(with_abundance=True))
        else:
            a.add_many(self)

        return a

    @deprecated(deprecated_in="3.5", removed_in="4.0",
                current_version=VERSION,
                details='Use downsample(num=...) instead.')
    def downsample_n(self, new_num):
        "Copy this object and downsample new object to num=``new_num``."
        return self.downsample(num=new_num)

    @deprecated(deprecated_in="3.5", removed_in="4.0",
                current_version=VERSION,
                details='Use scaled instead.')
    def downsample_max_hash(self, *others):
        """Copy this object and downsample new object to min of ``*others``.

        Here, ``*others`` is one or more MinHash objects.
        """
        max_hashes = [x.max_hash for x in others]
        new_max_hash = min(self.max_hash, *max_hashes)
        new_scaled = _get_scaled_for_max_hash(new_max_hash)

        return self.downsample_scaled(new_scaled)

    @deprecated(deprecated_in="3.5", removed_in="4.0",
                current_version=VERSION,
                details='Use downsample(scaled=...) instead.')
    def downsample_scaled(self, new_scaled):
        """Copy this object and downsample new object to scaled=``new_scaled``.
        """
        return self.downsample(scaled=new_scaled)

    @deprecated(deprecated_in="3.3", removed_in="4.0",
                current_version=VERSION,
                details='Use count_common or set methods instead.')
    def intersection(self, other, in_common=False):
        """Calculate the intersection between ``self`` and ``other``, and
        return ``(mins, size)`` where ``mins`` are the hashes in common, and
        ``size`` is the number of hashes.

        if ``in_common``, return the actual hashes. Otherwise, mins will be
        empty.
        """
        if not isinstance(other, MinHash):
            raise TypeError("Must be a MinHash!")

        if self.num != other.num:
            err = "must have same num: {} != {}".format(self.num, other.num)
            raise TypeError(err)

        if in_common:
            # TODO: copy from buffer to Python land instead,
            # this way involves more moving data around.
            combined_mh = self.copy_and_clear()
            combined_mh.merge(self)
            combined_mh.merge(other)

            size = len(combined_mh)
            common = set(self.get_mins())
            common.intersection_update(other.get_mins())
        else:
            size = self._methodcall(lib.kmerminhash_intersection, other._get_objptr())
            common = set()

        return common, max(size, 1)

    def flatten(self):
        """Return a new MinHash with track_abundance=False."""
        # create new object:
        a = MinHash(
            self.num, self.ksize, self.is_protein, self.dayhoff, self.hp,
            False, self.seed, self.max_hash
        )
        a.add_many(self)

        return a

    def jaccard(self, other, downsample=False):
        "Calculate Jaccard similarity of two MinHash objects."
        if self.num != other.num:
            err = "must have same num: {} != {}".format(self.num, other.num)
            raise TypeError(err)
        return self._methodcall(lib.kmerminhash_similarity, other._get_objptr(), True, downsample)

    @deprecated(deprecated_in="3.3", removed_in="4.0",
                current_version=VERSION,
                details="Use 'similarity' instead of compare.")
    def compare(self, other, downsample=False):
        "Calculate Jaccard similarity of two sketches."
        return self.jaccard(other, downsample=downsample)


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
        return self._methodcall(lib.kmerminhash_angular_similarity,
                                other._get_objptr())

    def is_compatible(self, other):
        return self._methodcall(lib.kmerminhash_is_compatible, other._get_objptr())

    def contained_by(self, other, downsample=False):
        """\
        Calculate how much of self is contained by other.
        """
        if not len(self):
            return 0.0

        return self.count_common(other, downsample) / len(self)

    @deprecated(deprecated_in="3.3", removed_in="4.0",
                current_version=VERSION,
                details="Use 'contained_by' with downsample=True instead.")
    def containment_ignore_maxhash(self, other):
        """Calculate contained_by, with downsampling.
        """
        return self.contained_by(other, downsample=True)

    def __iadd__(self, other):
        if not isinstance(other, MinHash):
            raise TypeError("Must be a MinHash!")
        self._methodcall(lib.kmerminhash_merge, other._get_objptr())
        return self

    merge = __iadd__

    def set_abundances(self, values, clear=True):
        """Set abundances for hashes from ``values``, where
        ``values[hash] = abund``
        """
        if self.track_abundance:
            hashes = []
            abunds = []

            for h, v in values.items():
                hashes.append(h)
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

    @deprecated(deprecated_in="3.5", removed_in="4.0",
                current_version=VERSION,
                details='Use the moltype property instead.')
    def is_molecule_type(self, molecule):
        """Check if this MinHash is a particular human-readable molecule type.

        Supports 'protein', 'dayhoff', 'hp', 'DNA'.
        @CTB deprecate for 4.0?
        """
        if molecule.lower() not in ('protein', 'dayhoff', 'hp', 'dna'):
            raise ValueError("unknown moltype in query, '{}'".format(molecule))
        return molecule == self.moltype

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
