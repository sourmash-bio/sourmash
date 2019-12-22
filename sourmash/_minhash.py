# -*- coding: UTF-8 -*-
from __future__ import unicode_literals, division

import math
import copy

from ._compat import string_types, range_type
from ._lowlevel import ffi, lib
from .utils import RustObject, rustcall, decode_str
from .exceptions import SourmashError

# default MurmurHash seed
MINHASH_DEFAULT_SEED = 42


def get_minhash_default_seed():
    return MINHASH_DEFAULT_SEED


# we use the 64-bit hash space of MurmurHash only
# this is 2 ** 64 - 1 in hexadecimal
MINHASH_MAX_HASH = 0xFFFFFFFFFFFFFFFF


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


def to_bytes(s):
    # Allow for strings, bytes or int
    # Single item of byte string = int

    if isinstance(s, bytes):
        return s

    if not isinstance(s, string_types + (bytes, int)):
        raise TypeError("Requires a string-like sequence")

    if isinstance(s, string_types):
        s = s.encode("utf-8")
    elif isinstance(s, int):
        s = bytes([s])

    return s


def hash_murmur(kmer, seed=MINHASH_DEFAULT_SEED):
    "hash_murmur(string, [,seed])\n\n"
    "Compute a hash for a string, optionally using a seed (an integer). "
    "The current default seed is returned by hash_seed()."

    return lib.hash_murmur(to_bytes(kmer), seed)


class MinHash(RustObject):
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
        if max_hash and scaled:
            raise ValueError("cannot set both max_hash and scaled")
        elif scaled:
            max_hash = get_max_hash_for_scaled(scaled)

        if max_hash and n:
            raise ValueError("cannot set both n and max_hash")

        if not n and not (max_hash or scaled):
            raise ValueError("cannot omit both n and scaled")

        if dayhoff or hp:
            is_protein = False

        self._objptr = lib.kmerminhash_new(
            n, ksize, is_protein, dayhoff, hp, seed, int(max_hash), track_abundance
        )
        self.__dealloc_func__ = lib.kmerminhash_free

        if mins:
            if track_abundance:
                self.set_abundances(mins)
            else:
                self.add_many(mins)

    def __copy__(self):
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

    def __getstate__(self):  # enable pickling
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
        (n, ksize, is_protein, dayhoff, hp, mins, _, track_abundance, max_hash, seed) = tup

        self.__del__()
        self._objptr = lib.kmerminhash_new(
            n, ksize, is_protein, dayhoff, hp, seed, max_hash, track_abundance
        )
        if track_abundance:
            self.set_abundances(mins)
        else:
            self.add_many(mins)

    def __reduce__(self):
        return (
            MinHash,
            (
                self.num,
                self.ksize,
                self.is_protein,
                self.dayhoff,
                self.hp,
                self.track_abundance,
                self.seed,
                self.max_hash,
                self.get_mins(with_abundance=self.track_abundance),
                0,
            ),
        )

    def __eq__(self, other):
        return self.__getstate__() == other.__getstate__()

    def copy_and_clear(self):
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
        self._methodcall(lib.kmerminhash_add_sequence, to_bytes(sequence), force)

    def add(self, kmer):
        "Add kmer into sketch."
        self.add_sequence(kmer)

    def add_many(self, hashes):
        "Add many hashes in at once."
        if isinstance(hashes, MinHash):
            self._methodcall(lib.kmerminhash_add_from, hashes._objptr)
        else:
            for hash in hashes:
                self._methodcall(lib.kmerminhash_add_hash, hash)

    def remove_many(self, hashes):
        "Add many hashes in at once."
        self._methodcall(lib.kmerminhash_remove_many, list(hashes), len(hashes))

    def update(self, other):
        "Update this estimator from all the hashes from the other."
        self.add_many(other)

    def __len__(self):
        return self._methodcall(lib.kmerminhash_get_mins_size)

    def get_mins(self, with_abundance=False):
        size = self._methodcall(lib.kmerminhash_get_mins_size)
        mins_ptr = self._methodcall(lib.kmerminhash_get_mins)

        if with_abundance and self.track_abundance:
            abunds_ptr = self._methodcall(lib.kmerminhash_get_abunds)
            result = dict(zip(ffi.unpack(mins_ptr, size), ffi.unpack(abunds_ptr, size)))
            lib.kmerminhash_slice_free(abunds_ptr, size)
        else:
            result = ffi.unpack(mins_ptr, size)

        lib.kmerminhash_slice_free(mins_ptr, size)
        return result

    def get_hashes(self):
        return self.get_mins()

    def subtract_mins(self, other):
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
        if self.max_hash:
            return get_scaled_for_max_hash(self.max_hash)
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
        return self._methodcall(lib.kmerminhash_add_hash, h)

    def translate_codon(self, codon):
        try:
            return rustcall(lib.sourmash_translate_codon,
                            to_bytes(codon)).decode('utf-8')
        except SourmashError as e:
            raise ValueError(e.message)

    def count_common(self, other):
        if not isinstance(other, MinHash):
            raise TypeError("Must be a MinHash!")
        return self._methodcall(lib.kmerminhash_count_common, other._get_objptr())

    def downsample_n(self, new_num):
        if self.num and self.num < new_num:
            raise ValueError("new sample n is higher than current sample n")

        a = MinHash(
            new_num, self.ksize, self.is_protein, self.dayhoff, self.hp, self.track_abundance, self.seed, 0
        )
        if self.track_abundance:
            a.set_abundances(self.get_mins(with_abundance=True))
        else:
            a.add_many(self)

        return a

    def downsample_max_hash(self, *others):
        max_hashes = [x.max_hash for x in others]
        new_max_hash = min(self.max_hash, *max_hashes)
        new_scaled = get_scaled_for_max_hash(new_max_hash)

        return self.downsample_scaled(new_scaled)

    def downsample_scaled(self, new_num):
        if self.num:
            raise ValueError("num != 0 - cannot downsample a standard MinHash")

        max_hash = self.max_hash
        if max_hash is None:
            raise ValueError("no max_hash available - cannot downsample")

        old_scaled = get_scaled_for_max_hash(self.max_hash)
        if old_scaled > new_num:
            raise ValueError(
                "new scaled {} is lower than current sample scaled {}".format(
                    new_num, old_scaled
                )
            )

        new_max_hash = get_max_hash_for_scaled(new_num)

        a = MinHash(
            0,
            self.ksize,
            self.is_protein,
            self.dayhoff,
            self.hp,
            self.track_abundance,
            self.seed,
            new_max_hash,
        )
        if self.track_abundance:
            a.set_abundances(self.get_mins(with_abundance=True))
        else:
            a.add_many(self)

        return a

    def intersection(self, other, in_common=False):
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
            common.intersection_update(combined_mh.get_mins())
        else:
            size = self._methodcall(lib.kmerminhash_intersection, other._get_objptr())
            common = set()

        return common, max(size, 1)

    def compare(self, other):
        if self.num != other.num:
            err = "must have same num: {} != {}".format(self.num, other.num)
            raise TypeError(err)
        return self._methodcall(lib.kmerminhash_compare, other._get_objptr())
    jaccard = compare

    def similarity(self, other, ignore_abundance=False):
        """Calculate similarity of two sketches.

        If the sketches are not abundance weighted, or ignore_abundance=True,
        compute Jaccard similarity.

        If the sketches are abundance weighted, calculate a distance metric
        based on the cosine similarity.

        Note, because the term frequencies (tf-idf weights) cannot be negative,
        the angle will never be < 0deg or > 90deg.

        See https://en.wikipedia.org/wiki/Cosine_similarity
        """

        # if either signature is flat, calculate Jaccard only.
        if not (self.track_abundance and other.track_abundance) or ignore_abundance:
            return self.jaccard(other)
        else:
            return self._methodcall(lib.kmerminhash_similarity,
                                    other._get_objptr(),
                                    ignore_abundance)

    def is_compatible(self, other):
        return self._methodcall(lib.kmerminhash_is_compatible, other._get_objptr())

    def contained_by(self, other):
        """\
        Calculate how much of self is contained by other.
        """
        if not len(self):
            return 0.0

        return self.count_common(other) / len(self)

    def containment_ignore_maxhash(self, other):
        if len(self) == 0:
            return 0.0

        if not isinstance(other, MinHash):
            raise TypeError("Must be a MinHash!")

        return self._methodcall(lib.kmerminhash_containment_ignore_maxhash, other._get_objptr())

    def __iadd__(self, other):
        if not isinstance(other, MinHash):
            raise TypeError("Must be a MinHash!")
        self._methodcall(lib.kmerminhash_merge, other._get_objptr())
        return self

    merge = __iadd__

    def set_abundances(self, values):
        if self.track_abundance:
            added = 0

            for k, v in sorted(values.items()):
                if not self.max_hash or k <= self.max_hash:
                    self._methodcall(lib.kmerminhash_mins_push, k)
                    self._methodcall(lib.kmerminhash_abunds_push, v)
                    added += 1
                    if self.num > 0 and added >= self.num:
                        break
        else:
            raise RuntimeError(
                "Use track_abundance=True when constructing "
                "the MinHash to use set_abundances."
            )

    def add_protein(self, sequence):
        ksize = self.ksize // 3
        if len(sequence) < ksize:
            return

        aa_kmers = (sequence[i:i + ksize] for i in range(0, len(sequence) - ksize + 1))
        if self.is_protein:
            for aa_kmer in aa_kmers:
                self._methodcall(
                    lib.kmerminhash_add_word, to_bytes(aa_kmer)
        )
        elif self.dayhoff:
            for aa_kmer in aa_kmers:
                dayhoff_kmer = ''
                for aa in aa_kmer:
                    data = rustcall(lib.sourmash_aa_to_dayhoff, to_bytes(aa))
                    dayhoff_letter = data.decode('utf-8')
                    dayhoff_kmer += dayhoff_letter
                self._methodcall(
                    lib.kmerminhash_add_word, to_bytes(dayhoff_kmer)
                )
        elif self.hp:
            for aa_kmer in aa_kmers:
                hp_kmer = ''
                for aa in aa_kmer:
                    data = rustcall(lib.sourmash_aa_to_hp, to_bytes(aa))
                    hp_letter = data.decode('utf-8')
                    hp_kmer += hp_letter
                self._methodcall(
                    lib.kmerminhash_add_word, to_bytes(hp_kmer)
                )
        else:
            raise ValueError("Invalid protein type")

    def is_molecule_type(self, molecule):
        if self.is_protein and molecule == 'protein':
            return True
        elif self.dayhoff and molecule == 'dayhoff':
            return True
        elif self.hp and molecule == 'hp':
            return True
        elif molecule.upper() == "DNA" and self.is_dna:
            return True

        return False
