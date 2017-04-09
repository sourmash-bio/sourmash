#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import re
import math
from ._minhash import MinHash
import os

# retrieve VERSION from sourmash_lib/VERSION.
thisdir = os.path.dirname(__file__)
version_file = open(os.path.join(thisdir, 'VERSION'))
VERSION = version_file.read().strip()

khmer_available = False
try:
    import khmer
    khmer_available = True
except ImportError:
    pass

DEFAULT_SEED=MinHash(1,1).seed

class Estimators(object):
    """
    A simple bottom n-sketch MinHash implementation.

    Usage::

       E = Estimators(n=1000, ksize=31)
       E.add_sequence(dna)
       ...
       similarity = E.jaccard(other_E)

    ``Estimator`` supports the pickle protocol.
    """

    def __init__(self, n=None, ksize=None, is_protein=False,
                 with_cardinality=False, track_abundance=False,
                 max_hash=0, seed=DEFAULT_SEED):
        """\
        Create a new MinHash estimator with size n and k-mer size ksize.

        is_protein - compute hashes from amino acid translation (False)
        with_cardinality - count total unique k-mers (False)
        track_abundance - track abundance of k-mers as well as presence (False)
        max_hash - only admit hash values under this number (not set)
        seed - hash function seed
        """
        if n is None:
            raise ValueError("n is required")
        if ksize is None:
            raise ValueError("ksize is required")

        self.num = n
        self.ksize = ksize
        self.is_protein = False
        self.max_hash = max_hash
        self.seed = seed
        if is_protein:
            self.is_protein = True

        self.hll = None
        if with_cardinality:
            if is_protein:
                raise Exception("Cannot do cardinality counting with protein")
            if not khmer_available:
                raise Exception("Error: to do cardinality counting, " + \
                                "we require the khmer package.")
            self.hll = khmer.HLLCounter(.01, ksize)
        self.track_abundance = track_abundance

        # initialize sketch to size n
        self.mh = MinHash(n, ksize,
                          is_protein=is_protein,
                          track_abundance=track_abundance,
                          max_hash=max_hash,
                          seed=seed)

    def is_molecule_type(self, molecule):
        if molecule == 'dna' and not self.mh.is_protein():
            return True
        if molecule == 'protein' and self.mh.is_protein():
            return True
        return False

    def __getstate__(self):             # enable pickling
        with_abundance = False
        if self.track_abundance:
            with_abundance = True

        return (self.num, self.ksize, self.is_protein,
                self.mh.get_mins(with_abundance=with_abundance),
                self.hll, self.track_abundance, self.max_hash,
                self.seed)

    def __setstate__(self, tup):
        (self.num, self.ksize, self.is_protein, mins, hll,
         self.track_abundance, self.max_hash, self.seed) = tup
        self.mh = MinHash(self.num, self.ksize,
                          is_protein=self.is_protein,
                          track_abundance=self.track_abundance,
                          max_hash=self.max_hash,
                          seed=self.seed)

        if not self.track_abundance:
            self.add_many(mins)
        else:
            self.mh.set_abundances(mins)

        self.hll = hll

    def __eq__(self, other):
        return self.__getstate__() == other.__getstate__()

    def add(self, kmer):
        "Add kmer into sketch."
        self.mh.add_sequence(kmer)

    def add_many(self, hashes):
        "Add many hashes in at once."
        for hash in hashes:
            self.mh.add_hash(hash)

    def update(self, other):
        "Update this estimator from all the hashes from the other."
        self.add_many(other.mh.get_mins())

    def get_hashes(self):
        "Get the list of hashes."
        return self.mh.get_mins()

    def add_sequence(self, seq, force=False):
        "Sanitize and add a sequence to the sketch."
        self.mh.add_sequence(seq, force)
        if self.hll is not None:
            # @CTB hack to deal with HLL behavior around non-ACGTN.
            if force:
                seq = re.sub('[^ACGT]', 'A', seq)
            self.hll.consume_string(seq)

    def jaccard(self, other):
        "Calculate Jaccard index of two sketches."
        return self.mh.compare(other.mh)

    def similarity_ignore_maxhash(self, other):
        a = set(self.mh.get_mins())
        b = set(other.mh.get_mins())

        overlap = a.intersection(b)
        return float(len(overlap)) / float(len(a))


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

        if not self.track_abundance or ignore_abundance:
            return self.jaccard(other)
        else:
            a = self.mh.get_mins(with_abundance=True)
            b = other.mh.get_mins(with_abundance=True)

            prod = dotproduct(a, b)
            prod = min(1.0, prod)

            distance = 2*math.acos(prod) / math.pi
            return 1.0 - distance

    def count_common(self, other):
        "Calculate number of common k-mers between two sketches."
        return self.mh.count_common(other.mh)


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


def test_dotproduct_1():
    a = {'x': 1}
    assert dotproduct(a, a, normalize=True) == 1.0

    a = {'x': 1}
    b = {'x': 1}
    assert dotproduct(a, b, normalize=True) == 1.0

    c = {'x': 1, 'y': 1}
    prod = dotproduct(c, c, normalize=True)
    assert round(prod, 2) == 1.0

    # check a.c => 45 degree angle
    a = {'x': 1}
    c = {'x': 1, 'y': 1}

    angle = 45
    rad = math.radians(angle)
    cosval = math.cos(rad)
    prod = dotproduct(a, c, normalize=True)
    assert round(prod, 2) == 0.71
    assert round(cosval, 2) == round(prod, 2)

    c = {'x': 1, 'y': 1}
    d = {'x': 1, 'y': 1}
    prod = dotproduct(c, d, normalize=True)
    assert round(prod, 2) == 1.0

    a = {'x': 1}
    e = {'y': 1}
    assert dotproduct(a, e, normalize=True) == 0.0


def test_dotproduct_zeroes():
    a = {'x': 1}
    b = {}

    assert dotproduct(a, b) == 0.0
    assert dotproduct(b, a) == 0.0
