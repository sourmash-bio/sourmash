#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import re
import math
from ._minhash import MinHash, dotproduct

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
                 track_abundance=False, max_hash=0, seed=DEFAULT_SEED):
        """\
        Create a new MinHash estimator with size n and k-mer size ksize.

        is_protein - compute hashes from amino acid translation (False)
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
        self.track_abundance = track_abundance

        # initialize sketch to size n
        self.mh = MinHash(n, ksize,
                          is_protein=is_protein,
                          track_abundance=track_abundance,
                          max_hash=max_hash,
                          seed=seed)

    def is_molecule_type(self, molecule):
        return self.mh.is_molecule_type(molecule)

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
        self.mh.add_sequence(kmer)

    def add_many(self, hashes):
        self.mh.add_many(hashes)

    def update(self, other):
        "Update this estimator from all the hashes from the other."
        self.mh.update(other.mh)

    def get_hashes(self):
        "Get the list of hashes."
        return self.mh.get_mins()

    def add_sequence(self, seq, force=False):
        "Sanitize and add a sequence to the sketch."
        self.mh.add_sequence(seq, force)

    def jaccard(self, other):
        "Calculate Jaccard index of two sketches."
        return self.mh.compare(other.mh)

    def similarity_ignore_maxhash(self, other):
        return self.mh.similarity_ignore_maxhash(other.mh)

    def similarity(self, other, ignore_abundance=False):
        return self.mh.similarity(other.mh, ignore_abundance)

    def count_common(self, other):
        "Calculate number of common k-mers between two sketches."
        return self.mh.count_common(other.mh)
