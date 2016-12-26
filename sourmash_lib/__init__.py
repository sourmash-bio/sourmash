#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import re

khmer_available = False
try:
    import khmer
    khmer_available = True
except ImportError:
    pass

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

    def __init__(self, n=None, ksize=None, protein=False,
                 with_cardinality=False, track_abundance=False):
        "Create a new MinHash estimator with size n and k-mer size ksize."
        from . import _minhash

        if n is None:
            raise ValueError("n is required")
        if ksize is None:
            raise ValueError("ksize is required")

        self.num = n
        self.ksize = ksize
        self.is_protein = False
        if protein:
            self.is_protein = True

        self.hll = None
        if with_cardinality:
            if protein:
                raise Exception("Cannot do cardinality counting with protein")
            if not khmer_available:
                raise Exception("Error: to do cardinality counting, " + \
                                "we require the khmer package.")
            self.hll = khmer.HLLCounter(.01, ksize)
        self.track_abundance = track_abundance

        # initialize sketch to size n
        self.mh = _minhash.MinHash(n, ksize, protein, track_abundance)

    def is_molecule_type(self, molecule):
        if molecule == 'dna' and not self.mh.is_protein():
            return True
        if molecule == 'protein' and self.mh.is_protein():
            return True
        return False

    def __getstate__(self):             # enable pickling
        return (self.num, self.ksize, self.is_protein, self.mh.get_mins(),
                self.hll, self.track_abundance)

    def __setstate__(self, tup):
        from . import _minhash

        (self.num, self.ksize, self.is_protein, mins, hll, self.track_abundance) = tup
        self.mh = _minhash.MinHash(self.num, self.ksize, self.is_protein, self.track_abundance)
        for m in mins:
            self.mh.add_hash(m)
        self.hll = hll

    def __eq__(self, other):
        return self.__getstate__() == other.__getstate__()

    def add(self, kmer):
        "Add kmer into sketch."
        self.mh.add_sequence(kmer)

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

    def similarity(self, other, ignore_abundance=False):
        if not self.track_abundance or ignore_abundance:
            return self.jaccard(other)
        else:
            a = self.mh.get_mins(with_abundance=True)
            b = other.mh.get_mins(with_abundance=True)

            common_abund = 0
            total_abund = 0
            for k, abundance in a.items():
                common_abund += b.get(k, 0)
                total_abund += abundance
            return common_abund / float(total_abund)

    def count_common(self, other):
        "Calculate number of common k-mers between two sketches."
        return self.mh.count_common(other.mh)
