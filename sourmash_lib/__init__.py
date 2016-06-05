#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import sys
import argparse
import itertools
import string

class Estimators(object):
    """
    A simple bottom n-sketch MinHash implementation.

    Usage::

       E = Estimators(n=1000, ksize=31)
       E.add_sequence(dna)
       ...
       E.jaccard(other_E)

    ``Estimator`` supports the pickle protocol.
    """

    def __init__(self, n=None, ksize=None, protein=False):
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

        # initialize sketch to size n
        self.mh = _minhash.MinHash(n, ksize, protein)

    def __getstate__(self):             # enable pickling
        return (self.num, self.ksize, self.is_protein, self.mh.get_mins())

    def __setstate__(self, tup):
        from . import _minhash

        (self.num, self.ksize, self.is_protein, mins) = tup
        self.mh = _minhash.MinHash(self.num, self.ksize, self.is_protein)
        for m in mins:
            self.mh.add_hash(m)

    def __eq__(self, other):
        print(self.__getstate__())
        print(other.__getstate__())
        return self.__getstate__() == other.__getstate__()

    def add(self, kmer):
        "Add kmer into sketch."
        self.mh.add_sequence(kmer)

    def add_sequence(self, seq, force=False):
        "Sanitize and add a sequence to the sketch."
        self.mh.add_sequence(seq, force)

    def jaccard(self, other):
        "Calculate Jaccard index of two sketches."
        return self.mh.compare(other.mh)
    similarity = jaccard

    def count_common(self, other):
        "Calculate number of common k-mers between two sketches."
        return self.mh.count_common(other.mh)


def test_jaccard_1():
    E1 = Estimators(n=5, ksize=20)
    E2 = Estimators(n=5, ksize=20)

    for i in [1, 2, 3, 4, 5]:
        E1.mh.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.mh.add_hash(i)

    assert round(E1.jaccard(E2), 2) == 4 / 5.0
    assert round(E2.jaccard(E1), 2) == 4 / 5.0


def test_jaccard_2_difflen():
    E1 = Estimators(n=5, ksize=20)
    E2 = Estimators(n=5, ksize=20)

    for i in [1, 2, 3, 4, 5]:
        E1.mh.add_hash(i)
    for i in [1, 2, 3, 4]:
        E2.mh.add_hash(i)

    assert round(E1.jaccard(E2), 2) == 4 / 5.0
    assert round(E2.jaccard(E1), 2) == 4 / 4.0


def test_common_1():
    E1 = Estimators(n=5, ksize=20)
    E2 = Estimators(n=5, ksize=20)

    for i in [1, 2, 3, 4, 5]:
        E1.mh.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.mh.add_hash(i)

    assert E1.count_common(E2) == 4
    assert E2.count_common(E1) == 4


def test_dna_mh():
    e1 = Estimators(n=5, ksize=4)
    e2 = Estimators(n=5, ksize=4)

    seq = 'ATGGCAGTGACGATGCCAG'
    e1.add_sequence(seq)
    for i in range(len(seq) - 3):
        e2.add(seq[i:i + 4])

    assert e1.mh.get_mins() == e2.mh.get_mins()
    print(e1.mh.get_mins())
    assert 726311917625663847 in e1.mh.get_mins()
    assert 3697418565283905118 in e1.mh.get_mins()


def test_protein_mh():
    e1 = Estimators(n=5, ksize=6, protein=True)
    e2 = Estimators(n=5, ksize=6, protein=True)

    seq = 'ATGGCAGTGACGATGCCG'
    e1.add_sequence(seq)

    for i in range(len(seq) - 5):
        kmer = seq[i:i + 6]
        e2.add(kmer)

    assert e1.mh.get_mins() == e2.mh.get_mins()
    assert 901193879228338100 in e1.mh.get_mins()


def test_pickle():
    import pickle
    from io import BytesIO

    e1 = Estimators(n=5, ksize=6, protein=False)
    seq = 'ATGGCAGTGACGATGCCG'
    e1.add_sequence(seq)

    fp = BytesIO()
    pickle.dump(e1, fp)

    fp2 = BytesIO(fp.getvalue())
    e2 = pickle.load(fp2)

    assert e1.mh.get_mins() == e2.mh.get_mins()
    assert e1.num == e2.num
    assert e1.ksize == e2.ksize
    assert e1.is_protein == e2.is_protein


def test_bad_construct_1():
    try:
        e1 = Estimators(ksize=6, protein=False)
        assert 0, "require n in constructor"
    except ValueError:
        pass


def test_bad_construct_2():
    try:
        e1 = Estimators(n=100, protein=False)
        assert 0, "require ksize in constructor"
    except ValueError:
        pass
