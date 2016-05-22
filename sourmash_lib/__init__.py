#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import sys
import argparse
import itertools
import string

from . import _minhash

class Estimators(object):
    """
    A simple bottom n-sketch MinHash implementation.
    """

    def __init__(self, n=None, ksize=None, protein=False):
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception

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
        (self.num, self.ksize, self.is_protein, mins) = tup
        self.mh = _minhash.MinHash(self.num, self.ksize, self.is_protein)
        for m in mins:
            self.mh.add_hash(m)

    def add(self, kmer):
        "Add kmer into sketch, keeping sketch sorted."
        self.mh.add_sequence(kmer)

    def add_sequence(self, seq):
        "Sanitize and add a sequence to the sketch."
        seq = seq.upper().replace('N', 'G')
        self.mh.add_sequence(seq)

    def jaccard(self, other):
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


def test_dna_mh():
    e1 = Estimators(n=5, ksize=4)
    e2 = Estimators(n=5, ksize=4)

    seq = 'ATGGCAGTGACGATGCCAG'
    e1.add_sequence(seq)
    for i in range(len(seq) - 3):
        e2.add(seq[i:i + 4])

    assert e1.mh.get_mins() == e2.mh.get_mins()
    print(e1.mh.get_mins())
    assert 58416682101672111 in e1.mh.get_mins()
    assert 726311917625663847 in e1.mh.get_mins()


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
