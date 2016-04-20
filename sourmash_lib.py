#! /usr/bin/env python
"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import khmer
import screed
import argparse
import itertools
import string
import _minhash

class Estimators(object):
    """
    A simple bottom n-sketch MinHash implementation.
    """

    def __init__(self, n=None, max_prime=1e10, ksize=None, protein=False):
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception

        self.num = n
        self.ksize = ksize
        self.is_protein = False
        if protein:
            self.is_protein = True

        # get a prime to use for hashing
        p = get_prime_lt_x(max_prime)
        self.p = p

        # initialize sketch to size n
        self.mh = _minhash.MinHash(n, ksize, p, protein)

    def __getstate__(self):             # enable pickling
        return (self.num, self.ksize, self.is_protein, self.p,
                self.mh.get_mins())

    def __setstate__(self, tup):
        (self.num, self.ksize, self.is_protein, self.p, mins) = tup
        self.mh = _minhash.MinHash(self.num, self.ksize, self.p,
                                   self.is_protein)
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
        _mins = self.mh.get_mins()
        truelen = len(_mins)
        
        return self.common(other) / float(truelen)
    similarity = jaccard

    def common(self, other):
        "Calculate number of common k-mers between two sketches."
        if self.ksize != other.ksize:
            raise Exception("different k-mer sizes - cannot compare")
        if self.p != other.p:
            raise Exception("different primse - cannot compare")

        common = 0
        self_mins = self.mh.get_mins()
        other_mins = other.mh.get_mins()
        for val in _yield_overlaps(self_mins, other_mins):
            common += 1
        return common


class CompositionSketch(object):
    def __init__(self, n=None, max_prime=1e10, ksize=None, prefixsize=None):
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception
        if prefixsize is None:
            raise Exception

        self.prefixsize = prefixsize
        self.ksize = ksize
        self.threshold = 0.01

        # get a prime to use for hashing
        p = get_prime_lt_x(max_prime)
        self.p = p

        # initialize 4**prefixsize MinHash sketches
        self.sketches = []
        for i in range(4**prefixsize):
            self.sketches.append(Estimators(n, self.p, ksize))

    def add(self, kmer):
        idx = khmer.forward_hash(kmer, self.prefixsize)
        E = self.sketches[idx]
        
        hash = khmer.hash_murmur3(kmer)
        E.add(hash)

    def add_sequence(self, seq):
        "Sanitize and add a sequence to the sketch."
        seq = seq.upper().replace('N', 'G')
        for kmer in kmers(seq, self.ksize):
            self.add(kmer)

    def jaccard(self, other):
        total = 0.
        count = 0
        for n, v in enumerate(self.sketches):
            try:
                total += v.jaccard(other.sketches[n])
                count += 1
            except ValueError:
                pass
        return total / float(count)

    def similarity(self, other):
        matches = 0
        count = 0
        for n, v in enumerate(self.sketches):
            try:
                f = v.jaccard(other.sketches[n])
                count += 1
                if f > self.threshold:
                    matches += 1
            except ValueError:
                pass
        return matches / float(count)


def _yield_overlaps(x1, x2):
    "yield common hash values while iterating over two sorted lists of hashes"
    i = 0
    j = 0
    try:
        while 1:
            while x1[i] < x2[j]:
                i += 1
            while x1[i] > x2[j]:
                j += 1
            if x1[i] == x2[j]:
                yield x1[i]
                i += 1
                j += 1
    except IndexError:
        return

# taken from khmer 2.0; original author Jason Pell.
def is_prime(number):
    """Check if a number is prime."""
    if number < 2:
        return False
    if number == 2:
        return True
    if number % 2 == 0:
        return False
    for _ in range(3, int(number ** 0.5) + 1, 2):
        if number % _ == 0:
            return False
    return True


def get_prime_lt_x(target):
    """Backward-find a prime smaller than (or equal to) target.

    Step backwards until a prime number (other than 2) has been
    found.

    Arguments: target -- the number to step backwards from
    """
    if target == 1 and number == 1:
        return 1

    i = int(target)
    if i % 2 == 0:
        i -= 1
    while i > 0:
        if is_prime(i):
            return i
        i -= 2

    if len(primes) != number:
        raise RuntimeError("unable to find a prime number < %d" % (target))


def test_jaccard_1():
    E1 = Estimators(n=5, ksize=20)
    E2 = Estimators(n=5, ksize=20)

    for i in [1, 2, 3, 4, 5]:
        E1.mh.add_hash(i)
    for i in [1, 2, 3, 4, 6]:
        E2.mh.add_hash(i)
    
    assert E1.jaccard(E2) == 4/5.0
    assert E2.jaccard(E1) == 4/5.0


def test_jaccard_2_difflen():
    E1 = Estimators(n=5, ksize=20)
    E2 = Estimators(n=5, ksize=20)

    for i in [1, 2, 3, 4, 5]:
        E1.mh.add_hash(i)
    for i in [1, 2, 3, 4]:
        E2.mh.add_hash(i)
    
    assert E1.jaccard(E2) == 4/5.0
    assert E2.jaccard(E1) == 4/4.0


def test_yield_overlaps():
    x1 = [1, 3, 5]
    x2 = [2, 4, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 0
               

def test_yield_overlaps_2():
    x1 = [1, 3, 5]
    x2 = [1, 2, 4, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 1
    assert len(list(_yield_overlaps(x2, x1))) == 1
               

def test_yield_overlaps_3():
    x1 = [1, 3, 6]
    x2 = [1, 2, 6]
    assert len(list(_yield_overlaps(x1, x2))) == 2
    assert len(list(_yield_overlaps(x2, x1))) == 2


def test_dna_mh():
    e1 = Estimators(n=5, ksize=4)
    e2 = Estimators(n=5, ksize=4)

    seq = 'ATGGCAGTGACGATGCCAG'
    e1.add_sequence(seq)
    for i in range(len(seq) - 3):
        e2.add(seq[i:i+4])

    assert e1.mh.get_mins() == e2.mh.get_mins()
    assert 1149966211 in e1.mh.get_mins()
    assert 530237262 in e1.mh.get_mins()


def test_protein_mh():
    e1 = Estimators(n=5, ksize=6, protein=True)
    e2 = Estimators(n=5, ksize=6, protein=True)

    seq = 'ATGGCAGTGACGATGCCG'
    e1.add_sequence(seq)

    for i in range(len(seq) - 5):
        kmer = seq[i:i+6]
        e2.add(kmer)

    assert e1.mh.get_mins() == e2.mh.get_mins()
    assert 1054538492 in e1.mh.get_mins()
