from __future__ import print_function
import khmer
import screed
import argparse
import itertools

def yield_overlaps(x1, x2):
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

def kmers(seq, ksize):
    for i in range(len(seq) - ksize + 1):
        yield seq[i:i+ksize]


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


class Estimators(object):
    """
    A simple bottom n-sketch MinHash implementation.
    """

    def __init__(self, n=None, max_prime=1e10, ksize=None):
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception
        
        _kh = khmer.Countgraph(ksize, 1, 1)
        self._kh = _kh
        self.ksize = ksize

        # get a prime to use for hashing
        p = get_prime_lt_x(max_prime)
        self.p = p

        # initialize sketch to size n
        self._mins = [p]*n
        
    def add(self, hash):
        "Add hash into _mins, keeping _mins in sorted order."
        _mins = self._mins
        h = hash % self.p
        
        if h >= _mins[-1]:
            return

        for i, v in enumerate(_mins):
            if h < v:
                _mins.insert(i, h)
                _mins.pop()
                return
            elif h == v:
                return
            # else: h > v, keep on going.


    def add_sequence(self, seq):
        seq = seq.upper().replace('N', 'G')
        for kmer in kmers(seq, self._kh.ksize()):
            h = khmer.hash_murmur3(kmer)
            self.add(h)

    def jaccard(self, other):
        if self.ksize != other.ksize:
            raise Exception("different k-mer sizes - cannot compare")

        common = 0
        for _ in yield_overlaps(self._mins, other._mins):
            common += 1

        return float(common) / float(len(self._mins))

    def common(self, other):
        if self.ksize != other.ksize:
            raise Exception("different k-mer sizes - cannot compare")

        common = 0
        for _ in yield_overlaps(self._mins, other._mins):
            common += 1
        return common

    def _truncate(self, n):
        self._mins = self._mins[:n]
    

class CompositionSketchEstimator(object):
    def __init__(self, n=None, max_prime=1e10, ksize=None):
        pass


def test_jaccard_1():
    E1 = Estimators(n=0)
    E2 = Estimators(n=0)

    E1._mins = [1, 2, 3, 4, 5]
    E2._mins = [1, 2, 3, 4, 6]
    assert E1.jaccard(E2) == 4/5.0
    assert E2.jaccard(E1) == 4/5.0


def test_jaccard_2_difflen():
    E1 = Estimators(n=0)
    E2 = Estimators(n=0)

    E1._mins = [1, 2, 3, 4, 5]
    E2._mins = [1, 2, 3, 4]
    assert E1.jaccard(E2) == 4/5.0
    assert E2.jaccard(E1) == 4/4.0


def test_yield_overlaps():
    x1 = [1, 3, 5]
    x2 = [2, 4, 6]
    assert len(list(yield_overlaps(x1, x2))) == 0
               

def test_yield_overlaps_2():
    x1 = [1, 3, 5]
    x2 = [1, 2, 4, 6]
    assert len(list(yield_overlaps(x1, x2))) == 1
    assert len(list(yield_overlaps(x2, x1))) == 1
               

def test_yield_overlaps_2():
    x1 = [1, 3, 6]
    x2 = [1, 2, 6]
    assert len(list(yield_overlaps(x1, x2))) == 2
    assert len(list(yield_overlaps(x2, x1))) == 2
