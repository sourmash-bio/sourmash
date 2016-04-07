#! /usr/bin/env python2
import khmer
import screed
import argparse
import itertools

K = 31
NUM_ESTIMATORS=100

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

class Estimators(object):
    """
    A simple bottom n-sketch MinHash implementation.
    """

    def __init__(self, n=100, max_prime=1e10, ksize=K):
        _kh = khmer.Countgraph(ksize, 1, 1)
        self._kh = _kh

        # get a prime to use for hashing
        p = khmer.get_n_primes_near_x(1, max_prime)[0]
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
        #hs = self._kh.get_kmer_hashes(seq)
        #for h in hs:
        #    self.add(h)
        for kmer in kmers(seq, self._kh.ksize()):
            h = khmer.hash_murmur3(kmer)
            self.add(h)

    def jaccard(self, other):
        common = 0
        for _ in yield_overlaps(self._mins, other._mins):
            common += 1

        return float(common) / float(len(self._mins))

    def common(self, other):
        common = 0
        for _ in yield_overlaps(self._mins, other._mins):
            common += 1
        return common
    

class CompositionSketchEstimator(object):
    def __init__(self, n=100, max_prime=1e10, ksize=K):
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

    
class WeightedEstimators(Estimators):
    def _init_kh(self):
        self._kh = khmer.Countgraph(ksize, 1e8, 4)

    def add_sequence(self, seq):
        seq = seq.upper().replace('N', 'G')
        hs = self._kh.get_kmer_hashes(seq)
        self._kh.consume(seq)           # track k-mer abundances
        for h in hs:
            self.add(h * self._kh.get(h)) # weight by k-mer abundances


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sequences1')
    parser.add_argument('sequences2')
    parser.add_argument('-k', '--ksize', type=int, default=K)
    parser.add_argument('-n', '--num_estimators', type=int,
                        default=NUM_ESTIMATORS)
    parser.add_argument('-w', '--weighted', action='store_true')
    args = parser.parse_args()

    if args.weighted:
        print('using weighted estimator')
        E = WeightedEstimators(n=args.num_estimators, ksize=args.ksize)
        E2 = WeightedEstimators(n=args.num_estimators, ksize=args.ksize)
    else:
        print('using unweighted estimator')
        E = Estimators(n=args.num_estimators, ksize=args.ksize)
        E2 = Estimators(n=args.num_estimators, ksize=args.ksize)
        

    print('reading both')
    n = 0
    for r1, r2 in itertools.izip(screed.open(args.sequences1),
                                 screed.open(args.sequences2)):
        E.add_sequence(r1.sequence)
        E2.add_sequence(r2.sequence)
        n += 1

        if n % 10000 == 0:
            jaccard = E.jaccard(E2)
            print(n, 'similarity', args.sequences1, args.sequences2, jaccard)
            #if n >= 100000:
            #    break


if __name__ == '__main__':
    main()
