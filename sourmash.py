#! /usr/bin/env python2
import khmer
import screed
import argparse

K = 19

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
        hs = self._kh.get_kmer_hashes(seq)
        for h in hs:
            self.add(h)

    def jaccard(self, other):
        # this is stupid and badly implemented but I'm pretty sure it works :)
        jaccard = 0
        d = {}
        n = len(self._mins)

        for i in range(n):
            key = self._mins[i]
            d[key] = 1
            assert key != self.p

        for i in range(n):
            key = other._mins[i]
            d[key] = 1 + d.get(key, 0)
            assert key != other.p

        for v in d.values():
            if v == 2:
                jaccard += 1

        return float(jaccard) / float(n)


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
    parser.add_argument('-w', '--weighted', action='store_true')
    args = parser.parse_args()

    if args.weighted:
        print 'using weighted estimator'
        E = WeightedEstimators()
        E2 = WeightedEstimators()
    else:
        print 'using unweighted estimator'
        E = Estimators()
        E2 = Estimators()
        

    print 'reading first'
    for record in screed.open(args.sequences1):
        E.add_sequence(record.sequence)

    print 'reading second'
    for record in screed.open(args.sequences2):
        E2.add_sequence(record.sequence)

    jaccard = E.jaccard(E2)

    print 'similarity', args.sequences1, args.sequences2, jaccard


if __name__ == '__main__':
    main()
