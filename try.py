#! /usr/bin/env python2
import khmer
import screed
import argparse

K = 31

class Estimators(object):
    def __init__(self, n=100, max_prime=1e10, ksize=K):
        _kh = khmer.Countgraph(ksize, 1, 1)
        self._kh = _kh
        
        _mins = dict()
        _primes = dict()
        for i in range(n):
            _mins[i] = None

        ps = khmer.get_n_primes_near_x(n, max_prime)
        for i, p in zip(range(n), ps):
            _primes[i] = p

        self._mins = _mins
        self._primes = _primes

    def add(self, hash):
        _mins = self._mins
        for i, p in self._primes.items():
            h = hash % p
            if _mins[i] != None:
                _mins[i] = min(h, _mins[i])
            else:
                _mins[i] = h

    def add_sequence(self, seq):
        hs = self._kh.get_kmer_hashes(seq.upper().replace('N', 'G'))
        for h in hs:
            self.add(h)
        

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('sequences1')
    parser.add_argument('sequences2')
    args = parser.parse_args()

    E = Estimators()
    E2 = Estimators()

    print 'first'
    for record in screed.open(args.sequences1):
        E.add_sequence(record.sequence)

    print 'second'
    for record in screed.open(args.sequences2):
        E2.add_sequence(record.sequence)

    jaccard = 0
    n = len(E._mins)
    for i in range(n):
        #print E._mins[i], E2._mins[i], E._mins[i] == E2._mins[i]
        if E._mins[i] == E2._mins[i]:
            jaccard += 1

    print 'similarity', args.sequences1, args.sequences2, \
          float(jaccard) / float(n)

if __name__ == '__main__':
    main()
