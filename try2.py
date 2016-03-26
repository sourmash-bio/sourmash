#! /usr/bin/env python2
import khmer
import screed
import argparse

K = 19

class Estimators(object):
    def __init__(self, n=100, max_prime=1e10, ksize=K):
        _kh = khmer.Countgraph(ksize, 1, 1)
        self._kh = _kh

        p = khmer.get_n_primes_near_x(1, max_prime)[0]
        self.p = p
        self._mins = [p]*n
        
    def add(self, hash):
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
    d = {}
    n = len(E._mins)

    for i in range(n):
        key = E._mins[i]
        d[key] = 1
        assert key != E2.p

    for i in range(n):
        key = E2._mins[i]
        d[key] = 1 + d.get(key, 0)
        assert key != E2.p

    for v in d.values():
        if v == 2:
            jaccard += 1

    print 'similarity', args.sequences1, args.sequences2, \
          float(jaccard) / float(n)

if __name__ == '__main__':
    main()
