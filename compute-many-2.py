#! /usr/bin/env python
import sourmash
import argparse
import screed
from cPickle import load

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dumpfile')
    args = parser.parse_args()

    emins = load(open(args.dumpfile))
    estimators = []
    for (filename, mins) in emins:
        E = sourmash.Estimators()
        E._mins = mins
        estimators.append((filename, E))

    for f, E in estimators:
        for f2, E2 in estimators:
            print f, f2, E.jaccard(E2)
    

if __name__ == '__main__':
    main()
    
