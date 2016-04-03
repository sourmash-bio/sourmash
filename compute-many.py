#! /usr/bin/env python
import sourmash
import argparse
import screed
from cPickle import dump

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()

    estimators = []
    for filename in args.filenames:
        E = sourmash.Estimators()
        for n, record in enumerate(screed.open(filename)):
            if n % 10000 == 0:
                print '...', filename, n
                if n > 0:
                    break
            E.add_sequence(record.sequence)

        estimators.append((filename, E._mins))

    fp = open('picklejar', 'w')
    dump(estimators, fp)
    fp.close()
    

if __name__ == '__main__':
    main()
    
