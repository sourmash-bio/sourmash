#! /usr/bin/env python
import sourmash
import argparse
import screed
from pickle import load
import sig

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('dumpfile')
    args = parser.parse_args()

    labels = dict([ i.strip().split(' ', 2) for i in open('labels.txt') ])

    emins = load(open(args.dumpfile, 'rb'))
    estimators = []
    for (filename, mins) in emins:
        E = sourmash.Estimators()
        E._mins = mins
        estimators.append((filename, E))

    for (filename, E) in estimators:
        S = sig.SourmashSignature('titus@idyll.org',
                                  E,
                                  filename=filename,
                                  name=labels.get(filename, ''))

        sigfile = filename + '.sig'
        print('opening', sigfile)
        fp = open(sigfile, 'w')
        fp.write(S.save())
        fp.close()


main()
