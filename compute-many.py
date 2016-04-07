#! /usr/bin/env python
import sourmash
import argparse
import screed
from pickle import dump, load

LEN_CUTOFF=80

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-p', '--picklejar', default='picklejar')
    args = parser.parse_args()

    already_done = set()
    estimators = []

    try:
        emins = load(open(args.picklejar, 'rb'))
        for (filename, mins) in emins:
            already_done.add(filename)
    except IOError:
        pass

    for filename in args.filenames:
        if filename in already_done:
            print('skipping', filename)
            continue

        #for n, record in enumerate(screed.open(filename)):
        #    if len(record.sequence) < LEN_CUTOFF:
        #        print 'skipping - too short reads', filename
                
            
        E = sourmash.Estimators()
        for n, record in enumerate(screed.open(filename)):
            if n % 10000 == 0:
                print('...', filename, n)
            E.add_sequence(record.sequence)

        estimators.append((filename, E._mins))

    fp = open(args.picklejar, 'wb') #, encoding='utf-8')
    dump(estimators, fp)
    fp.close()
    

if __name__ == '__main__':
    main()
    
