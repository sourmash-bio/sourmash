#! /usr/bin/env python
import sourmash
import argparse
import screed
import sig
import os

LEN_CUTOFF=80

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    args = parser.parse_args()

    for filename in args.filenames:
        sigfile = filename + '.sig'
        if os.path.exists(sigfile):
            print('skipping', filename, '- already done')
            continue

        E = sourmash.Estimators()
        for n, record in enumerate(screed.open(filename)):
            if n % 10000 == 0:
                print('...', filename, n)
            E.add_sequence(record.sequence)

        signature = sig.SourmashSignature('titus@idyll.org',
                                          E,
                                          filename=filename)
        siglist = [signature]
        data = sig.save_signatures(siglist)
        fp = open(sigfile, 'w')
        fp.write(data)
        fp.close()
    

if __name__ == '__main__':
    main()
    
