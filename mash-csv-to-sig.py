#! /usr/bin/env python
import sys
import argparse
import csv
import sourmash_lib
from sourmash_lib import signature

def main():
    p = argparse.ArgumentParser()
    p.add_argument('mash_csvfile')
    p.add_argument('-o', '--output', type=argparse.FileType('wt'),
                   default=sys.stdout)
    args = p.parse_args()

    with open(args.mash_csvfile, 'r') as fp:
        reader = csv.reader(fp, delimiter=',')
        siglist = []
        for row in reader:
            hashfn = row[0]
            hashseed = int(row[1])

            assert hashfn == 'murmur64'
            assert hashseed == 42

            _, _, ksize, name, hashes = row
            ksize = int(ksize)

            hashes = hashes.strip()
            hashes = list(map(int, hashes.split(' ' )))

            e = sourmash_lib.Estimators(len(hashes), ksize)
            for h in hashes:
                e.mh.add_hash(h)
            sig = signature.SourmashSignature('titus@idyll.org',
                                              e, filename=name)
            siglist.append(sig)
            print('loaded signature:', name,
                  sig.md5sum()[:8], file=sys.stderr)

        print('saving %d signatures to YAML' % (len(siglist),),
              file=sys.stderr)
        yaml = signature.save_signatures(siglist, args.output)


if __name__ == '__main__':
    main()
