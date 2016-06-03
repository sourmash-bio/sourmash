#! /usr/bin/env python
import argparse
import csv
import sourmash_lib
from sourmash_lib import signature

def main():
    p = argparse.ArgumentParser()
    p.add_argument('mash_csvfile')
    args = p.parse_args()

    with open(args.mash_csvfile, 'r') as fp:
        reader = csv.reader(fp, delimiter=',')
        for row in reader:
            print(row[0])
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
            outfile = sig.md5sum()[:8] + '.sig'

            print('wrote signature for %s to %s' % (name, outfile))
            with open(outfile, 'wt') as outfp:
                yaml = signature.save_signatures([sig])
                outfp.write(yaml)


if __name__ == '__main__':
    main()
