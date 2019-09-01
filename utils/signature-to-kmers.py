#! /usr/bin/env python3
"""
Given a signature file and a collection of sequences, output all of the
k-mers and sequences that match a hashval in the signature file.

NOTE: for now, only works for DNA.
"""
import sys
import argparse
import sourmash
from sourmash import MinHash
from sourmash import sourmash_args
from sourmash._minhash import hash_murmur
import screed
import csv
from sourmash.logging import notify, error


NOTIFY_EVERY_BP=1e7


def get_kmers_for_hashvals(sequence, hashvals, ksize):
    for start in range(0, len(sequence) - ksize + 1):
        kmer = sequence[start:start + ksize]
        hashval = hash_murmur(kmer)
        if hashval in hashvals:
            yield kmer, hashval


def main():
    p = argparse.ArgumentParser()
    p.add_argument('query') 					# signature file
    p.add_argument('seqfiles', nargs='+')		# sequence files from which to look for matches
    p.add_argument('--output-sequences', type=str, default=None,
                   help='save matching sequences to this file.')
    p.add_argument('--output-kmers', type=str, default=None,
                   help='save matching kmers to this file.')
    args = p.parse_args()

    # set up the outputs.
    seqout_fp = None
    if args.output_sequences:
        seqout_fp = open(args.output_sequences, 'wt')

    kmerout_fp = None
    if args.output_kmers:
        kmerout_fp = open(args.output_kmers, 'wt')
        kmerout_w = csv.writer(kmerout_fp)
        kmerout_w.writerow(['kmer', 'hashval'])

    if not (seqout_fp or kmerout_fp):
        error("No output options given!")
        return(-1)

    # first, load the signature and extract the hashvals
    sigobj = sourmash.load_one_signature(args.query)
    query_hashvals = set(sigobj.minhash.get_mins())

    # create an empty minhash object of the same type...
    empty_mh = sigobj.minhash.copy_and_clear()

    # track found kmers
    found_kmers = {}

    # now, iterate over the input sequences and output those that overlap
    # with hashes!
    n_seq = 0
    n = 0 # bp loaded
    m = 0 # bp in found sequences
    p = 0 # number of k-mers found
    watermark = NOTIFY_EVERY_BP
    for filename in args.seqfiles:
        for record in screed.open(filename):
            n += len(record.sequence)
            n_seq += 1
            while n >= watermark:
                sys.stderr.write('... {} {} {}\r'.format(n_seq, watermark, filename))
                watermark += NOTIFY_EVERY_BP

            empty_mh.add_sequence(record.sequence, force=True)
            if empty_mh:
                hashes = set(empty_mh.get_mins())
                is_match = hashes.intersection(query_hashvals)

                if is_match:
                    if kmerout_fp:
                        ksize = empty_mh.ksize
                        sequence = record.sequence.upper()

                        # now do the hard work of finding the matching k-mers!
                        for kmer, hashval in get_kmers_for_hashvals(sequence,
                                                                    is_match,
                                                                    ksize):
                            found_kmers[kmer] = hashval

                        p = len(found_kmers)

                    # write out sequence
                    if seqout_fp:
                        seqout_fp.write('>{}\n{}\n'.format(record.name,
                                        record.sequence))
                        m += len(record.sequence)

                empty_mh = empty_mh.copy_and_clear()

    if seqout_fp:
        notify('read {} bp, wrote {} bp in matching sequences', n, m)

    if kmerout_fp:
        for kmer, hashval in found_kmers.items():
            kmerout_w.writerow([kmer, str(hashval)])
        notify('read {} bp, found {} kmers matching hashvals', n,
               len(found_kmers))


if __name__ == '__main__':
    sys.exit(main())
