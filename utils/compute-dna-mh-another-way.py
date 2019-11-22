#! /usr/bin/env python
"""
Use the MurmurHash library mmh3 and separate Python code to calculate
a MinHash signature for input DNA sequence, as a way to do an
external check on our C++ implementation.

The output of this is used in test_sourmash.py to verify our C++ code.
"""

__complementTranslation = { "A": "T", "C": "G", "G": "C", "T": "A", "N": "N" }
def complement(s):
    """
    Return complement of 's'.
    """
    c = "".join(__complementTranslation[n] for n in s)
    return c


def reverse(s):
    """
    Return reverse of 's'.
    """
    r = "".join(reversed(s))
    return r


def kmers(seq, k):
    for start in range(len(seq) - k + 1):
        yield seq[start:start + k]

###

K = 21

import sys, screed
import mmh3
import sourmash
print('imported sourmash:', sourmash, file=sys.stderr)
import sourmash.signature

record = next(iter(screed.open(sys.argv[1])))
print('loaded', record.name, file=sys.stderr)
revcomp = reverse(complement((record.sequence)))

mh = sourmash.MinHash(ksize=K, n=500, is_protein=False)

#
# compute the actual hashes to insert by breaking down the sequence
# into k-mers and applying MurmurHash to each one; here, the only
# interesting thing that is done by add_hash is to keep only the
# (numerically) lowest n=500 hashes.
#
# this method of hash computation is exactly how sourmash does it
# internally, and should be approximately the same as what mash does.
#

for fwd_kmer in kmers(record.sequence, K):
    rev_kmer = reverse(complement(fwd_kmer))
    if fwd_kmer < rev_kmer:
        kmer = fwd_kmer
    else:
        kmer = rev_kmer

    hash = mmh3.hash64(kmer, seed=42)[0]

    # convert to unsigned int if negative
    if hash < 0:
        hash += 2**64

    mh.add_hash(hash)

s = sourmash.signature.SourmashSignature('', mh, name=record.name)
print(sourmash.signature.save_signatures([s]))
