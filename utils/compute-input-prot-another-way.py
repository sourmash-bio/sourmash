#! /usr/bin/env python
"""
Use the MurmurHash library mmh3 and separate Python code to calculate
a MinHash signature for input protein sequence, as a way to do an
external check on our C++ implementation.

The output of this is used in test_sourmash.py to verify our C++ code.
"""

dna_to_aa={'TTT':'F','TTC':'F', 'TTA':'L','TTG':'L',
                'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
                'TAT':'Y','TAC':'Y', 'TAA':'*','TAG':'*','TGA':'*',
                'TGT':'C','TGC':'C', 'TGG':'W',
                'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
                'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
                'CAT':'H','CAC':'H', 'CAA':'Q','CAG':'Q',
                'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
                'ATT':'I','ATC':'I','ATA':'I', 'ATG':'M',
                'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
                'AAT':'N','AAC':'N', 'AAA':'K','AAG':'K',
                'AGT':'S','AGC':'S', 'AGA':'R','AGG':'R',
                'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
                'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
                'GAT':'D','GAC':'D', 'GAA':'E','GAG':'E',
                'GGT':'G','GGC':'G','GGA':'G','GGG':'G'}


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


def peptides(seq, start):
    for i in range(start, len(seq), 3):
        yield dna_to_aa.get(seq[i:i+3], "X")


def translate(seq):
    for i in range(3):
        pep = peptides(seq, i)
        yield "".join(pep)

    revcomp = reverse(complement((seq)))
    for i in range(3):
        pep = peptides(revcomp, i)
        yield "".join(pep)

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

mh = sourmash.MinHash(ksize=K, n=500, is_protein=True)
prot_ksize = int(K / 3)

for kmer in kmers(record.sequence, prot_ksize):
    hash = mmh3.hash64(kmer, seed=42)[0]

    # convert to unsigned int if negative
    if hash < 0:
        hash += 2**64

    mh.add_hash(hash)

s = sourmash.signature.SourmashSignature('', mh, name=record.name)
print(sourmash.signature.save_signatures([s]))
