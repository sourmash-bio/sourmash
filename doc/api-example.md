# `sourmash` Python examples

## A first example: two k-mers


Define two sequences:

```
>>> seq1 = "ATGGCA"
>>> seq2 = "AGAGCA"

```

Create two minhashes using 3-mers, and add the sequences:

```
>>> import sourmash
>>> E1 = sourmash.MinHash(n=20, ksize=3)
>>> E1.add_sequence(seq1)

```

```
>>> E2 = sourmash.MinHash(n=20, ksize=3)
>>> E2.add_sequence(seq2)

```

One of the 3-mers (out of 7) overlaps, so Jaccard index is 1/7:

```
>>> round(E1.jaccard(E2), 2)
0.14

```

and of course the minhashes match themselves:

```
>>> E1.jaccard(E1)
1.0

```

We can add sequences and query at any time --

```
>>> E1.add_sequence(seq2)
>>> x = E1.jaccard(E2)
>>> round(x, 3)
0.571

```

## Consuming files


Suppose we want to create MinHash sketches from genomes --

```
>>> import glob, pprint
>>> genomes = glob.glob('data/GCF*.fna.gz')
>>> genomes = list(sorted(genomes))
>>> pprint.pprint(genomes)
['data/GCF_000005845.2_ASM584v2_genomic.fna.gz',
 'data/GCF_000006945.1_ASM694v1_genomic.fna.gz',
 'data/GCF_000783305.1_ASM78330v1_genomic.fna.gz']

```

We have to read them in (here using screed), but then they can be fed
into `add_sequence` directly; here we set `force=True` in `add_sequence`
to skip over k-mers containing characters other than ACTG, rather than
raising an exception.

(Note, just for speed reasons, we'll truncate the sequences to 50kb in length.)

```
>>> import screed
>>> minhashes = []
>>> for g in genomes:
...     E = sourmash.MinHash(n=500, ksize=31)
...     for record in screed.open(g):
...         E.add_sequence(record.sequence[:50000], True)
...     minhashes.append(E)

```

And now the minhashes can be compared against each other:

```
>>> import sys
>>> for i, e in enumerate(minhashes):
...    _ = sys.stdout.write(genomes[i][:20] + ' ')
...    for j, e2 in enumerate(minhashes):
...       x = e.jaccard(minhashes[j])
...       _ = sys.stdout.write(str(round(x, 3)) + ' ')
...    _= sys.stdout.write('\n')
data/GCF_000005845.2 1.0 0.0 0.0 
data/GCF_000006945.1 0.0 1.0 0.0 
data/GCF_000783305.1 0.0 0.0 1.0 

```

Note that the comparisons are quite quick; most of the time is spent in
making the minhashes, which can be saved and loaded easily.

## Saving and loading signature files

Signature files encapsulate MinHashes in JSON, and provide a way to
add some metadata to MinHashes.

```
>>> from sourmash import SourmashSignature, save_signatures
>>> sig1 = SourmashSignature(minhashes[0], name=genomes[0][:20])
>>> with open('/tmp/genome1.sig', 'wt') as fp:
...   save_signatures([sig1], fp)

```

Here, `/tmp/genome1.sig` is a JSON file that can now be loaded and
compared -- first, load:

```
>>> from sourmash import load_one_signature
>>> sigfp = open('/tmp/genome1.sig', 'rt')
>>> loaded_sig = load_one_signature(sigfp)
```

then compare:

```
>>> loaded_sig.jaccard(sig1)
1.0
>>> sig1.jaccard(loaded_sig)
1.0

```

## Manipulating signatures and their hashes.

It is relatively straightforward to work directly with hashes.

First, load two signatures:

```
>>> sigfile1 = 'tests/test-data/genome-s10.fa.gz.sig'
>>> sig1 = load_one_signature(sigfile1, ksize=21, select_moltype='DNA')
>>> sigfile2 = 'tests/test-data/genome-s11.fa.gz.sig'
>>> sig2 = load_one_signature(sigfile2, ksize=21, select_moltype='DNA')

```

Then, get the hashes, and (e.g.) compute the union:

```
>>> hashes1 = set(sig1.minhash.get_mins())
>>> hashes2 = set(sig2.minhash.get_mins())
>>> hash_union = hashes1.union(hashes2)
>>> print('{} hashes in union of {} and {}'.format(len(hash_union), len(hashes1), len(hashes2)))
1000 hashes in union of 500 and 500

```

## Working with `scaled` signatures and SBTs.

Suppose we calculate signatures with a scaled value of 10,000, and
create a search tree from them like so:

```
sourmash compute --scaled 10000 -f data/GCF*.fna.gz
sourmash index -k 31 test.sbt GCF*.sig
```

How do we load the SBT and search it with a DNA sequence,
from within Python?

The SBT filename is `test.sbt`, as above:
```
>>> SBT_filename = test.sbt
```

and with it we can load the SBT:
```
>>> tree = sourmash_lib.load_sbt_index(SBT_filename)
```

Now, load a DNA sequence:

```
>>> query_seq = next(iter(screed.open(filename))).sequence
>>> print('got {} DNA characters to query'.format(len(query_seq)))
```

and create a signature:
```
>>> minhash = sourmash_lib.MinHash(ksize=KSIZE, n=0, scaled=10000)
>>> minhash.add_sequence(query_seq)

>>> query_sig = sourmash_lib.SourmashSignature('', minhash, name='my favorite query')
```

Now do a search --

```
>>> threshold = 0.1
                                           
>>> for found_sig, similarity in sourmash.search_sbt_index(tree, query_sig, threshold):
...    print((query_sig.name(), found_sig.name(), similarity,))
```
