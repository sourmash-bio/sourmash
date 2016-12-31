============================
``sourmash`` Python examples
============================

A very simple example: two k-mers
---------------------------------

Define two sequences:

>>> seq1 = "ATGGCA"
>>> seq2 = "AGAGCA"

Create two estimators using 3-mers, and add the sequences:

>>> import sourmash_lib
>>> E1 = sourmash_lib.Estimators(n=20, ksize=3)
>>> E1.add_sequence(seq1)

>>> E2 = sourmash_lib.Estimators(n=20, ksize=3)
>>> E2.add_sequence(seq2)

One of the 3-mers (out of 4) overlaps, so Jaccard index is 1/4:

>>> E1.jaccard(E2)
0.25

and of course the estimators match themselves:

>>> E1.jaccard(E1)
1.0

We can add sequences and query at any time --

>>> E1.add_sequence(seq2)
>>> x = E1.jaccard(E2)
>>> round(x, 3)
0.571

Consuming files
---------------

Suppose we want to create MinHash sketches from genomes --

>>> import glob, pprint
>>> genomes = glob.glob('data/GCF*.fna.gz')
>>> genomes = list(sorted(genomes))
>>> pprint.pprint(genomes)
['data/GCF_000005845.2_ASM584v2_genomic.fna.gz',
 'data/GCF_000006945.1_ASM694v1_genomic.fna.gz',
 'data/GCF_000783305.1_ASM78330v1_genomic.fna.gz']

We have to read them in (here using screed), but then they can be fed
into 'add_sequence' directly; here we set 'force=True' in ``add_sequence``
to ignore non-ACTGN characters.

(Note, just for speed reasons, we'll truncate the sequences to 50kb in length.)
  
>>> import screed
>>> estimators = []
>>> for g in genomes:
...     E = sourmash_lib.Estimators(n=500, ksize=31)
...     for record in screed.open(g):
...         E.add_sequence(record.sequence[:50000], True)
...     estimators.append(E)

And now the estimators can be compared against each other:

>>> import sys
>>> for i, e in enumerate(estimators):
...    _ = sys.stdout.write(genomes[i][:20] + ' ')
...    for j, e2 in enumerate(estimators):
...       x = e.jaccard(estimators[j])
...       _ = sys.stdout.write(str(round(x, 3)) + ' ')
...    _= sys.stdout.write('\n')
data/GCF_000005845.2 1.0 0.0 0.0 
data/GCF_000006945.1 0.0 1.0 0.0 
data/GCF_000783305.1 0.0 0.0 1.0 

Note that the comparisons are quite quick; most of the time is spent in
making the estimators, which can be saved and loaded easily.

Saving and loading signature files
----------------------------------

>>> from sourmash_lib import signature
>>> sig1 = signature.SourmashSignature('titus@idyll.org', estimators[0],
...                                    name=genomes[0][:20])
>>> with open('/tmp/genome1.sig', 'wt') as fp:
...   signature.save_signatures([sig1], fp)

Here, ``/tmp/genome1.sig`` is a YAML file that can now be loaded and
compared -- first, load:

>>> sigfp = open('/tmp/genome1.sig', 'rt')
>>> siglist = list(signature.load_signatures(sigfp))
>>> loaded_sig = siglist[0]

then compare:

>>> loaded_sig.estimator.jaccard(sig1.estimator)
1.0
>>> sig1.estimator.jaccard(loaded_sig.estimator)
1.0
