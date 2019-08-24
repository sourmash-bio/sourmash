# `sourmash` API examples

## A first example: two k-mers


Define two sequences:

```
>>> seq1 = "ATGGCA"
>>> seq2 = "AGAGCA"

```

Create two MinHashes using 3-mers, and add the sequences:

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

and of course the MinHashes match themselves:

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
>>> loaded_sig = load_one_signature('/tmp/genome1.sig')

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

## sourmash MinHash objects and manipulations

sourmash supports two basic kinds of signatures, MinHash and modulo hash
signatures. MinHash signatures are equivalent to mash signatures;
they are limited in size, and very effective for comparing genomes and
other data sets that are of similar size. The key parameter for MinHash
signatures is `num`, which specifies the maximum number of hashes to
be collected for a given input data set.

```
>>> signum = sourmash.MinHash(n=500, ksize=31)

```
Because of this parameter, below we'll call them 'num' signatures.

Modulo hash (or 'scaled') signatures are specific to sourmash and they
enable an expanded range of metagenome analyses, with the downside
that they can become arbitrarily large.  The key parameter for modulo
hash signatures is `scaled`, which specifies the average sampling rate
for hashes for a given input data set.  A scaled factor  of 1000 means that,
on average, 1 in 1000 k-mers will be turned into a hash for later
comparisons; this is a sort of compression factor, in that a 5 Mbp
genome will yield approximately 5000 hash values with a scaled factor
of 1000 (5000 x 1000 = 5,000,000).

```
>>> sigscaled = sourmash.MinHash(n=0, ksize=31, scaled=1000)

```
Note also that with a scaled factor of 1, the signature will contain **all**
of the k-mers.


----

You can differentiate between num signatures and scaled signatures by
looking at the `num` and `scaled` attributes on a MinHash object:

```
>>> signum.num
500
>>> signum.scaled
0
>>> sigscaled.num
0
>>> sigscaled.scaled
1000

```

The MinHash class is otherwise identical between the two types of signatures.

Note that you cannot compute Jaccard similarity or containment for
MinHash objects with different num or scaled values (or different ksizes):

```
>>> signum2 = sourmash.MinHash(n=1000, ksize=31)
>>> signum.similarity(signum2)
Traceback (most recent call last):
  ...
TypeError: must have same num: 500 != 1000

```

You can make signatures compatible by downsampling; see the next
sections.

### A brief introduction to MinHash object methods and attributes

MinHash objects have the following methods and attributes:

* `ksize`, `num`, and `scaled` - the basic parameters used to create a MinHash object.
* `get_mins()` - retrieve all of the hashes contained in this object.
* `add_sequence(seq)` - hash sequence and add hash values.
* `add(hash)` and `add_many(hashvals)` - add hash values directly.
* `similarity(other)` - calculate Jaccard similarity with the other MinHash object.
* `contained_by(other)` - calculate the Jaccard containment of self by other.
* `copy_and_clear()` - make an empty copy of a MinHash object with the same parameters.
* `__len__()` - return the number of actual hash values. Note you can also do `len(mh)`, where `mh` is a MinHash object.

### Downsampling signatures

Num and scaled signatures can always be downsampled without referring
back to the original data.

Let's start by loading 50kb of genomic sequence in to memory:
```
>>> genomes = glob.glob('data/GCF*.fna.gz')
>>> genomes = list(sorted(genomes))
>>> genome = genomes[0]
>>> record = next(iter(screed.open(genome)))
>>> sequence = record.sequence[:50000]

```

Now, suppose we make a num signature limited to 1000 hashes:

```
>>> larger = sourmash.MinHash(n=1000, ksize=31)
>>> larger.add_sequence(sequence)
>>> len(larger)
1000

```

We can downsample this to 500 by extracting the hashes and using
`add_many` to add them to a new MinHash like so:

```
>>> hashvals = larger.get_mins()
>>> smaller = sourmash.MinHash(n=500, ksize=31)
>>> smaller.add_many(hashvals)
>>> len(smaller)
500

```

Also note that there's a convenience function that does the same thing,
faster!
```
>>> smaller2 = larger.downsample_n(500)
>>> smaller2 == smaller
True

```

The same can be done with scaled MinHashes:

```
>>> large_scaled = sourmash.MinHash(n=0, ksize=31, scaled=100)
>>> large_scaled.add_sequence(sequence)
>>> len(large_scaled)
459
>>> small_scaled = sourmash.MinHash(n=0, ksize=31, scaled=500)
>>> small_scaled.add_many(large_scaled.get_mins())
>>> len(small_scaled)
69

```

And, again, there's a convenience function that you can use:
```
>>> small_scaled2 = large_scaled.downsample_scaled(500)
>>> small_scaled == small_scaled2
True

```

### Converting between `num` and `scaled` signatures

(Beware, these are confusing techniques for working with hashes that
are easy to get wrong! We suggest
[posting questions in the issue tracker](https://github.com/dib-lab/sourmash/issues)
as you go, if you are interested in exploring this area!)

The hashing function used is identical between num and scaled signatures,
so the hash values themselves are compatible - it's the comparison between
collections of them that doesn't work.

But, in some circumstances, num signatures can be extracted from
scaled signatures, and vice versa.  We haven't yet implemented a
Python API for this in sourmash, but you can hack it together yourself
quite easily, and a conversion utility is implemented through the command
line in `sourmash signature downsample`.


To extract a num MinHash object from a scaled MinHash, first create or load
your MinHash, and then extract the hash values:
```
>>> num_mh = sourmash.MinHash(n=1000, ksize=31)
>>> num_mh.add_sequence(sequence)
>>> hashvals = num_mh.get_mins()

```

Now, create the new scaled MinHash object and add the hashes to it:

```
>>> scaled_mh = sourmash.MinHash(n=0, ksize=31, scaled=10000)
>>> scaled_mh.add_many(hashvals)

```

and you are done!

The same works in reverse, of course:
```
>>> scaled_mh = sourmash.MinHash(n=0, ksize=31, scaled=50)
>>> scaled_mh.add_sequence(sequence)
>>> hashvals = scaled_mh.get_mins()
>>> num_mh = sourmash.MinHash(n=500, ksize=31)
>>> num_mh.add_many(hashvals)

```

So... when can you do this extraction reliably?

You can extract num MinHashes from scaled MinHashes whenever the
maximum hash value in the num MinHash is greater than or equal to the
`max_hash` attribute of the scaled MinHash.

You can extract scaled MinHashes to num MinHashes whenever there are
more hash values in the scaled MinHash than num.

Yoda sayeth: *When understand these two sentences you can, use this code may
you.*

(You can also take a look at the logic in `sourmash signature
downsample` if you are interested.)

## Working with fast search trees (Sequence Bloom Trees, or SBTs)

Suppose we have a number of signatures calculated with `--scaled`, like so:

```
sourmash compute --scaled 10000 data/GCF*.fna.gz
```

and now we want to create a Sequence Bloom Tree (SBT) so that we can
search them efficiently.  You can do this with `sourmash index`, but
you can also access the Python API directly.

### Creating a search tree

Let's start by using 'glob' to grab some example signatures from the
test data in the sourmash repository:

```
>>> import glob
>>> input_filenames = glob.glob('tests/test-data/doctest-data/GCF*.sig')

```

Now, create a tree:

```
>>> import sourmash
>>> tree = sourmash.create_sbt_index()

```

Load each signature, and add it to the tree:

```
>>> from sourmash.sbtmh import SigLeaf
>>> for filename in input_filenames:
...     sig = sourmash.load_one_signature(filename, ksize=31)
...     leaf = SigLeaf(sig.md5sum(), sig)
...     tree.add_node(leaf)

```
(note, you'll need to make sure that all of the signatures are compatible
with each other! The `sourmash index` command does all of the necessary
checks.)

Now, save the tree:

```
>>> filename = tree.save('/tmp/test.sbt.json')

```

### Loading and search SBTs

How do we load the SBT and search it with a DNA sequence,
from within Python?

The SBT filename is `/tmp/test.sbt.json`, as above:
```
>>> SBT_filename = '/tmp/test.sbt.json'

```

and with it we can load the SBT:
```
>>> tree = sourmash.load_sbt_index(SBT_filename)

```

Now, load a DNA sequence:

```
>>> filename = 'data/GCF_000005845.2_ASM584v2_genomic.fna.gz'
>>> query_seq = next(iter(screed.open(filename))).sequence
>>> print('got {} DNA characters to query'.format(len(query_seq)))
got 4641652 DNA characters to query

```

and create a scaled signature:
```
>>> minhash = sourmash.MinHash(ksize=31, n=0, scaled=10000)
>>> minhash.add_sequence(query_seq)

>>> query_sig = sourmash.SourmashSignature(minhash, name='my favorite query')

```

Now do a search --

```
>>> threshold = 0.1
                                           
>>> for found_sig, similarity in sourmash.search_sbt_index(tree, query_sig, threshold):
...    print(query_sig.name())
...    print(found_sig.name())
...    print(similarity)
my favorite query
NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
1.0

```

et voila!
