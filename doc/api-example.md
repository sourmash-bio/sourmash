# `sourmash` Python API examples

All of sourmash's functionality is available via its [Python API](api.md). Below are both basic and advanced examples that use the API to accomplish common tasks.

```{contents}
   :depth: 2
```

## A first example: two k-mers

Define two sequences:

```
>>> seq1 = "ATGGCA"
>>> seq2 = "AGAGCA"

```

Create two MinHashes using 3-mers, and add the sequences:

```
>>> import sourmash
>>> mh1 = sourmash.MinHash(n=0, ksize=3, scaled=1)
>>> mh1.add_sequence(seq1)

```

```
>>> mh2 = sourmash.MinHash(n=0, ksize=3, scaled=1)
>>> mh2.add_sequence(seq2)

```

One of the 3-mers (out of 7) overlaps, so Jaccard index is 1/7:

```
>>> round(mh1.jaccard(mh2), 2)
0.14

```

and of course the MinHashes match themselves:

```
>>> mh1.jaccard(mh1)
1.0

```

We can add sequences to the MinHash objects and query at any time --

```
>>> mh1.add_sequence(seq2)
>>> x = mh1.jaccard(mh2)
>>> round(x, 3)
0.571

```

## Introduction: k-mers, molecule types, and hashing.

### DNA k-mers

The basis of sourmash is k-mers. Suppose we have a 35 base DNA sequence, and
we break it into four 31-mers:
```
>>> K = 31
>>> dnaseq = "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGC"
>>> for i in range(0, len(dnaseq) - K + 1):
...    kmer = dnaseq[i:i+K]
...    print(i, kmer)
0 ATGCGAGTGTTGAAGTTCGGCGGTACATCAG
1 TGCGAGTGTTGAAGTTCGGCGGTACATCAGT
2 GCGAGTGTTGAAGTTCGGCGGTACATCAGTG
3 CGAGTGTTGAAGTTCGGCGGTACATCAGTGG
4 GAGTGTTGAAGTTCGGCGGTACATCAGTGGC

```

sourmash uses a hash function (by default MurmurHash) that converts
each k-mer into 64-bit numbers. These numbers form the basis of everything
else sourmash does; the k-mer strings are not used internally at all.

The easiest way to access the hash
function is via the `seq_to_hashes` method on `MinHash` objects, which
returns a list:
```
>>> from sourmash import MinHash
>>> mh = MinHash(n=0, ksize=K, scaled=1)
>>> for i in range(0, len(dnaseq) - K + 1):
...    kmer = dnaseq[i:i+K]
...    print(i, kmer, mh.seq_to_hashes(kmer))
0 ATGCGAGTGTTGAAGTTCGGCGGTACATCAG [7488148386897425535]
1 TGCGAGTGTTGAAGTTCGGCGGTACATCAGT [3674733966066518639]
2 GCGAGTGTTGAAGTTCGGCGGTACATCAGTG [2135725670290847794]
3 CGAGTGTTGAAGTTCGGCGGTACATCAGTGG [14521729668397845245]
4 GAGTGTTGAAGTTCGGCGGTACATCAGTGGC [15919051675656106963]

```

Note that this is the same as using the MurmurHash hash function with a seed
of 42 and taking the first 64 bits.

Because DNA is double-stranded and has no inherent directionality, but
computers represent DNA with only one strand, it's important for
sourmash to represent both strands. sourmash does this by building a
canonical representation for each k-mer so that reverse-complement
sequences match to their forward sequence.

Underneath, sourmash DNA hashing does this by taking each k-mer,
building the reverse complement, choosing the lexicographically lesser of
the two, and then hashes it - for example, for the first and second
k-mers above, you get:

```
>>> from sourmash.minhash import hash_murmur
>>> from screed import rc
>>> kmer_1 = "ATGCGAGTGTTGAAGTTCGGCGGTACATCAG"
>>> hash_murmur(kmer_1)
7488148386897425535
>>> hash_murmur(kmer_1) == mh.seq_to_hashes(kmer_1)[0]
True
>>> kmer_2 = "TGCGAGTGTTGAAGTTCGGCGGTACATCAGT"
>>> hash_murmur(kmer_2)
6388498842523164783
>>> kmer_2_rc = rc(kmer_2)
>>> kmer_2_rc
'ACTGATGTACCGCCGAACTTCAACACTCGCA'
>>> hash_murmur(kmer_2_rc) == mh.seq_to_hashes(kmer_2)[0]
True

```
where the second k-mer's reverse complement starts with 'A' and is therefore
chosen for hashing by sourmash. This method was chosen to be compatible
with [mash](https://mash.readthedocs.io/.

### Protein-based encodings

By default, `MinHash` objects work with DNA. However, sourmash
supports amino acid, Dayhoff, and hydrophobic-polar (hp) encodings as
well. The Dayhoff and hp encodings support degenerate
matching that is less stringent than exact matches.

The simplest way to use a protein `MinHash` object is to create one and
call `add_protein` on it --

```
>>> K = 7
>>> mh = MinHash(0, K, is_protein=True, scaled=1)
>>> protseq = "MVKVYAPAS"
>>> mh.add_protein(protseq)

```

This creates three 7-mers from the sequence and hashes them:
```
>>> list(sorted(mh.hashes))
[5834377656419371297, 8846570680426381265, 10273850291677879123]

```

As with DNA k-mers, above, you can also use `seq_to_hashes` to generate
the hashes for protein k-mers, if you add the `is_protein=True` flag:
```
>>> for i in range(0, len(protseq) - K + 1):
...    kmer = protseq[i:i+K]
...    print(i, kmer, mh.seq_to_hashes(kmer, is_protein=True))
0 MVKVYAP [5834377656419371297]
1 VKVYAPA [10273850291677879123]
2 KVYAPAS [8846570680426381265]

```

In this case, the k-mers are always hashed in the forward direction
(because protein sequence always has an orientation, unlike DNA).

sourmash also supports the
[Dayhoff](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff#Table_of_Dayhoff_encoding_of_amino_acids)
and hydrophobic-polar encodings; here, amino acids are first mapped to
their encodings and then hashed. So, for example, the amino acid sequence
`CADHIF*` is mapped to `abcdef*` in the Dayhoff encoding:

```
>>> mh = MinHash(0, K, dayhoff=True, scaled=1)
>>> h1 = mh.seq_to_hashes('CADHIF*', is_protein=True)[0]
>>> h1
12061624913065022412
>>> h1 == hash_murmur('abcdef*')
True

```

### Translating DNA into protein-based encodings.

If you use `add_sequence(...)` to add DNA sequence to a protein encoding,
or call `seq_to_hashes(...)` on a protein encoding without `is_protein=True`,
sourmash will *translate* the sequences in all possible reading frames
and hash the translated amino acids. The k-mer size for the `MinHash`
is used as the k-mer size of the amino acids, i.e. 7 aa is 21 DNA bases.

```
>>> dnaseq = "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGC"
>>> len(dnaseq)
35
>>> K = 7
>>> mh = MinHash(n=0, ksize=K, is_protein=True, scaled=1)
>>> mh.add_sequence(dnaseq)
>>> len(mh)
30

```
Here, 30 hashes are added to the `MinHash` object because we are starting
with a 35 base DNA sequence, and using 21-mers of DNA (7-mer of protein),
which give us 15 distinct 21-mers in the three forward frames, and
15 distinct 21-mers in the three reverse-complement frames, for a total
of 30.

If a dayhoff or hp `MinHash` object is used, then `add_sequence(...)`
will first translate each sequence into protein space in all six frames,
and then further encode the sequences into Dayhoff or hp encodings.

### Invalid characters in DNA and protein sequences

sourmash detects invalid DNA characters (non-ACTG) and raises an
exception on `add_sequence`, unless `force=True`, in which case
k-mers containing invalid characters are ignored.
```
>>> dnaseq = "nTGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGC"
>>> K = 31
>>> mh = MinHash(n=0, ksize=K, scaled=1)
>>> mh.add_sequence(dnaseq)
Traceback (most recent call last):
    ...
ValueError: invalid DNA character in input k-mer: NTGCGAGTGTTGAAGTTCGGCGGTACATCAG
>>> mh.add_sequence(dnaseq, force=True)
>>> len(mh)
4

```

For protein sequences, sourmash does not currently do any invalid
character detection; k-mers are hashed as they are, and can only be
matched by an identical k-mer (with the same invalid character).
(Please [file an issue](https://github.com/dib-lab/sourmash/issues) if
you'd like us to change this!)
```
>>> K = 7
>>> mh = MinHash(n=0, ksize=K, is_protein=True, scaled=1)
>>> protseq = "XVKVYAPAS"
>>> mh.add_protein(protseq)
>>> len(mh)
3

```

For the Dayhoff and hp encodings on top of amino acids, invalid amino
acids (any letter for which no encoded character exists) are replaced with
'X'.
```
>>> K = 7
>>> mh = MinHash(n=0, ksize=K, dayhoff=True, scaled=1)
>>> protseq1 = ".VKVYAPAS"
>>> hashes1 = mh.seq_to_hashes(protseq1, is_protein=True)
>>> protseq2 = "XVKVYAPAS"
>>> hashes2 = mh.seq_to_hashes(protseq2, is_protein=True)
>>> hashes1 == hashes2
True

```

### Extracting both k-mers and hashes for a sequence

As of sourmash 4.2.2, `MinHash` objects provide a method called
`kmers_and_hashes` that will return the k-mers and their corresponding
hashes for an input sequence --

```
>>> mh = MinHash(n=0, ksize=31, scaled=1)
>>> dnaseq = "ATGCGAGTGTTGAAGTTCGGCGGTACATCAGTGGC"
>>> for kmer, hashval in mh.kmers_and_hashes(dnaseq):
...    print(kmer, hashval)
ATGCGAGTGTTGAAGTTCGGCGGTACATCAG 7488148386897425535
TGCGAGTGTTGAAGTTCGGCGGTACATCAGT 3674733966066518639
GCGAGTGTTGAAGTTCGGCGGTACATCAGTG 2135725670290847794
CGAGTGTTGAAGTTCGGCGGTACATCAGTGG 14521729668397845245
GAGTGTTGAAGTTCGGCGGTACATCAGTGGC 15919051675656106963

```

This works for protein `MinHash` objects as well, of course, although
you have to provide the `is_protein` flag, since `MinHash` objects assume
input sequence is DNA otherwise --

```
>>> K = 7
>>> mh = MinHash(n=0, ksize=K, is_protein=True, scaled=1)
>>> protseq = "XVKVYAPAS"
>>> for (kmer, hashval) in mh.kmers_and_hashes(protseq, is_protein=True):
...     print(kmer, hashval)
XVKVYAP 3140823561012061964
VKVYAPA 10273850291677879123
KVYAPAS 8846570680426381265

```

For translated `MinHash`, the k-mers and hashes corresponding to all
six reading frames are returned.

```
>>> dnaseq = "ATGCGAGTGTTGAAGTTCGGCGGTA"
>>> K = 7
>>> mh = MinHash(n=0, ksize=K, is_protein=True, scaled=1)
>>> for (kmer, hashval) in mh.kmers_and_hashes(dnaseq):
...    print(kmer, hashval)
ATGCGAGTGTTGAAGTTCGGC 16652503548557650904
CGAGTGTTGAAGTTCGGCGGT 9978056796243419534
TACCGCCGAACTTCAACACTC 2748622134668949083
CGCCGAACTTCAACACTCGCA 4263227699724621735
TGCGAGTGTTGAAGTTCGGCG 14299765336094039482
GAGTGTTGAAGTTCGGCGGTA 18155608748862746902
ACCGCCGAACTTCAACACTCG 14490181201772650983
GCCGAACTTCAACACTCGCAT 17205086974168937105
GCGAGTGTTGAAGTTCGGCGG 13354527969598897281
CCGCCGAACTTCAACACTCGC 16506504121672505595

```

In all cases, the k-mers are taken from the sequence itself, so the
k-mers will match to the input sequence, even when there are multiple
k-mers that hash to the same value (e.g. in the case of reverse
complements, or DNA k-mers that are translated to the same amino acid
sequence).

Note that sourmash also provides a `translate_codon` function if you
need to get the specific amino acids -

```
>>> from sourmash.minhash import translate_codon
>>> kmer = 'ATGCGAGT'
>>> for start in range(0, len(kmer) - 3 + 1, 3):
...    codon = kmer[start:start+3]
...    print(codon, translate_codon(codon))
ATG M
CGA R

```

### Summary

In sum,
* `MinHash.add_sequence(...)` converts DNA sequence into DNA or protein k-mers, and then hashes them and stores them.
* `MinHash.add_protein(...)` converts protein sequence into protein k-mers, and then hashes them and stores them.
* `MinHash.seq_to_hashes(...)` will give you the hash values without adding them to the `MinHash` object.
* `MinHash.kmers_and_hashes(...)` will provide tuples of `(kmer, hashval)` for an input sequence.
* The `dayhoff` and `hp` encodings can be calculated on amino acid k-mers as well, using `MinHash` objects.

Note that this is the code that is used by the command-line
functionality in `sourmash sketch`, so the results at the command-line
will match the results from the Python API.

## Set operations on hashes

All of the hashes in a `MinHash` object are available via the `hashes`
property:

```
>>> mh1 = sourmash.MinHash(n=0, ksize=3, scaled=1)
>>> seq1 = "ATGGCA"
>>> mh1.add_sequence(seq1)
>>> seq2 = "AGAGCA"
>>> mh1.add_sequence(seq2)
>>> list(mh1.hashes)
[1274996984489324440, 2529443451610975987, 3115010115530738562, 5059920851104263793, 5740495330885152257, 8652222673649005300, 18398176440806921933]

```

and you can easily do your own set operations with `.hashes` - e.g.
the following calculates the Jaccard similarity (intersection over union) of two 
```
>>> s1 = set(mh1.hashes)
>>> s2 = set(mh2.hashes)
>>> round(len(s1 & s2) / len(s1 | s2), 3)
0.571

```
However, the MinHash class also supports a number of basic operations - the following operations work directly on the hashes:
```
>>> combined = mh1 + mh2
>>> combined += mh1
>>> combined.remove_many(mh1.hashes)
>>> combined.add_many(mh2.hashes)

```

You can create an empty copy of a MinHash object with `copy_and_clear`:
```
>>> new_mh = mh1.copy_and_clear()

```

and you can also access the various parameters of a MinHash object directly as properties --
```
>>> mh1.ksize
3
>>> mh1.scaled
1
>>> mh1.num
0
>>> mh1.is_dna
True
>>> mh1.is_protein
False
>>> mh1.dayhoff
False
>>> mh1.hp
False
>>> mh1.moltype
'DNA'

```
see the "Advanced" section, below, for a more complete discussion of MinHash objects.

## Creating MinHash sketches programmatically, from genome files

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

(Note, just for speed reasons, we're truncating the sequences to 50kb in length.)

```
>>> import screed
>>> minhashes = []
>>> for g in genomes:
...     mh = sourmash.MinHash(n=500, ksize=31)
...     for record in screed.open(g):
...         mh.add_sequence(record.sequence[:50000], True)
...     minhashes.append(mh)

```

And now the resulting MinHash objects can be compared against each other:

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
building the minhashes.

## Plotting dendrograms and matrices

If you're interested in building comparison matrices and dendrograms,
please see the notebook
[Building plots from `sourmash compare` output](plotting-compare.md).

## Saving and loading signature files

Signature files encapsulate MinHashes in JSON, and provide a way to
wrap MinHash objects with some metadata (the name and filename). To save signatures, use `save_signatures` with a list of signatures and a Python file pointer:

```
>>> from sourmash import SourmashSignature, save_signatures
>>> from tempfile import mkdtemp
>>> sig1 = SourmashSignature(minhashes[0], name=genomes[0][:20])
>>> tempdir = mkdtemp(suffix = "temp")
>>> with open(tempdir + '/genome1.sig', 'wt') as fp:
...   save_signatures([sig1], fp)

```

Here, `genome1.sig` is a JSON file that can now be loaded and
compared -- first, load it using `load_one_signature`:

```
>>> from sourmash import load_one_signature
>>> loaded_sig = load_one_signature(tempdir + '/genome1.sig')

```

then compare:

```
>>> loaded_sig.jaccard(sig1)
1.0
>>> sig1.jaccard(loaded_sig)
1.0

```

There are two primary signature loading functions - `load_one_signature`, used above, which loads exactly one signature or else raises an exception; and the  powerful and more generic `load_file_as_signatures`, which takes in a filename or directory containing a collection of signatures and returns the individual signatures -- for example, you can load all of the signatures under the `tempdir` created above like so,

```
>>> loaded_sigs = list(sourmash.load_file_as_signatures(tempdir))

```

Both `load_file_as_signatures` and `load_one_signature` take molecule type and k-mer size selectors, e.g.
```
>>> loaded_sigs = load_one_signature(tempdir + '/genome1.sig', select_moltype='DNA', ksize=31)

```
will load precisely one signature containing a DNA MinHash created at k-mer size of 31.

## Going from signatures back to MinHash objects and their hashes -

Once you load a signature, you can go back to its MinHash object with
`.minhash`; e.g.

First, load two signatures:

```
>>> sigfile1 = 'tests/test-data/genome-s10.fa.gz.sig'
>>> sig1 = load_one_signature(sigfile1, ksize=21, select_moltype='DNA')
>>> sigfile2 = 'tests/test-data/genome-s11.fa.gz.sig'
>>> sig2 = load_one_signature(sigfile2, ksize=21, select_moltype='DNA')

```

Then, get the hashes, and (e.g.) calculate the union:

```
>>> hashes1 = set(sig1.minhash.hashes.keys())
>>> hashes2 = set(sig2.minhash.hashes.keys())
>>> hash_union = hashes1.union(hashes2)
>>> print(f'{len(hash_union)} hashes in union of {len(hashes1)} and {len(hashes2)}')
1000 hashes in union of 500 and 500

```

## Advanced features of sourmash MinHash objects - `scaled` and `num`

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
enable containment operations that are useful for metagenome analyses. The tradeoff is that unlike num MinHashes, they can become arbitrarily large.  The key parameter for modulo hash signatures is `scaled`, which specifies the average sampling rate
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

You cannot calculate Jaccard similarity or containment for
MinHash objects with different num or scaled values (or different ksizes):

```
>>> signum2 = sourmash.MinHash(n=1000, ksize=31)
>>> signum.jaccard(signum2)
Traceback (most recent call last):
  ...
TypeError: must have same num: 500 != 1000

```

However, you can make signatures compatible by downsampling; see the next
sections.

### A brief introduction to MinHash object methods and attributes

MinHash objects have the following methods and attributes:

* `ksize`, `num`, and `scaled` - the basic parameters used to create a MinHash object.
* `hashes` - retrieve all of the hashes contained in this object.
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
>>> hashvals = larger.hashes.keys()
>>> smaller = sourmash.MinHash(n=500, ksize=31)
>>> smaller.add_many(hashvals)
>>> len(smaller)
500

```

Also note that there's a convenience function that does the same thing,
faster!
```
>>> smaller2 = larger.downsample(num=500)
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
>>> small_scaled.add_many(large_scaled.hashes.keys())
>>> len(small_scaled)
69

```

And, again, there's a convenience function that you can use:
```
>>> small_scaled2 = large_scaled.downsample(scaled=500)
>>> small_scaled == small_scaled2
True

```

### Converting between `num` and `scaled` signatures

(Beware, these are confusing techniques for working with hashes that
are easy to get wrong! We suggest
[posting questions in the issue tracker](https://github.com/sourmash-bio/sourmash/issues)
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
>>> hashvals = num_mh.hashes.keys()

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
>>> hashvals = scaled_mh.hashes.keys()
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

## Working with indexed collections of signatures

If you want to search large collections of signatures, sourmash provides
two different indexing strategies, together with a generic `Index` class
that supports a common API for searching the collections.

The first indexing strategy is a Sequence Bloom Tree, which is
designed to support fast and efficient containment operations on large
collections of signatures.  SBTs are an _on disk_ search structure, so
they are a low-memory way to search collections.

To use SBTs from the command line, we first
need to create some `scaled` signatures:

```
sourmash sketch dna -p scaled=10000 data/GCF*.fna.gz --outdir data/
```

and then build a Sequence Bloom Tree (SBT) index with `sourmash
index`, like so:

```
sourmash index foo.sbt.zip data/GCF*.sig -k 31
```

Here, sourmash is storing the entire SBT in a single portable Zip file.

### Creating an on-disk SBT in Python

Let's start by using 'glob' to grab some example signatures from the
test data in the sourmash repository:

```
>>> import glob
>>> input_filenames = glob.glob('tests/test-data/doctest-data/GCF*.sig')

```

Now, create an SBT:

```
>>> import sourmash.sbtmh
>>> tree = sourmash.sbtmh.create_sbt_index()

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
checks, but the Python API doesn't.)

Now, save the tree:

```
>>> filename = tree.save(tempdir + '/test.sbt.zip')

```

### Loading and searching SBTs

How do we load the SBT and search it with a DNA sequence,
from within Python?

The SBT filename is `test.sbt.zip`, as above:
```
>>> SBT_filename = tempdir + '/test.sbt.zip'

```

and with it we can load the SBT:
```
>>> tree = sourmash.load_file_as_index(SBT_filename)

```

Now, load a DNA sequence:

```
>>> filename = 'data/GCF_000005845.2_ASM584v2_genomic.fna.gz'
>>> query_seq = next(iter(screed.open(filename))).sequence
>>> print(f'got {len(query_seq)} DNA characters to query')
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
>>> for similarity, found_sig, filename in tree.search(query_sig, threshold=0.1):
...    print(query_sig)
...    print(found_sig)
...    print(similarity)
my favorite query
NC_000913.3 Escherichia coli str. K-12 substr. MG1655, complete genome
1.0

```

et voila!

### In-memory databases: the LCA or "reverse index" database.

The LCA database lets you work with large collections of signatures in
memory.

The LCA database was initially designed to support individual hash
queries for taxonomic operations - hence its name, which stands for
"Lowest Common Ancestor." However, it supports all of the standard
`Index` operations, just like the SBT. 

First, let's create an LCA database programmatically.

```
>>> from sourmash.lca import LCA_Database
>>> db = LCA_Database(ksize=31, scaled=10000, moltype='DNA')

```

Now, let's load in all of the signatures from the test directory:

```
>>> for sig in sourmash.load_file_as_signatures('tests/test-data/doctest-data', ksize=31):
...    hashes_inserted = db.insert(sig)
...    print(f"Inserted {hashes_inserted} hashes into db.")
Inserted 493 hashes into db.
Inserted 490 hashes into db.
Inserted 525 hashes into db.

```

and now you have an `Index` class that supports all the generic index operations (below). You can save an LCA Database to disk with `db.save(filename)`, and load it with `sourmash.load_file_as_index`, below.

### The `Index` class API.

The `Index` class supports a generic API for SBTs, LCAs, and other collections of signatures.

To load an SBT or an LCA database from a file, use `sourmash.load_file_as_index`:
```
>>> sbt_db = sourmash.load_file_as_index('tests/test-data/prot/protein.sbt.zip')
>>> lca_db = sourmash.load_file_as_index('tests/test-data/prot/protein.lca.json.gz')

```

`Index` objects provide `search`, `insert`, `load`, `save`, and `__len__`. The signatures can be accessed directly via the  `.signatures()` method, which returns an iterable.  Last but not least, `Index.select(ksize=..., moltype=...)` will return a view on the Index object that contains only signatures with the desired k-mer size/molecule type.
