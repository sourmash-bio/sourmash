# Using sourmash: a practical guide

```{contents}
   :depth: 2
```

So! You've installed sourmash, run a few of the tutorials and commands,
and now you actually want to *use* it.  This guide is here to answer some
of your questions, and explain why we can't answer others.

(If you have additional questions, please [file an issue!](https://github.com/sourmash-bio/sourmash/issues))

## What k-mer size(s) should I use?

You can build signatures at a variety of k-mer sizes all at once, and
(unless you are working with very large metagenomes) the resulting
signature files will still be quite small.  

We suggest including k=31 and k=51.  k=51 gives you the most stringent
matches, and has very few false positives. k=31 may be more sensitive
at the genus level.

Why 31 and 51, specifically? To a large extent these numbers were
picked out of a hat, based on our reading of papers like the
[Metapalette paper (Koslicki and Falush, 2016](http://msystems.asm.org/content/1/3/e00020-16). You
could go with k=49 or k=53 and probably get very similar results to
k=51.  The general rule is that longer k-mer sizes are less prone to
false positives. But you can pick your own parameters.

One additional wrinkle is that we provide a number of
[precalculated databases](databases.md) at k=21, k=31, and k=51.
It is often convenient to calculate signatures at these sizes so that
you can use these databases.

You'll notice that all of the above numbers are odd.  That is to avoid
occasional minor complications from palindromes in numerical
calculations, where the forward and reverse complements of a k-mer are
identical.  This cannot happen if k is odd. It is not enforced by sourmash,
however, and it probably doesn't really matter.

(When we have blog posts or publications providing more formal
guidance, we'll link to them here!)

## What resolution should my signatures be / how should I create them?

sourmash supports two ways of choosing the resolution or size of
your signatures: using `num` to specify the maximum number of hashes,
or `scaled` to specify the compression ratio.  Which should you use?

We suggest calculating all your signatures using `-p scaled=1000`.
This will give you a compression ratio of 1000-to-1 while making it
possible to detect regions of similarity in the 10kb range.

For comparison with more traditional MinHash approaches like `mash`,
if you have a 5 Mbp genome and use `-p scaled=1000`, you will extract
approximately 5000 hashes. So a scaled of 1000 is equivalent to using
`-p num=5000` with mash on a 5 Mbp genome.

The difference between using `num` and `scaled` is in metagenome
analysis: fixing the number of hashes with `num` limits your ability to
detect rare organisms, or alternatively results in very large
signatures (e.g. if you use n larger than 10000).  `scaled` will scale
your resolution with the diversity of the metagenome.

You can read more about this in this blog post from the mash folk,
[Mash Screen: What's in my sequencing run?](https://genomeinformatics.github.io/mash-screen/) What
we do with sourmash and `scaled` is similar to the 'modulo hash'
mentioned in that blog post.

(Again, when we have formal guidance on this based on benchmarks, we'll
link to it here.)

## What kind of input data does sourmash work on?

sourmash has been used most extensively with Illumina read data sets
and assembled genomes, transcriptomes, and metagenomes.  The high error
rate of PacBio and Nanopore sequencing is problematic for k-mer based
approaches and we have not yet explored how to tune parameters for
this kind of sequencing.

On a more practical note, `sourmash sketch` will autodetect FASTA and
FASTQ formats, whether they are uncompressed, gzipped, or bzip2-ed.
Nothing special needs to be done.

## How should I prepare my data?

Raw Illumina read data sets should be k-mer abundance trimmed to get rid of
the bulk of erroneous kmers. We suggest a command like the following,
using [trim-low-abund from the khmer project](https://peerj.com/preprints/890/) --

```
trim-low-abund.py -C 3 -Z 18 -V -M 2e9 <all of your input read files>
```

This is safe to use on genomes, metagenomes, and transcriptomes.  If you
are working with large genomes or diverse metagenomes, you may need to
increase the `-M` parameter to use more memory.

See [the khmer docs for trim-low-abund.py](https://khmer.readthedocs.io/en/v2.1.2/user/scripts.html#trim-low-abund-py) and [the semi-streaming preprint](https://peerj.com/preprints/890/) for more information.

For high coverage genomic data, you can do very stringent trimming with
an absolute cutoff, e.g.

```
trim-low-abund.py -C 10 -M 2e9 <all of your input read files>
```

will eliminate all k-mers that appear fewer than 10 times in your data
set.  This kind of trimming will dramatically reduce your sensitivity
when working with metagenomes and transcriptomes, however, where there
are always real low-abundance k-mers present.

## Could you just give us the !#%#!$ command line?

Sorry, yes! See below.

### Calculating signatures for read files:

```
trim-low-abund -C 3 -Z 18 -V -M 2e9 input-reads-1.fq input-reads-2.fq ...
sourmash sketch dna -p scaled=1000,k=21,k=31,k=51 input-reads*.fq.abundtrim \
    --merge SOMENAME -o SOMENAME-reads.sig
```

The first command trims off low-abundance k-mers from high-coverage
reads; the second takes all the trimmed read files, subsamples k-mers
from them at 1000:1, and outputs a single merged signature named
'SOMENAME' into the file `SOMENAME-reads.sig`.

### Creating signatures for individual genome files:

```
sourmash sketch dna -p scaled=1000,k=21,k=31,k=51 *.fna.gz --name-from-first
```

This command creates signatures for all `*.fna.gz` files, and names
each signature based on the first FASTA header in each file (that's
what the option `--name-from-first` does). The signatures will be placed
in `*.fna.gz.sig`.

### Creating signatures from a collection of genomes in a single file:

```
sourmash sketch dna -p scaled=1000,k=21,k=31,k=51 file.fa --singleton
```

This creates signatures for all individual FASTA sequences in `file.fa`,
names them based on their FASTA headers, and places them all in a single
`.sig` file, `file.fa.sig`.  (This behavior is triggered by the option
`--singleton`, which tells sourmash to treat each individual sequence in
the file as an independent sequence.)

## How do I store and search collections of signatures?

sourmash supports a variety of signature loading and storage options for
flexibility.  If you have only a few hundred signatures, here are some
options -

* you can put all your signature files in a directory and search them all
  using the path to the directory.
* you can use `sourmash sig cat` to concatenate multiple signatures into a
  single file.
* you can compress any signature file using `gzip` and sourmash will
  load them.

If you have more than a few hundred genome signatures that you
regularly search, it might be worth creating an indexed database of
them that will support faster searches.

sourmash supports two types of indexed databases: Sequence Bloom
Trees, or SBTs; and reverse indices, or LCAs.  (You can read more
detail about their implementation and design considerations
[in Chapter 2 of Dr. Luiz Irber's thesis, "Efficient indexing of collections of signatures"](https://github.com/luizirber/phd/releases/download/2020.09.28/thesis.pdf).)

### Sequence Bloom Tree (SBT) indexed databases

Sequence Bloom Trees (SBTs) (see
[Solomon and Kingsford, 2016](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4804353/))
are on disk databases that support low-memory query of 10s-100s of
thousands of signatures.  They can be created using `sourmash index`.

SBTs are the lowest-memory way to run search or gather on a collection
of signatures. The tradeoff is that they may be quite large on disk,
because SBTs also contain intermediate nodes in the tree.  The default
way to store SBTs is in a Zip file, named `.sbt.zip`, that can be
built and searched directly from the command line.

### Reverse indexed (LCA) databases

Reverse indexed or LCA databases are *in-memory* databases that, once
loaded from disk, support fast search and gather across 10s of thousands
of signatures.  They can be created using `sourmash lca index` ([docs](command-line.md#sourmash-lca-index-build-an-lca-database))

LCA databases are currently stored in JSON files (that can be gzipped).
As these files get larger, the time required to load them from disk
can be substantial.

LCA databases are also currently (sourmash 2.0-4.0) the only databases
that support the inclusion of taxonomic information in the database,
and there is an associated collection of commands
[under `sourmash lca`](command.md#sourmash-lca-subcommands-for-taxonomic-classification).
However, they can also be used as regular indexed databases for search
and gather as above.

(These are called "LCA databases" because they originally were created
to support "lowest common ancestor" taxonomic analyses, e.g. like
Kraken; their functionality has evolved a lot since, but their name
hasn't changed to match!)
