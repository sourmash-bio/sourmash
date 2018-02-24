# Using sourmash: a practical guide

So! You've installed sourmash, run a few of the tutorials and commands,
and now you actually want to *use* it.  This guide is here to answer some
of your questions, and explain why we can't answer others.

(If you have additional questions, please [file an issue!](https://github.com/dib-lab/sourmash/issues))

## What k-mer size(s) should I use?

You can build signatures at a variety of k-mer sizes all at once, and
(unless you are working with very large metagenomes) the resulting
signature files will still be quite small.  

We suggest including k=31 or k=51.  k=51 gives you the most stringent
matches, and has very few false positives. k=31 may be more sensitive
at the genus level.

Why 31 and 51, specifically? To a large extent these numbers were
picked out of a hat, based on our reading of papers like the
[Metapalette paper (Koslicki and Falush, 2016](http://msystems.asm.org/content/1/3/e00020-16). You
could go with k=49 or k=53 and probably get very similar results to
k=51.  The general rule is that longer k-mer sizes are less prone to
false positives. But you can pick your own parameters.

One additional wrinkle is that we provide a number of
[precomputed databases](databases.html) at k=21, k=31, and k=51.
It is often convenient to calculate signatures at these sizes so that
you can use these databases.

You'll notice that all of the above numbers are odd.  That is to avoid
occasional minor complications from palindromes in numerical
calculations, where the forward and reverse complements of a k-mer are
identical.  This cannot happen if k is odd. It is not enforced by sourmash,
however, and it probably doesn't really matter.

(When we have blog posts or publications providing more formal
guidance, we'll link to them here!)

## What scaled values should I use?

Right now it's a bit hand-wavy, but:

We suggest calculating all your signatures using `--scaled
1000`.  This will give you a compression ratio of 1000-to-1 while
making it possible to detect regions of similarity in the 10kb range.

For comparison with more traditional MinHash approaches like `mash`,
if you have a 5 Mbp genome and use `--scaled 1000`, you will extract
approximately 5000 hashes. So a scaled of 1000 is equivalent to using
`-n 5000` with mash.

Again, when we have formal guidance on this, we'll link to it here.

## What kind of input data does sourmash work on?

sourmash has been used most extensively with Illumina read data sets
and assembled genomes, transcriptomes, metagenomes.  The high error
rate of PacBio and Nanopore sequencing is problematic for k-mer based
approaches and we have not yet explored how to tune parameters for
this kind of sequencing.

On a more practical note, `sourmash compute` should autodetect FASTA,
FASTQ, whether they are uncompressed, gzipped, or bzip2-ed.  Nothing
special needs to be done.

## How should I prepare my data?

Raw Illumina read data sets should be k-mer abundance trimmed to get rid of
the bulk of erroneous kmers. We suggest a command like the following,
(using trim-low-abund from the khmer project @CTB link) --

```
trim-low-abund.py -C 3 -V -M 2e9 <all of your input read files>
```

This is safe to use on genomes, metagenomes, and transcriptomes.  If you
are working with large genomes or diverse metagenomes, you may need to
increase the `-M` parameter to use more memory.

See (link to preprint and khmer docs @CTB) for more information.

For high coverage genomic data, you can do very stringent trimming with
an absolute cutoff, e.g.

```
trim-low-abund.py -C 10 -M 2e9 <all of your input read files>
```

this will eliminate all k-mers that have lower than 10 abundance across
your data sets.  This kind of trimming will dramatically reduce your
sensitivity when working with metagenomes and transcriptomes, however, where
there are always real low-abundance k-mers present.
