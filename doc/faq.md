# Frequently Asked Questions

## How is sourmash different from mash?

mash is an awesome piece of software that inspired sourmash - hence the
name we chose for sourmash! But are they the same? No, they are not!

mash is based on MinHash sketching, a technique for randomly selecting
a very small set of k-mers to represent a potentially very large
sequence. MinHash sketches can be used to estimate Jaccard similarity
(a distance metric that lets you find closely related sequences)
with very high accuracy, under one condition - that the two sequences
be of similar size.

So, MinHash sketching is great when comparing bacterial genomes, which are
all around 2-10 MB in size.

But... MinHash sketching doesn't work when comparing _metagenomes_ to
genomes, because metagenomes are usually _much_ larger than genomes.

So we built sourmash around a different kind of sketch - FracMinHash -
which can estimate Jaccard similarity between sets of very different
sizes. FracMinHash sketches also support overlap and containment analysis,
and are convenient in a variety of other ways - see our paper on
FracMinHash, [Lightweight compositional analysis of metagenomes with
FracMinHash and minimum metagenome covers](@@).

There are some drawbacks to FracMinHash sketches - read on!

## What are the drawbacks to FracMinHash and sourmash?

There are two drawbacks to FracMinHash relative to MinHash. Neither
one is a showstopper, but you should know about them!

One drawback is that FracMinHash sketches can be a _lot_ bigger than
MinHash sketches - technically speaking, MinHash sketches are bounded
in size (you pick the number of hashes to keep!) while FracMinHash
sketches can get arbitrarily large. In practice this means that when
you're sketching large metagenomes, you end up with large sketches.
How large depends on the _cardinality_ of the metagenome k-mers - if a
metagenome has 1 billion unique k-mers, a FracMinHash sketch with scaled
parameter of 1000 will contain a million hashes.

The other drawback is that FracMinHash sketches _don't work well_ for
very small sequences. Our default parameter choice for DNA
(scaled=1000) works well for finding 10 kb or larger matches between
sequences - some simple Poisson matching math suggests that about 99.98%
of 10kb overlaps will be found with scaled=1000.

## What k-mer size(s) should I use with sourmash?

The short answer is: for DNA, use k=31.

Slightly longer answer: when we look at the k-mer distribution
across all of the bacterial genomes in GTDB, we find that 99% or
more of 31-mers are _genome_, _species_ or _genus_ specific.

If you go lower (say, k=21), then you get a few percent of k-mers
that match above the genus level - family or above.

If you go higher (k=51), a higher percentage of k-mers are genome-specific.

For the core sourmash operations - search, gather, and compare -
we believe (with evidence!) that (a) the differences between k=21,
k=31, and k=51 are negligible; and that (b) k=31 works fine for
most day-to-day use of sourmash.

We also provide [Genbank and GTDB databases](@@) for k=21, k=31, and
k=51, too.

## How do k-mer-based analyses compare with read mapping?

