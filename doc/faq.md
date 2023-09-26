# Frequently Asked Questions

```{contents} Contents
:depth: 3
```

## How is sourmash different from mash?

[mash](https://mash.readthedocs.org/) is an awesome piece of software
that inspired sourmash - hence the name we chose for sourmash! But are
they the same? No, they are not!

mash is based on
[MinHash sketching](https://en.wikipedia.org/wiki/MinHash), a
technique for randomly selecting a very small set of k-mers to
represent a potentially very large sequence. MinHash sketches can be
used to estimate Jaccard similarity (a distance metric that lets you
find closely related sequences) with very high accuracy, under one
condition - that the two sequences be of similar size.

So, MinHash sketching is great when comparing bacterial genomes, which are
all around 2-10 MB in size.

But... MinHash sketching doesn't work when comparing _metagenomes_ to
genomes, because metagenomes are usually _much_ larger than genomes.

So we built sourmash around a different kind of sketch - FracMinHash -
which can estimate Jaccard similarity between sets of very different
sizes. FracMinHash sketches also support overlap and containment
analysis, and are convenient in a variety of other ways - see our
paper on FracMinHash,
[Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2).

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
sequences - some simple Poisson matching math suggests that about
99.98% of 10kb overlaps will be found with scaled=1000. @CTB verify I
think this is 5kb.

## What k-mer size(s) should I use with sourmash?

The short answer is: for DNA, use k=31.

Slightly longer answer: when we look at the k-mer distribution
across all of the bacterial genomes in GTDB, we find that 99% or
more of 31-mers are _genome_, _species_ or _genus_ specific.

If you go lower (say, k=21), then you get a few percent of k-mers
that match above the genus level - family or above.

If you go higher (k=51), a higher percentage of k-mers are genome-specific.

For the core sourmash operations - search, gather, and compare - we
believe (with evidence!) that (a) the differences between k=21, k=31,
and k=51 are negligible; and that (b) k=31 works fine for most
day-to-day use of sourmash.

We also provide [Genbank and GTDB databases](databases.md) for k=21,
k=31, and k=51.

For some background on k-mer specificity, we recommend this paper:
[MetaPalette: a k-mer Painting Approach for Metagenomic Taxonomic Profiling and Quantification of Novel Strain Variation](https://journals.asm.org/doi/10.1128/msystems.00020-16),
Koslicki & Falush, 2016.

## How do k-mer-based analyses compare with read mapping?

tl;dr very well! But it's a bit one sided: if k-mers match, reads will
map, but not necessarily vice versa. So read mapping rates are almost always
higher than k-mer matching rates.

### Mapping reads to reference vs k-mer "detection" or containment

Let's start by looking at a simple example: suppose you have Illumina
shotgun sequencing of a new isolate, and you want to compare it to a
reference genome for a member of the same species.  You calculate
k-mer containment (using e.g. `sourmash search --containment
ref.sig.gz isolate.sig.gz` @@verify) of the reference genome in the isolate shotgun
sequence to be 65%. You then map the reads from the isolate data to the
reference genome. What should you expect to see?

What you should see is that 65% or more of the reference genome is
covered by at least one read. This is known as the mapping-based
"detection" of the reference genome in the read data set, and k-mer
detection typically *underestimates* mapping-based detection.

(If you want to know _how many of the reads will map_, you need to use
a different number that is output by `sourmash gather` - although, for
single genomes with only a few repeats, the percentage of reads that
map should match the k-mer detection number you calculate with
containment. Read the section below on metagenome analysis for more
information on read mapping rates!)

### Why do k-mers underestimate read mapping?

K-mers underestimate read mapping because k-mers rely on exact matches,
while read mapping tolerates mismatches. So unless you are looking at
entirely error free data with no strain variation, mismatches will
sometimes prevent k-mers from matching.

For some math: consider a k-mer of size 21. It will only match to a
specific 21 base sequence exactly. If you are matching two isolate genomes
that align completely but have a 95% average nucleotide identity, then
one out of every 20 bases will be different (on average), and, on average,
every 21-mer will be different! (In practice you'd expect about 1/3 of the
k-mers to match if the variation is random, since there will be many 21-base
windows that don't have any mismatches.) However, alignment-based approaches
like read mapping will happily align across single-base mismatches, and so
even at 95% ANI you would expect many reads to align.

In practice, the numbers will differ, but the intuition remains!

### How do read mapping rates for metagenomes compare with k-mer statistics?

Shotgun metagenome sequencing adds several complicating factors when
comparing metagenome reads to a reference database. First, each genome
will generally be present at a different abundance in the metagenome.
Second, metagenomes rarely contain the exact strain genome present in
a reference database; often the genomes in the metagenome differ both
in terms of SNPs _and_ accessory elements from what's in the reference
database. And third, metagenomes often contain complicated mixtures of
strains - even if one strain is dominant.

The `sourmash gather` method is built for analyzing metagenomes against
reference databases, and it does so by finding the shortest list of
reference genomes that "cover" all of the k-mers in the metagenome.
This list is arranged by how many k-mers from the metagenome are covered
by that entry in the output: the first match is the biggest, and the
second match is the second biggest, and so on. We call this a "minimum
metagenome cover" and it is described in the Irber et al., 2022, paper below.

When we construct the minimum metagenome cover, it correlates well with
mapping (per Irber et al., 2022), with one caveat: you need to assign
reads to the reference genomes in the rank order output by gather.
This is needed to properly assign reads that map to multiple genomes
to the "best" match - the one that will capture the most reads.

`sourmash gather` is still a research method but it seems to work
pretty well - in particular, it is both highly sensitive _and_ highly
specific in taxonomic benchmarking. Please ask questions as you have them!

### Further reading and links on k-mers and mapping:

* The paper
  [Biogeographic Distribution of Five Antarctic Cyanobacteria Using Large-Scale k-mer Searching with sourmash branchwater](https://www.biorxiv.org/content/10.1101/2022.10.27.514113v1),
  Lumian et al., 2022, shows that k-mer detection almost always
  underestimates mapping, and k-mer abundance analysis is always more
  conservative than mapping-based analyses.

* The paper
  [Deriving confidence intervals for mutation rates across a wide range of evolutionary distances using FracMinHash](https://genome.cshlp.org/content/33/7/1061),
  Rahman Hera et al., 2023, shows how to translate between average
  nucleotide identity (ANI) and k-mer statistics.

* The paper
[Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2),
Irber et al., 2022, describes how `sourmash gather` assigns k-mers
from metagenomes to a set of reference genomes, and shows that read
mapping correlates pretty well with k-mer overlap. Note that it is
focused on *systematic decomposition of metagenomes against reference
databases*, so it tackles the question of analyzing an entire
metagenome against all available references, not just a single
matching genome.
