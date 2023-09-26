# Classifying signatures: `search`, `gather`, and `lca` methods.

```{contents} Contents
:depth: 3
```

sourmash provides several different techniques for doing
classification and breakdown of signatures.

## Searching for similar samples with `search`.

The `sourmash search` command is most useful when you are looking for
high similarity matches to other signatures; this is the most basic use
case for MinHash searching.  The command takes a query signature and one
or more search signatures, and finds all the matches it can above a particular
threshold.

By default `search` will find matches with high [*Jaccard
similarity*](https://en.wikipedia.org/wiki/Jaccard_index), which will
consider all of the k-mers in the union of the two samples.
Practically, this means that you will only find matches if there is
both high overlap between the samples *and* relatively few k-mers that
are disjoint between the samples.  This is effective for finding genomes
or transcriptomes that are similar but rarely works well for samples
of vastly different sizes.

One useful modification to `search` is to calculate containment with
`--containment` instead of the (default) similarity; this will find
matches where the query is contained within the subject, but the
subject may have many other k-mers in it. For example, if you are using
a plasmid as a query, you would use `--containment` to find genomes
that contained that plasmid.  `gather` (discussed below) uses containment
analysis only.

See [the main sourmash tutorial](tutorial-basic.md#make-and-search-a-database-quickly)
for information on using `search` with and without `--containment`.

## Analyzing metagenomic samples with `gather`

Neither search option (similarity or containment) is effective when
comparing or searching with metagenomes, which typically contain a
mixture of many different genomes.  While you might use containment to
see if a query genome is present in one or more metagenomes, a common
question to ask is the reverse: **what genomes are in my metagenome?**
An alternative phrasing is this: **what reference genomes should I map
my metagenomic reads to?**

The main approach we provide in sourmash is `sourmash gather`. This
constructs the shortest possible list of reference genomes that cover
all of the known k-mers in a metagenome. We call this a *minimum
metagenome cover*.

From an algorithmic perspective, `gather` generates a minimum set
cover for a query metagenome, using the reference database you give
it.  The minimum set cover is calculated using a greedy approximation
algorithm.  Essentially, `gather` takes a query metagenome and
searches the database for the most highly contained genome; it then
subtracts that match from the metagenome, and repeats.  At the end it
reports how much of the metagenome remains unknown.  The
[basic sourmash tutorial](tutorial-basic.md#whats-in-my-metagenome)
has some sample output from using gather with GenBank.  See Appendix A
at the bottom of this page for more technical details.

The `gather` method is described in
[Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers, Irber et al., 2022](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2).
Our benchmarking in that paper and also in
[Evaluation of taxonomic classification and profiling methods for long-read shotgun metagenomic sequencing datasets, Portik et al., 2022](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05103-0)
suggests that it is a very sensitive and specific method for
analyzing metagenomes.

## Taxonomic profiling with sourmash

sourmash supports two basic kinds of taxonomic profiling, under the
`lca` and `tax` modules. **As of 2023, we strongly recommend the
`tax`-based profiling approach.**

But first, let's back up! By default, there is no structured taxonomic
information available in sourmash signatures or collections.  What
this means is that you will have to provide your own mapping from a
match to some taxonomic hierarchy.  Both the `lca` and `tax` modules
support identifier-based taxonomic mappings, in which identifiers
from the signature names can be linked to the standard seven NCBI/GTDB
taxonomy ranks - superkingdom, phylum, class, order, family, genus, and
species. These are typically provided in a spreadsheet _separately_ from
the signature collection. The `tax` module also supports `lins` taxonomies,
for which we have a tutorial.

There are several advantages that this approach affords sourmash. One
is that sourmash is not tied closely to a specific taxonomy - you can
use either GTDB or NCBI as you wish. Another advantage is that you can
create your own custom taxonomic ranks and even use them with private
databases of genomes to classify your own metagenomes.

The main disadvantage of sourmash's approach to taxonomy is that
sourmash doesn't classify individual metagenomic reads to either a
genome or a taxon. (Note that we're not sure this can be done robustly
in practice - neither short nor long reads typically contain enough
information to uniquely identify a single genome, especially if there
are many genomes from the same species present in the database.)  If
you want to do this, we suggest running `sourmash gather` first, and
then mapping the reads to the matching genomes; then you can determine
which read maps to which genome. This is the approach taken by
[the genome-grist pipeline](https://dib-lab.github.io/genome-grist/).

<!-- link to tutorials and examples -->

### Using `sourmash tax` to do taxonomic analysis

We recommend using the `tax` module to do taxonomic classification of
genomes (with `tax genome`) and metagenomes (with `tax metagenome`).
The `tax` module commands operate downstream of `sourmash gather`,
which builds a minimum set cover of the query against the database -
intuitively speaking, this is the shortest possible list of genomes
that the query would map to.  Then, both `tax genome` and `tax
metagenome` take the CSV output of `sourmash gather` and produce
taxonomic profiles.  (You can read more about minimum set covers
in
[Lightweight compositional analysis of metagenomes with FracMinHash and minimum metagenome covers, Irber et al., 2022](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2).)

The `tax metagenome` approach was benchmarked in
[Evaluation of taxonomic classification and profiling methods for long-read shotgun metagenomic sequencing datasets, Portik et al., 2022](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-022-05103-0)
and appears to be both very accurate and very sensitive, unless you're
using Nanopore data or other data types that have a high sequencing
error rate.

It's important to note that taxonomy based on multiple k-mers is very,
very specific and if you get a match, it's pretty reliable. On the
converse, however, k-mer identification is very brittle with respect
to evolutionary divergence, so if you don't get a match it may only
mean that the specific species or genus you're searching for isn't in
the database.

### Using `sourmash lca` to do taxonomic analysis

The `sourmash lca` module supports taxonomic classification using
single hashes, corresponding to single k-mers, in an approach inspired
by Kraken. Briefly, you first build an LCA database using `lca index`,
which takes a taxonomy spreadsheet and a collection of sketches.  Then,
you can use `lca classify` to classify single-genome sketches or 
`lca summarize` to classify metagenomes.

The `lca` approach is not published anywhere, but we're happy to discuss
it in more detail; just [post to the issue tracker](https://github.com/sourmash-bio/sourmash/issues).

While we do not recommend the `lca` approach for general taxonomic
classification purposes (see below!), it remains useful for certain
kinds of diagnostic evaluation of sequences, so we are leaving the
functionality in sourmash.

### `sourmash tax` vs `sourmash lca`

Why do we recommend using the `tax` module over `lca`? `sourmash lca`
was the first implementation in sourmash, and over the years we've
found that it is prone to false positives: that is, individual k-mers
are very sensitive but are often misassigned to higher taxonomic ranks
than they need to be, either because of contamination in the reference
database or because the taxonomy is not based directly on genome
similarity.  Instead of using single k-mers, `sourmash gather` estimates
the best matching genome based on combinations of k-mers, which is much
more specific than the LCA approach; only then is a taxonomy assigned
using `sourmash tax`.

The bottom line is that in our experience, `sourmash tax` is as
sensitive as `lca`, and a lot more specific. Please let us know if you
discover differently!

## Abundance weighting

By default, sourmash tracks k-mer presence, *not* their abundance. The
proportions and fractions reported also ignore abundance. So, if
`sourmash gather` reports that a genome is 5% of a metagenome, it is
reporting Jaccard containment of that genome in the metagenome, and it
is ignoring information like the number of reads in the metagenome
that come from that genome.  Similarly, when `sourmash compare`
compares genome or metagenome signatures, it's reporting Jaccard
similarity *without* abundance.

However, it is possible to take into account abundance information by
computing signatures with `-p abund`. The abundance
information will be used if it's present in the signature, and it can
be ignored with `--ignore-abundance` in any signature comparison.

There are two ways that abundance weighting can be used. One is in
containment queries for metagenomes, e.g. with `sourmash
gather`, and the other is in comparisons of abundance-weighted signatures,
e.g. with `sourmash search` and `sourmash compare`.  Below, we refer to the
first as "abundance projection" and the second as "angular similarity".

### Projecting abundances in `sourmash gather`:

`sourmash gather` can report approximate abundance information for
containment queries against genome databases.  This will give you
numbers that (approximately) match what you get from counting mapped
reads.

If you create your input signatures with `-p abund`,
`sourmash gather` will use that information
to calculate an abundance-weighted result.  This will weight
each match to a hash value by the multiplicity of the hash value in
the query signature.  You can turn off this behavior with
`--ignore-abundance`.

For example, if you have a metagenome composed of two equal sized genomes
A and B, with A present at 10 times the abundance of B, `gather` on
abundance-weighted signatures will report that approximately 91% of the
metagenome is A and approximately 9% is B. (If you use `--ignore-abundance`,
then `gather` will report approximately 50:50, since the genomes are equal
sized.)

You can also get count-like information from the CSV output of `sourmash
gather`; the column `median_abund` contains the median abundance of the k-mers
in the match to the given genome.

Please see Appendix B, below, for some actual numbers and output.

**Buyer beware:** There are substantial challenges in doing this kind
of analysis on real metagenomic samples, relating to genome representation
and strain overlap; see [this issue](https://github.com/sourmash-bio/sourmash/issues/461) for a discussion.

### Computing signature similarity with angular similarity.

If signatures that have abundance information are compared with
`sourmash search` or `sourmash compare`, the default comparison is
done with
[angular similarity](https://en.wikipedia.org/wiki/Cosine_similarity#Angular_distance_and_similarity). This
is a distance metric based on cosine similarity, and it is suitable
for use in clustering.

For more information on the value of this kind of comparison for
metagenomics, please see the simka paper,
[Multiple comparative metagenomics using multiset k-mer counting](https://peerj.com/articles/cs-94/),
Benoit et al., 2016.

**Implementation note:** Angular similarity searches cannot be done on
SBT or LCA databases currently; you have to provide lists of signature
files to `sourmash search` and `sourmash compare`.  sourmash will
provide a warning if you run `sourmash search` on an LCA or SBT with
an abundance-weighted query, and automatically apply `--ignore-abundance`.

### Estimating ANI from FracMinHash comparisons.

As of v4.4, `sourmash` can estimate Average Nucleotide Identity (ANI)
between two FracMinHash ("scaled") sketches. `sourmash compare` can now
produce a matrix of ANI values estimated from Jaccard, Containment,
or Max Containment by specifying `--ani` (optionally along with search type,
e.g. `--containment`). `sourmash search`, `sourmash prefetch`, and
`sourmash gather` will now output ANI estimates to output CSVs.

Note that while ANI can be estimated from either the Jaccard Index or
the Containment Index, ANI from Containment is preferable (more accurate).
For `sourmash search`, `sourmash prefetch`, and `sourmash gather`, you can
optionally return confidence intervals around containment-derived ANI estimates,
which take into account the impact of the scaling factor (via `--estimate-ani-ci`).

For details on ANI estimation, please see our preprint "Debiasing FracMinHash and
deriving confidence intervals for mutation rates across a wide range of evolutionary
distances," [here](https://www.biorxiv.org/content/10.1101/2022.01.11.475870v2),
Hera et al., 2022.

## What commands should I use?

It's not always easy to figure that out, we know! We're thinking about
better tutorials and documentation constantly.

We suggest the following approach:

* build some signatures and do some searches, to get some basic familiarity
  with sourmash;

* explore the available databases;

* then ask questions [via the issue tracker](https://github.com/sourmash-bio/sourmash/issues) and we will do our best to help you out!

This helps us figure out what people are actually interested in doing, and
any help we provide via the issue tracker will eventually be added into the
documentation.

## Appendix A: how `sourmash gather` works.

The sourmash gather algorithm works as follows:

* find the best match in the database, based on containment;
* subtract that match from the query;
* repeat.
* when the number of shared hashes between the _remaining_ query and the
  best match drops below `threshold_bp/scaled` (or is zero), break out of
  the loop.

The output below is the CSV output for a fictional metagenome.

The first column, `f_unique_to_query`, is the fraction of the database
match that is _unique_ with respect to the original query. It will
always decrease as you get more matches. The sum of
`f_unique_to_query` across all rows is what is reported in by gather
as the fraction of query k-mers hit by the recovered matches
(unweighted) and should never be greater than 1!

The second column, `f_match_orig`, is how much of the match is in the
_original_ query.  For this fictional metagenome, each match is
entirely contained in the original query. This is the number you would
get by running `sourmash search --containment <match> <metagenome>`.

The third column, `f_match`, is how much of the match is in the remaining
query metagenome, after all of the previous matches have been removed.

The fourth column, `f_orig_query`, is how much of the original query
belongs to the match. This is the number you'd get by running
`sourmash search --containment <metagenome> <match>`.

```
f_unique_to_query      f_match_orig  f_match                f_orig_query
0.3321964529331514     1.0           1.0                    0.3321964529331514
0.13096862210095497    1.0           1.0                    0.13096862210095497
0.11527967257844475    1.0           0.898936170212766      0.12824010914051842
0.10709413369713507    1.0           1.0                    0.10709413369713507
0.10368349249658936    1.0           0.3134020618556701     0.33083219645293316
```
Where there are overlapping matches (e.g. to multiple
*E. coli* species in a gut metagenome) the overlaps will be represented
multiple times in this column.

A few quick notes for the algorithmic folk out there --

* the key innovation for gather is that it looks for **groups** of
  k-mers in the databases, and picks the best group (based on
  containment). It does not treat k-mers individually.
* because of this, gather does not saturate as databases grow in size,
  and in fact should only become more sensitive and specific as we
  increase database size. (Although of course it may get a lot
  slower...)

## Appendix B: sourmash gather and signatures with abundance information

Below is a discussion of a synthetic set of test cases using three
randomly generated (fake) genomes of the same size, with two even read
data set abundances of 2x each, and a third read data set with 20x.

### Data set prep

First, we make some synthetic data sets:

* r1.fa with 2x coverage of genome s10
* r2.fa with 20x coverage of genome s10.
* r3.fa with 2x coverage of genome s11.

then we make signature s10-s11 with r1 and r3, i.e. 1:1 abundance, and
make signature s10x10-s11 with r2 and r3, i.e. 10:1 abundance.

### A first experiment: 1:1 abundance.

When we project r1+r3, 1:1 abundance, onto s10, s11, and s12 genomes
with gather:

```
sourmash gather r1+r3 genome-s10.sig genome-s11.sig genome-s12.sig
```

we get:

```
overlap     p_query p_match avg_abund
---------   ------- ------- ---------
394.0 kbp     49.6%   78.5%       1.8    genome-s10.fa.gz
376.0 kbp     50.4%   80.0%       1.9    genome-s11.fa.gz
```

* approximately 50% of each query matching (first column, `p_query`)
* approximately 80% of subject genome's contents being matched (this is due to the low coverage of 2 used to build queries; `p_match`)
* approximately 2.0 abundance (third column, `avg_abund`)
* no match to genome s12.

### A second experiment: 10:1 abundance.

When we project r2+r3, 10:1 abundance, onto s10, s11, and s12 genomes
with gather:

```
sourmash gather r2+r3 genome-s10.sig genome-s11.sig genome-s12.sig
```

we get:

```
overlap     p_query p_match avg_abund
---------   ------- ------- ---------
0.5 Mbp       91.0%  100.0%      14.5    tests/test-data/genome-s10.fa.gz
376.0 kbp      9.0%   80.0%       1.9    tests/test-data/genome-s11.fa.gz
```

* approximately 91% of s10 matching
* approximately 9% of s11 matching
* approximately 100% of the high coverage genome being matched, with only 80% of the low coverage genome
* approximately 2.0 abundance (third column, avg_abund) for s11, and (very) approximately 20x abundance for genome s10.

The cause of the poor approximation for genome s10 is unclear; it
could be due to low coverage, or small genome size, or something
else. More investigation needed.
