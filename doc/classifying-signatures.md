# Classifying signatures: `search`, `gather`, and `lca` methods.

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

See [the main sourmash
tutorial](http://sourmash.readthedocs.io/en/latest/tutorials.html#make-and-search-a-database-quickly)
for information on using `search` with and without `--containment`.

## Breaking down metagenomic samples with `gather` and `lca`

Neither search option (similarity or containment) is effective when
comparing or searching with metagenomes, which typically have a
mixture of many different genomes.  While you might use containment to
see if a query genome is present in one or more metagenomes, a common
question to ask is the reverse: **what genomes are in my metagenome?**

We have implemented two approaches in sourmash to do this.

One approach uses taxonomic information from e.g. GenBank to classify
individual k-mers, and then infers taxonomic distributions of
metagenome contents from the presence of these individual
k-mers. (This is the approach pioneered by
[Kraken](https://ccb.jhu.edu/software/kraken/) and used by many other tools.)
`sourmash lca` can be used to classify individual genome bins with
`classify`, or summarize metagenome taxonomy with `summarize`.  The
[sourmash lca tutorial](http://sourmash.readthedocs.io/en/latest/tutorials-lca.html)
shows how to use the `lca classify` and `lca summarize` commands, and also
provides guidance on building your own database.

The other approach, `gather`, breaks a metagenome down into individual
genomes based on greedy partitioning. Essentially, it takes a query
metagenome and searches the database for the most highly contained
genome; it then subtracts that match from the metagenome, and repeats.
At the end it reports how much of the metagenome remains unknown.  The
[basic sourmash
tutorial](http://sourmash.readthedocs.io/en/latest/tutorials.html#what-s-in-my-metagenome)
has some sample output from using gather with GenBank.  See the appendix at
the bottom of this page for more technical details.

Some benchmarking on CAMI suggests that `gather` is a very accurate
method for doing strain-level resolution of genomes. More on
that as we move forward!

## To do taxonomy, or not to do taxonomy?

By default, there is no structured taxonomic information available in
sourmash signatures or SBT databases of signatures.  Generally what
this means is that you will have to provide your own mapping from a
match to some taxonomic hierarchy.  This is generally the case when
you are working with lots of genomes that have no taxonomic
information.

The `lca` subcommands, however, work with LCA databases, which contain
taxonomic information by construction.  This is one of the main
differences between the `sourmash lca` subcommands and the basic
`sourmash search` functionality.  So the `lca` subcommands will generally
output structured taxonomic information, and these are what you should look
to if you are interested in doing classification.

The command `lca gather` applies the `gather` algorithm to search an
LCA database; it reports taxonomy.

It's important to note that taxonomy based on k-mers is very, very
specific and if you get a match, it's pretty reliable. On the
converse, however, k-mer identification is very brittle with respect
to evolutionary divergence, so if you don't get a match it may only mean
that the particular species isn't known.

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
computing signatures with `--track-abundance`. The abundance
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

If you compute your input signatures with `--track-abundance`, both
`sourmash gather` and `sourmash lca gather` will use that information
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

**Buyer beware:** There are substantial challenges in doing this kind
of analysis on real metagenomic samples, relating to genome representation
and strain overlap; see [this issue](https://github.com/dib-lab/sourmash/issues/461) for a discussion.

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
Benoit et al., 2016. Initial comparisons of metagenome similarity
approximations computed with sourmash to the output of simka suggest a
significant correlation.

**Implementation note:** Angular similarity searches cannot be done on
SBT or LCA databases currently; you have to provide lists of signature
files to `sourmash search` and `sourmash compare`.  sourmash will
provide a warning if you run `sourmash search` on an LCA or SBT with
an abundance-weighted query, and automatically apply `--ignore-abundance`.

## What commands should I use?

It's not always easy to figure that out, we know! We're thinking about
better tutorials and documentation constantly.

We suggest the following approach:

* build some signatures and do some searches, to get some basic familiarity
  with sourmash;

* explore the available databases;

* then ask questions [via the issue tracker](https://github.com/dib-lab/sourmash/issues) and we will do our best to help you out!

This helps us figure out what people are actually interested in doing, and
any help we provide via the issue tracker will eventually be added into the
documentation.

## Appendix: how `sourmash gather` works.

The sourmash gather algorithm works as follows:

* find the best match in the database, based on containment;
* subtract that match from the query;
* repeat.

The output below is the CSV output for a fictional metagenome.

The first column, `f_unique_to_query`, is the fraction of the database
match that is _unique_ with respect to the original query. It will
always decrease as you get more matches.

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

A few quick notes for the algorithmic folk out there --

* the key innovation for gather is that it looks for **groups** of
  k-mers in the databases, and picks the best group (based on
  containment). It does not treat k-mers individually.
* because of this, gather does not saturate as databases grow in size,
  and in fact should only become more sensitive and specific as we
  increase database size. (Although of course it may get a lot
  slower...)
