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
that contained that plasmid.

See [the main sourmash
tutorial](http://sourmash.readthedocs.io/en/latest/tutorials.html#make-and-search-a-database-quickly)
for information on using `search` with and without `--containment`.

## Breaking down metagenomic samples with `gather` and `lca`

Neither search option (similarity or containment) is effective when
comparing or searching with metagenomes, which typically have a
mixture of many different genomes.  While you might use containment to
see if a query genome is present in one or more metagenomes, a common
question to ask is the reverse: **what genomes are in my metagenome?**

We have implemented two algorithms in sourmash to do this.

One algorithm uses taxonomic information from e.g. GenBank to classify
individual k-mers, and then infers taxonomic distributions of
metagenome contents from the presence of these individual
k-mers. (This is the approach pioneered by
[Kraken](https://ccb.jhu.edu/software/kraken/) and many other tools.)
`sourmash lca` can be used to classify individual genome bins with
`classify`, or summarize metagenome taxonomy with `summarize`.  The
[sourmash lca tutorial](http://sourmash.readthedocs.io/en/latest/tutorials-lca.html)
shows how to use the `lca classify` and `summarize` commands, and also
provides guidance on building your own database.

The other approach, `gather`, breaks a metagenome down into individual
genomes based on greedy partitioning. Essentially, it takes a query
metagenome and searches the database for the most highly contained
genome; it then subtracts that match from the metagenome, and repeats.
At the end it reports how much of the metagenome remains unknown.  The
[basic sourmash
tutorial](http://sourmash.readthedocs.io/en/latest/tutorials.html#what-s-in-my-metagenome)
has some sample output from using gather with GenBank.

Our preliminary benchmarking suggests that `gather` is the most accurate
method available for doing strain-level resolution of genomes. More on that
as we move forward!

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

If you compute your input signatures with `--track-abundance`, both
`sourmash gather` and `sourmash lca gather` will use that information
to calculate an abundance-weighted result.  Briefly, this will weight
each match to a hash value by the multiplicity of the hash value in
the query signature.  You can turn off this behavior with
`--ignore-abundance`.

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
