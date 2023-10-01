# Welcome to sourmash!

sourmash is a command-line tool and Python/Rust library for
**metagenome analysis** and **genome comparison** with k-mers.  It
supports the compositional analysis of metagenomes, rapid search of
large sequence databases, and flexible taxonomic analysis with both
NCBI and GTDB taxonomies. sourmash works well with sequences 30kb or
larger, including bacterial and viral genomes.

You might try sourmash if you want to -

* identify which reference genomes to map your metagenomic reads to
* search all Genbank microbial genomes with a sequence query
* cluster many genomes by similarity
* taxonomically classify genomes or metagenomes against NCBI and/or GTDB;
* search thousands of metagenomes with a query genome or sequence

Underneath, sourmash uses [FracMinHash sketches](https://www.biorxiv.org/content/10.1101/2022.01.11.475838) for fast and
lightweight sequence comparison; FracMinHash builds on
[MinHash sketching](https://en.wikipedia.org/wiki/MinHash) to support both Jaccard similarity
_and_ containment analyses with k-mers.  This significantly expands
the range of operations that can be done quickly and in low
memory. sourmash also implements a number of new and powerful analysis
techniques, including minimum metagenome covers and alignment-free ANI
estimation.

sourmash is inspired by [mash](https://mash.readthedocs.io), and
supports most mash analyses. sourmash also implements an expanded set
of functionality for metagenome and taxonomic analysis.

sourmash development was initiated with a grant from the Moore
Foundation under the Data Driven Discovery program, and has been
supported by further funding from the NIH and NSF. Please see
[funding acknowledgements](funding.md) for details!

## Using sourmash

### Tutorials and examples

These tutorials are command line tutorials that should work on Mac OS
X and Linux. They require about 5 GB of disk space and 5 GB of RAM.

* [The first sourmash tutorial - making signatures, comparing, and searching](tutorial-basic.md)

* [Using sourmash LCA to do taxonomic classification](tutorials-lca.md)

* [Analyzing the genomic and taxonomic composition of an environmental genome using GTDB and sample-specific MAGs with sourmash](tutorial-lemonade.md)

* [Some sourmash command line examples!](sourmash-examples.md)

### How-To Guides

* Installing sourmash

* [Classifying genome sketches](classifying-signatures.md)

* [Working with private collections of genome sketches.](sourmash-collections.md)

* [Using the `LCA_Database` API.](using-LCA-database-API.ipynb)

* [Building plots from `sourmash compare` output](plotting-compare.md).

* [A short guide to using sourmash output with R](other-languages.md).

### How sourmash works under the hood

* [An introduction to k-mers for genome comparison and analysis](kmers-and-minhash.md)
* [Support, versioning, and migration between versions](support.md)

### Reference material

* [UNIX command-line documentation](command-line.md)
* [Genbank and GTDB databases and taxonomy files](databases.md)
* [Python examples using the API](api-example.md)
* [Publications about sourmash](publications.md)
* [A guide to the internals of sourmash](sourmash-internals.md)
* [Funding acknowledgements](funding.md)

## Developing and extending sourmash

* [Releasing a new version of sourmash](release.md)
