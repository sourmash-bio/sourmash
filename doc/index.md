# Welcome to sourmash!

sourmash is a command-line tool and Python/Rust library for
**metagenome analysis** and **genome comparison** using k-mers.  It
supports the compositional analysis of metagenomes, rapid search of
large sequence databases, and flexible taxonomic profiling with both
NCBI and GTDB taxonomies
([see our prepared databases for more information](databases.md)). sourmash
works well with sequences 30kb or larger, including bacterial and
viral genomes.

You might try sourmash if you want to -

* identify which reference genomes to use for metagenomic read mapping;
* search all Genbank microbial genomes with a sequence query;
* cluster hundreds or thousands of genomes by similarity;
* taxonomically classify genomes or metagenomes against NCBI and/or GTDB;
* search thousands of metagenomes with a query genome or sequence;

<span style="color:red"><b>New!</b></span> The sourmash project also supports
[querying all 1 million publicly available metagenomes in the Sequence Read
Archive](https://branchwater.sourmash.bio/). Give it a try!

Our **vision**: sourmash strives to support biologists in analyzing
modern sequencing data at high resolution and with full context,
including all public reference genomes and metagenomes.

## Mission statement

This project's mission is to provide practical tools and approaches for
analyzing extremely large sequencing data sets, with an emphasis on
high resolution results. Our designs follow these guiding principles:

* genomic and metagenomic analyses should be able to make use of all
  available reference genomes.
* metagenomic analyses should support assembly independent approaches,
  to avoid biases stemming from low coverage or high strain
  variability.
* private and public databases should be equally well supported.
* a variety of data structures and algorithms are necessary to support
  a wide set of use cases, including efficient command-line analysis,
  real-time queries, and massive-scale batch analyses.
* our tools should be well behaved members of the bioinformatics
  analysis tool ecosystem, and use common installation approaches,
  standard formats, and semantic versioning.
* our tools should be robustly tested, well documented, and supported.
* we discuss scientific and computational tradeoffs and make specific
  recommendations where possible, admitting uncertainty as needed.

## How does sourmash work?

Underneath, sourmash uses [FracMinHash sketches](https://www.biorxiv.org/content/10.1101/2022.01.11.475838) for fast and
lightweight sequence comparison; FracMinHash builds on
[MinHash sketching](https://en.wikipedia.org/wiki/MinHash) to support both Jaccard similarity
_and_ containment analyses with k-mers.  This significantly expands
the range of operations that can be done quickly and in low
memory. sourmash also implements a number of new and powerful techniques
for analysis, including minimum metagenome covers and alignment-free ANI
estimation.

sourmash is inspired by [mash](https://mash.readthedocs.io), and
supports most mash analyses. sourmash also implements an expanded set
of functionality for metagenome and taxonomic analysis.

While sourmash is currently single-threaded, the [branchwater plugin for sourmash](https://github.com/sourmash-bio/sourmash_plugin_branchwater) provides faster and lower-memory multithreaded implementations of several important sourmash features - sketching, searching, and gather (metagenome decomposition). It does so by implementing higher-level functions in Rust on top of the core Rust library of sourmash. As a result it provides some of the same functionality as sourmash, but 10-100x faster and in 10x lower memory. Note that this code is functional and tested, but does not have all of the features of sourmash. Code and features will be integrated back into sourmash as they mature.

sourmash development was initiated with a grant from the Moore
Foundation under the Data Driven Discovery program, and has been
supported by further funding from the NIH and NSF. Please see
[funding acknowledgements](funding.md) for details!

## Using sourmash

<!-- Please remember to add things to sidebar.md ;) -->

### Tutorials and examples

These tutorials are command line tutorials that should work on Mac OS
X and Linux. They require about 5 GB of disk space and 5 GB of RAM.

* [Installing sourmash with conda](tutorial-install.md)

* [The first sourmash tutorial - making signatures, comparing, and searching](tutorial-basic.md)

* [Using sourmash LCA to do taxonomic classification](tutorials-lca.md)

* [Analyzing the genomic and taxonomic composition of an environmental genome using GTDB and sample-specific MAGs with sourmash](tutorial-lemonade.md)

* [Some sourmash command line examples!](sourmash-examples.ipynb)

### How-To Guides

* [Classifying genome and metagenome sketches](classifying-signatures.md)

* [Working with private collections of genome sketches](sourmash-collections.ipynb)

* [Using the `LCA_Database` API](using-LCA-database-API.ipynb)

* [Building plots from `sourmash compare` output](plotting-compare.ipynb).

* [A short guide to using sourmash output with R](other-languages.md).

## Frequently Asked Questions

* [Frequently asked questions](faq.md)

### How sourmash works under the hood

* [An introduction to k-mers for genome comparison and analysis](kmers-and-minhash.ipynb)
* [Support, versioning, and migration between versions](support.md)

### Reference material

* [Full table of contents for all docs](toc.md)
* [UNIX command-line documentation](command-line.md)
* [Genbank and GTDB databases and taxonomy files](databases.md)
* [Python examples using the API](api-example.md)
* [Publications about sourmash](publications.md)
* [A guide to the internal design and structure of sourmash](sourmash-internals.md)
* [Funding acknowledgements](funding.md)

## Developing and extending sourmash

* [Releasing a new version of sourmash](release.md)

<!--

This toctree sets the sidebar menu, but is otherwise hidden so that it
doesn't show up redundantly at the bottom of the index page.

-->

```{toctree}
---
maxdepth: 1
hidden: true
---

sidebar
command-line
databases
api-example
```
