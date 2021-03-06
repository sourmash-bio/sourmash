# sourmash v2.0 release notes

We are pleased to announce release 2.0 of sourmash!  This release is
the first full release of sourmash since 1.0, and provides much
expanded core functionality as well as many more utility functions.

The full documentation for this version is available at
https://sourmash.readthedocs.io/en/v2.0/. Documentation for the latest
released version is at http://sourmash.readthedocs.io/en/stable/.

**This release breaks compatibility with sourmash 1.0: the sourmash
Python API, command-line, and signature file formats have all changed
substantially.** Please post questions about migrating to sourmash 2.0
in the
[sourmash issue tracker](https://github.com/dib-lab/sourmash/issues/new).

## Major new features since 1.0

This is a list of substantial new features and functionality in sourmash 2.0.

* Added Sequence Bloom Tree search to enable similarity and containment queries on very large collections of signatures in low memory; see `sourmash index`, `sourmash search`, and `sourmash gather` in [the command line documentation](../command-line.md).
* Added "LCA databases" for fast searching of large databases in not-so-low memory; see [`sourmash lca index` in command-line docs](../command-line.md#sourmash-lca-subcommands-for-taxonomic-classification).
* Created [precomputed databases](../databases.md) for most of GenBank genomes.
* Added taxonomic reporting functionality in the `sourmash lca` submodule - [see command-line docs](../command-line.md#sourmash-lca-subcommands-for-taxonomic-classification).
* Added signature manipulation utilities in the `sourmash signature` submodule - [see command-line docs](../command-line.md#sourmash-signature-subcommands-for-signature-manipulation)
* Introduced new modulo hash or "scaled" signatures for containment analysis; see [Using sourmash: a practical guide](../using-sourmash-a-guide.md#what-resolution-should-my-signatures-be-how-should-i-create-them) and [more details in the Python API examples](../api-example.md#advanced-features-of-sourmash-minhash-objects-scaled-and-num).
* Switched to using JSON instead of YAML for signatures.
* Many performance optimizations!
* Many more tests!
* A much cleaner and more robust [Python API](../api-example.md).
* Installation via bioconda is now recommended (and actively maintained :)
* Support for building signatures from BAM files.

## Other features of note

* Renamed Python library to `sourmash` from `sourmash_lib`.
* Signatures now default to a license of CC0.
* Implemented `sourmash watch` for classifying streams of data (see [blog post](http://ivory.idyll.org/blog/2017-sourmash-sra-microbial-wgs.html)).
* Implemented `sourmash categorize` for classifying many signatures (see [blog post](http://ivory.idyll.org/blog/2017-sourmash-sra-microbial-wgs.html)).
* CSV output enabled from most commands.
