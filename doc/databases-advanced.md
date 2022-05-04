# sourmash databases - advanced usage information.

sourmash uses a variety of different mechanisms and formats for storing, organizing, and searching signatures. Some of these mechanisms, "collections", just store the signatures; others ("indexed" databases) provide indices on the signatures for fast content-based search. _Most_ of the mechanisms now use manifests that permit fast selection and loading of signatures based on metadata. Below we refer to "databases" generically as any on-disk storage mechanism for sourmash signatures.

Which database type is best to use depends on what you're doing - which is what this document is about! In general, however, sourmash should be fast enough that database choice will only impacts performance when searching 1000s of signatures, or doing many 1000s of searches.

The recommended file extensions below are conventions used to signal the output format when using `-o` with `sourmash sketch` and the `sourmash sig` subcommands; so, for example, `sourmash sketch dna *.fa -o xyz.zip` will output signatures in the .zip format.

sourmash will automatically detect and load the database, based on the database _content_ in most cases.

Unless noted otherwise, the below database formats are supported in all release since sourmash v3.5.

## How are signatures actually stored?

sourmash signatures are typically serialized into JSON for on-disk storage, with rare exceptions (SQLite and LCA databases). The internal sourmash code automatically detects and properly handles compressed (gzipped) JSON data.

## Storing JSON in `.sig` and `.sig.gz` files: the original format.

Multiple signatures can be stored in a single JSON file. However, this file will be loaded in its entirety by sourmash, even if you only select one for later analysis.

This is the least efficient way to store multiple signatures, because all of the JSON must be loaded before any signature can be selected or searched. But it is the oldest format and so a lot of our documentation describes it!

## Storing signatures in `.zip` files: the **recommended** format.

**This is our recommended format for storing collections of signatures. It is supported as of sourmash v4.1.**

Multiple signatures can be stored in a single .zip file. The best way to construct that zip file is from within sourmash, by specifying `-o filename.zip` when outputting signatures. Zip files created from within sourmash will automatically have manifests; this enables rapid subselection and direct loading of signatures via e.g. picklists.

Zip files are not indexed by content, so they can be slow for searching. But they are small, and provide a good compromise between disk size (small), flexibility (can store any mixture of signatures), and speed (good for `gather`, not good for `search`).

Zip file collections can contain any number of signatures, of any type (`num` or `scaled`, DNA/protein/dayhoff/hp).

You can create your own Zip files by simply zipping any number of `.sig` or `.sig.gz` files into a .zip file, and sourmash will read this. However, since this zip file will not have a manifest, it will not be fast for certain operations that rely on manifests for speed, such as picklists and `sourmash sig summarize`. So we recommend using sourmash to create zip file collections with manifests.

### Storing signatures in SQLite databases

As of sourmash 4.4, we support storing signatures directly in a [SQLite](https://www.sqlite.org/index.html) database (`-o .sqldb`). This is a fast, low-memory, on-disk format that is suitable for use with `search` and can support multiple simultaneous queries. However, the resulting file is also rather large, so we do not distribute databases in this format.

SQLite databases are implemented as an [inverted index](https://en.wikipedia.org/wiki/Inverted_index), with hashes stored directly in a table.

SQLite databases are limited to scaled signatures, and can only contain sketches with the same scaled value across the entire database. They *can* store multiple molecule types.

While SQLite databases are a new format, they seem promising, especially when disk space is not a concern and/or when memory is limited. We particularly recommend them for use as LCA databases (see next section) where they are a considerable improvement over the legacy JSON format.

### Other Indexed collections - SBTs and LCAs.

We provide two other indexed collection formats, Sequence Bloom Trees (SBTs) and LCA databases.

SBTs implement our version of [Sequence Bloom Trees](http://www.cs.cmu.edu/~ckingsf/software/bloomtree/), a fast tree-based index that support rapid `search` for matches; they are particularly effective when searching for *best* matches across large databases. They are relatively low memory and typically about twice the size of .zip files on disk. They can be constructed with `sourmash index`.

LCA databases are [inverted indices](https://en.wikipedia.org/wiki/Inverted_index) that support individual hash lookup. They provide fast `search` and `gather`, and also support all of the `sourmash lca` subcommands for hash-based taxonomic analysis. There are two LCA database formats, JSON and SQLite; JSON is small on disk but JSON LCA databases consume a lot of memory when loaded, while SQLite LCA databases are large on disk but low-memory and fast. JSON LCA databases do not support multiprocess queries. LCA databases can be constructed with `sourmash lca index`.

Both SBTs and LCA databases can only store homogenous collections of signature types - all signatures must have the same molecule type and scaled or num value. Furthermore, LCA databases can only store scaled signatures.

We recommend SBT and LCA databases for use only in specific situations - e.g. SBTs are great for single-genome "best match" search for SBTs, and `sourmash lca` commands require LCA databases.

### Manifests

Manifests are catalogs of signature metadata - name, molecule type, k-mer size, and other information - that can be used to select specific signatures for searching or processing. Typically when using manifests the actual signatures themselves are not loaded until they are needed.

As of sourmash 4.4 manifests can be *directly* loaded from the command line as standalone collections. This lets manifests serve as a catalog of signatures stored in many different locations.

Standalone manifests are preferable to both directory storage and pathlists (below), because they support fast selection and direct lazy loading. They are the most effective solution for managing custom collections of thousands to millions of signatures.

Manifests can be created with `sourmash sig manifest` and `sourmash sig check`. For complex situations, we recommend using custom Python scripts to manage them - for example, see [sigs-to-manifest.py in database-examples](https://github.com/sourmash-bio/database-examples/blob/main/sigs-to-manifest.py).

Sourmash supports two manifest file formats - CSV and SQLite. SQLite manifests are much faster than CSV manifests in exchange for extra disk space.

### Directories

Directory hierarchies of signatures are read natively by sourmash, and can be created or extended by specifying `-o dirname/` (with a trailing slash).

To read from a directory, specify the directory name on the sourmash command line. When reading from directories, the entire directory hierarchy is traversed and all `.sig` and `.sig.gz` files are loaded as signatures. If `--force` is specified, _all_ files will be read, and failures will be ignored.

When directories are specified as outputs, the signatures will be saved by their complete md5sum underneath the directory.

We don't particularly recommend storing signatures in directory hierarchies, since most of their use cases are now covered by other approaches.

### Pathlists

Pathlists are text files containing paths to one or more sourmash databases; any type of sourmash-readable collection can be listed.

The paths in pathlists can be relative or absolute within the file system. If they are relative, they must resolve with respect to the current working directory of the sourmash command.

We don't recommend using pathlists any more, since the original use cases are now supported with picklists, but they are still supported!

Pathlists are not output by any sourmash commands.

## Storing taxonomies

sourmash supports taxonomic information output via the `sourmash lca` and `sourmash tax` subcommands. Both sets of commands rely on the same 7 taxonomic ranks: superkingdom, phylum, class, order, family, genus, and species (with limited support for a 'strain' rank). And both sets of subcommands take lineage spreadsheets that link specific identifiers to taxonomic lineages.

Lineage spreadsheets can be provided in two on-disk formats, CSV and SQLite.

CSV is the original format, and consists of separate columns for identifier and each taxonomic rank.

SQLite taxonomy databases are typically built from CSV using `sourmash tax prepare`. They contain a single table, `sourmash_taxonomy`, with columns for `ident` and each taxonomic rank. Only the `sourmash tax` command supports SQLite taxonomy databases.

## Appendix: SQLite complexities

The SQLite implementation of signature storage, metadata manifests, and LCA databases is all bundled into a single SQLite database. Beacuse of this, sourmash must examine the database tables to decide what kind of sourmash structure the database is - the logic is roughly this:

* does the database store both sketch information and taxonomy information? It's an LCA database!
* if it has sketch information but no taxonomy information, it's just a regular index.
* if it only has manifest information, it's a manifest!
* if it only has taxonomy information, it's a taxonomy!

This is complicated by several other details -

* we can treat SQLite databases with sketch information as read-only manifests, but because the sketch information is tightly coupled to the manifest table, we cannot insert new manifest entries;
* we can treat SQLite databases with sketch information as read/write taxonomy files, since the taxonomy information is not tightly coupled to the sketches;

Last but not least, the hashes in SQLite are stored as signed 64-bit integers and must be converted to unsigned 64-bit numbers internally by sourmash; negative numbers in the SQLite table represent unsigned ints that are larger than 2**63 - 1. Please see [this blog post](http://ivory.idyll.org/blog/2022-storing-ulong-in-sqlite-sourmash.html) for more information.

The SQLite schema itself is not very complicated and can be used for lineage and manifest querying by other scripts. However, we recommend doing hash value querying/search via the Python code.
