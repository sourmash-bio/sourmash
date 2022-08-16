# A guide to the internals of sourmash

```{contents} Contents
:depth: 3
```

sourmash was created in 2015, and has been repeatedly reorganized,
refactored, and optimized to support ever larger databases, faster
queries, and new use cases. We've also regularly added new
functionality and features.  So sourmash can be pretty complicated
internally, and our user-facing documentation only covers a fraction
of its potential!

This document is a brain dump intended for expert users and sourmash
developers who want to understand how, why, and when to use various
sourmash features. It is unlikely ever to be comprehensive, so the
information you are interested in may not yet exist in this document,
but we are always happy to add to it -
[just ask in an issue!](https://github.com/sourmash-bio/sourmash/issues)

## Signatures and sketches

sourmash operates on sketches. Each sketch is a collection of hashes.
Each sketch is contained in a signature.

Internally, sketches (class `MinHash`) contain the following information:
* a set of hashes;
* an optional abundance for each hash (when `track_abund` is True);
* a seed;
* a k-mer size;
* a molecule type;
* either a `num` (for MinHash) or a `scaled` value (for FracMinHash);

Signature objects (class `SourmashSignature`) contain a sketch (property `.minhash`) as well as additional information:
* an optional `name`
* an optional `filename`
* a license (currently must be CC0);
* an `md5sum(...)` method that returns a hash of the sketch.

For now, we tend to refer to signatures and sketches interchangeably,
because they are almost entirely 1:1 in the code base (but see [sourmash#616](https://github.com/sourmash-bio/sourmash/issues/616)).

The default signature interchange/serialization format is JSON, optionally
gzipped. This format is written and read by Rust code.

In general, a lot of effort in sourmash is spent managing collections of
signatures _before_ actually doing comparisons with them; see manifests,
and `Index` objects, below.

### Making sketches

Sketches are produced by hashing k-mers with murmurhash and then
keeping either the lowest `num` hashes (for MinHashes sketches) or
keeping all hashes below `2**64 / scaled` (for FracMinHash sketches).

The default MinHash sketches use parameters so that they are
compatible with mash sketches.

See [utils/compute-dna-mh-another-way.py](https://github.com/sourmash-bio/sourmash/blob/latest/utils/compute-dna-mh-another-way.py) for details on how k-mers are
hashed.

Note that if hashes are produced some other way (with a different hash
function) or from some source other than DNA, sourmash can still work
with them; only `sourmash sketch` does input checking.

### Compatibility checking

The point of the signatures and sketches is to enable certain kinds of
rapid comparisons - Jaccard similarity and overlap analysis,
specifically. However, these comparisons can only be done between
compatible sketches.

Here, "compatible" means -
* the same MurmurHash seed (default 42);
* the same k-mer size/ksize (see k-mer sizes, below);
* the same molecule type (see molecule types, below);
* the same `num` or `scaled` (although see [this downsampling discussion](api-example.md#downsampling-signatures), and the next two sections);

sourmash uses selectors (`Index.select(...)`) to select sketches with
compatible ksizes, molecule types, and sketch types.

### Scaled (FracMinHash) sketches support similarity and containment

Per our discussion in [Irber et al., 2022](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2), FracMinHash sketches can always be compared
by downsampling to the max of the two scaled values.  (This is not always
true of sketches stored in indexed collections, e.g. SBTs; see [sourmash#1799](https://github.com/sourmash-bio/sourmash/issues/1799).)

In practice, sourmash does all necessary downsampling dynamically, but
returns the original sketches. This means that (for example) you can
do a low-resolution/high-scaled large scale search by specifying a
high scaled value, and then do a higher resolution comparison with
only the highly similar matches the results to do a more refined (see
below, Speeding up `gather` and `search`.)

### Num (MinHash) sketches support Jaccard similarity

"Regular" MinHash (or "num MinHash") sketches are implemented the same
way as in mash.  However, they are less well supported in sourmash,
because they don't offer the same opportunities for metagenome
analysis.  (See also [sourmash#1354](https://github.com/sourmash-bio/sourmash/issues/1354).)

Num MinHash sketches can always be compared by downsampling to a
common `num` value. This may need to be done manually using `sourmash
sig downsample`, however.

## K-mer sizes

There is no explicit restriction on k-mer sizes built into sourmash.

For highly specific genome and metagenome comparisons, we typically
use k=21, k=31, or k=51. For a longer discussion, see [Assembly Free Analysis with K-mers](https://github.com/mblstamps/stamps2022/blob/main/kmers_and_sourmash/2022-stamps-assembly-free%20class.pdf) from STAMPS 2022.

## Molecule types - DNA, protein, Dayhoff, and hydrophobic-polar

sourmash supports four different sequence encodings, which we refer to
as "molecule: DNA (`--dna`), protein (`--protein`), Dayhoff,
(`--dayhoff`), and hydrophobic-polar (`--hp`).

All FracMinHash sketches have exactly one molecule type, and can only
be compared to the same molecule type (and ksize).

DNA moltype sketches can be constructed from DNA input sequences using
`sourmash sketch dna`.

Protein, Dayhoff, and HP moltype sketches can be constructed from
protein input sequences using `sourmash sketch protein`, or from DNA
input sequences using `sourmash sketch translate`; `translate` will
translate in all six reading frames (see also
[orpheum](https://github.com/czbiohub/orpheum) from
[Botvinnik et al., 2021](https://www.biorxiv.org/content/10.1101/2021.07.09.450799v1)).
By default protein sketches will be created; dayhoff sketches can be
created by including `dayhoff` in the param string, e.g. `sourmash
sketch protein -p dayhoff`, and hydrophobic-polar sketches can be
built with `hp` in the param string.

## Manifests

sourmash makes extensive use of signature manifests to support rapid
selection and lazy loading of signatures based on signature metadata
(name, ksize, moltype, etc.) See
[Blog post: Scaling sourmash to millions of samples](http://ivory.idyll.org/blog/2021-sourmash-scaling-to-millions.html)
for some of the motivation.

Manifests are an internal format that is not meant to be particularly
human readable, but the CSV format can be loaded into a spreadsheet
program if you're curious :).

If in a zipfile (`.zip`) or SBT zipfile (`.sbt.zip`), manifests must
be named `SOURMASH-MANIFEST.csv`. They can also be stored directly on
disk in CSV/gzipped CSV, or in a sqlite database; see
`sourmash sig manifest`, `sourmash sig check`, and `sourmash sig collect`
for manifest creation, management, and export utilities.

Where signatures are stored individually in `Index` collections,
e.g. in a zipfile, manifests may be stored alongside them; for other
subclasses of `Index` such as the inverted indices, manifests are
generated dynamically by the class itself.

Currently (sourmash 4.x) manifests do not contain information about the
hash seed or sketch license. This will be fixed in the future - see [sourmash#1849](https://github.com/sourmash-bio/sourmash/issues/1849).

Manifests are very flexible and, especially when stored in a sqlite
database, can be extremely performant for organizing hundreds of
thousands to millions of sketches.  Please see `StandaloneManifestIndex`
for a lazy-loading `Index` class that supports such massive-scale
organization.

## Index implementations

The `Index` class and its various subclasses (in `sourmash.index`) are
containers that provide an API for organizing, selecting, and
searching (potentially) large numbers of signatures.

`sourmash sig summarize` is a good way to determine what type of `Index`
class is used to handle a collection.

Loading and saving of `Index` objects is handled separately from the
class: loading can be done in Python via the
`sourmash.load_file_as_index(...)` method, while creation and/or
updating of `Index` objects is done via
`sourmash.sourmash_args.SaveSignaturesToLocation(...)`.  These are the
same APIs used by the command-line functionality.

There are quite a few different `Index` subclasses and they all have
distinct features.  We have a high-level guide to which collection
type to use [here](command-line.md#loading-many-signatures).

Conceptually, `Index` classes are either organized around storing
individual signatures often with metadata that permits loading,
selecting, and/or searching them more efficiently
(e.g. `ZipFileLinearIndex` and `SBTs`); or they store signatures
as inverted indices (`LCA_Database` and `SqliteIndex`) that permit
certain kinds of fast queries.

Unless otherwise noted, the `Index` classes below can be loaded
concurrently in "read only" mode - that is, you should build the
collection _once_, and then use it from multiple processes. We
currently do not test for or support concurrent read/write. Note also
that (generally speaking) memory footprints will be additive, so
loading the same `LCA_Database` twice will consume twice the memory.
(If you're interested in concurrency, we suggest using the sqlite
containers - see `SqliteIndex`.)

### In-memory storage and search.

The simplest way to handle collections of signatures is to load them
into memory, but it is also the least performant and most memory
intensive mechanism!

`LinearIndex` and `MultiIndex` both support sketches loaded from
JSON files; both will load the sketches once and then keep them in
memory.  `LinearIndex` does not use manifests while `MultiIndex` builds
a manifest as it loads the sketches.

Note that `MultiIndex` is the class used to load signatures from
pathlists, directory hierarchies, and so on; because it stores
sketches in memory, this can incur a significant memory penalty (see
[sourmash#1899](https://github.com/sourmash-bio/sourmash/issues/1899)).  Therefore where possible we suggest building a standalone
manifest (`StandaloneManifestIndex`) to do lazy loading from the disk
instead; you can use `sourmash sig collect` to do this.

### Zipfile collections

`ZipFileLinearIndex` stores signature files in a zip file with an
accompanying manifest.  This is the most versatile and compressed
option for working with large collections of sketches - it supports
rapid selection and loading of specific sketches from disk, and can
store and search any mixture of sketch types (ksize, molecule type,
scaled values, etc.)

By default, `ZipFileLinearIndex` stores one signature (equiv. one
sketch) in each member file in the zip archive. Each signature is
stored uncompressed. The accompanying manifest stores the full member
file path in `internal_location`, so that sketches can be retrieved
directly.

Searching a `ZipFileLinearIndex` is done linearly, as implied by the
name. This is fine for `gather` but if you are doing repeated queries
with `search` you may want to use an SBT or LCA database instead; see
below.

In the future we expect to parallelize searching `ZipFileLinearIndex`
files in Rust; see [sourmash#1752](https://github.com/sourmash-bio/sourmash/issues/1752).

`ZipFileLinearIndex` does support zip files without manifests as well
as multiple signatures in a single file; this was originally intended
to support simply zipping entire directory hierarchies into a zipfile.
However, this slows down performance and is not recommended.  If you
have an existing zipfile (or really any collection of signatures) and
you want to turn them into a proper `ZipFileLinearIndex`, you can use
`sig cat <collection(s)> -o combined.zip` to create a
`ZipFileLinearIndex` file named `combined.zip` that will have a
manifest and signatures broken out into individual files.

### Sequence Bloom Trees (SBTs)

Sequence Bloom Trees (SBTs; see
[the Kingsford Lab page for details](http://www.cs.cmu.edu/~ckingsf/software/bloomtree/))
provide a faster (but more memory intensive) on-disk storage and
search mechanism.  In brief, SBTs implement a binary tree organization
of an arbitrary number of signatures; each internal node is a Bloom
filter containing all of the hashes for the nodes below them. This
permits potentially rapid elimination of irrelevant nodes on search.

SBTs are restricted to storing and searching sketches with the
same/single k-mer size and molecule type, as well as either a single
num value or a single scaled value.

We suggest using SBTs when you are doing multiple Jaccard search or
containment searches with genomes via `sourmash search`.

### Lowest common ancestor (LCA) databases

The `LCA_Database` index class stores signatures in an inverted index,
where a Python dictionary is used to link individual hashes back to
signatures and/or taxonomic lineages. This supports the individualized
hash analyses used in the `lca` submodule.

LCA databases only support a single ksize, moltype, and scaled. They
can only be used with FracMinHash (scaled) sketches.

The default `LCA_Database` class is serialized via JSON, and loads
everything into memory when requested. The load time incurs a
significant latency penalty when used from the command line, as well
as having a potentially large memory footprint; this makes it
difficult to use the default `LCA_Database` for very large databases,
e.g. genbank bacteria.

The newer `LCA_SqliteDatabase` (based on `SqliteIndex`, described
below) also supports LCA-style queries, and is stored on disk, is fast to
load, and uses very little memory. The tradeoff is that the underlying
sqlite database can be quite large.  `LCA_SqliteDatabase` should also
support rapid concurrent access (see [sourmash#909](https://github.com/sourmash-bio/sourmash/issues/909)).

Both types of LCA database can be constructed with `sourmash lca index`.

### SqliteIndex

The `SqliteIndex` storage class uses sqlite3 to store hashes and
sketch information for search and retrieval; see
[this blog post](http://ivory.idyll.org/blog/2022-storing-ulong-in-sqlite-sourmash.html)
for background information and details. These are fast, low-memory,
on-disk databases, with the tradeoff that they can be quite large.
This is probably currently the best solution for concurrent access to
sketches via e.g. a Web server (see also [sourmash#909](https://github.com/sourmash-bio/sourmash/issues/909)).

`SqliteIndex` can only contain FracMinHash sketches and can only store
sketches with the same scaled parameter. However, it can store 
multiple ksizes and moltypes as long as the same scaled is used.

`SqliteIndex` objects can be constructed using `sourmash sig cat
... -o filename.sqldb`.

### Standalone manifests

The `StandaloneManifestIndex` class loads standalone manifests generated
by `sourmash sig collect`. They support rapid selection and lazy loading
on potentially extremely large collections of signatures.

The underlying mechanism uses the `internal_location` field of
manifests to point to the container file. When particular sketches are
requested, the container file is loaded into an `Index` object with
`sourmash.load_file_as_index` and the `md5` values of the requested
sketches are used as a picklist to retrieve the desired signatures.

Thus, while standalone manifests can point at any kind of container,
including JSON files or LCA databases, they are most efficient when
`internal_location` points at a file with either a single sketch in
it, or a manifest that supports direct loading of sketches. Therefore,
we suggest using standalone manifest indices.

Note that searching a standalone manifest is currently done through a
linear iteration, and does not use any features of indexed containers
such as SBTs or LCAs.  This is fine for `gather` with the default
approach, but is probably suboptimal for a `search`.

### Pathlists and `--from-file`

All (or most) sourmash commands natively support taking in lists of
signature collections via pathlists, `--from-file`, or paths to
directories. This is useful for situations where you have thousands of
signature files and don't want to provide them explicitly on the
command line; you can simply put a list of the files in a text file,
and pass it in directly (or use `--from-file` to pass it in).

Both pathlists and files passed to `--from-file` contain a list of
paths to be loaded; relatives paths will be interpreted relative to
the current working directory of sourmash.  Pathlists should be
universally available on sourmash commands.  When `--from-file` is
available for a command, sourmash will behave as if the file paths in
the file were provided on the command line.

We suggest avoiding pathlists. Instead, we suggest using `--from-file`
or a standalone manifest index (generated with `sourmash sig
collect`). This is because the signatures from pathlists are loaded
into memory (see `MultiIndex`, above) it is generally a bad idea to
use them - they may be slow to load and may consume a lot of
memory. They also do not support good loading error messages;
see [sourmash#1414](https://github.com/sourmash-bio/sourmash/issues/1414).

### Extensions for outputting index classes

Most commands that support saving signatures will save them in a
variety of formats, based on the extension provided (see
[sourmash#1890](https://github.com/sourmash-bio/sourmash/issues/1890)
for exceptions). The supported extensions are -

* `.zip` for `ZipFileLinearIndex`
* `.sqldb` for `SqliteIndex`
* `.sig` or `.sig.gz` for JSON/gzipped JSON
* `dirname/` to save in a directory hierarchy

The default signature save format is JSON, if the extension is not
recognized.

## Speeding up `gather` and `search`

There are two primary search commands in sourmash: `gather` and
`search`.

`gather` calculates a minimum metagenome cover as discussed in [Irber et al., 2022](https://www.biorxiv.org/content/10.1101/2022.01.11.475838v2). It
is mostly intended for querying a database with a metagenome, although
it can be used with genome queries, as well. It depends on overlap
analyses and can only be used with FracMinHash sketches.

`search` does a straight Jaccard similarity search on MinHash and
FracMinHash sketches (or, with `--containment`, a containment search
on FracMinHash sketches). It is typically used to find matches to a
query genome sketch in a large database of sketches.

The `prefetch` command does a containment search and is intended for
power users; it is a standalone implementation of the prefetch
algorithm discussed below for `gather`.  It only works with FracMinHash
sketches.

Note that all of these commands work with any and all `Index`
collection/container types, and will return the same results however
the collections are organized - see the "online behavior" section,
below. In practice this means that you can provide additional
collections of signatures via the command line without building a
combined index of all your signatures. It also means that the only
reason to choose different collections/containers is for
optimization - you should select the containers that help you achieve
the desired performance characteristics for your search
(i.e. the right memory/time/disk space tradeoffs).

### Running `search` many times on the same database

`search` typically is used to search a large database of sketches for
all similarity or containment matches above a threshold. Depending on
the query and the database, certain kinds of database indices may make
search much faster, especially when only a few matches are expected.

If you are doing many searches against the same database, indexing the
database as an SBT (with `sourmash index`) or as a `SqliteIndex`/sqldb
database is likely to provide a significant speed increase, albeit
with increased memory usage (SBT) or increased disk space (sqldb).

Conversely, `ZipFileLinearIndex` and the default `LCA_Database` are likely
to be poor choices for many searches - the former only supports linear
searches, and the latter needs to be loaded from disk and deserialized each
time.

### Running `gather` once

`gather` is typically used to search a metagenome against a large
database of sketches, as part of finding a minimum set cover. This can
be quite slow! Our current implementation (as of sourmash 4.1.0, [pull request sourmash#1370](https://github.com/sourmash-bio/sourmash/pull/1370)) does a single pass across the database to find all matches
with an overlap above the provided threshold, and then organizes
the matches for rapid min-set-cov analysis. This single pass across the
database is called a "prefetch", and it is also implemented in the
`prefetch` subcommand.

With this single pass approach, benchmarks - [sourmash#2014](https://github.com/sourmash-bio/sourmash/issues/2014) - show that a
linearly searchable database is performant enough to be used with
`gather`.  We therefore suggest using a `ZipFileLinearIndex` container
with gather, or in cases where low-memory concurrency is desired, a
`SqliteIndex` container.

### Using `prefetch` and `gather` together

If you want to use `prefetch` independently of `gather`, you can use
the prefetch output as a picklist passed into gather - see
[picklists](#picklists), below.  This can be useful when you want to
experiment with different threshold parameters for `gather` - first,
do a very sensitive/low-threshold search with `prefetch` and save the
results to a CSV file with `-o`,

Repeated gathers and searches.

Using prefetch explicitly.

### Using a higher scaled value

With FracMinHash sketches, you can downsample the query to make both
`search` and `gather` _much_ faster.  A good rule of thumb is to use a
scaled value that is about 5x smaller than the minimum overlap to
detect; so, if you want to be able to detect 50kb of similarity, you
can use a scaled value of 10,000. Conversely, the default scaled value
of 1,000 (for DNA sketches) should robustly detect overlaps of 5kb.

You can supply `--scaled` to `gather` and `prefetch` to dynamically
downsample the query FracMinHash. For `search` you will need to use
`sourmash sig downsample` to generate a downsampled sketch.

### Running `gather` many times - `multigather`

In situations where loading the search database is slow (e.g.
`LCA_Database` or zipfiles with very large manifests), the `sourmash
multigather` command supports many queries against many databases.

(We don't particularly suggest using `multigather`; we would prefer
to make search databases faster. But it's there! :)

## Taxonomy and assigning lineages

All sourmash taxonomy handling is done within the `lca` and `tax`
subcommands (CLI) and submodules (Python).

In the case of the `lca` subcommands, the taxonomic information is
incorporated into the LCA database construction (see the `lca index`
command), while the `tax` subcommands load taxonomic information
on demand from taxonomy databases (CSVs or databases).

sourmash anchors all taxonomy to identifiers, and uses the signature
name to do so - this is the name as set by the `--name` parameter to
`sourmash sketch`, and output by `sourmash sig describe` as the
`signature:` field.

### Identifier handling

sourmash prefers identifiers to be the first space-separated token in
the signature name.  This token can contain any alphanumeric letters
other than space (@@ctb check me), and should contain at most one
period.  The version of the identifier will be the component after
the period.

So, for example, for a signature name of

```
CP001941.1 Aciduliprofundum boonei T469, complete genome
```
the identifier would be `CP001941.1` and the version would be 1.
There are no other constraints placed on the identifier, and
versions are not handled in any special way other than as below.

The `lca index` and `tax` commands both support some modified
identifier handling in sourmash 3.x and 4.x, but in the future, we
plan to deprecate these as they mostly cause confusion and internal
complexity.

The two modifiers are:

* `--keep-full-identifiers` will use the entire signature
name instead of just the first space-separated token. It is by default
off (set to False).

* `--keep-identifier-versions` turns on keeping the full identifier,
including what is after the first period. It is by default off (set to
False), stripping identifiers of their version on load. When it is on (True), identifiers are not stripped of their version on load.

### Taxonomies, or lineage spreadsheets

sourmash supports arbitrary (free) taxonomies, and new taxonomic
lineages can be created and used internally as long as they
are provided in the appropriate spreadsheet format.

You can also mix and match taxonomies as you need; for example, it is
entirely legitimate in sourmash-land to combine the GTDB taxonomy for
bacterial and archaeal sequence classification, with the NCBI taxonomy
for eukaryotic and viral sequence classification.  (You probably don't
want to mix and match within superkingdoms, though!)

As of sourmash v4, lineage spreadsheets should contain columns for
superkingdom, phylum, class, order, family, genus, and species.  Some
commands may also support a 'strain' column, although this is
inconsistently handled within sourmash internally.

For spreadsheet organization, `lca index` expects the columns to be
present in order from superkingdom on down, while the `tax`
subcommands use CSV column headers instead.  We are planning to
consolidate around the `tax` subcommand handling in the future (see [sourmash#2198](https://github.com/sourmash-bio/sourmash/issues/2198)).

@@link to example spreadsheets, talk about how to extract.

### `LCA_SqliteDatabase` - a special case

The `LCA_SqliteDatabase` index class can serve multiple purposes: as
an index of sketches (for regular search and gather); as a taxonomy
database for use with the `tax` subcommands; and as an LCA database
for use with the `lca` subcommands.

When used as a taxonomy database, an `LCA_SqliteDatabase` file
contains the same SQL tables as a sqlite taxonomy database.

When used as an LCA database, an `LCA_SqliteDatabase` dynamically loads
the taxonomic lineages from the sqlite database and applies them to the
individual hashes, permitting the same kind of hash-to-lineage query
capability as the `LCA_Database`.

## Picklists

Picklists are a generic mechanism used to select a (potentially small)
subset of signatures for search/display.

The general idea of picklists is that you create a list of signatures
you're interested in - by name, or identifier, or md5sum - and then
supply that list in a csvfile on the command line via `--picklist`.
For example, `--picklist list.csv:colname:ident` would load the
values in the column named `colname` in the file `list.csv` as identifiers
to be used to restrict the search.

The support picklist column types are `name`, `ident`
(space-delimieted identifier), `identprefix` (identifier with version
removed), `md5`, `md5prefix8`, and `md5short`.  Generally the `md5`
and derived values are used to reference signatures found some other
way with sourmash, while the identifiers are more broadly useful.

There are also four special column types that can be used without a column
name: `gather`, `prefetch`, `search`, and `manifest`. These take the
CSV output of the respective sourmash commands as inputs for picklists,
so that you can use prefetch to generate a picklist and then use that
picklist with `--picklist prefetch_out.csv.gz::prefetch`.

### Differing internal behavior 

Picklists behave differently with different `Index` classes.

For indexed databases like SBT, LCA, and `SqliteIndex`, the search is
done _first_, and then only those results that match the picklist are
selected.

For linear search databases like `ZipFileLinearIndex` or standalone
manifests, picklists are _first_ used to subselect the desired
signatures, and only those signatures are searched.

This means that picklists can dramatically speed up searches on some
`Index` types, but won't affect performance on others. But
the results will be the same.

### Taxonomy / lineage spreadsheets as picklists

Note that lineage CSV spreadsheets, as consumed by `sourmash tax` commands
and as output by `sourmash tax grep`, can be used as `ident` picklists.

## Similarity matrices with `sourmash compare`

## ANI

## Online and streaming

n+1 problem

## Formats natively understood by sourmash

sourmash should always autodetect the format of a collection or
database, in most cases based on its content (and not its
filename). Please file a bug report if this doesn't work for you!

`sourmash sig summarize` is a good way to examine the properties of a
signature collection.

### Reading and writing gzipped CSV files

(As of sourmash v4.5)

When a CSV filename is specified (e.g. `sourmash gather ... -o
mygather.csv`), you can always provide a name that ends with `.gz` to
produce a gzip-compressed file instead. This can save quite a bit of
space for prefetch results and manifests in particular!

All sourmash commands that take in a CSV (via manifest, or picklist,
or taxonomy) will autodetect a gzipped CSV based on content (the file
does not need to end with `.gz`). The one exception is manifests,
where the CSV needs to end with `.gz` to be loaded as a gzipped CSV;
see
[sourmash#2214](https://github.com/sourmash-bio/sourmash/issues/2214)
for an issue to fix this.
