# Support, Versioning, and Migration

```{contents}
   :depth: 3
```

## Asking questions and filing bugs

We do our best to support sourmash users! Users have found important
bugs, and some of our best features have come from user
requests. Please help us improve sourmash for everyone by asking
questions as you have them!

Please ask questions and file bug descriptions [on the GitHub issue tracker for sourmash, sourmash-bio/sourmash/issues][0].

You can also ask questions of Titus on Twitter at [@ctitusbrown][1].

[0]:https://github.com/sourmash-bio/sourmash/issues
[1]:https://twitter.com/ctitusbrown/

## Versioning and stability of features and APIs

We do our best to guarantee stability of features and APIs within
major versions - because of this, upgrading from (e.g.) sourmash v3.4 to
sourmash v3.5 should be a simple matter of installing the new version.

We also recommend using _version pinning_ for software and workflows
that depend on sourmash, e.g. specifying `sourmash >=3,<4` for
software that is tested with sourmash 3.x. Read on for details!

Upgrading major versions (to sourmash 4.0, for example) will often involve
more work; see the [next section](#upgrading-versions) for more
our suggested process.

### Semantic versioning

Our goal is to support the use of sourmash in pipelines and
applications by communicating clearly about bug fixes, feature
additions, and feature changes in sourmash.  Versions are tagged in a
`vMAJOR.MINOR.PATCH` format, following the [Semantic Versioning]
convention.  From their definition:

"Given a version number MAJOR.MINOR.PATCH, increment the:

* MAJOR version when you make incompatible API changes,
* MINOR version when you add functionality in a backwards compatible manner, and
* PATCH version when you make backwards compatible bug fixes."

So, for example,

* Major releases, like v4.0.0, may break backwards compatibility at
  the command line as well as top-level Python/Rust APIs.
* Minor releases, like v4.1.0, will remain backwards compatible but
  may introduce significant new features.
* Patch releases, like v4.1.1, are for minor bug fixes; full backwards
  compatibility is retained.

We do sometimes (rarely!) alter behavior in minor versions by fixing
bugs; this will be documented in release notes.

### Version pinning

For software and workflows that depend on sourmash, we recommend
pinning versions to the current _major_ release of sourmash.

For example, with Python toolchains such as pip, you should be able to use:

```
sourmash>=3,<4
```
to pin the version requirement to any sourmash v3.x release.

For conda, the same syntax should work.

### Command line stability

We intend that all command-line commands, command-line options, input
formats, and output formats will be fully backwards compatible within
major versions. That is, you should never see old behavior change when
you upgrade within a major sourmash release (barring bug fixes!). Moreover,
if you rely on a feature introduced in v3.3.0, that feature will not break
in v3.4.0, but will also not be backported to version 3.2.0.

### Python API

We intend to guarantee the Python API at the top level, i.e.
functions and classes available from the `sourmash` top-level module
will be stable within major versions.

The latest minor release (e.g. v3.5) before a new major release (v4.0)
will contain deprecations for all top-level API changes at the time of
the first major release.  See below for our suggested migration
procedure.

### Python version support

sourmash v3.x supports Python 2.7 as well as Python 3.x, through Python 3.8.

sourmash v4.0 dropped support for versions of Python before Python 3.7,
and our intent is that it will support as-yet unreleased versions of Python 3.x
(e.g. 3.10) moving forward.

For future versions of sourmash, we plan to follow the
[Numpy NEP 29](https://numpy.org/neps/nep-0029-deprecation_policy.html)
proposal for Python version support. For example, this
would mean that we would drop support for Python 3.7 on December 26,
2021.

### Rust API

The Rust API is not yet at 1.0 and should not be regarded as stable.

## Upgrading major versions

If you depend on sourmash, we recommend using the following process:

* pin sourmash to the major version you developed against, e.g. `sourmash >=3,<4`.
* when ready to upgrade sourmash, upgrade to the latest minor release within that major version (e.g. sourmash 3.5.x).
* scan for deprecations that affect you, check [the release notes](https://github.com/sourmash-bio/sourmash/releases), 
and fix any major issues noted.
* upgrade to the next major version (e.g. sourmash 4.0) and run your integration tests or workflow.
* fix outstanding issues.

In particular, we recommend upgrading major versions of sourmash in
isolation, without adding any new features to your software.

### Migrating from sourmash v3.x to sourmash v4.x.

If you want to upgrade workflows and scripts from prior releases of
sourmash to sourmash v4.0, we suggest doing this in two stages.

First, upgrade to the latest version of sourmash 3.5.x (currently
[v3.5.1](https://github.com/sourmash-bio/sourmash/releases/tag/v3.5.1)),
which is compatible with all files and command lines used in previous
versions of sourmash (v2.x and v3.x). After upgrading to 3.5.x, scan
the sourmash output for deprecation warnings and fix those.

Next, upgrade to the latest version of 4.x, which will introduce some
backwards incompatibilities based upon the deprecation warnings.

The major changes are detailed below; please see the
[full release notes for 4.0](release-notes/sourmash-4.0.md) for all
the details and links to the code changes.

### Sourmash command line

If you use sourmash from the command line, there are a few major changes in 4.0 that you should know about.

First, **`sourmash compute` is deprecated in favor of [`sourmash sketch`](sourmash-sketch.md)**, which provides quite a bit more flexibility in creating signatures.

Second, **`sourmash index` will now save databases in the Zip format (`.sbt.zip`) instead of the old JSON+subdirectory format** (see [updated docs](command-line.md#sourmash-index-build-an-sbt-index-of-signatures)). You can revert to the old behavior by explicitly specifying the `.sbt.json` filename for output when running `sourmash index`.

Third, all sourmash commands that operate on signatures should now be able to directly read from lists of signatures in signature files, SBT databases, LCA databases, directories, and files containing lists of filenames (see [updated docs](command-line.md#advanced-command-line-usage)).

Fourth, if you use `sourmash lca` commands, **`sourmash lca gather` has been removed**. In addition, there are some **changes in how `summarize` works**: it now uses abundances by default, and no longer combines all signatures before summarizing. Specify `--ignore-abundance` and combine your signatures using `sourmash sig merge` to recover the old behavior. Note also that `lca summarize` now includes a new column, `filename`, in the CSV output.

Finally, **k-mer sizes have changed for amino acid sequences** in v4. If you use protein, Dayhoff, or HP signatures, we now interpret k-mer sizes differently on the command line. Briefly, k-mer sizes for protein/dayhoff/hp signatures are now the size of the k-mer in amino acid space, *not* the space of the k-mer in DNA space (as previously used). In practice this means that you need to divide all your old k-mer sizes by 3 when working with k-mers in amino acid space!

Note also that while `sourmash compute` still behaves the same way in v4.x as it did in sourmash 3.5.x, `sourmash sketch translate` and `sourmash sketch protein` both use the *new* approach to amino acid k-mer sizes, as do all of the the command line options for searching, manipulation, and display. Again, in practice this means that you need to divide all your old k-mer sizes by 3 if they apply to amino acid k-mers.

There are several minor changes where error messages should occur appropriately:
* `--traverse-directory` is no longer needed on the command line for `sourmash index` or other functions; directory traversal happens automatically.
* the command lines for `sourmash index` and `sourmash lca index` no longer require signature files to be specified, which can break existing command lines. To fix this, reorder arguments so that any signatures are specified at the end of the command line.

### Python API

First, all k-mer sizes for `protein`, `dayhoff`, and `hp` signatures have changed in the Python layer to be "correct", i.e., to be the size of the protein k-mer. Previously they were 3\*k, i.e. based on the size of the DNA k-mer from which the protein sequence would have been created.

Second, the `MinHash` class API has changed significantly!
* `get_mins()` has been deprecated in favor of `.hashes`, which is a dictionary that contains abundances.
* `merge` now just modifies `MinHash` objects in-place, and no longer returns the merged object; use `__iadd__` (`+=`) for the old behavior, or `__add__` (`+`) to create a new merged object.
* `max_hash` has been deprecated in favor of `scaled`.
* instead of `downsample_scaled(s)` use `downsample(scaled=s)`
* instead of `downsample_n(m)` use `downsample(num=m)`
* `is_molecule_type` has been replaced with a property, `moltype` -- instead of `is_molecule_type(t)` use `moltype == t`.


Third, `SourmashSignature` objects no longer have a `name()` method but instead a `name` property, which can be assigned to. This property is now `None` when no name has been assigned. Note that `str(sig)` should now be used to retrieve a display name, and should replace all previous uses of `sig.name()`.

Fourth, a few top-level functions have been deprecated: `load_signatures(...)`, `load_one_signature(...)`, `create_sbt_index(...)`, and `load_sbt_index(...)`.
* `load_signatures(...)`, `load_one_signature(...)` should be replaced with `load_file_as_signatures(...)`. Note there is currently no top-level way to load signatures from strings. For now, if you need that functionality, you can use `sourmash.signature.load_signatures(...)` and `sourmash.signature.load_one_signature(...)`, but please be aware that these are not considered part of the public API that is under semantic versioning, so they may change in the next minor point release; this is tracked in  https://github.com/sourmash-bio/sourmash/issues/1312.
* `load_sbt_index(...)` have been deprecated.  Please use `load_file_as_index(...)` instead.
* `create_sbt_index(...)` has been deprecated. There is currently no replacement, although you can use it directly from `sourmash.sbtmh` if necessary.

Fifth, directory traversal now happens by default when loading signatures, so remove `traverse=True` arguments to several functions in `sourmash_args` - `load_dbs_and_sigs`, `load_file_as_index`, `and load_file_as_signatures`.

Please post questions and concerns to the
[sourmash issue tracker](https://github.com/sourmash-bio/sourmash/issues)
and we'll be happy to help!
