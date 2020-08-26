# Support, Versioning, and Migration

## Asking questions and filing bugs

We do our best to support sourmash users! Users have found important
bugs, and some of our best features have come from user
requests. Please help us improve sourmash for everyone by asking
questions as you have them!

Please ask questions and file bug descriptions [on the GitHub issue tracker for sourmash, dib-lab/sourmash/issues][0].

You can also ask questions of Titus on Twitter at [@ctitusbrown][1].

[0]:https://github.com/dib-lab/sourmash/issues
[1]:https://twitter.com/ctitusbrown/

## Versioning

Our goal is to support the use of sourmash in pipelines and
applications by communicating clearly about bug fixes, feature
additions, and feature changes.  Versions are tagged in a
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

### Rust API

The Rust API is not yet at 1.0 and should not be regarded as stable.

### How to "pin" sourmash versions

If you are relying on sourmash in a pipeline or application, we
suggest specifying your version requirements at the major release,
e.g. in conda you would specify `sourmash>=3,<4` to rely on sourmash
v3.x features.

Release notes for minor and patch versions are available on the
[GitHub releases page](https://github.com/dib-lab/sourmash/releases).

[Semantic Versioning]: https://semver.org/

### Python version support

sourmash v3.x supports Python 2.7 as well as Python 3.x, through Python 3.8.

sourmash v4.0 dropped support for version of Python before Python 3.7,
and our intent is that it will support as-yet unreleased versions of Python 3.x
(e.g. 3.9) moving forward.

For future versions of sourmash, we plan to follow the
[Numpy NEP 29](https://numpy.org/neps/nep-0029-deprecation_policy.html)
proposal for Python version support in the future. For example, this
would mean that we would drop support for Python 3.7 on December 26,
2021.

## Migrating from sourmash v3.x to sourmash 4.x.

Prior to the release of sourmash v4, we are adding deprecation
warnings and/or future warnings to all APIs and modules in sourmash
v3.x that are being removed in v4.0. If you are using the Python API,
we suggest you use the following procedure to migrate:

* first, install the latest version of sourmash v3, which should be v3.5.0 or later.
* then, turn on `DeprecationWarning`s in your code per [the warnings module documentation](https://docs.python.org/3/library/warnings.html#overriding-the-default-filter).
* now, run python with the argument `-W error` to turn warnings into errors.
* fix all errors!
* finally, upgrade to sourmash v4.0.
