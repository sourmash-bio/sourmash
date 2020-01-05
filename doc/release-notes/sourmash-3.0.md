# sourmash v3.0 release notes

We are pleased to announce release 3.0 of sourmash!  This release includes internal changes like the new Rust backend (replacing C++) and lays the groundwork for future improvements.

The full documentation for this version is available at
https://sourmash.readthedocs.io/en/v3.0.0/. Documentation for the latest
released version is at http://sourmash.readthedocs.io/en/stable/.

**This release is compatible with sourmash 2.0: the sourmash
Python API, command-line, signature and databases file formats are all the same.** We are releasing 3.0 to indicate the build system and internal implementation changed.

Please post questions about migrating to sourmash 3.0
in the [sourmash issue tracker](https://github.com/dib-lab/sourmash/issues/new).

## Highlighted changes since 2.0

This is a list of substantial new features and functionality since sourmash 2.0. For more details check the [releases page on GitHub](https://github.com/dib-lab/sourmash/releases).

Features:

- Replacing C++ with Rust ([#424](https://github.com/dib-lab/sourmash/pull/424))
- Create an `Index` abstract base class ([#556](https://github.com/dib-lab/sourmash/pull/556))
- Dayhoff and HP encoding for proteins ([#689](https://github.com/dib-lab/sourmash/pull/689)) ([#758](https://github.com/dib-lab/sourmash/pull/758))
- Add `sourmash signature filter` to do abundance filtering. ([#748](https://github.com/dib-lab/sourmash/pull/748))
- Parallelized compare function with multiprocessing ([#709](https://github.com/dib-lab/sourmash/pull/709))
- add compute signatures for 10x bam file ([#713](https://github.com/dib-lab/sourmash/pull/713))

Improvements:

- improve error handling in `sourmash lca index`. ([#798](https://github.com/dib-lab/sourmash/pull/798))
- Include more base deps: numpy, scipy and matplotlib ([#770](https://github.com/dib-lab/sourmash/pull/770))
- use `bam2fasta` package to simplify `sourmash compute` ([#768](https://github.com/dib-lab/sourmash/pull/768))
- add a `--abundances-from` flag to sourmash signature intersect, to preserve abundances ([#747](https://github.com/dib-lab/sourmash/pull/747))
- Compare outputs can be saved to an output dir ([#715](https://github.com/dib-lab/sourmash/pull/715))
