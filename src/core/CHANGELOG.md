# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.11.0] - 2021-07-07

Added:

- Add HyperLogLog implementation (#1223)

Changed:

- Update `MinHash.set_abundances` to remove hash if 0 abund; handle negative abundances. (#1575)
- Improving `MinHash.remove_many(...)` performance (#1571)
- Improved intersection and union calculations (#1475)
- Bump MSRV to 1.42 (and other dep fixes) (#1461)
- Rework the `find` functionality for `Index` classes (#1392)
- Rationalize `SourmashSignature.name` and `str(sig)` (#1179)

Fixed:

- Fix needless borrows as suggested by clippy (#1636)
- Fix Rust 1.59 lints (#1600)
- Clean up clippy lints from 1.52 (#1505)
- Fix clippy lints introduced in 1.51 (#1407)
- CI/Rust: update and fix cbindgen config (#1473)
- pin needletail version to keep MSRV at 1.37 (#1393)
- Update proptest requirement from 0.9.6 to 1.0.0 (#1344)
- Fix clippy lints introduced in 1.50 and update nix configs (#1332)
- Update finch requirement from 0.3.0 to 0.4.1 (#1290)
- update rand for test, and activate "js" feature for getrandom (#1275)
- Fix new clippy warnings from Rust 1.49 (#1267)
- CI: small build fixes (#1252)

Removed:

- Remove 10x support in compute (#1229)

## [0.10.0] - 2020-10-08

Added:

- Add `clear` option to set_abundances(...) method (#1046)

Changed:

- Replace mx by scaled (#1139)

Fixed:

- Fix Rust panic error in signature creation (#1172)
- Update typed-builder requirement from 0.6.0 to 0.7.0 (#1121)
- update CI for latest branch name change (#1150)
- Update typed-builder requirement from 0.6.0 to 0.7.0 (#1121)

## [0.9.0] - 2020-07-13

Added:

- Cache md5sum calculation (#1058)
- Expose more of the API for wasm (signature and ComputeParameters) (#1058)
- Getters and setters for ComputeParameters (#1058)

Changed: 

- Migrate from failure to thiserror (#1058)
- Bump MSRV to 1.37 (#1058)

Fixed: 

- Use the derive feature in serde instead of serde_derive (#1058)
- Use nohash-hasher crate instead of previous NoHashHasher from finch.
- Update typed-builder to 0.6.0 (#1058)
- stricter niffler versions and add new gz feature to it (#1070)

## [0.8.0] - 2020-06-26

Added:

- compute-optimized MinHash (for small scaled or large cardinalities) (#1045)

## [0.7.0] - 2020-05-12

Changed:

- Hide internal representation in core (#986)

Fixed: 

- update FFI and cbindgen (#986)

## [0.6.0] - 2020-04-28

Added:

- Nodegraph implementation based on khmer.Nodegraph (#799)

## [0.5.0] - 2020-02-08

Added:

- add_hash_with_abundance method in core library (#892)

Changed:

- More refactoring of MinHash comparison code (#882)
- Replace mins_push and abunds_push with set_abundances (#887)

Fixed:

- add_hash with num doesn't set abundances properly (#891)

## [0.4.0] - 2020-01-26

Added:

- Compute improvements: Parameter sets for defining signatures, add_protein implemented (#845)
- add_many for faster insertion of multiple hashes (#826)

Changed:

- Compare/similarity now have a downsample argument (#856)

Fixed:

- Improve sketching performance with lookup tables for complement and DNA validation (#861) (#865)
- Use tarpaulin instead of grcov (#862)
- set up publishing workflow for NPM and crates.io (#824)

## [0.3.0] - 2020-01-05

Added:

- Similarity with abundance method for MinHash (#808)
- Experimental support for indices in Rust (#773)
- Experimental SBT with MQF internal nodes in Rust (#772)

Changed:

- Make the sourmash crate library-only (#812)

Fixed:

- Use once_cell instead of lazy_static and lazy-init (#815)
- Fix mem leak in get_mins (#807)
- Fixes for WASI and WASM compilation (#771) (#723)

[unreleased]: https://github.com/sourmash-bio/sourmash/compare/r0.11.0...HEAD
[0.11.0]: https://github.com/sourmash-bio/sourmash/compare/r0.10.0...r0.11.0
[0.10.0]: https://github.com/sourmash-bio/sourmash/compare/r0.9.0...r0.10.0
[0.9.0]: https://github.com/sourmash-bio/sourmash/compare/r0.9.0...r0.10.0
[0.8.0]: https://github.com/sourmash-bio/sourmash/compare/r0.8.0...r0.9.0
[0.7.0]: https://github.com/sourmash-bio/sourmash/compare/r0.7.0...r0.8.0
[0.6.0]: https://github.com/sourmash-bio/sourmash/compare/r0.6.0...r0.7.0
[0.5.0]: https://github.com/sourmash-bio/sourmash/compare/r0.5.0...r0.6.0
[0.4.0]: https://github.com/sourmash-bio/sourmash/compare/r0.4.0...r0.5.0
[0.3.0]: https://github.com/sourmash-bio/sourmash/compare/r0.3.0...r0.4.0
