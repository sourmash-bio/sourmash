# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.13.0] - 2024-02-23

MSRV: 1.65

Changes/additions:

* Calculate all gather stats in rust; use for rocksdb gather (#2943)
* adjust protein ksize for record/manifest (#3019)
* Allow changing storage location for a collection in RevIndex (#3015)
* make core Manifest booleans python compatible (core) (#3007)

Updates:
* Bump roaring from 0.10.2 to 0.10.3 (#3014)
* Bump histogram from 0.9.0 to 0.9.1 (#3002)
* Bump chrono from 0.4.33 to 0.4.34 (#3000)
* Bump web-sys from 0.3.67 to 0.3.68 (#2998)
* Bump num-iter from 0.1.43 to 0.1.44 (#2997)
* Bump wasm-bindgen-test from 0.3.40 to 0.3.41 (#2996)

## [0.12.1] - 2024-02-10

MSRV: 1.65

Changes/additions:

* bump rust core version to r0.12.1 (#2988)
* Clean up and refactor `KmerMinHash::merge` in core (#2973)
* core: add scaled selection to manifest; add helper functions for collection and sig/sketch usage (#2948)
* core: enable downsample within select (#2931)

Updates:

* Deps: update typed-builder and histogram, bump MSRV to 1.65 (#2858)
* Bump tempfile from 3.9.0 to 3.10.0 (#2979)
* Bump rkyv from 0.7.43 to 0.7.44 (#2978)
* Bump memmap2 from 0.9.3 to 0.9.4 (#2958)
* Bump chrono from 0.4.31 to 0.4.33 (#2957)
* Bump serde from 1.0.195 to 1.0.196 (#2956)
* Bump serde_json from 1.0.111 to 1.0.113 (#2955)
* Bump web-sys from 0.3.66 to 0.3.67 (#2939)
* Bump wasm-bindgen-test from 0.3.39 to 0.3.40 (#2938)
* Bump rayon from 1.8.0 to 1.8.1 (#2937)
* Bump ouroboros from 0.18.2 to 0.18.3 (#2936)
* Bump histogram from 0.8.4 to 0.9.0 (#2935)
* Bump wasm-bindgen from 0.2.89 to 0.2.90 (#2925)
* Bump histogram from 0.8.3 to 0.8.4 (#2923)
* Bump serde_json from 1.0.110 to 1.0.111 (#2902)
* Bump serde from 1.0.194 to 1.0.195 (#2901)
* Bump serde_json from 1.0.108 to 1.0.110 (#2896)
* Bump ouroboros from 0.18.1 to 0.18.2 (#2894)
* Bump tempfile from 3.8.1 to 3.9.0 (#2893)
* Bump memmap2 from 0.9.2 to 0.9.3 (#2889)
* Bump memmap2 from 0.9.0 to 0.9.2 (#2882)
* Bump rkyv from 0.7.42 to 0.7.43 (#2880)
* Bump ouroboros from 0.18.0 to 0.18.1 (#2875)
* Bump once_cell from 1.18.0 to 1.19.0 (#2874)
* Bump rkyv from 0.7.40 to 0.7.42 (#2863)
* Bump csv from 1.2.0 to 1.3.0 (#2862)
* Bump roaring from 0.10.1 to 0.10.2 (#2865)
* Bump web-sys from 0.3.65 to 0.3.66 (#2864)
* Bump byteorder from 1.4.3 to 1.5.0 (#2866)
* Bump proptest from 1.3.1 to 1.4.0 (#2837)

## [0.12.0] - 2023-11-26

MSRV: 1.64

Added:

- Initial implementation for `Manifest`, `Selection`, and `Picklist` following
  the Python API. (#2230)
- `Collection` is a new abstraction for working with a set of signatures. A
  collection needs a `Storage` for holding the signatures (on-disk, in-memory,
  or remotely), and a `Manifest` to describe the metadata for each signature. (#2230)
- Expose CSV parsing and RocksDB errors. (#2230)
- New module `sourmash::index::revindex::disk_revindex` with the on-disk
  RevIndex implementation based on RocksDB. (#2230)
- Add `iter` and `iter_mut` methods for `Signature`. (#2230)
- Add `load_sig` and `save_sig` methods to `Storage` trait for higher-level data
  manipulation and caching. (#2230)
- Add `spec` method to `Storage` to allow constructing a concrete `Storage` from
  a string description. (#2230)
- Add `InnerStorage` for synchronizing parallel access to `Storage`
  implementations. (#2230)
- Add `MemStorage` for keeping signatures in-memory (mostly for debugging and
  testing). (#2230)
- Add new `branchwater` feature (enabled by default), which can be disabled by
  downstream projects to limit bringing heavy dependencies like rocksdb. (#2230)
- Add new `rkyv` feature (disabled by default), making `MinHash` serializable
  with the `rkyv` crate. (#2230)
- Add semver checks for CI (so we bump versions accordingly, or avoid breaking
  changes). (#2230)
- Add cargo deny config. (#2724)
- Benchmarks for seq_to_hashes in protein mode. (#1944)
- Oxidize ZipStorage. (#1909)
- Move greyhound-core into sourmash. (#1238)
- add `MinHash.kmers_and_hashes(...)` and `sourmash sig kmers`. (#1695)
- Produce list of hashes from a sequence. (#1653)

Changed:

- Rename `HashFunctions` variants to follow camel-case, so `Murmur64Protein`
  instead of `murmur64_protein`. (#2230)
- `LinearIndex` is now implemented as a thin layer on top of `Collection`. (#2230)
- Move `GatherResult` to `sourmash::index` module. (#2230)
- Move `sourmash::index::revindex` to `sourmash::index::mem_revindex` (this is
  the Greyhound version of revindex, in-memory only). It was also refactored
  internally to build a version of a `LinearIndex` that will be merged in the
  future with `sourmash::index::LinearIndex`. (#2230)
- Move `select` method from `Index` trait into a separate `Select` trait,
  and implement it for `Signature` based on the new `Selection` API. (#2230)
- Move `SigStore` into `sourmash::storage` module, and remove the generic. Now
  it always stores `Signature`. Also implement `Select` for it. (#2230)
- Disable `musllinux` wheels (need to figure out how to build rocksdb for it). (#2230)
- Reorganize traits for easier wasm and native compilation. (#1836)
- Adjust dayhoff and hp encodings to tolerate stop codons in the protein sequence. (#1673)

Fixed:

- Reduce features combinations on Rust checks (takes much less time to run). (#2230)
- Build: MSRV check for 1.64. (#2680)
- maturin: move deprecated definition from Cargo.toml to pyproject.toml. (#2597)
- Fix broken crates.io badge. (#2556)
- Fix unnecessary typecasts in Rust. (#2366)
- Fix `Signature.minhash` API during `sourmash sketch`. (#2329)
- Return Err for angular_similarity when abundance tracking is off. (#2327)
- Update various descriptions to talk about k-mers, not just DNA. (#2137)
- Fix downsample_scaled in `core`. (#2108)
- speed up `SeqToHashes` `translate`. (#1946)
- Speed-up `SeqToHashes()`. (#1938)
- Fix containment calculation for nodegraphs. (#1862)
- Fix panic bug in `sourmash sketch` dna with bad input and `--check-sequence`. (#1702)
- Fix Rust panic in `MinHash.seq_to_hashes`. (#1701)
- Beta lints. (#2841 #2630 #2596 #2298 #1791 #1786 #1760)

Removed:

- Remove BIGSI and SBT code. (#2732)

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
