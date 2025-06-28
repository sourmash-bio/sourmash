# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [unreleased]

## [0.21.0] - 2025-06-28

MSRV: 1.74

Changes/additions:

* make `RevIndex.len()` and `RevIndex.signatures()` use picklist, if set (#3647)
* try fixing inline variables in rust `println!` (#3640)

Updates:

* Build(deps): Bump cfg-if from 1.0.0 to 1.0.1 (#3689)
* Build(deps): Bump roaring from 0.10.12 to 0.11.0 (#3702)
* Build(deps): Bump getset from 0.1.5 to 0.1.6 (#3700)
* Build(deps): Bump proptest from 1.6.0 to 1.7.0 (#3674)
* Build(deps): Bump camino from 1.1.9 to 1.1.10 (#3669)
* Build(deps): Bump criterion from 0.5.1 to 0.6.0 (#3655)
* Build(deps): Bump tempfile from 3.19.1 to 3.20.0 (#3639)

## [0.20.0] - 2025-05-06

MSRV: 1.74

Changes/additions:

* make the RocksDB handle directly accessible to external code (#3468)
* implement manifest retrieval from Rust via FFI for `RevIndex` (#3630)
* fully support skip-mers at the Python level; provide documentation (#3627)
* remove support for python 3.10 (#3606)
* fix linear gather in Rust (#3605)
* impl full mem-based `RevIndex` in Rust, and add Python support for mem- and disk-based `RevIndex` (#3545)
* Minhash deserialize hashfunction errorhandling (#3560)
* fix beta clippy errors (#3548)

Updates:

* Bump rand from 0.9.0 to 0.9.1 (#3620)
* Bump roaring from 0.10.10 to 0.10.12 (#3608)
* Bump log from 0.4.26 to 0.4.27 (#3587)
* Bump tempfile from 3.19.0 to 3.19.1 (#3588)
* Bump tempfile from 3.18.0 to 3.19.0 (#3582)
* Bump serde from 1.0.218 to 1.0.219 (#3576)
* Bump tempfile from 3.17.1 to 3.18.0 (#3575)
* Bump histogram from 0.11.2 to 0.11.3 (#3574)
* Bump serde_json from 1.0.139 to 1.0.140 (#3566)
* Bump getset from 0.1.4 to 0.1.5 (#3567)
* Bump needletail from 0.6.1 to 0.6.3 (#3553)
* Bump serde_json from 1.0.138 to 1.0.139 (#3552)
* Bump serde from 1.0.217 to 1.0.218 (#3550)
* Bump log from 0.4.25 to 0.4.26 (#3549)
* Bump tempfile from 3.16.0 to 3.17.1 (#3539)

## [0.19.0] - 2025-02-12

MSRV: 1.74

Changes/additions:

* update MSRV to 1.74, niffler to 3.0.0 (#3530)
* update to rocksdb 0.23 (#3456)
* remove finch conversion, support zstd and lzma in wasm (#3521)

Updates:

* Bump serde_json from 1.0.133 to 1.0.138 (#3453) (#3490) (#3500) (#3518)
* Bump tempfile from 3.14.0 to 3.16.0 (#3472) (#3519)
* Bump liblzma from 0.3.5 to 0.3.6 (#3526)
* Bump rand from 0.8.5 to 0.9.0 (#3512)
* Bump log from 0.4.22 to 0.4.25 (#3501)
* Bump histogram from 0.11.1 to 0.11.2 (#3498)
* Bump getset from 0.1.3 to 0.1.4 (#3499)
* Bump roaring from 0.10.9 to 0.10.10 (#3489)
* Bump ouroboros from 0.18.4 to 0.18.5 (#3491)
* Bump itertools from 0.13.0 to 0.14.0 (#3471)
* Bump serde from 1.0.216 to 1.0.217 (#3464)


## [0.18.0] - 2024-12-20

MSRV: 1.66

Changes/additions:

* add skipmer capacity to sourmash python layer via ffi (#3446)
* add skipmers; switch to reading frame approach for translation, skipmers (#3395)
* Refactor: Use to_writer/from_reader across the codebase (#3443)
* adjust `Signature::name()` to return `Option<String>` instead of `filename()` and `md5sum()` (#3434)
* propagate zipfile errors (#3431)

Updates:

* Bump proptest from 1.5.0 to 1.6.0 (#3437)
* Bump roaring from 0.10.8 to 0.10.9 (#3438)
* Bump serde from 1.0.215 to 1.0.216 (#3436)
* Bump statrs from 0.17.1 to 0.18.0 (#3426)
* Bump roaring from 0.10.7 to 0.10.8 (#3423)
* Bump needletail from 0.6.0 to 0.6.1 (#3427)
* Bump web-sys from 0.3.72 to 0.3.74 (#3411)
* Bump js-sys from 0.3.72 to 0.3.74 (#3412)
* Bump roaring from 0.10.6 to 0.10.7 (#3413)
* Bump serde_json from 1.0.132 to 1.0.133 (#3402)
* Bump serde from 1.0.214 to 1.0.215 (#3403)

## [0.17.2] - 2024-11-15

MSRV: 1.66

Changes/additions:

* enforce a single scaled on a `CollectionSet` (#3397)
* change `sig_from_record` to use scaled from `Record` to downsample (#3387)

Updates:

* Upgrade rocksdb to 0.22.0, bump MSRV to 1.66  (#3383)
* Bump thiserror from 1.0.68 to 2.0.3 (#3389)
* Bump csv from 1.3.0 to 1.3.1 (#3390)
* Bump tempfile from 3.13.0 to 3.14.0 (#3391)

## [0.17.1] - 2024-11-11

Changes/additions:
* fix: Avoid re-calculating md5sum on clone and conversion to KmerMinHashBTree (#3385)
* build: simplify Rust release (#3392)

## [0.17.0] - 2024-11-05

Changes/additions:
* standardize on u32 for scaled, and introduce `ScaledType` (#3364)
* panic when `FSStorage::load_sig` encounters more than one `Signature` in a JSON record (#3333)

Updates:

* Bump needletail from 0.5.1 to 0.6.0 (#3376)
* Bump histogram from 0.11.0 to 0.11.1 (#3377)
* Bump serde from 1.0.210 to 1.0.214 (#3368)
* Bump serde_json from 1.0.128 to 1.0.132 (#3358)
* Fix clippy lints from 1.83 beta (#3357)

## [0.16.0] - 2024-10-15

MSRV: 1.65

Changes/additions:

* refactor `calculate_gather_stats` to disallow repeated downsampling (#3352)
* improve downsampling behavior on `KmerMinHash`; fix `RevIndex::gather` bug around `scaled`. (#3342)
* derive Hash for `HashFunctions` (#3344)

Updates:

* Bump web-sys from 0.3.70 to 0.3.72 (#3354)
* Bump tempfile from 3.12.0 to 3.13.0 (#3340)


## [0.15.2] - 2024-09-25

MSRV: 1.65

Changes/additions:
* add `Manifest::intersect_manifest` to Rust core (#3305)
* propagate error from `RocksDB::open` on bad directory (#3306, #3307)

Updates:

* Bump getset from 0.1.2 to 0.1.3 (#3328)
* Bump memmap2 from 0.9.4 to 0.9.5 (#3326)
* Bump codspeed-criterion-compat from 2.6.0 to 2.7.2 (#3324)
* Bump serde_json from 1.0.127 to 1.0.128 (#3316)
* Bump serde from 1.0.209 to 1.0.210 (#3318)
* Bump serde from 1.0.208 to 1.0.209 (#3310)
* Bump serde_json from 1.0.125 to 1.0.127 (#3309)

## [0.15.1] - 2024-08-20

MSRV: 1.65

Changes/additions:

* Misc Rust updates to core (#3297)
* Implement resumability for revindex (#3275)
* Resolve issue for high precision MLE estimation (#3296)
* Added union method to HLL (#3293)

Updates:

* Bump camino from 1.1.7 to 1.1.9 (#3301)
* Bump web-sys from 0.3.69 to 0.3.70 (#3299)
* Bump serde_json from 1.0.120 to 1.0.125 (#3288) (#3280) (#3267) (#3302)
* Bump serde from 1.0.204 to 1.0.208 (#3289) (#3298)
* Bump tempfile from 3.10.1 to 3.12.0 (#3279) (#3287)

## [0.15.0] - 2024-07-27

MSRV: 1.65

Changes/additions:

* RocksDB storage and self-contained RevIndex with internal storage #3250
* Enable codspeed for Rust perf tracking (#3231)

Updates

* Bump roaring from 0.10.5 to 0.10.6 (#3245)
* Bump serde from 1.0.203 to 1.0.204 (#3244)
* Bump counter from 0.5.7 to 0.6.0 (#3235)
* Bump log from 0.4.21 to 0.4.22 (#3236)
* Bump serde_json from 1.0.117 to 1.0.120 (#3234)
* Bump proptest from 1.4.0 to 1.5.0 (#3222)

## [0.14.1] - 2024-06-19

MSRV: 1.65

Changes/additions:

* adjust how ANI is calculated in the revindex code. (#3218)

Updates:

* Bump histogram from 0.10.2 to 0.11.0 (#3216)
* Bump histogram from 0.10.1 to 0.10.2 (#3207)
* Bump statrs from 0.16.1 to 0.17.1 (#3205)
* Bump roaring from 0.10.4 to 0.10.5 (#3206)
* Bump primal-check from 0.3.3 to 0.3.4 (#3208)
* Bump niffler from 2.5.0 to 2.6.0 (#3204)

## [0.14.0] - 2024-06-10

MSRV: 1.65

Changes/additions:

* fix cargo fmt for updated `disk_revindex.rs` code (#3197)
* fix RocksDB-based gather & other rust-based infelicities revealed by plugins (#3193)
* use correct denominator in f_unique_to_query (#3138)
* fix clippy warnings about max_value (#3146)
* allow get/set record.filename (#3121)

Updates:

* Bump statrs from 0.16.0 to 0.16.1 (#3186)
* Bump serde from 1.0.202 to 1.0.203 (#3175)
* Bump ouroboros from 0.18.3 to 0.18.4 (#3176)
* Bump itertools from 0.12.1 to 0.13.0 (#3166)
* Bump camino from 1.1.6 to 1.1.7 (#3169)
* Bump serde from 1.0.201 to 1.0.202 (#3168)
* Bump serde_json from 1.0.116 to 1.0.117 (#3159)
* Bump serde from 1.0.200 to 1.0.201 (#3160)
* Bump roaring from 0.10.3 to 0.10.4 (#3142)
* Bump histogram from 0.10.0 to 0.10.1 (#3141)
* Bump num-iter from 0.1.44 to 0.1.45 (#3140)
* Bump serde from 1.0.199 to 1.0.200 (#3144)
* Bump serde from 1.0.198 to 1.0.199 (#3130)
* Bump serde_json from 1.0.115 to 1.0.116 (#3124)
* Bump serde from 1.0.197 to 1.0.198 (#3122)
* Bump histogram from 0.9.1 to 0.10.0 (#3109)
* Bump enum_dispatch from 0.3.12 to 0.3.13 (#3102)
* Bump serde_json from 1.0.114 to 1.0.115 (#3101)
* Bump rayon from 1.9.0 to 1.10.0 (#3098)

## [0.13.1] - 2024-03-23

MSRV: 1.65

Changes/additions:

* Implement file parsing for webassembly (#3047)
* fix `calculate_gather_stats` `threshold=0` bug (#3052)
* fix clippy beta issues (#3088)

Updates:

* Bump wasm-bindgen-test from 0.3.41 to 0.3.42 (#3063)
* Bump web-sys from 0.3.68 to 0.3.69 (#3061)
* Bump log from 0.4.20 to 0.4.21 (#3062)
* Bump rayon from 1.8.1 to 1.9.0 (#3058)
* Bump tempfile from 3.10.0 to 3.10.1 (#3059)
* Bump serde_json from 1.0.113 to 1.0.114 (#3044)
* Bump serde from 1.0.196 to 1.0.197 (#3045)
* Bump itertools from 0.12.0 to 0.12.1 (#3043)

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

[unreleased]: https://github.com/sourmash-bio/sourmash/compare/r0.19.0...HEAD
[0.19.0]: https://github.com/sourmash-bio/sourmash/compare/r0.18.0...r0.19.0
[0.18.0]: https://github.com/sourmash-bio/sourmash/compare/r0.17.2...r0.18.0
[0.17.2]: https://github.com/sourmash-bio/sourmash/compare/r0.17.1...r0.17.2
[0.17.1]: https://github.com/sourmash-bio/sourmash/compare/r0.17.0...r0.17.1
[0.17.0]: https://github.com/sourmash-bio/sourmash/compare/r0.16.0...r0.17.0
[0.16.0]: https://github.com/sourmash-bio/sourmash/compare/r0.15.1...r0.16.0
[0.15.1]: https://github.com/sourmash-bio/sourmash/compare/r0.15.0...r0.15.1
[0.15.0]: https://github.com/sourmash-bio/sourmash/compare/r0.14.1...r0.15.0
[0.14.1]: https://github.com/sourmash-bio/sourmash/compare/r0.14.0...r0.14.1
[0.14.0]: https://github.com/sourmash-bio/sourmash/compare/r0.13.1...r0.14.0
[0.13.1]: https://github.com/sourmash-bio/sourmash/compare/r0.13.0...r0.13.1
[0.13.0]: https://github.com/sourmash-bio/sourmash/compare/r0.12.1...r0.13.0
[0.12.1]: https://github.com/sourmash-bio/sourmash/compare/r0.12.0...r0.12.1
[0.12.0]: https://github.com/sourmash-bio/sourmash/compare/r0.11.0...r0.12.0
[0.11.0]: https://github.com/sourmash-bio/sourmash/compare/r0.10.0...r0.11.0
[0.10.0]: https://github.com/sourmash-bio/sourmash/compare/r0.9.0...r0.10.0
[0.9.0]: https://github.com/sourmash-bio/sourmash/compare/r0.9.0...r0.10.0
[0.8.0]: https://github.com/sourmash-bio/sourmash/compare/r0.8.0...r0.9.0
[0.7.0]: https://github.com/sourmash-bio/sourmash/compare/r0.7.0...r0.8.0
[0.6.0]: https://github.com/sourmash-bio/sourmash/compare/r0.6.0...r0.7.0
[0.5.0]: https://github.com/sourmash-bio/sourmash/compare/r0.5.0...r0.6.0
[0.4.0]: https://github.com/sourmash-bio/sourmash/compare/r0.4.0...r0.5.0
[0.3.0]: https://github.com/sourmash-bio/sourmash/compare/r0.3.0...r0.4.0
