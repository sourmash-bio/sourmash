//! Quickly search, compare, and analyze genomic and metagenomic data sets.
//!
//! This is the core library supporting sourmash, a command-line tool and Python
//! library for computing [MinHash sketches][0] from genomic sequences,
//! comparing them to each other, and plotting the results.
//! This allows you to estimate sequence similarity between even very
//! large data sets quickly and accurately.
//!
//! More information for sourmash as an application is available in [the
//! documentation](https://sourmash.readthedocs.io/en/latest/)
//!
//! [0]: https://en.wikipedia.org/wiki/MinHash
//!
//! This crate is organized using the following concepts:
//!
//! - A **dataset** contains genomic sequencing data. It can be assembled or raw
//!   reads from a sequencing experiment, or more generally anything
//!   representable with FASTA/FASTQ formats in a nucleotide ("ACGTN") or
//!   amino acid ("ACDEFGHIKLMNPQRSTVWY*") alphabet.
//!   Encodings are defined in the [`encodings`] submodule, as well as functions to
//!   work with them.
//!
//! - A **sketch** is a representation of a **dataset**. It is usually sublinear
//!   in size, meaning it uses less space than the original dataset, while still
//!   allowing limited operations (set membership, similarity or cardinality estimation).
//!   Examples include [MinHash][0], [Bloom Filter][1] or [HyperLogLog][2] sketches.
//!   Sketches are implemented in the [`sketch`] submodule.
//!
//! - A **signature** is a collection of **sketches** derived from the same
//!   **dataset**. Signatures can be compared as long as they have compatible
//!   sketches between them.
//!   Signatures are implemented in the [`signature`] submodule.
//!
//! - An **index** is a collection of **signatures**, possibly containing (or
//!   implemented as) optimized data structures allowing faster search.
//!   Indices are implemented in the [`index`] submodule.
//!
//!  [1]: https://en.wikipedia.org/wiki/Bloom_filter
//!  [2]: https://en.wikipedia.org/wiki/HyperLogLog

// TODO: remove this line and update all the appropriate type names for 1.0
#![allow(clippy::upper_case_acronyms)]

pub mod errors;
pub use errors::SourmashError as Error;

pub mod cmd;

pub mod index;

pub mod signature;
pub mod sketch;

pub mod encodings;

#[cfg(feature = "from-finch")]
pub mod from;

use cfg_if::cfg_if;
use murmurhash3::murmurhash3_x64_128;

cfg_if! {
    if #[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))] {
        pub mod wasm;
    } else {
        pub mod ffi;
    }
}

type HashIntoType = u64;

pub fn _hash_murmur(kmer: &[u8], seed: u64) -> u64 {
    murmurhash3_x64_128(kmer, seed).0
}
