#![deny(missing_docs)]

//! High-level access to CLI commands and utility functions
//!
//! This module mirrors parts of the [sourmash commands in Python][0],
//! and defines structures to make it easier to integrate large parameters sets
//! between the Rust and Python code.
//!
//! [0]: https://github.com/dib-lab/sourmash/blob/latest/src/sourmash/commands.py
//!
//! While Python is the primary consumer, this also works as a convenient function
//! for [building signatures from a set of sketcking parameters][1].
//!
//! [1]: https://github.com/luizirber/2020-11-02-greyhound/blob/ef457b7517b259246080d2bb22634bfd22683e6f/frontend/src/native_worker.rs#L41

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

use getset::{CopyGetters, Getters, Setters};
use typed_builder::TypedBuilder;

use crate::encodings::HashFunctions;
use crate::signature::Signature;
use crate::sketch::minhash::{max_hash_for_scaled, KmerMinHashBTree};
use crate::sketch::Sketch;

impl Signature {
    /// Build a new [`Signature`] from a set of [`ComputeParameters`].
    ///
    /// Since [`ComputeParameters`] has many default values, this is a convenient
    /// method for initializing sketches in a [`Signature`].
    pub fn from_params(params: &ComputeParameters) -> Signature {
        let template = build_template(params);

        Signature::builder()
            .hash_function("0.murmur64")
            .name(params.merge.clone())
            .filename(None)
            .signatures(template)
            .build()
    }
}

/// Parameters that can be used to construct sketches.
///
/// Initially derived from the parameters for [`sourmash compute`][0],
/// with the intention of making it easier to generate a set of sketches with
/// combinations of these parameters.
///
/// Some options are combinational (generate new sketches), while other are
/// exclusive (if set, all sketches will have this behavior). Setting too many
/// combinational options can lead to a large number of sketches being
/// generated.
///
/// ## Combinational
///
/// - ksizes
/// - dna
/// - dayhoff
/// - hp
/// - protein
///
/// ## Exclusive
///
/// - scaled
/// - num_hashes
/// - singleton
/// - name_from_first
/// - seed
/// - input_is_protein
/// - track_abundance
/// - merge
/// - license
///
/// ## Behavior when adding sequences
///
/// - check_sequence
/// - force
/// - randomize
/// - processes
/// - output
///
/// [0]: https://github.com/dib-lab/sourmash/blob/6b5806cf528583b864e1969739f65508c980ebd3/src/sourmash/cli/compute.py#L49L127
#[allow(dead_code)]
#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(TypedBuilder, CopyGetters, Getters, Setters)]
pub struct ComputeParameters {
    /// List of k-mer sizes to generate
    #[getset(get = "pub", set = "pub")]
    #[builder(default = vec![21, 31, 51])]
    ksizes: Vec<u32>,

    /// Complain if input sequence is invalid
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    check_sequence: bool,

    /// Build nucleotide sketches
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = true)]
    dna: bool,

    /// Build [Dayhoff-encoded](https://en.wikipedia.org/wiki/Margaret_Oakley_Dayhoff#Table_of_Dayhoff_encoding_of_amino_acids) sketches
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    dayhoff: bool,

    /// Build hydrophobic-polar-encoded amino acid signatures
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    hp: bool,

    /// Compute a sketch for each sequence record individually
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    singleton: bool,

    /// Choose number of hashes as 1 in FRACTION of input k-mers
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 0u64)]
    scaled: u64,

    /// Recompute signature even if the file exists
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    force: bool,

    /// Output computed signatures to this directory
    #[getset(get = "pub", set = "pub")]
    #[builder(default = None)]
    output: Option<String>, // TODO: check

    /// Number of hashes to use in each sketch
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 500u32)]
    num_hashes: u32,

    /// Build a protein signature (by default a nucleotide signature is used)
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    protein: bool,

    /// Name the signature using the name for first sequence in file
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    name_from_first: bool,

    /// seed used by MurmurHash; default = 42
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 42u64)]
    seed: u64,

    /// Consume protein sequences - no translation needed
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    input_is_protein: bool,

    /// Merge all input files into one signature file with the specified name
    #[getset(get = "pub", set = "pub")]
    #[builder(default = None)]
    merge: Option<String>,

    /// Track k-mer abundances in the generated signature
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    track_abundance: bool,

    /// Shuffle the list of input filenames
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = false)]
    randomize: bool,

    /// Signature license. Currently only CC0 is supported
    #[getset(get = "pub", set = "pub")]
    #[builder(default = "CC0".into())]
    license: String,

    /// Number of processes to use when adding sequences to sketches
    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default = 2usize)]
    processes: usize,
}

impl Default for ComputeParameters {
    fn default() -> Self {
        Self::builder().build()
    }
}

/// Build a collection of sketches from a set of [`ComputeParameters`] options.
pub fn build_template(params: &ComputeParameters) -> Vec<Sketch> {
    let max_hash = max_hash_for_scaled(params.scaled);

    params
        .ksizes
        .iter()
        .flat_map(|k| {
            let mut ksigs = vec![];

            if params.protein {
                ksigs.push(Sketch::LargeMinHash(
                    KmerMinHashBTree::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_protein)
                        .max_hash(max_hash)
                        .seed(params.seed)
                        .abunds(if params.track_abundance {
                            Some(Default::default())
                        } else {
                            None
                        })
                        .build(),
                ));
            }

            if params.dayhoff {
                ksigs.push(Sketch::LargeMinHash(
                    KmerMinHashBTree::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_dayhoff)
                        .max_hash(max_hash)
                        .seed(params.seed)
                        .abunds(if params.track_abundance {
                            Some(Default::default())
                        } else {
                            None
                        })
                        .build(),
                ));
            }

            if params.hp {
                ksigs.push(Sketch::LargeMinHash(
                    KmerMinHashBTree::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_hp)
                        .max_hash(max_hash)
                        .seed(params.seed)
                        .abunds(if params.track_abundance {
                            Some(Default::default())
                        } else {
                            None
                        })
                        .build(),
                ));
            }

            if params.dna {
                ksigs.push(Sketch::LargeMinHash(
                    KmerMinHashBTree::builder()
                        .num(params.num_hashes)
                        .ksize(*k)
                        .hash_function(HashFunctions::murmur64_DNA)
                        .max_hash(max_hash)
                        .seed(params.seed)
                        .abunds(if params.track_abundance {
                            Some(Default::default())
                        } else {
                            None
                        })
                        .build(),
                ));
            }

            ksigs
        })
        .collect()
}
