#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

use getset::{CopyGetters, Getters, Setters};
use typed_builder::TypedBuilder;

use crate::index::MHBT;
use crate::signature::Signature;
use crate::sketch::minhash::{max_hash_for_scaled, HashFunctions, KmerMinHashBTree};
use crate::sketch::Sketch;
use crate::Error;

pub fn prepare(index_path: &str) -> Result<(), Error> {
    let mut index = MHBT::from_path(index_path)?;

    // TODO equivalent to fill_internal in python
    //unimplemented!();

    index.save_file(index_path, None)?;

    Ok(())
}

impl Signature {
    pub fn from_params(params: &ComputeParameters) -> Signature {
        let template = build_template(&params);

        Signature::builder()
            .hash_function("0.murmur64")
            .name(params.merge.clone())
            .filename(None)
            .signatures(template)
            .build()
    }
}

#[allow(dead_code)]
#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(TypedBuilder, CopyGetters, Getters, Setters)]
pub struct ComputeParameters {
    #[getset(get = "pub", set = "pub")]
    #[builder(default_code = "vec![21, 31, 51]")]
    ksizes: Vec<u32>,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    check_sequence: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "true")]
    dna: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    dayhoff: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    hp: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    singleton: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "0usize")]
    count_valid_reads: usize,

    #[getset(get = "pub", set = "pub")]
    #[builder(default_code = "None")]
    barcodes_file: Option<String>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "1500usize")]
    line_count: usize,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "None")]
    rename_10x_barcodes: Option<bool>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "None")]
    write_barcode_meta_csv: Option<bool>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "None")]
    save_fastas: Option<bool>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "0u64")]
    scaled: u64,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    force: bool,

    #[getset(get = "pub", set = "pub")]
    #[builder(default_code = "None")]
    output: Option<String>, // TODO: check

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "500u32")]
    num_hashes: u32,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    protein: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    name_from_first: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "42u64")]
    seed: u64,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    input_is_protein: bool,

    #[getset(get = "pub", set = "pub")]
    #[builder(default_code = "None")]
    merge: Option<String>,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    track_abundance: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    randomize: bool,

    #[getset(get = "pub", set = "pub")]
    #[builder(default_code = "default_license()")]
    license: String,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "false")]
    input_is_10x: bool,

    #[getset(get_copy = "pub", set = "pub")]
    #[builder(default_code = "2usize")]
    processes: usize,
}

fn default_license() -> String {
    "CC0".to_string()
}

impl Default for ComputeParameters {
    fn default() -> Self {
        Self::builder().build()
    }
}

pub fn build_template(params: &ComputeParameters) -> Vec<Sketch> {
    let max_hash = max_hash_for_scaled(params.scaled).unwrap_or(0);

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
