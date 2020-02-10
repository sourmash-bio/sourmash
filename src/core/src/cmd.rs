#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

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

#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
pub struct ComputeParameters {
    pub(crate) ksizes: Vec<u32>,
    pub(crate) check_sequence: bool,
    pub(crate) dna: bool,
    pub(crate) dayhoff: bool,
    pub(crate) hp: bool,
    pub(crate) singleton: bool,
    pub(crate) count_valid_reads: usize,
    pub(crate) barcodes_file: Option<String>, // TODO: check
    pub(crate) line_count: usize,
    pub(crate) rename_10x_barcodes: Option<bool>, // TODO: check
    pub(crate) write_barcode_meta_csv: Option<bool>, // TODO: check
    pub(crate) save_fastas: Option<bool>,         // TODO: check
    pub(crate) scaled: u64,
    pub(crate) force: bool,
    pub(crate) output: Option<String>, // TODO: check
    pub(crate) num_hashes: u32,
    pub(crate) protein: bool,
    pub(crate) name_from_first: bool,
    pub(crate) seed: u64,
    pub(crate) input_is_protein: bool,
    pub(crate) merge: Option<String>,
    pub(crate) track_abundance: bool,
    pub(crate) randomize: bool,
    pub(crate) license: String,
    pub(crate) input_is_10x: bool,
    pub(crate) processes: usize,
}

impl Default for ComputeParameters {
    fn default() -> Self {
        ComputeParameters {
            ksizes: vec![21, 31, 51],
            check_sequence: false,
            dna: true,
            dayhoff: false,
            hp: false,
            singleton: false,
            count_valid_reads: 0,
            barcodes_file: None,
            line_count: 1500,
            rename_10x_barcodes: None,
            write_barcode_meta_csv: None,
            save_fastas: None,
            scaled: 0,
            force: false,
            output: None,
            num_hashes: 500,
            protein: false,
            name_from_first: false,
            seed: 42,
            input_is_protein: false,
            merge: None,
            track_abundance: false,
            randomize: false,
            license: "CC0".into(),
            input_is_10x: false,
            processes: 2,
        }
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
