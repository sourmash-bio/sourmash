//! # Indexing structures for fast similarity search
//!
//! An index organizes signatures to allow for fast similarity search.
//! Some indices also support containment searches.

pub mod linear;

#[cfg(not(target_arch = "wasm32"))]
#[cfg(feature = "branchwater")]
pub mod revindex;

pub mod search;

use std::path::Path;

use getset::{CopyGetters, Getters, Setters};

use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::encodings::Idx;
use crate::index::search::{search_minhashes, search_minhashes_containment};
use crate::prelude::*;
use crate::Result;
// use crate::sketch::hyperloglog::HyperLogLog;
// use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use crate::signature::SigsTrait;
use crate::sketch::minhash::KmerMinHash;

#[derive(TypedBuilder, CopyGetters, Getters, Setters, Serialize, Deserialize, Debug, PartialEq)]
pub struct GatherResult {
    #[getset(get_copy = "pub")]
    intersect_bp: usize,

    #[getset(get_copy = "pub")]
    f_orig_query: f64,

    #[getset(get_copy = "pub")]
    f_match: f64,

    f_unique_to_query: f64,
    f_unique_weighted: f64,
    average_abund: usize,
    median_abund: usize,
    std_abund: usize,

    #[getset(get = "pub")]
    filename: String,

    #[getset(get = "pub")]
    name: String,

    #[getset(get = "pub")]
    md5: String,

    #[serde(skip)]
    match_: Signature,

    f_match_orig: f64,
    unique_intersect_bp: usize,
    gather_result_rank: usize,
    remaining_bp: usize,
}

impl GatherResult {
    pub fn get_match(&self) -> Signature {
        self.match_.clone()
    }
}

type SigCounter = counter::Counter<Idx>;

pub trait Index<'a> {
    type Item: Comparable<Self::Item>;
    //type SignatureIterator: Iterator<Item = Self::Item>;

    fn find<F>(&self, search_fn: F, sig: &Self::Item, threshold: f64) -> Result<Vec<&Self::Item>>
    where
        F: Fn(&dyn Comparable<Self::Item>, &Self::Item, f64) -> bool,
    {
        Ok(self
            .signature_refs()
            .into_iter()
            .flat_map(|node| {
                if search_fn(&node, sig, threshold) {
                    Some(node)
                } else {
                    None
                }
            })
            .collect())
    }

    fn search(
        &self,
        sig: &Self::Item,
        threshold: f64,
        containment: bool,
    ) -> Result<Vec<&Self::Item>> {
        if containment {
            self.find(search_minhashes_containment, sig, threshold)
        } else {
            self.find(search_minhashes, sig, threshold)
        }
    }

    //fn gather(&self, sig: &Self::Item, threshold: f64) -> Result<Vec<&Self::Item>>;

    fn insert(&mut self, node: Self::Item) -> Result<()>;

    fn batch_insert(&mut self, nodes: Vec<Self::Item>) -> Result<()> {
        for node in nodes {
            self.insert(node)?;
        }

        Ok(())
    }

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<()>;

    fn load<P: AsRef<Path>>(path: P) -> Result<()>;

    fn signatures(&self) -> Vec<Self::Item>;

    fn signature_refs(&self) -> Vec<&Self::Item>;

    fn len(&self) -> usize {
        self.signature_refs().len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /*
    fn iter_signatures(&self) -> Self::SignatureIterator;
    */
}

impl<'a, N, L> Comparable<L> for &'a N
where
    N: Comparable<L>,
{
    fn similarity(&self, other: &L) -> f64 {
        (*self).similarity(other)
    }

    fn containment(&self, other: &L) -> f64 {
        (*self).containment(other)
    }
}

#[derive(TypedBuilder, CopyGetters, Getters, Setters, Serialize, Deserialize, Debug)]
pub struct FastGatherResult {
    #[serde(skip)]
    orig_query: KmerMinHash, // can we make this a reference?
    #[serde(skip)]
    query: KmerMinHash,

    #[serde(skip)]
    match_: Signature,

    remaining_hashes: usize,

    gather_result_rank: usize,

    total_orig_query_abund: u64,

    #[serde(skip)] // do we need this at all?
    selection: Selection, // borrow / use static??
}

impl FastGatherResult {
    pub fn get_match(&self) -> Signature {
        self.match_.clone()
    }
}

// let match_mh = fgres.match_.select(&selection).minhash();
//how does the sig.sketches().swap_remove(0) in prepare_query compare with sig.minhash()?
// do we need to do prepare_query here?
// let match_mh = prepare_query(match_sig.into(), &selection)
//     .expect("Couldn't find a compatible MinHash");

pub fn calculate_gather_stats(fgres: FastGatherResult) -> GatherResult {
    // Calculate stats
    let name = fgres.match_.name();
    let md5 = fgres.match_.md5sum();
    let match_mh = fgres
        .match_
        .minhash()
        .expect("Couldn't find a compatible MinHash");
    let match_size = match_mh.size();
    let f_match = match_size as f64 / match_mh.size() as f64;
    let unique_intersect_bp = match_mh.scaled() as usize * match_size;
    let (intersect_orig, _) = match_mh.intersection_size(&fgres.orig_query).unwrap(); //?;
    let intersect_bp = (match_mh.scaled() * intersect_orig) as usize;

    let f_unique_to_query = match_size as f64 / fgres.orig_query.size() as f64;
    let f_orig_query = intersect_orig as f64 / fgres.orig_query.size() as f64;
    let f_match_orig = intersect_orig as f64 / match_mh.size() as f64;

    let filename = fgres.match_.filename();

    // todo
    let f_unique_weighted = 0.0;
    let average_abund: usize = 0;
    let median_abund: usize = 0;
    let std_abund: usize = 0;

    // weight common by query abundances
    // let (common_abunds, total_common_weighted) = match_mh.weighted_intersect_size(match_mh, abunds_from=fgres.orig_query);
    // //Calculate abund-related metrics
    // let f_unique_weighted = total_common_weighted as f64 / fgres.total_orig_query_abund as f64;
    // // mean, median, std of abundances
    // let average_abund: f64 = mean(&common_abunds.iter().map(|&x| x as f64).collect::<Vec<f64>>());
    // let median_abund: f64 = median(&common_abunds.iter().map(|&x| x as f64).collect::<Vec<f64>>());
    // let std_abund: f64 = statistics::std_dev(&common_abunds.iter().map(|&x| x as f64).collect::<Vec<f64>>(), None)?;

    // remaining_bp is the number of base pairs in the query that are not in the match (or any prior match)
    let remaining_bp = fgres.remaining_hashes * match_mh.scaled() as usize;
    let result = GatherResult::builder()
        .intersect_bp(intersect_bp)
        .f_orig_query(f_orig_query)
        .f_match(f_match)
        .f_unique_to_query(f_unique_to_query)
        .f_unique_weighted(f_unique_weighted)
        .average_abund(average_abund)
        .median_abund(median_abund)
        .std_abund(std_abund)
        .filename(filename)
        .name(name)
        .md5(md5)
        .match_(fgres.match_.into())
        .f_match_orig(f_match_orig)
        .unique_intersect_bp(unique_intersect_bp)
        .gather_result_rank(fgres.gather_result_rank)
        .remaining_bp(remaining_bp)
        .build();
    result
}
