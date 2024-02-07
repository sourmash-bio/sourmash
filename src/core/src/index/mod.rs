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
use crate::signature::SigsTrait;
use crate::sketch::minhash::KmerMinHash;
use crate::storage::SigStore;
use crate::Result;

// Todo: evaluate get vs get_copy for added fields
#[derive(TypedBuilder, CopyGetters, Getters, Setters, Serialize, Deserialize, Debug, PartialEq)]
pub struct GatherResult {
    #[getset(get_copy = "pub")]
    intersect_bp: usize,

    #[getset(get_copy = "pub")]
    f_orig_query: f64,

    #[getset(get_copy = "pub")]
    f_match: f64,

    #[getset(get_copy = "pub")]
    f_unique_to_query: f64,

    #[getset(get_copy = "pub")]
    f_unique_weighted: f64,

    #[getset(get_copy = "pub")]
    average_abund: usize,

    #[getset(get_copy = "pub")]
    median_abund: usize,

    #[getset(get_copy = "pub")]
    std_abund: usize,

    #[getset(get = "pub")]
    filename: String,

    #[getset(get = "pub")]
    name: String,

    #[getset(get = "pub")]
    md5: String,

    #[serde(skip)]
    match_: Signature,

    #[getset(get_copy = "pub")]
    f_match_orig: f64,

    #[getset(get_copy = "pub")]
    unique_intersect_bp: usize,

    #[getset(get_copy = "pub")]
    gather_result_rank: usize,

    #[getset(get_copy = "pub")]
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
    orig_query: KmerMinHash, // make static/ref?

    #[serde(skip)]
    query: KmerMinHash,

    #[serde(skip)]
    match_: SigStore, // sigstore ok? or signature?

    match_size: usize,

    remaining_hashes: usize,

    gather_result_rank: usize,
    // total_orig_query_abund: u64,
}

impl FastGatherResult {
    pub fn get_match(&self) -> Signature {
        self.match_.clone().into()
    }
}

pub fn calculate_gather_stats(fgres: FastGatherResult) -> Result<GatherResult> {
    // Calculate stats
    let name = fgres.match_.name();
    let md5 = fgres.match_.md5sum();
    let match_mh = fgres
        .match_
        .minhash()
        .expect("Couldn't find a compatible MinHash");

    // basics
    let filename = fgres.match_.filename();
    //bp remaining in subtracted query
    let remaining_bp = fgres.remaining_hashes * match_mh.scaled() as usize;

    // stats for this match vs original query
    // downsample here?? this seems like unnecesary work for each round -- how to get around?
    let orig_query_ds = &fgres.orig_query.downsample_scaled(match_mh.scaled())?;

    let (intersect_orig, _) = match_mh.intersection_size(&orig_query_ds).unwrap(); //?;
    let intersect_bp = (match_mh.scaled() * intersect_orig) as usize;
    let f_orig_query = intersect_orig as f64 / fgres.orig_query.size() as f64;
    let f_match_orig = intersect_orig as f64 / match_mh.size() as f64;

    // stats for this match vs current (subtracted) query
    let f_match = fgres.match_size as f64 / match_mh.size() as f64;
    let unique_intersect_bp = match_mh.scaled() as usize * fgres.match_size;
    let f_unique_to_query = fgres.match_size as f64 / fgres.query.size() as f64;

    // set up non-abundance weighted values
    let mut f_unique_weighted = f_unique_to_query;
    let mut average_abund = 1.0;
    let mut median_abund = 1.0;
    let mut std_abund = 0.0;

    // If abundance, calculate abund-related metrics (vs current query)
    if fgres.query.track_abundance() {
        // need current downsampled query here to get f_unique_weighted
        let (abunds, matched_abund) = match match_mh.inflated_abundances(&fgres.query) {
            Ok((abunds, matched_abund)) => (abunds, matched_abund),
            Err(e) => {
                return Err(e);
            }
        };
        let match_mh_total_abund = match_mh.sum_abunds();
        f_unique_weighted = match_mh_total_abund as f64 / matched_abund as f64;

        average_abund = matched_abund as f64 / abunds.len() as f64;

        // todo
        median_abund = 1.0;
        std_abund = 0.0;
        // let median_abund: f64 = median(&common_abunds.iter().map(|&x| x as f64).collect::<Vec<f64>>());
        // let std_abund: f64 = statistics::std_dev(&common_abunds.iter().map(|&x| x as f64).collect::<Vec<f64>>(), None)?;
    }

    // do want the f_weighted / abundance info for the original query?
    // f_weighted = match_mh_total_abund as f64 / fgres.total_orig_query_abund as f64;

    let result = GatherResult::builder()
        .intersect_bp(intersect_bp)
        .f_orig_query(f_orig_query)
        .f_match(f_match)
        .f_unique_to_query(f_unique_to_query)
        .f_unique_weighted(f_unique_weighted)
        .average_abund(average_abund as usize)
        .median_abund(median_abund as usize)
        .std_abund(std_abund as usize)
        .filename(filename)
        .name(name)
        .md5(md5)
        .match_(fgres.match_.into())
        .f_match_orig(f_match_orig)
        .unique_intersect_bp(unique_intersect_bp)
        .gather_result_rank(fgres.gather_result_rank)
        .remaining_bp(remaining_bp)
        .build();
    Ok(result)
}

#[cfg(test)]
mod test_calculate_gather_stats {
    use super::*;
    use crate::cmd::ComputeParameters;
    use crate::encodings::HashFunctions;
    use crate::signature::Signature;
    use crate::sketch::minhash::KmerMinHash;
    use crate::sketch::Sketch;
    // TODO: use f64::EPSILON when we bump MSRV
    const EPSILON: f64 = 0.01;

    #[test]
    fn test_calculate_gather_stats() {
        let scaled = 10;
        let params = ComputeParameters::builder()
            .ksizes(vec![3])
            .scaled(scaled)
            .build();

        let mut match_sig = Signature::from_params(&params);
        // create two minhash
        let mut match_mh = KmerMinHash::new(scaled, 3, HashFunctions::Murmur64Dna, 42, true, 0);
        match_mh.add_hash(1);
        match_mh.add_hash(2);
        match_mh.add_hash(3);
        match_mh.add_hash(5);

        match_sig.reset_sketches();
        match_sig.push(Sketch::MinHash(match_mh));
        match_sig.set_filename("match-filename");
        match_sig.set_name("match-name");

        eprintln!("num_sketches: {:?}", match_sig.size());
        eprintln!("match_md5: {:?}", match_sig.md5sum());

        // Setup orig_query minhash with abundances and non-matching hash
        let mut orig_query = KmerMinHash::new(scaled, 3, HashFunctions::Murmur64Dna, 42, true, 0);
        orig_query.add_hash_with_abundance(1, 1);
        orig_query.add_hash_with_abundance(2, 3);
        orig_query.add_hash_with_abundance(4, 6); // Non-matching hash
        orig_query.add_hash_with_abundance(5, 9);
        orig_query.add_hash_with_abundance(6, 1);

        let mut query = orig_query.clone();
        let rm_hashes = vec![1];
        query.remove_many(rm_hashes.as_slice()).unwrap(); // remove hash 1

        // build FastGatherResult
        let fgres = FastGatherResult::builder()
            .orig_query(orig_query)
            .query(query)
            .match_(match_sig.into())
            .match_size(2) // 2  -- only 2 hashes match, one was previously consumed
            .remaining_hashes(30) // arbitrary
            .gather_result_rank(5) // arbitrary
            // .total_orig_query_abund(20) // sum of orig_query abundances
            .build();

        let result = calculate_gather_stats(fgres).unwrap();
        // first, print all results
        eprintln!("result: {:?}", result);
        assert_eq!(result.filename(), "match-filename");
        assert_eq!(result.name(), "match-name");
        assert_eq!(result.md5(), "d70f195edbc052d647a288e3d46b3b2e");
        assert_eq!(result.gather_result_rank, 5);
        assert_eq!(result.remaining_bp, 300);

        // results from match vs subtracted query
        assert_eq!(result.f_match, 0.5);
        assert_eq!(result.unique_intersect_bp, 20);
        assert_eq!(result.f_unique_to_query, 0.5);
        assert!((result.f_unique_weighted - 0.333).abs() < EPSILON);
        assert_eq!(result.average_abund, 6);
        // todo: implement median and std
        assert_eq!(result.median_abund, 1);
        assert_eq!(result.std_abund, 0);

        // results from match vs orig_query
        assert_eq!(result.intersect_bp, 30);
        assert_eq!(result.f_orig_query, 0.6);
        assert_eq!(result.f_match_orig, 0.75);
    }
}
