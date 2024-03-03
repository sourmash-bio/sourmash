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
use stats::{median, stddev};
use typed_builder::TypedBuilder;

use crate::ani_utils::{ani_ci_from_containment, ani_from_containment};
use crate::encodings::Idx;
use crate::index::search::{search_minhashes, search_minhashes_containment};
use crate::prelude::*;
use crate::selection::Selection;
use crate::signature::SigsTrait;
use crate::sketch::minhash::KmerMinHash;
use crate::storage::SigStore;
use crate::Result;

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
    average_abund: f64,

    #[getset(get_copy = "pub")]
    median_abund: f64,

    #[getset(get_copy = "pub")]
    std_abund: f64,

    #[getset(get = "pub")]
    filename: String,

    #[getset(get = "pub")]
    name: String,

    #[getset(get = "pub")]
    md5: String,

    #[serde(skip)]
    match_: SigStore,

    #[getset(get_copy = "pub")]
    f_match_orig: f64,

    #[getset(get_copy = "pub")]
    unique_intersect_bp: usize,

    #[getset(get_copy = "pub")]
    gather_result_rank: usize,

    #[getset(get_copy = "pub")]
    remaining_bp: usize,

    #[getset(get_copy = "pub")]
    n_unique_weighted_found: usize,

    #[getset(get_copy = "pub")]
    total_weighted_hashes: usize,

    #[getset(get_copy = "pub")]
    sum_weighted_found: usize,

    #[getset(get_copy = "pub")]
    query_containment_ani: f64,

    #[getset(get_copy = "pub")]
    #[serde(skip_serializing_if = "Option::is_none")]
    query_containment_ani_ci_low: Option<f64>,

    #[getset(get_copy = "pub")]
    #[serde(skip_serializing_if = "Option::is_none")]
    query_containment_ani_ci_high: Option<f64>,

    #[getset(get_copy = "pub")]
    match_containment_ani: f64,

    #[getset(get_copy = "pub")]
    #[serde(skip_serializing_if = "Option::is_none")]
    match_containment_ani_ci_low: Option<f64>,

    #[getset(get_copy = "pub")]
    #[serde(skip_serializing_if = "Option::is_none")]
    match_containment_ani_ci_high: Option<f64>,

    #[getset(get_copy = "pub")]
    average_containment_ani: f64,

    #[getset(get_copy = "pub")]
    max_containment_ani: f64,
}

impl GatherResult {
    pub fn get_match(&self) -> Signature {
        self.match_.clone().into()
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

// note all mh should be selected/downsampled prior to being passed in here.
#[allow(clippy::too_many_arguments)]
pub fn calculate_gather_stats(
    orig_query: &KmerMinHash,
    query: KmerMinHash,
    match_sig: SigStore,
    match_size: usize,
    gather_result_rank: usize,
    sum_weighted_found: usize,
    total_weighted_hashes: usize,
    calc_abund_stats: bool,
    calc_ani_ci: bool,
    confidence: Option<f64>,
) -> Result<GatherResult> {
    // get match_mh
    let match_mh = match_sig.minhash().unwrap();
    //bp remaining in subtracted query
    let remaining_bp = (query.size() - match_size) * query.scaled() as usize;

    // stats for this match vs original query
    let (intersect_orig, _) = match_mh.intersection_size(orig_query).unwrap();
    let intersect_bp = (match_mh.scaled() * intersect_orig) as usize;
    let f_orig_query = intersect_orig as f64 / orig_query.size() as f64;
    let f_match_orig = intersect_orig as f64 / match_mh.size() as f64;

    // stats for this match vs current (subtracted) query
    let f_match = match_size as f64 / match_mh.size() as f64;
    let unique_intersect_bp = match_mh.scaled() as usize * match_size;
    let f_unique_to_query = match_size as f64 / query.size() as f64;

    // // get ANI values
    let ksize = match_mh.ksize() as f64;
    let query_containment_ani = ani_from_containment(f_unique_to_query, ksize);
    let match_containment_ani = ani_from_containment(f_match, ksize);
    let mut query_containment_ani_ci_low = None;
    let mut query_containment_ani_ci_high = None;
    let mut match_containment_ani_ci_low = None;
    let mut match_containment_ani_ci_high = None;

    if calc_ani_ci {
        let n_unique_kmers = match_mh.n_unique_kmers();
        let (qani_low, qani_high) = ani_ci_from_containment(
            f_unique_to_query,
            ksize,
            match_mh.scaled(),
            n_unique_kmers,
            confidence,
        )?;
        query_containment_ani_ci_low = Some(qani_low);
        query_containment_ani_ci_high = Some(qani_high);

        let (mani_low, mani_high) = ani_ci_from_containment(
            f_match,
            ksize,
            match_mh.scaled(),
            n_unique_kmers,
            confidence,
        )?;
        match_containment_ani_ci_low = Some(mani_low);
        match_containment_ani_ci_high = Some(mani_high);
    }

    let average_containment_ani = (query_containment_ani + match_containment_ani) / 2.0;
    let max_containment_ani = f64::max(query_containment_ani, match_containment_ani);

    // set up non-abundance weighted values
    let mut f_unique_weighted = f_unique_to_query;
    let mut average_abund = 1.0;
    let mut median_abund = 1.0;
    let mut std_abund = 0.0;
    // should these default to the unweighted numbers?
    let mut n_unique_weighted_found = 0;
    let mut sum_total_weighted_found = 0;

    // If abundance, calculate abund-related metrics (vs current query)
    if calc_abund_stats {
        // take abunds from subtracted query
        let (abunds, unique_weighted_found) = match match_mh.inflated_abundances(&query) {
            Ok((abunds, unique_weighted_found)) => (abunds, unique_weighted_found),
            Err(e) => {
                return Err(e);
            }
        };

        n_unique_weighted_found = unique_weighted_found as usize;
        sum_total_weighted_found = sum_weighted_found + n_unique_weighted_found;
        f_unique_weighted = n_unique_weighted_found as f64 / total_weighted_hashes as f64;

        average_abund = n_unique_weighted_found as f64 / abunds.len() as f64;

        // todo: try to avoid clone for these?
        median_abund = median(abunds.iter().cloned()).unwrap();
        std_abund = stddev(abunds.iter().cloned());
    }

    let result = GatherResult::builder()
        .intersect_bp(intersect_bp)
        .f_orig_query(f_orig_query)
        .f_match(f_match)
        .f_unique_to_query(f_unique_to_query)
        .f_unique_weighted(f_unique_weighted)
        .average_abund(average_abund)
        .median_abund(median_abund)
        .std_abund(std_abund)
        .filename(match_sig.filename())
        .name(match_sig.name())
        .md5(match_sig.md5sum())
        .match_(match_sig)
        .f_match_orig(f_match_orig)
        .unique_intersect_bp(unique_intersect_bp)
        .gather_result_rank(gather_result_rank)
        .remaining_bp(remaining_bp)
        .n_unique_weighted_found(n_unique_weighted_found)
        .query_containment_ani(query_containment_ani)
        .query_containment_ani_ci_low(query_containment_ani_ci_low)
        .query_containment_ani_ci_high(query_containment_ani_ci_high)
        .match_containment_ani_ci_low(match_containment_ani_ci_low)
        .match_containment_ani_ci_high(match_containment_ani_ci_high)
        .match_containment_ani(match_containment_ani)
        .average_containment_ani(average_containment_ani)
        .max_containment_ani(max_containment_ani)
        .sum_weighted_found(sum_total_weighted_found)
        .total_weighted_hashes(total_weighted_hashes)
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
    // use std::f64::EPSILON;
    // TODO: use f64::EPSILON when we bump MSRV
    const EPSILON: f64 = 0.01;

    #[test]
    fn test_calculate_gather_stats() {
        let scaled = 10;
        let params = ComputeParameters::builder()
            .ksizes(vec![31])
            .scaled(scaled)
            .build();

        let mut match_sig = Signature::from_params(&params);
        // create two minhash
        let mut match_mh = KmerMinHash::new(scaled, 31, HashFunctions::Murmur64Dna, 42, true, 0);
        match_mh.add_hash_with_abundance(1, 5);
        match_mh.add_hash_with_abundance(3, 3);
        match_mh.add_hash_with_abundance(5, 2);
        match_mh.add_hash_with_abundance(8, 2);
        match_mh.add_hash_with_abundance(11, 2); // Non-matching hash

        match_sig.reset_sketches();
        match_sig.push(Sketch::MinHash(match_mh.clone()));
        match_sig.set_filename("match-filename");
        match_sig.set_name("match-name");

        eprintln!("num_sketches: {:?}", match_sig.size());
        eprintln!("match_md5: {:?}", match_sig.md5sum());

        // Setup orig_query minhash with abundances and non-matching hash
        let mut orig_query = KmerMinHash::new(scaled, 31, HashFunctions::Murmur64Dna, 42, true, 0);
        orig_query.add_hash_with_abundance(1, 3);
        orig_query.add_hash_with_abundance(3, 2);
        orig_query.add_hash_with_abundance(5, 1);
        orig_query.add_hash_with_abundance(6, 1); // Non-matching hash
        orig_query.add_hash_with_abundance(8, 1);
        orig_query.add_hash_with_abundance(10, 1); // Non-matching hash

        let query = orig_query.clone();
        let total_weighted_hashes = orig_query.sum_abunds();

        let match_size = 4;
        let gather_result_rank = 0;
        let calc_abund_stats = true;
        let calc_ani_ci = false;
        let result = calculate_gather_stats(
            &orig_query,
            query,
            match_sig.into(),
            match_size,
            gather_result_rank,
            0,
            total_weighted_hashes.try_into().unwrap(),
            calc_abund_stats,
            calc_ani_ci,
            None,
        )
        .unwrap();
        // first, print all results
        assert_eq!(result.filename(), "match-filename");
        assert_eq!(result.name(), "match-name");
        assert_eq!(result.md5(), "f54b271a62fb7e2856e7b8a33e741b6e");
        assert_eq!(result.gather_result_rank, 0);
        assert_eq!(result.remaining_bp, 20);

        // results from match vs current query
        assert_eq!(result.f_match, 0.8);
        assert_eq!(result.unique_intersect_bp, 40);
        assert_eq!(result.f_unique_to_query, 4.0 / 6.0);
        eprintln!("{}", result.f_unique_weighted);
        assert_eq!(result.f_unique_weighted, 7. / 9.);
        assert_eq!(result.average_abund, 1.75);
        assert_eq!(result.median_abund, 1.5);
        assert_eq!(result.std_abund, 0.82915619758885);

        // results from match vs orig_query
        assert_eq!(result.intersect_bp, 40);
        assert_eq!(result.f_orig_query, 4.0 / 6.0);
        assert_eq!(result.f_match_orig, 4.0 / 5.0);

        assert!((result.average_containment_ani - 0.98991665567826).abs() < EPSILON);
        assert!((result.match_containment_ani - 0.9928276657672302).abs() < EPSILON);
        assert!((result.query_containment_ani - 0.9870056455892898).abs() < EPSILON);
        assert!((result.max_containment_ani - 0.9928276657672302).abs() < EPSILON);

        assert_eq!(result.total_weighted_hashes, 9);
        assert_eq!(result.n_unique_weighted_found, 7);
        assert_eq!(result.sum_weighted_found, 7);
    }
}
