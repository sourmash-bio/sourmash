use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};

use camino::Utf8PathBuf as PathBuf;
use log::info;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::collection::CollectionSet;
use crate::encodings::Idx;
use crate::index::{GatherResult, Index, Selection, SigCounter};
use crate::selection::Select;
use crate::signature::SigsTrait;
use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use crate::sketch::Sketch;
use crate::storage::SigStore;
use crate::Result;

/// Supports parallel search without a particular index.
pub struct LinearIndex {
    collection: CollectionSet,
    template: Sketch,
}

impl LinearIndex {
    pub fn from_collection(collection: CollectionSet) -> Self {
        let sig = collection.sig_for_dataset(0).unwrap();
        let template = sig.sketches().swap_remove(0);
        Self {
            collection,
            template,
        }
    }

    pub fn sig_for_dataset(&self, dataset_id: Idx) -> Result<SigStore> {
        self.collection.sig_for_dataset(dataset_id)
    }

    pub fn collection(&self) -> &CollectionSet {
        &self.collection
    }

    pub fn template(&self) -> &Sketch {
        &self.template
    }

    pub fn location(&self) -> Option<String> {
        unimplemented!()
    }

    pub fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
        let processed_sigs = AtomicUsize::new(0);

        let template = self.template();

        #[cfg(feature = "parallel")]
        let sig_iter = self.collection.par_iter();

        #[cfg(not(feature = "parallel"))]
        let sig_iter = self.collection.iter();

        let counters = sig_iter.filter_map(|(dataset_id, record)| {
            let filename = record.internal_location();

            let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {i} reference sigs");
            }

            let search_sig = self
                .collection
                .sig_for_dataset(dataset_id)
                .unwrap_or_else(|_| panic!("error loading {filename:?}"));

            let mut search_mh = None;
            if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
                search_mh = Some(mh);
            };
            let search_mh = search_mh.expect("Couldn't find a compatible MinHash");

            let (large_mh, small_mh) = if query.size() > search_mh.size() {
                (query, search_mh)
            } else {
                (search_mh, query)
            };

            let (size, _) = small_mh
                .intersection_size(large_mh)
                .unwrap_or_else(|_| panic!("error computing intersection for {filename:?}"));

            if size == 0 {
                None
            } else {
                let mut counter: SigCounter = Default::default();
                counter[&(dataset_id as Idx)] += size as usize;
                Some(counter)
            }
        });

        let reduce_counters = |mut a: SigCounter, b: SigCounter| {
            a.extend(&b);
            a
        };

        #[cfg(feature = "parallel")]
        let counter = counters.reduce(SigCounter::new, reduce_counters);

        #[cfg(not(feature = "parallel"))]
        let counter = counters.fold(SigCounter::new(), reduce_counters);

        counter
    }

    pub fn search(
        &self,
        counter: SigCounter,
        similarity: bool,
        threshold: usize,
    ) -> Result<Vec<String>> {
        let mut matches = vec![];
        if similarity {
            unimplemented!("TODO: threshold correction")
        }

        for (dataset_id, size) in counter.most_common() {
            if size >= threshold {
                matches.push(
                    self.collection
                        .record_for_dataset(dataset_id)?
                        .internal_location()
                        .to_string(),
                );
            } else {
                break;
            };
        }
        Ok(matches)
    }

    pub fn gather_round(
        &self,
        dataset_id: Idx,
        match_size: usize,
        query: &KmerMinHash,
        round: usize,
        orig_query: &KmerMinHash,
    ) -> Result<GatherResult> {
        let match_path = self
            .collection
            .record_for_dataset(dataset_id)?
            .internal_location()
            .into();
        let match_sig = self.collection.sig_for_dataset(dataset_id)?;
        let result = self.stats_for_match(
            match_sig,
            query,
            match_size,
            match_path,
            round as u32,
            orig_query,
        )?;
        Ok(result)
    }

    fn stats_for_match(
        &self,
        match_sig: SigStore,
        query: &KmerMinHash,
        match_size: usize,
        match_path: PathBuf,
        gather_result_rank: u32,
        orig_query: &KmerMinHash,
    ) -> Result<GatherResult> {
        let template = self.template();

        let mut match_mh = None;
        if let Some(Sketch::MinHash(mh)) = match_sig.select_sketch(template) {
            match_mh = Some(mh);
        }
        let match_mh = match_mh.expect("Couldn't find a compatible MinHash");

        // Calculate stats
        let f_orig_query = match_size as f64 / query.size() as f64;
        let f_match = match_size as f64 / match_mh.size() as f64;
        let filename = match_path.into_string();
        let name = match_sig.name();
        let unique_intersect_bp = (match_mh.scaled() as usize * match_size) as u64;

        let (intersect_hashes, _) = match_mh.intersection_size(query)?;
        let intersect_bp: u64 = match_mh.scaled() as u64 * intersect_hashes;

        let f_unique_to_query = intersect_hashes as f64 / orig_query.size() as f64;
        let match_ = match_sig;

        // TODO: all of these
        let f_unique_weighted = 0.;
        let average_abund = 0.;
        let median_abund = 0.;
        let std_abund = 0.;
        let md5 = "".into();
        let f_match_orig = 0.;
        let remaining_bp = 0;
        let total_weighted_hashes = 0;
        let n_unique_weighted_found = 0;
        let query_containment_ani = 0.0;
        let match_containment_ani = 0.0;
        let max_containment_ani = 0.0;
        let average_containment_ani = 0.0;
        let query_containment_ani_ci_low = None;
        let query_containment_ani_ci_high = None;
        let match_containment_ani_ci_low = None;
        let match_containment_ani_ci_high = None;
        let sum_weighted_found = 0;

        Ok(GatherResult {
            intersect_bp,
            f_orig_query,
            f_match,
            f_unique_to_query,
            f_unique_weighted,
            average_abund,
            median_abund,
            std_abund,
            filename,
            name,
            md5,
            match_,
            f_match_orig,
            unique_intersect_bp,
            gather_result_rank,
            remaining_bp,
            sum_weighted_found,
            total_weighted_hashes,
            n_unique_weighted_found,
            query_containment_ani,
            query_containment_ani_ci_low,
            query_containment_ani_ci_high,
            match_containment_ani,
            match_containment_ani_ci_low,
            match_containment_ani_ci_high,
            max_containment_ani,
            average_containment_ani,
        })
    }

    pub fn gather(
        &self,
        mut counter: SigCounter,
        threshold: usize,
        orig_query: &KmerMinHash,
    ) -> std::result::Result<Vec<GatherResult>, Box<dyn std::error::Error>> {
        let mut query = KmerMinHashBTree::from(orig_query.clone());
        let mut match_size = usize::MAX;
        let mut matches = vec![];
        let template = self.template();

        // iterate over matches, progressively removing intersections from the
        // counters.
        while match_size > threshold && !counter.is_empty() {
            let (dataset_id, size) = counter.most_common()[0];
            if threshold == 0 && size == 0 {
                break;
            }

            match_size = if size >= threshold {
                size
            } else {
                break;
            };

            let query_mh = KmerMinHash::from(query.clone());
            let result =
                self.gather_round(dataset_id, match_size, &query_mh, matches.len(), orig_query)?;

            // Prepare counter for finding the next match by decrementing
            // all hashes found in the current match in other datasets
            // TODO: maybe par_iter?
            let mut to_remove: HashSet<Idx> = Default::default();
            to_remove.insert(dataset_id);

            // retrieve the match
            let dataset_sig = self.collection.sig_for_dataset(dataset_id)?;
            let mut match_mh = None;
            if let Some(Sketch::MinHash(mh)) = dataset_sig.select_sketch(template) {
                match_mh = Some(mh);
            }
            let match_mh = match_mh.expect("Couldn't find a compatible MinHash");

            let (isect_hashes, _) = match_mh.intersection(&query_mh)?;
            let mut isect_mh = match_mh.clone();
            isect_mh.clear();
            let _ = isect_mh.add_many(&isect_hashes);

            query.remove_many(isect_mh.iter_mins().copied())?;

            // CTB: could redo this entire loop using a CounterGather-style
            // struct, with peek/consume, I 'spose.
            for (dataset, value) in counter.iter_mut() {
                let dataset_sig = self.collection.sig_for_dataset(*dataset)?;
                let mut match_mh = None;
                if let Some(Sketch::MinHash(mh)) = dataset_sig.select_sketch(template) {
                    match_mh = Some(mh);
                }
                let match_mh = match_mh.expect("Couldn't find a compatible MinHash");

                // take the intersection of this match with the best
                // intersection & remove from counter.
                let (matched_hashes, _) = isect_mh.intersection(match_mh)?;
                let mut this_isect_mh = match_mh.clone();
                this_isect_mh.clear();
                this_isect_mh.add_many(&matched_hashes)?;

                if this_isect_mh.size() > *value {
                    to_remove.insert(*dataset);
                } else {
                    *value -= this_isect_mh.size();
                };
            }

            to_remove.iter().for_each(|dataset_id| {
                counter.remove(dataset_id);
            });
            matches.push(result);
        }
        Ok(matches)
    }

    pub fn signatures_iter(&self) -> impl Iterator<Item = SigStore> + '_ {
        (0..self.collection.len()).map(move |dataset_id| {
            self.collection
                .sig_for_dataset(dataset_id as Idx)
                .expect("error loading sig")
        })
    }
}

impl Select for LinearIndex {
    fn select(self, selection: &Selection) -> Result<Self> {
        let Self {
            collection,
            template,
        } = self;
        let collection = collection.into_inner().select(selection)?.try_into()?;

        Ok(Self {
            collection,
            template,
        })
    }
}

impl Index<'_> for LinearIndex {
    type Item = SigStore;

    fn insert(&mut self, _node: Self::Item) -> Result<()> {
        unimplemented!()
    }

    fn save<P: AsRef<std::path::Path>>(&self, _path: P) -> Result<()> {
        unimplemented!()
    }

    fn load<P: AsRef<std::path::Path>>(_path: P) -> Result<()> {
        unimplemented!()
    }

    fn len(&self) -> usize {
        self.collection.len()
    }

    fn signatures(&self) -> Vec<Self::Item> {
        self.collection()
            .iter()
            .map(|(i, p)| {
                self.collection()
                    .sig_for_dataset(i as Idx)
                    .unwrap_or_else(|_| panic!("Error processing {}", p.internal_location()))
            })
            .collect()
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        unimplemented!()
    }
}
