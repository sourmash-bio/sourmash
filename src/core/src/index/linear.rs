use std::collections::HashSet;
use std::sync::atomic::{AtomicUsize, Ordering};

use camino::Utf8PathBuf as PathBuf;
use log::info;

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::collection::CollectionSet;
use crate::encodings::Idx;
use crate::index::{GatherResult, Index, Selection, SigCounter};
use crate::manifest::Manifest;
use crate::selection::Select;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
use crate::storage::{InnerStorage, SigStore, Storage};
use crate::Result;

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
        if let Some(_storage) = &self.storage() {
            // storage.path()
            unimplemented!()
        } else {
            None
        }
    }

    pub fn storage(&self) -> Option<InnerStorage> {
        Some(self.collection.storage.clone())
    }

    pub fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
        let processed_sigs = AtomicUsize::new(0);

        let search_sigs: Vec<_> = self
            .collection
            .manifest
            .internal_locations()
            .map(PathBuf::from)
            .collect();

        let template = self.template();

        #[cfg(feature = "parallel")]
        let sig_iter = search_sigs.par_iter();

        #[cfg(not(feature = "parallel"))]
        let sig_iter = search_sigs.iter();

        let counters = sig_iter.enumerate().filter_map(|(dataset_id, filename)| {
            let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} reference sigs", i);
            }

            let search_sig = if let Some(storage) = &self.storage() {
                let sig_data = storage
                    .load(filename.as_str())
                    .unwrap_or_else(|_| panic!("error loading {:?}", filename));

                Signature::from_reader(sig_data.as_slice())
            } else {
                Signature::from_path(filename)
            }
            .unwrap_or_else(|_| panic!("Error processing {:?}", filename))
            .swap_remove(0);

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
                .unwrap_or_else(|_| panic!("error computing intersection for {:?}", filename));

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
                    self.collection.manifest[dataset_id as usize]
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
    ) -> Result<GatherResult> {
        let match_path = if self.collection.manifest.is_empty() {
            ""
        } else {
            self.collection.manifest[dataset_id as usize]
                .internal_location()
                .as_str()
        }
        .into();
        let match_sig = self.collection.sig_for_dataset(dataset_id)?;
        let result = self.stats_for_match(&match_sig, query, match_size, match_path, round)?;
        Ok(result)
    }

    fn stats_for_match(
        &self,
        match_sig: &Signature,
        query: &KmerMinHash,
        match_size: usize,
        match_path: PathBuf,
        gather_result_rank: usize,
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
        let unique_intersect_bp = match_mh.scaled() as usize * match_size;

        let (intersect_orig, _) = match_mh.intersection_size(query)?;
        let intersect_bp = (match_mh.scaled() * intersect_orig) as usize;

        let f_unique_to_query = intersect_orig as f64 / query.size() as f64;
        let match_ = match_sig.clone();

        // TODO: all of these
        let f_unique_weighted = 0.;
        let average_abund = 0;
        let median_abund = 0;
        let std_abund = 0;
        let md5 = "".into();
        let f_match_orig = 0.;
        let remaining_bp = 0;

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
        })
    }

    pub fn gather(
        &self,
        mut counter: SigCounter,
        threshold: usize,
        query: &KmerMinHash,
    ) -> std::result::Result<Vec<GatherResult>, Box<dyn std::error::Error>> {
        let mut match_size = usize::max_value();
        let mut matches = vec![];
        let template = self.template();

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

            let result = self.gather_round(dataset_id, match_size, query, matches.len())?;

            // Prepare counter for finding the next match by decrementing
            // all hashes found in the current match in other datasets
            // TODO: maybe par_iter?
            let mut to_remove: HashSet<Idx> = Default::default();
            to_remove.insert(dataset_id);

            for (dataset, value) in counter.iter_mut() {
                let dataset_sig = self.collection.sig_for_dataset(*dataset)?;
                let mut match_mh = None;
                if let Some(Sketch::MinHash(mh)) = dataset_sig.select_sketch(template) {
                    match_mh = Some(mh);
                }
                let match_mh = match_mh.expect("Couldn't find a compatible MinHash");

                let (intersection, _) = query.intersection_size(match_mh)?;
                if intersection as usize > *value {
                    to_remove.insert(*dataset);
                } else {
                    *value -= intersection as usize;
                };
            }
            to_remove.iter().for_each(|dataset_id| {
                counter.remove(dataset_id);
            });
            matches.push(result);
        }
        Ok(matches)
    }

    pub fn manifest(&self) -> Manifest {
        self.collection.manifest.clone()
    }

    pub fn set_manifest(&mut self, new_manifest: Manifest) -> Result<()> {
        self.collection.manifest = new_manifest;
        Ok(())
    }

    pub fn signatures_iter(&self) -> impl Iterator<Item = SigStore> + '_ {
        // FIXME temp solution, must find better one!
        (0..self.collection.manifest.len()).map(move |dataset_id| {
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

impl<'a> Index<'a> for LinearIndex {
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
        self.collection.manifest.len()
    }

    fn signatures(&self) -> Vec<Self::Item> {
        self.collection()
            .manifest
            .internal_locations()
            .map(PathBuf::from)
            .map(|p| {
                self.collection()
                    .storage
                    .load_sig(p.as_str())
                    .unwrap_or_else(|_| panic!("Error processing {:?}", p))
            })
            .collect()
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        unimplemented!()
    }
}
