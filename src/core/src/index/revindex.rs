use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use getset::{CopyGetters, Getters, Setters};
use log::{debug, info};
use nohash_hasher::BuildNoHashHasher;
use serde::{Deserialize, Serialize};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

use crate::encodings::{Color, Colors, Idx};
use crate::index::{Index, Selection, SigStore};
use crate::manifest::Manifest;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
use crate::storage::{Storage, ZipStorage};
use crate::Error;
use crate::HashIntoType;

type SigCounter = counter::Counter<Idx>;

#[derive(Serialize, Deserialize)]
struct HashToColor(HashMap<HashIntoType, Color, BuildNoHashHasher<HashIntoType>>);

impl HashToColor {
    fn new() -> Self {
        HashToColor(HashMap::<
            HashIntoType,
            Color,
            BuildNoHashHasher<HashIntoType>,
        >::with_hasher(BuildNoHashHasher::default()))
    }

    fn get(&self, hash: &HashIntoType) -> Option<&Color> {
        self.0.get(hash)
    }

    fn retain(&mut self, hashes: &HashSet<HashIntoType>) {
        self.0.retain(|hash, _| hashes.contains(hash))
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn add_to(&mut self, colors: &mut Colors, dataset_id: usize, matched_hashes: Vec<u64>) {
        let mut color = None;

        matched_hashes.into_iter().for_each(|hash| {
            color = Some(colors.update(color, &[dataset_id as Idx]).unwrap());
            self.0.insert(hash, color.unwrap());
        });
    }

    fn reduce_hashes_colors(
        a: (HashToColor, Colors),
        b: (HashToColor, Colors),
    ) -> (HashToColor, Colors) {
        let ((small_hashes, small_colors), (mut large_hashes, mut large_colors)) =
            if a.0.len() > b.0.len() {
                (b, a)
            } else {
                (a, b)
            };

        small_hashes.0.into_iter().for_each(|(hash, color)| {
            large_hashes
                .0
                .entry(hash)
                .and_modify(|entry| {
                    // Hash is already present.
                    // Update the current color by adding the indices from
                    // small_colors.
                    let ids = small_colors.indices(&color);
                    let new_color = large_colors.update(Some(*entry), ids).unwrap();
                    *entry = new_color;
                })
                .or_insert_with(|| {
                    // In this case, the hash was not present yet.
                    // we need to create the same color from small_colors
                    // into large_colors.
                    let ids = small_colors.indices(&color);
                    let new_color = large_colors.update(None, ids).unwrap();
                    assert_eq!(new_color, color);
                    new_color
                });
        });

        (large_hashes, large_colors)
    }
}

// Use rkyv for serialization?
// https://davidkoloski.me/rkyv/
#[derive(Serialize, Deserialize)]
pub struct RevIndex {
    linear: LinearRevIndex,
    hash_to_color: HashToColor,
    colors: Colors,
}

#[derive(Serialize, Deserialize)]
pub struct LinearRevIndex {
    sig_files: Manifest,

    #[serde(skip)]
    ref_sigs: Option<Vec<SigStore>>,

    template: Sketch,

    #[serde(skip)]
    storage: Option<Arc<ZipStorage>>,
}

impl LinearRevIndex {
    pub fn new(
        sig_files: Option<Manifest>,
        template: &Sketch,
        keep_sigs: bool,
        ref_sigs: Option<Vec<Signature>>,
        storage: Option<ZipStorage>,
    ) -> Self {
        if ref_sigs.is_none() && sig_files.is_none() {
            todo!("throw error, one need to be set");
        }

        let ref_sigs = if let Some(ref_sigs) = ref_sigs {
            Some(ref_sigs.into_iter().map(|m| m.into()).collect())
        } else if keep_sigs {
            let search_sigs: Vec<_> = sig_files
                .as_ref()
                .unwrap()
                .internal_locations()
                .map(PathBuf::from)
                .collect();

            #[cfg(feature = "parallel")]
            let sigs_iter = search_sigs.par_iter();

            #[cfg(not(feature = "parallel"))]
            let sigs_iter = search_sigs.iter();

            Some(
                sigs_iter
                    .map(|ref_path| {
                        if let Some(storage) = &storage {
                            let sig_data = storage
                                .load(ref_path.to_str().unwrap_or_else(|| {
                                    panic!("error converting path {:?}", ref_path)
                                }))
                                .unwrap_or_else(|_| panic!("error loading {:?}", ref_path));
                            Signature::from_reader(sig_data.as_slice())
                                .unwrap_or_else(|_| panic!("Error processing {:?}", ref_path))
                                .swap_remove(0)
                                .into()
                        } else {
                            Signature::from_path(&ref_path)
                                .unwrap_or_else(|_| panic!("Error processing {:?}", ref_path))
                                .swap_remove(0)
                                .into()
                        }
                    })
                    .collect(),
            )
        } else {
            None
        };

        let storage = storage.map(Arc::new);

        LinearRevIndex {
            sig_files: sig_files.unwrap(),
            template: template.clone(),
            ref_sigs,
            storage,
        }
    }

    fn index(
        self,
        threshold: usize,
        merged_query: Option<KmerMinHash>,
        queries: Option<&[KmerMinHash]>,
    ) -> RevIndex {
        let processed_sigs = AtomicUsize::new(0);

        let search_sigs: Vec<_> = self
            .sig_files
            .internal_locations()
            .map(PathBuf::from)
            .collect();

        #[cfg(feature = "parallel")]
        let sig_iter = search_sigs.par_iter();

        #[cfg(not(feature = "parallel"))]
        let sig_iter = search_sigs.iter();

        let filtered_sigs = sig_iter.enumerate().filter_map(|(dataset_id, filename)| {
            let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} reference sigs", i);
            }

            let search_sig = if let Some(storage) = &self.storage {
                let sig_data = storage
                    .load(
                        filename
                            .to_str()
                            .unwrap_or_else(|| panic!("error converting path {:?}", filename)),
                    )
                    .unwrap_or_else(|_| panic!("error loading {:?}", filename));

                Signature::from_reader(sig_data.as_slice())
            } else {
                Signature::from_path(&filename)
            }
            .unwrap_or_else(|_| panic!("Error processing {:?}", filename))
            .swap_remove(0);

            RevIndex::map_hashes_colors(
                dataset_id,
                &search_sig,
                queries,
                &merged_query,
                threshold,
                &self.template,
            )
        });

        #[cfg(feature = "parallel")]
        let (hash_to_color, colors) = filtered_sigs.reduce(
            || (HashToColor::new(), Colors::default()),
            HashToColor::reduce_hashes_colors,
        );

        #[cfg(not(feature = "parallel"))]
        let (hash_to_color, colors) = filtered_sigs.fold(
            (HashToColor::new(), Colors::default()),
            HashToColor::reduce_hashes_colors,
        );

        RevIndex {
            hash_to_color,
            colors,
            linear: self,
        }
    }

    pub fn location(&self) -> Option<String> {
        if let Some(storage) = &self.storage {
            storage.path()
        } else {
            None
        }
    }

    pub fn storage(&self) -> Option<Arc<ZipStorage>> {
        self.storage.clone()
    }

    pub fn select(mut self, selection: &Selection) -> Result<Self, Error> {
        let manifest = self.sig_files.select_to_manifest(selection)?;
        self.sig_files = manifest;

        Ok(self)
        /*
        # if we have a manifest, run 'select' on the manifest.
        manifest = self.manifest
        traverse_yield_all = self.traverse_yield_all

        if manifest is not None:
            manifest = manifest.select_to_manifest(**kwargs)
            return ZipFileLinearIndex(self.storage,
                                      selection_dict=None,
                                      traverse_yield_all=traverse_yield_all,
                                      manifest=manifest,
                                      use_manifest=True)
        else:
            # no manifest? just pass along all the selection kwargs to
            # the new ZipFileLinearIndex.

            assert manifest is None
            if self.selection_dict:
                # combine selects...
                d = dict(self.selection_dict)
                for k, v in kwargs.items():
                    if k in d:
                        if d[k] is not None and d[k] != v:
                            raise ValueError(f"incompatible select on '{k}'")
                    d[k] = v
                kwargs = d

            return ZipFileLinearIndex(self.storage,
                                      selection_dict=kwargs,
                                      traverse_yield_all=traverse_yield_all,
                                      manifest=None,
                                      use_manifest=False)
        */
    }

    pub fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
        let processed_sigs = AtomicUsize::new(0);

        // TODO: Some(ref_sigs) case

        let search_sigs: Vec<_> = self
            .sig_files
            .internal_locations()
            .map(PathBuf::from)
            .collect();

        #[cfg(feature = "parallel")]
        let sig_iter = search_sigs.par_iter();

        #[cfg(not(feature = "parallel"))]
        let sig_iter = search_sigs.iter();

        let counters = sig_iter.enumerate().filter_map(|(dataset_id, filename)| {
            let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} reference sigs", i);
            }

            let search_sig = if let Some(storage) = &self.storage {
                let sig_data = storage
                    .load(
                        filename
                            .to_str()
                            .unwrap_or_else(|| panic!("error converting path {:?}", filename)),
                    )
                    .unwrap_or_else(|_| panic!("error loading {:?}", filename));

                Signature::from_reader(sig_data.as_slice())
            } else {
                Signature::from_path(&filename)
            }
            .unwrap_or_else(|_| panic!("Error processing {:?}", filename))
            .swap_remove(0);

            let mut search_mh = None;
            if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(&self.template) {
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
                counter[&(dataset_id as u64)] += size as usize;
                Some(counter)
            }
        });

        let reduce_counters = |mut a: SigCounter, b: SigCounter| {
            a.extend(&b);
            a
        };

        #[cfg(feature = "parallel")]
        let counter = counters.reduce(|| SigCounter::new(), reduce_counters);

        #[cfg(not(feature = "parallel"))]
        let counter = counters.fold(SigCounter::new(), reduce_counters);

        counter
    }

    pub fn search(
        &self,
        counter: SigCounter,
        similarity: bool,
        threshold: usize,
    ) -> Result<Vec<String>, Box<dyn std::error::Error>> {
        let mut matches = vec![];
        if similarity {
            unimplemented!("TODO: threshold correction")
        }

        for (dataset_id, size) in counter.most_common() {
            if size >= threshold {
                matches.push(
                    self.sig_files[dataset_id as usize]
                        .internal_location()
                        .to_str()
                        .unwrap()
                        .into(),
                );
            } else {
                break;
            };
        }
        Ok(matches)
    }

    fn gather_round(
        &self,
        dataset_id: u64,
        match_size: usize,
        query: &KmerMinHash,
        round: usize,
    ) -> Result<GatherResult, Error> {
        let match_path = if self.sig_files.is_empty() {
            PathBuf::new()
        } else {
            self.sig_files[dataset_id as usize].internal_location()
        };
        let match_sig = self.sig_for_dataset(dataset_id as usize)?;
        let result = self.stats_for_match(&match_sig, query, match_size, match_path, round)?;
        Ok(result)
    }

    fn sig_for_dataset(&self, dataset_id: usize) -> Result<SigStore, Error> {
        let match_path = if self.sig_files.is_empty() {
            PathBuf::new()
        } else {
            self.sig_files[dataset_id as usize].internal_location()
        };

        let match_sig = if let Some(refsigs) = &self.ref_sigs {
            refsigs[dataset_id as usize].clone()
        } else {
            let mut sig = if let Some(storage) = &self.storage {
                let sig_data = storage
                    .load(
                        match_path
                            .to_str()
                            .unwrap_or_else(|| panic!("error converting path {:?}", match_path)),
                    )
                    .unwrap_or_else(|_| panic!("error loading {:?}", match_path));
                Signature::from_reader(sig_data.as_slice())?
            } else {
                Signature::from_path(&match_path)?
            };
            // TODO: remove swap_remove
            sig.swap_remove(0).into()
        };
        Ok(match_sig)
    }

    fn stats_for_match(
        &self,
        match_sig: &Signature,
        query: &KmerMinHash,
        match_size: usize,
        match_path: PathBuf,
        gather_result_rank: usize,
    ) -> Result<GatherResult, Error> {
        let mut match_mh = None;
        if let Some(Sketch::MinHash(mh)) = match_sig.select_sketch(&self.template) {
            match_mh = Some(mh);
        }
        let match_mh = match_mh.expect("Couldn't find a compatible MinHash");

        // Calculate stats
        let f_orig_query = match_size as f64 / query.size() as f64;
        let f_match = match_size as f64 / match_mh.size() as f64;
        let filename = match_path.to_str().unwrap().into();
        let name = match_sig.name();
        let unique_intersect_bp = match_mh.scaled() as usize * match_size;

        let (intersect_orig, _) = match_mh.intersection_size(query)?;
        let intersect_bp = (match_mh.scaled() as u64 * intersect_orig) as usize;

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
    ) -> Result<Vec<GatherResult>, Box<dyn std::error::Error>> {
        let mut match_size = usize::max_value();
        let mut matches = vec![];

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
            let mut to_remove: HashSet<u64> = Default::default();
            to_remove.insert(dataset_id);

            for (dataset, value) in counter.iter_mut() {
                let dataset_sig = self.sig_for_dataset(*dataset as usize)?;
                let mut match_mh = None;
                if let Some(Sketch::MinHash(mh)) = dataset_sig.select_sketch(&self.template) {
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
        self.sig_files.clone()
    }

    pub fn set_manifest(&mut self, new_manifest: Manifest) -> Result<(), Error> {
        self.sig_files = new_manifest;
        Ok(())
    }

    pub fn signatures_iter(&self) -> impl Iterator<Item = SigStore> + '_ {
        if let Some(_sigs) = &self.ref_sigs {
            //sigs.iter().cloned()
            todo!("this works, but need to match return types")
        } else {
            // FIXME temp solution, must find better one!
            (0..self.sig_files.len())
                .map(move |dataset_id| self.sig_for_dataset(dataset_id).expect("error loading sig"))
        }
    }
}

impl<'a> Index<'a> for LinearRevIndex {
    type Item = SigStore;

    fn insert(&mut self, _node: Self::Item) -> Result<(), Error> {
        unimplemented!()
    }

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn len(&self) -> usize {
        if let Some(refs) = &self.ref_sigs {
            refs.len()
        } else {
            self.sig_files.len()
        }
    }

    fn signatures(&self) -> Vec<Self::Item> {
        if let Some(ref sigs) = self.ref_sigs {
            sigs.to_vec()
        } else {
            unimplemented!()
        }
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        unimplemented!()
    }
}

impl RevIndex {
    pub fn load<P: AsRef<Path>>(
        index_path: P,
        queries: Option<&[KmerMinHash]>,
    ) -> Result<RevIndex, Box<dyn std::error::Error>> {
        let (rdr, _) = niffler::from_path(index_path)?;
        let revindex = if let Some(qs) = queries {
            // TODO: avoid loading full revindex if query != None
            /*
            struct PartialRevIndex<T> {
                hashes_to_keep: Option<HashSet<HashIntoType>>,
                marker: PhantomData<fn() -> T>,
            }

            impl<T> PartialRevIndex<T> {
                pub fn new(hashes_to_keep: HashSet<u64>) -> Self {
                    PartialRevIndex {
                        hashes_to_keep: Some(hashes_to_keep),
                        marker: PhantomData,
                    }
                }
            }
            */

            let mut hashes: HashSet<u64> = HashSet::new();
            for q in qs {
                hashes.extend(q.iter_mins());
            }

            //let mut revindex: RevIndex = PartialRevIndex::new(hashes).deserialize(&rdr).unwrap();

            let mut revindex: RevIndex = serde_json::from_reader(rdr)?;
            revindex.hash_to_color.retain(&hashes);
            revindex
        } else {
            // Load the full revindex
            serde_json::from_reader(rdr)?
        };

        Ok(revindex)
    }

    pub fn new(
        search_sigs: &[PathBuf],
        template: &Sketch,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
        keep_sigs: bool,
    ) -> RevIndex {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        let linear = LinearRevIndex::new(Some(search_sigs.into()), template, keep_sigs, None, None);
        linear.index(threshold, merged_query, queries)
    }

    pub fn from_zipstorage(
        storage: ZipStorage,
        template: &Sketch,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
        keep_sigs: bool,
    ) -> Result<RevIndex, Error> {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        // Load manifest from zipstorage
        let manifest = Manifest::from_reader(storage.load("SOURMASH-MANIFEST.csv")?.as_slice())?;
        let search_sigs: Vec<_> = manifest.internal_locations().map(PathBuf::from).collect();

        let linear = LinearRevIndex::new(
            Some(search_sigs.as_slice().into()),
            template,
            keep_sigs,
            None,
            Some(storage),
        );

        Ok(linear.index(threshold, merged_query, queries))
    }

    fn merge_queries(qs: &[KmerMinHash], threshold: usize) -> Option<KmerMinHash> {
        if threshold == 0 {
            let mut merged = qs[0].clone();
            for query in &qs[1..] {
                merged.merge(query).unwrap();
            }
            Some(merged)
        } else {
            None
        }
    }

    pub fn new_with_sigs(
        search_sigs: Vec<Signature>,
        template: &Sketch,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
    ) -> RevIndex {
        // If threshold is zero, let's merge all queries and save time later
        let merged_query = queries.and_then(|qs| Self::merge_queries(qs, threshold));

        let linear = LinearRevIndex::new(
            Default::default(),
            template,
            false,
            search_sigs.into(),
            None,
        );

        linear.index(threshold, merged_query, queries)
    }

    fn map_hashes_colors(
        dataset_id: usize,
        search_sig: &Signature,
        queries: Option<&[KmerMinHash]>,
        merged_query: &Option<KmerMinHash>,
        threshold: usize,
        template: &Sketch,
    ) -> Option<(HashToColor, Colors)> {
        let mut search_mh = None;
        if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
            search_mh = Some(mh);
        }

        let search_mh = search_mh.expect("Couldn't find a compatible MinHash");
        let mut hash_to_color = HashToColor::new();
        let mut colors = Colors::default();

        if let Some(qs) = queries {
            if let Some(ref merged) = merged_query {
                let (matched_hashes, intersection) = merged.intersection(search_mh).unwrap();
                if !matched_hashes.is_empty() || intersection > threshold as u64 {
                    hash_to_color.add_to(&mut colors, dataset_id, matched_hashes);
                }
            } else {
                for query in qs {
                    let (matched_hashes, intersection) = query.intersection(search_mh).unwrap();
                    if !matched_hashes.is_empty() || intersection > threshold as u64 {
                        hash_to_color.add_to(&mut colors, dataset_id, matched_hashes);
                    }
                }
            }
        } else {
            let matched = search_mh.mins();
            let size = matched.len() as u64;
            if !matched.is_empty() || size > threshold as u64 {
                hash_to_color.add_to(&mut colors, dataset_id, matched);
            }
        };

        if hash_to_color.is_empty() {
            None
        } else {
            Some((hash_to_color, colors))
        }
    }

    pub fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
        query
            .iter_mins()
            .filter_map(|hash| self.hash_to_color.get(hash))
            .flat_map(|color| self.colors.indices(color))
            .cloned()
            .collect()
    }

    pub fn search(
        &self,
        counter: SigCounter,
        similarity: bool,
        threshold: usize,
    ) -> Result<Vec<String>, Box<dyn std::error::Error>> {
        self.linear.search(counter, similarity, threshold)
    }

    pub fn gather(
        &self,
        mut counter: SigCounter,
        threshold: usize,
        query: &KmerMinHash,
    ) -> Result<Vec<GatherResult>, Box<dyn std::error::Error>> {
        let mut match_size = usize::max_value();
        let mut matches = vec![];

        while match_size > threshold && !counter.is_empty() {
            let (dataset_id, size) = counter.most_common()[0];
            match_size = if size >= threshold { size } else { break };
            let result = self
                .linear
                .gather_round(dataset_id, match_size, query, matches.len())?;
            if let Some(Sketch::MinHash(match_mh)) =
                result.match_.select_sketch(&self.linear.template)
            {
                // Prepare counter for finding the next match by decrementing
                // all hashes found in the current match in other datasets
                for hash in match_mh.iter_mins() {
                    if let Some(color) = self.hash_to_color.get(hash) {
                        counter.subtract(self.colors.indices(color).cloned());
                    }
                }
                counter.remove(&dataset_id);
                matches.push(result);
            } else {
                unimplemented!()
            }
        }
        Ok(matches)
    }

    pub fn template(&self) -> Sketch {
        self.linear.template.clone()
    }

    // TODO: mh should be a sketch, or even a sig...
    pub(crate) fn find_signatures(
        &self,
        mh: &KmerMinHash,
        threshold: f64,
        containment: bool,
        _ignore_scaled: bool,
    ) -> Result<Vec<(f64, Signature, String)>, Error> {
        /*
        let template_mh = None;
        if let Sketch::MinHash(mh) = self.template {
            template_mh = Some(mh);
        };
        // TODO: throw error
        let template_mh = template_mh.unwrap();

        let tmp_mh;
        let mh = if template_mh.scaled() > mh.scaled() {
            // TODO: proper error here
            tmp_mh = mh.downsample_scaled(self.scaled)?;
            &tmp_mh
        } else {
            mh
        };

                if self.scaled < mh.scaled() && !ignore_scaled {
                        return Err(LcaDBError::ScaledMismatchError {
                                db: self.scaled,
                                query: mh.scaled(),
                        }
                        .into());
                }
        */

        // TODO: proper threshold calculation
        let threshold: usize = (threshold * (mh.size() as f64)) as _;

        let counter = self.counter_for_query(mh);

        debug!(
            "number of matching signatures for hashes: {}",
            counter.len()
        );

        let mut results = vec![];
        for (dataset_id, size) in counter.most_common() {
            let match_size = if size >= threshold { size } else { break };

            let match_path = if self.linear.sig_files.is_empty() {
                PathBuf::new()
            } else {
                self.linear.sig_files[dataset_id as usize].internal_location()
            };

            let ref_match;
            let match_sig = if let Some(refsigs) = &self.linear.ref_sigs {
                &refsigs[dataset_id as usize]
            } else {
                let mut sig = if let Some(storage) = &self.linear.storage {
                    let sig_data =
                        storage
                            .load(match_path.to_str().unwrap_or_else(|| {
                                panic!("error converting path {:?}", match_path)
                            }))
                            .unwrap_or_else(|_| panic!("error loading {:?}", match_path));
                    Signature::from_reader(sig_data.as_slice())?
                } else {
                    Signature::from_path(&match_path)?
                };
                // TODO: remove swap_remove
                ref_match = sig.swap_remove(0);
                &ref_match
            };

            let mut match_mh = None;
            if let Some(Sketch::MinHash(mh)) = match_sig.select_sketch(&self.linear.template) {
                match_mh = Some(mh);
            }
            let match_mh = match_mh.unwrap();

            if size >= threshold {
                let score = if containment {
                    size as f64 / mh.size() as f64
                } else {
                    size as f64 / (mh.size() + match_size - size) as f64
                };
                let filename = match_path.to_str().unwrap().into();
                let mut sig = match_sig.clone();
                sig.reset_sketches();
                sig.push(Sketch::MinHash(match_mh.clone()));
                results.push((score, sig, filename));
            } else {
                break;
            };
        }
        Ok(results)
    }
}

#[derive(CopyGetters, Getters, Setters, Serialize, Deserialize, Debug, PartialEq)]
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

    md5: String,
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

impl<'a> Index<'a> for RevIndex {
    type Item = Signature;

    fn insert(&mut self, _node: Self::Item) -> Result<(), Error> {
        unimplemented!()
    }

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn len(&self) -> usize {
        if let Some(refs) = &self.linear.ref_sigs {
            refs.len()
        } else {
            self.linear.sig_files.len()
        }
    }

    fn signatures(&self) -> Vec<Self::Item> {
        if let Some(ref sigs) = self.linear.ref_sigs {
            sigs.iter().map(|s| s.clone().into()).collect()
        } else {
            unimplemented!()
        }
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        unimplemented!()
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use crate::sketch::minhash::max_hash_for_scaled;

    #[test]
    fn revindex_new() {
        let max_hash = max_hash_for_scaled(10000);
        let template = Sketch::MinHash(
            KmerMinHash::builder()
                .num(0u32)
                .ksize(31)
                .max_hash(max_hash)
                .build(),
        );
        let search_sigs = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig".into(),
        ];
        let index = RevIndex::new(&search_sigs, &template, 0, None, false);
        assert_eq!(index.colors.len(), 3);
    }

    #[test]
    fn revindex_many() {
        let max_hash = max_hash_for_scaled(10000);
        let template = Sketch::MinHash(
            KmerMinHash::builder()
                .num(0u32)
                .ksize(31)
                .max_hash(max_hash)
                .build(),
        );
        let search_sigs = [
            "../../tests/test-data/gather/GCF_000006945.2_ASM694v2_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000007545.1_ASM754v1_genomic.fna.gz.sig".into(),
            "../../tests/test-data/gather/GCF_000008105.1_ASM810v1_genomic.fna.gz.sig".into(),
        ];

        let index = RevIndex::new(&search_sigs, &template, 0, None, false);
        /*
         dbg!(&index.colors.colors);
         0: 86
         1: 132
         2: 91
         (0, 1): 53
         (0, 2): 90
         (1, 2): 26
         (0, 1, 2): 261
         union: 739
        */
        //assert_eq!(index.colors.len(), 3);
        assert_eq!(index.colors.len(), 7);
    }

    #[test]
    fn revindex_from_zipstorage() {
        let max_hash = max_hash_for_scaled(100);
        let template = Sketch::MinHash(
            KmerMinHash::builder()
                .num(0u32)
                .ksize(57)
                .hash_function(crate::encodings::HashFunctions::murmur64_protein)
                .max_hash(max_hash)
                .build(),
        );
        let storage = ZipStorage::from_file("../../tests/test-data/prot/protein.zip")
            .expect("error loading zipfile");
        let index = RevIndex::from_zipstorage(storage, &template, 0, None, false)
            .expect("error building from ziptorage");

        assert_eq!(index.colors.len(), 3);

        let query_sig = Signature::from_path(
            "../../tests/test-data/prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig",
        )
        .expect("Error processing query")
        .swap_remove(0);
        let mut query_mh = None;
        if let Some(Sketch::MinHash(mh)) = query_sig.select_sketch(&template) {
            query_mh = Some(mh);
        }
        let query_mh = query_mh.expect("Couldn't find a compatible MinHash");

        let counter_rev = index.counter_for_query(query_mh);
        let counter_lin = index.linear.counter_for_query(query_mh);

        let results_rev = index.search(counter_rev, false, 0).unwrap();
        let results_linear = index.linear.search(counter_lin, false, 0).unwrap();
        assert_eq!(results_rev, results_linear);

        let counter_rev = index.counter_for_query(query_mh);
        let counter_lin = index.linear.counter_for_query(query_mh);

        let results_rev = index.gather(counter_rev, 0, query_mh).unwrap();
        let results_linear = index.linear.gather(counter_lin, 0, query_mh).unwrap();
        assert_eq!(results_rev.len(), 1);
        assert_eq!(results_rev, results_linear);
    }
}
