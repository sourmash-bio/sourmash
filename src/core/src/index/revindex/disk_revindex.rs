use std::cmp::max;
use std::collections::HashSet;
use std::hash::{BuildHasher, BuildHasherDefault};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::{Arc, RwLock};

use byteorder::{LittleEndian, WriteBytesExt};
use log::{info, trace};
use rayon::prelude::*;
use rocksdb::MergeOperands;

use crate::collection::{Collection, CollectionSet};
use crate::encodings::{Color, Idx};
use crate::index::revindex::{
    self as module, stats_for_cf, CounterGather, DatasetPicklist, Datasets, DbStats, QueryColors,
    RevIndexOps,
};
use crate::index::{calculate_gather_stats, GatherResult, SigCounter};
use crate::manifest::Manifest;
use crate::prelude::*;
use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use crate::sketch::Sketch;
use crate::storage::{
    rocksdb::{cf_descriptors, db_options, ALL_CFS, DB, HASHES, METADATA},
    InnerStorage, RocksDBStorage, Storage,
};
use crate::Result;

const DB_VERSION: u8 = 1;

// DB metadata saved in the METADATA column family
const MANIFEST: &str = "manifest";
const STORAGE_SPEC: &str = "storage_spec";
const VERSION: &str = "version";
const PROCESSED: &str = "processed";

fn compute_color(idxs: &Datasets) -> Color {
    let s = BuildHasherDefault::<twox_hash::Xxh3Hash128>::default();
    s.hash_one(idxs)
}

#[derive(Clone)]
pub struct DiskRevIndex {
    location: String,
    db: Arc<DB>,
    collection: Arc<CollectionSet>,
    processed: Arc<RwLock<Datasets>>,
}

pub(crate) fn merge_datasets(
    _: &[u8],
    existing_val: Option<&[u8]>,
    operands: &MergeOperands,
) -> Option<Vec<u8>> {
    let mut datasets = match existing_val {
        Some(val) => Datasets::from_slice(val).expect("cannot unpack slice"),
        None => Default::default(),
    };

    for op in operands {
        let new_vals = Datasets::from_slice(op).unwrap();
        datasets.union(new_vals);
    }
    // TODO: optimization! if nothing changed, skip as_bytes()
    datasets.as_bytes()
}

/* TODO: need the repair_cf variant, not available in rocksdb-rust yet
pub fn repair(path: &Path) {
    let opts = db_options();

    DB::repair(&opts, path).unwrap()
}
*/

impl DiskRevIndex {
    pub fn create(path: &Path, collection: CollectionSet) -> Result<module::RevIndex> {
        let mut opts = db_options();
        opts.create_if_missing(true);
        opts.create_missing_column_families(true);

        // prepare column family descriptors
        let cfs = cf_descriptors();

        let db = Arc::new(DB::open_cf_descriptors(&opts, path, cfs).unwrap());

        let processed_sigs = AtomicUsize::new(0);

        let collection = Arc::new(collection);
        let processed = Arc::new(RwLock::new(Self::load_processed(
            db.clone(),
            collection.clone(),
            true,
        )?));

        let index = Self {
            location: String::from(path.to_str().expect("cannot extract path")),
            db,
            collection,
            processed: processed.clone(),
        };

        index.collection.par_iter().for_each(|(dataset_id, _)| {
            // check if this dataset_id was processed already
            // call map_hashes_colors only if not already processed
            if !processed.read().unwrap().contains(&dataset_id) {
                let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                if i % 1000 == 0 {
                    info!("Processed {} reference sigs", i);
                }

                index.map_hashes_colors(dataset_id as Idx);

                // if cached in a new field in the RevIndex,
                // then update the cache too

                processed.write().unwrap().extend([dataset_id]);
            }
        });

        index.save_collection().expect("Error saving collection");

        info!("Compact SSTs");
        index.compact();
        info!(
            "Done! Processed {} reference sigs",
            processed_sigs.into_inner()
        );

        Ok(module::RevIndex::Disk(index))
    }

    pub fn open<P: AsRef<Path>>(
        path: P,
        read_only: bool,
        storage_spec: Option<&str>,
    ) -> Result<module::RevIndex> {
        let mut opts = db_options();
        opts.create_if_missing(true);
        opts.create_missing_column_families(true);

        // prepare column family descriptors
        let cfs = cf_descriptors();

        let db = if read_only {
            Arc::new(DB::open_cf_descriptors_read_only(
                &opts,
                path.as_ref(),
                cfs,
                false,
            )?)
        } else {
            Arc::new(DB::open_cf_descriptors(&opts, path.as_ref(), cfs)?)
        };

        let collection = Arc::new(Self::load_collection_from_rocksdb(
            db.clone(),
            storage_spec,
        )?);

        let processed = Arc::new(RwLock::new(Self::load_processed(
            db.clone(),
            collection.clone(),
            false,
        )?));

        Ok(module::RevIndex::Disk(Self {
            location: String::from(path.as_ref().to_str().expect("cannot extract path")),
            db,
            collection,
            processed,
        }))
    }

    /// Access to the DB is unsafe because RocksDB allows writing with
    /// a regular Arc<DB> handle; it doesn't have to be mut.
    pub unsafe fn db(&self) -> Arc<DB> {
        self.db.clone()
    }

    fn load_processed(
        db: Arc<DB>,
        collection: Arc<CollectionSet>,
        assume_empty: bool,
    ) -> Result<Datasets> {
        let cf_metadata = db.cf_handle(METADATA).unwrap();
        if let Some(rdr) = db.get_pinned_cf(&cf_metadata, PROCESSED)? {
            // convert rdr to Datasets
            Datasets::from_slice(&rdr)
        } else if assume_empty {
            Ok(Datasets::default())
        } else {
            let all_datasets: Vec<_> = (0..collection.manifest().len()).map(|v| v as Idx).collect();
            Ok(Datasets::new(&all_datasets))
        }
    }

    fn load_collection_from_rocksdb(
        db: Arc<DB>,
        storage_spec: Option<&str>,
    ) -> Result<CollectionSet> {
        let cf_metadata = db.cf_handle(METADATA).unwrap();

        let rdr = db.get_cf(&cf_metadata, VERSION)?.unwrap();
        assert_eq!(rdr[0], DB_VERSION);

        let rdr = db.get_cf(&cf_metadata, MANIFEST)?.unwrap();
        let manifest = Manifest::from_reader(&rdr[..])?;

        let spec = match storage_spec {
            Some(spec) => spec.into(),
            None => {
                let db_spec = db.get_cf(&cf_metadata, STORAGE_SPEC)?;
                String::from_utf8(db_spec.unwrap()).map_err(|e| e.utf8_error())?
            }
        };

        let storage = if spec == "rocksdb://" {
            InnerStorage::new(RocksDBStorage::from_db(db.clone()))
        } else {
            InnerStorage::from_spec(spec)?
        };

        Collection::new(manifest, storage).try_into()
    }

    fn save_collection(&self) -> Result<()> {
        let cf_metadata = self.db.cf_handle(METADATA).unwrap();

        // save DB version
        // TODO: probably should go together with a more general
        //       saving procedure used in create/update
        self.db.put_cf(&cf_metadata, VERSION, [DB_VERSION])?;

        // write manifest
        let mut wtr = vec![];
        {
            self.collection.manifest().to_writer(&mut wtr)?;
        }
        self.db.put_cf(&cf_metadata, MANIFEST, &wtr[..])?;

        // write storage spec
        let spec = self.collection.storage().spec();

        // TODO: check if spec if memstorage, would probably have to
        // save into rocksdb in that case!

        self.db.put_cf(&cf_metadata, STORAGE_SPEC, spec)?;

        Ok(())
    }

    fn map_hashes_colors(&self, dataset_id: Idx) {
        let search_sig = self
            .collection
            .sig_for_dataset(dataset_id)
            .expect("Couldn't find a compatible Signature");
        let search_mh = &search_sig.sketches()[0];

        let colors = Datasets::new(&[dataset_id]).as_bytes().unwrap();

        let cf_hashes = self.db.cf_handle(HASHES).unwrap();

        let hashes = match search_mh {
            Sketch::MinHash(mh) => mh.mins(),
            Sketch::LargeMinHash(mh) => mh.mins(),
            _ => unimplemented!(),
        };

        let mut hash_bytes = [0u8; 8];
        for hash in hashes {
            (&mut hash_bytes[..])
                .write_u64::<LittleEndian>(hash)
                .expect("error writing bytes");
            self.db
                .merge_cf(&cf_hashes, &hash_bytes[..], colors.as_slice())
                .expect("error merging");
        }

        // finished processing this dataset,
        // do a merge_cf in the PROCESSED key in metadata
        // to account for that.
        let cf_metadata = self.db.cf_handle(METADATA).unwrap();
        self.db
            .merge_cf(&cf_metadata, PROCESSED, colors.as_slice())
            .expect("error merging");
    }
}

impl RevIndexOps for DiskRevIndex {
    fn location(&self) -> &str {
        self.location.as_str()
    }

    fn counter_for_query(
        &self,
        query: &KmerMinHash,
        picklist: Option<DatasetPicklist>,
    ) -> SigCounter {
        info!("Collecting hashes");
        let cf_hashes = self.db.cf_handle(HASHES).unwrap();
        let hashes_iter = query.iter_mins().map(|hash| {
            let mut v = vec![0_u8; 8];
            (&mut v[..])
                .write_u64::<LittleEndian>(*hash)
                .expect("error writing bytes");
            (&cf_hashes, v)
        });

        info!("Multi get");
        self.db
            .multi_get_cf(hashes_iter)
            .into_iter()
            .filter_map(|r| r.ok().unwrap_or(None))
            .flat_map(|raw_datasets| {
                let new_vals = Datasets::from_slice(&raw_datasets).unwrap();

                // filter against picklist if need be.
                if let Some(pl) = &picklist {
                    let new_vals: HashSet<_> = new_vals
                        .into_iter()
                        .filter(|&i| pl.dataset_ids.contains(&i))
                        .collect();
                    Box::new(new_vals.into_iter())
                } else {
                    new_vals.into_iter()
                }
            })
            .collect()
    }

    fn prepare_gather_counters(
        &self,
        query: &KmerMinHash,
        picklist: Option<DatasetPicklist>,
    ) -> CounterGather {
        let cf_hashes = self.db.cf_handle(HASHES).unwrap();
        let hashes_iter = query.iter_mins().map(|hash| {
            let mut v = vec![0_u8; 8];
            (&mut v[..])
                .write_u64::<LittleEndian>(*hash)
                .expect("error writing bytes");
            (&cf_hashes, v)
        });

        /*
         build a HashToColors for query,
         and a QueryColors (Color -> Datasets) mapping.
         Loading Datasets from rocksdb for every hash takes too long.
        */
        let mut query_colors: QueryColors = Default::default();
        let mut counter: SigCounter = Default::default();

        info!("Building hash_to_color and query_colors");
        let hash_to_color = query
            .iter_mins()
            .zip(self.db.multi_get_cf(hashes_iter))
            .filter_map(|(k, r)| {
                let raw: Option<Vec<u8>> = r.ok().unwrap_or(None);

                if let Some(r) = raw {
                    let mut new_vals = Datasets::from_slice(&r).unwrap();

                    // filter by picklist?
                    if let Some(pl) = &picklist {
                        let val_set: Vec<Idx> = new_vals
                            .into_iter()
                            .filter(|&i| pl.dataset_ids.contains(&i))
                            .collect();
                        new_vals = Datasets::new(&val_set[..]);
                    }

                    if new_vals.len() > 0 {
                        let color = compute_color(&new_vals);
                        query_colors
                            .entry(color)
                            .or_insert_with(|| new_vals.clone());
                        counter.update(new_vals);
                        Some((*k, color))
                    } else {
                        None
                    }
                } else {
                    None
                }
            })
            .collect();

        CounterGather {
            counter,
            query_colors,
            hash_to_color,
        }
    }

    fn gather(
        &self,
        mut cg: CounterGather,
        threshold: usize,
        orig_query: &KmerMinHash,
        selection: Option<Selection>,
    ) -> Result<Vec<GatherResult>> {
        let match_size = usize::MAX;
        let mut matches = vec![];
        let mut query = KmerMinHashBTree::from(orig_query.clone());
        let mut sum_weighted_found = 0;
        let _selection = selection.unwrap_or_else(|| self.collection.selection());
        let total_weighted_hashes = orig_query.sum_abunds();

        // or set this with user --track-abundance?
        let calc_abund_stats = orig_query.track_abundance();

        // todo: let user pass these options in
        let calc_ani_ci = false;
        let ani_confidence_interval_fraction = None;

        while match_size > threshold && !cg.is_empty() {
            trace!("counter len: {}", cg.len());
            trace!("match size: {}", match_size);

            let result = cg.peek(threshold);
            if result.is_none() {
                break;
            }
            let (dataset_id, match_size) = result.unwrap();

            let match_sig = self.collection.sig_for_dataset(dataset_id)?;
            let match_mh = match_sig.minhash().unwrap().clone();

            // make downsampled minhashes
            let max_scaled = max(match_mh.scaled(), query.scaled());

            let match_mh = match_mh
                .downsample_scaled(max_scaled)
                .expect("cannot downsample match");

            // repeatedly downsample query, then extract to KmerMinHash
            // => calculate_gather_stats
            query = query
                .downsample_scaled(max_scaled)
                .expect("cannot downsample query");
            let query_mh = KmerMinHash::from(query.clone());

            // just calculate essentials here
            let gather_result_rank = matches.len() as u32;

            // grab the specific intersection:
            // Calculate stats
            let (gather_result, isect) = calculate_gather_stats(
                orig_query,
                query_mh,
                match_sig,
                match_size,
                gather_result_rank,
                sum_weighted_found,
                total_weighted_hashes,
                calc_abund_stats,
                calc_ani_ci,
                ani_confidence_interval_fraction,
            )
            .expect("could not calculate gather stats");

            // use intersection from calc_gather_stats to make a KmerMinHash.
            let mut isect_mh = match_mh.clone();
            isect_mh.clear();
            isect_mh.add_many(&isect.0)?;

            // keep track of the sum weighted found
            sum_weighted_found = gather_result.sum_weighted_found();
            matches.push(gather_result);

            trace!("Preparing counter for next round");
            // Prepare counter for finding the next match by decrementing
            // all hashes found in the current match in other datasets
            // TODO: not used at the moment, so just skip.
            query.remove_many(isect_mh.iter_mins().copied())?; // is there a better way?

            cg.consume(&isect_mh);
        }
        Ok(matches)
    }

    fn update(mut self, collection: CollectionSet) -> Result<module::RevIndex> {
        // TODO: verify new collection manifest is a superset of current one,
        //       and the initial chunk is the same
        self.collection.check_superset(&collection)?;
        info!("sigs in the original index: {}", self.collection.len());

        self.collection = Arc::new(collection);
        info!(
            "sigs in the new index once finished: {}",
            self.collection.len()
        );

        let processed = self.processed.clone();
        info!(
            "sigs left to process: {}",
            self.collection.len() - processed.read().unwrap().len()
        );

        // process the remainder
        let processed_sigs = AtomicUsize::new(0);

        self.collection.par_iter().for_each(|(dataset_id, _)| {
            // check if this dataset_id was processed already
            // call map_hashes_colors only if not already processed
            if !processed.read().unwrap().contains(&dataset_id) {
                let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                if i % 1000 == 0 {
                    info!("Processed {} reference sigs", i);
                }

                self.map_hashes_colors(dataset_id as Idx);

                // if cached in a new field in the RevIndex,
                // then update the cache too

                processed.write().unwrap().extend([dataset_id]);
            }
        });

        self.save_collection().expect("Error saving collection");

        info!("Compact SSTs");
        self.compact();

        info!(
            "Processed additional {} reference sigs",
            processed_sigs.into_inner()
        );

        Ok(module::RevIndex::Disk(self))
    }

    fn check(&self, quick: bool) -> DbStats {
        stats_for_cf(self.db.clone(), HASHES, true, quick)
    }

    fn compact(&self) {
        for cf_name in ALL_CFS {
            let cf = self.db.cf_handle(cf_name).unwrap();
            self.db.compact_range_cf(&cf, None::<&[u8]>, None::<&[u8]>)
        }
    }

    fn flush(&self) -> Result<()> {
        self.db.flush_wal(true)?;

        for cf_name in [HASHES, METADATA] {
            let cf = self.db.cf_handle(cf_name).unwrap();
            self.db.flush_cf(&cf)?;
        }

        Ok(())
    }

    fn collection(&self) -> &CollectionSet {
        &self.collection
    }

    fn internalize_storage(&mut self) -> Result<()> {
        // check if collection is already internal, if so return
        if self.collection.storage().spec() == "rocksdb://" {
            return Ok(());
        }

        // build new rocksdb storage from db
        let new_storage = RocksDBStorage::from_db(self.db.clone());

        // use manifest to copy from current storage to new one
        self.collection()
            .par_iter()
            .try_for_each(|(_, record)| -> Result<()> {
                let path = record.internal_location().as_str();
                let sig_data = self.collection.storage().load(path).unwrap();
                new_storage.save(path, &sig_data)?;
                Ok(())
            })?;

        // Replace storage for collection.
        // Using unchecked version because we just used the manifest
        // above to make sure the storage is still consistent
        unsafe {
            if let Some(v) = Arc::get_mut(&mut self.collection) {
                v.set_storage_unchecked(InnerStorage::new(new_storage))
            }
        }

        // write storage spec
        let cf_metadata = self.db.cf_handle(METADATA).unwrap();
        let spec = "rocksdb://";
        self.db.put_cf(&cf_metadata, STORAGE_SPEC, spec)?;

        Ok(())
    }

    fn convert(&self, _output_db: module::RevIndex) -> Result<()> {
        todo!()
        /*
        if let RevIndex::Color(db) = output_db {
            let other_db = db.db;

            let cf_hashes = self.db.cf_handle(HASHES).unwrap();

            info!("start converting colors");
            let mut color_bytes = [0u8; 8];
            let iter = self
                .db
                .iterator_cf(&cf_hashes, rocksdb::IteratorMode::Start);
            for (key, value) in iter {
                let datasets = Datasets::from_slice(&value).unwrap();
                let new_idx: Vec<_> = datasets.into_iter().collect();
                let new_color = Colors::update(other_db.clone(), None, new_idx.as_slice()).unwrap();

                (&mut color_bytes[..])
                    .write_u64::<LittleEndian>(new_color)
                    .expect("error writing bytes");
                other_db
                    .put_cf(&cf_hashes, &key[..], &color_bytes[..])
                    .unwrap();
            }
            info!("finished converting colors");

            info!("copying sigs to output");
            let cf_sigs = self.db.cf_handle(SIGS).unwrap();
            let iter = self.db.iterator_cf(&cf_sigs, rocksdb::IteratorMode::Start);
            for (key, value) in iter {
                other_db.put_cf(&cf_sigs, &key[..], &value[..]).unwrap();
            }
            info!("finished copying sigs to output");

            Ok(())
        } else {
            todo!()
        }
        */
    }
    fn find_signatures(
        &self,
        query_mh: &KmerMinHash,
        threshold: f64,
        picklist: Option<DatasetPicklist>,
    ) -> Result<Vec<(f64, Signature, String)>> {
        // do search
        let counter = self.counter_for_query(query_mh, picklist);

        // retrieve/convert matches. I don't think there's a simple way to
        // truncate this without going through all the matches, so it's
        // potentially (much) more expensive than prefetch.
        let filename = self.location();
        let results: Vec<(f64, Signature, String)> = counter
            .most_common()
            .into_iter()
            .filter_map(|(dataset_id, _size)| {
                let sig: Signature = self
                    .collection()
                    .sig_for_dataset(dataset_id)
                    .expect("dataset not found")
                    .into();

                let match_mh = sig.minhash().expect("cannot retrieve match");

                let f_match = if match_mh.scaled() != query_mh.scaled() {
                    let match_ds = match_mh
                        .clone()
                        .downsample_scaled(query_mh.scaled())
                        .expect("cannot downsample");
                    query_mh
                        .jaccard(&match_ds)
                        .expect("cannot calculate Jaccard")
                } else {
                    query_mh
                        .jaccard(match_mh)
                        .expect("cannot calculate Jaccard")
                };

                if f_match >= threshold {
                    Some((f_match, sig, filename.to_owned()))
                } else {
                    None
                }
            })
            .collect();

        Ok(results)
    }
}
