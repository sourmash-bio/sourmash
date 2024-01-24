use std::hash::{BuildHasher, BuildHasherDefault, Hash, Hasher};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use byteorder::{LittleEndian, WriteBytesExt};
use log::{info, trace};
use rayon::prelude::*;
use rocksdb::{ColumnFamilyDescriptor, MergeOperands, Options};

use crate::collection::{Collection, CollectionSet};
use crate::encodings::{Color, Idx};
use crate::index::revindex::{
    self as module, prepare_query, stats_for_cf, Datasets, DbStats, HashToColor, QueryColors,
    RevIndexOps, DB, HASHES, MANIFEST, METADATA, STORAGE_SPEC, VERSION,
};
use crate::index::{GatherResult, SigCounter};
use crate::manifest::Manifest;
use crate::prelude::*;
use crate::signature::SigsTrait;
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
use crate::storage::{InnerStorage, Storage};
use crate::Result;

use statistics::median;
use mean::mean;

const DB_VERSION: u8 = 1;

fn compute_color(idxs: &Datasets) -> Color {
    let s = BuildHasherDefault::<twox_hash::Xxh3Hash128>::default();
    let mut hasher = s.build_hasher();
    idxs.hash(&mut hasher);
    hasher.finish()
}

#[derive(Clone)]
pub struct RevIndex {
    db: Arc<DB>,
    collection: Arc<CollectionSet>,
}

fn merge_datasets(
    _: &[u8],
    existing_val: Option<&[u8]>,
    operands: &MergeOperands,
) -> Option<Vec<u8>> {
    let mut datasets = existing_val
        .and_then(Datasets::from_slice)
        .unwrap_or_default();

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

impl RevIndex {
    pub fn create(path: &Path, collection: CollectionSet) -> Result<module::RevIndex> {
        let mut opts = module::RevIndex::db_options();
        opts.create_if_missing(true);
        opts.create_missing_column_families(true);
        opts.prepare_for_bulk_load();

        // prepare column family descriptors
        let cfs = cf_descriptors();

        let db = Arc::new(DB::open_cf_descriptors(&opts, path, cfs).unwrap());

        let processed_sigs = AtomicUsize::new(0);

        let index = Self {
            db,
            collection: Arc::new(collection),
        };

        index.collection.par_iter().for_each(|(dataset_id, _)| {
            let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
            if i % 1000 == 0 {
                info!("Processed {} reference sigs", i);
            }

            index.map_hashes_colors(dataset_id as Idx);
        });

        index.save_collection().expect("Error saving collection");

        info!("Compact SSTs");
        index.compact();
        info!("Processed {} reference sigs", processed_sigs.into_inner());

        Ok(module::RevIndex::Plain(index))
    }

    pub fn open<P: AsRef<Path>>(path: P, read_only: bool) -> Result<module::RevIndex> {
        let mut opts = module::RevIndex::db_options();
        if !read_only {
            opts.prepare_for_bulk_load();
        }

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

        let collection = Arc::new(Self::load_collection_from_rocksdb(db.clone())?);

        Ok(module::RevIndex::Plain(Self { db, collection }))
    }

    fn load_collection_from_rocksdb(db: Arc<DB>) -> Result<CollectionSet> {
        let cf_metadata = db.cf_handle(METADATA).unwrap();

        let rdr = db.get_cf(&cf_metadata, VERSION)?.unwrap();
        assert_eq!(rdr[0], DB_VERSION);

        let rdr = db.get_cf(&cf_metadata, MANIFEST)?.unwrap();
        let manifest = Manifest::from_reader(&rdr[..])?;

        let spec = String::from_utf8(db.get_cf(&cf_metadata, STORAGE_SPEC)?.unwrap())
            .expect("invalid utf-8");

        let storage = if spec == "rocksdb://" {
            todo!("init storage from db")
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
        self.db.put_cf(&cf_metadata, VERSION, &[DB_VERSION])?;

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
    }
}

impl RevIndexOps for RevIndex {
    fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
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
                new_vals.into_iter()
            })
            .collect()
    }

    fn prepare_gather_counters(
        &self,
        query: &KmerMinHash,
    ) -> (SigCounter, QueryColors, HashToColor) {
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

        info!("Building hash_to_colors and query_colors");
        let hash_to_colors = query
            .iter_mins()
            .zip(self.db.multi_get_cf(hashes_iter))
            .filter_map(|(k, r)| {
                let raw = r.ok().unwrap_or(None);
                raw.map(|raw| {
                    let new_vals = Datasets::from_slice(&raw).unwrap();
                    let color = compute_color(&new_vals);
                    query_colors
                        .entry(color)
                        .or_insert_with(|| new_vals.clone());
                    counter.update(new_vals);
                    (*k, color)
                })
            })
            .collect();

        (counter, query_colors, hash_to_colors)
    }

    fn matches_from_counter(&self, counter: SigCounter, threshold: usize) -> Vec<(String, usize)> {
        info!("get matches from counter");
        counter
            .most_common()
            .into_iter()
            .filter_map(|(dataset_id, size)| {
                if size >= threshold {
                    let row = &self
                        .collection
                        .record_for_dataset(dataset_id)
                        .expect("dataset not found");
                    Some((row.name().into(), size))
                } else {
                    None
                }
            })
            .collect()
    }

    fn gather(
        &self,
        mut counter: SigCounter,
        query_colors: QueryColors,
        hash_to_color: HashToColor,
        threshold: usize,
        orig_query: &KmerMinHash,
        selection: Option<Selection>,
    ) -> Result<Vec<GatherResult>> {
        let mut match_size = usize::max_value();
        let mut matches = vec![];
        //let mut query: KmerMinHashBTree = orig_query.clone().into();
        let selection = selection.unwrap_or_else(|| self.collection.selection());
        let mut remaining_bp = orig_query.scaled() as usize * orig_query.size();

        while match_size > threshold && !counter.is_empty() {
            trace!("counter len: {}", counter.len());
            trace!("match size: {}", match_size);

            let (dataset_id, size) = counter.k_most_common_ordered(1)[0];
            match_size = if size >= threshold { size } else { break };

            let match_sig = self.collection.sig_for_dataset(dataset_id)?;

            // Calculate stats
            let name = match_sig.name();
            let gather_result_rank = matches.len();
            let match_ = match_sig.clone();
            let md5 = match_sig.md5sum();

            let match_mh = prepare_query(match_sig.into(), &selection)
                .expect("Couldn't find a compatible MinHash");
            let f_match = match_size as f64 / match_mh.size() as f64;
            let unique_intersect_bp = match_mh.scaled() as usize * match_size;
            let (intersect_orig, _) = match_mh.intersection_size(orig_query)?;
            let intersect_bp = (match_mh.scaled() * intersect_orig) as usize;

            let f_unique_to_query = match_size as f64 / orig_query.size() as f64;
            let f_orig_query = intersect_orig as f64 / orig_query.size() as f64;
            let f_match_orig = intersect_orig as f64 / match_mh.size() as f64;

            let filename = match_sig.filename();
            // weight common by query abundances
            (common_abunds, total_common_weighted) = query.weighted_intersect_size(match_mh, abunds_from=orig_query);
            let total_orig_query_abund = orig_query.abunds().iter().sum::<u64>();
            //Calculate abund-related metrics
            let f_unique_weighted = total_common_weighted as f64 / total_orig_query_abund as f64;
            // mean, median, std of abundances
            let average_abund: f64 = mean(&common_abunds.iter().map(|&x| x as f64).collect::<Vec<f64>>());
            let median_abund: f64 = median(&common_abunds.iter().map(|&x| x as f64).collect::<Vec<f64>>());
            let std_abund: f64 = statistics::std_dev(&common_abunds.iter().map(|&x| x as f64).collect::<Vec<f64>>(), None)?;

            // remaining_bp is the number of base pairs in the query that are not in the match (or any prior match)
            remaining_bp = remaining_bp - unique_intersect_bp;

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
                .match_(match_.into())
                .f_match_orig(f_match_orig)
                .unique_intersect_bp(unique_intersect_bp)
                .gather_result_rank(gather_result_rank)
                .remaining_bp(remaining_bp)
                .build();
            matches.push(result);

            trace!("Preparing counter for next round");
            // Prepare counter for finding the next match by decrementing
            // all hashes found in the current match in other datasets
            // TODO: not used at the moment, so just skip.
            //query.remove_many(match_mh.to_vec().as_slice())?;

            // TODO: Use HashesToColors here instead. If not initialized,
            //       build it.
            match_mh
                .iter_mins()
                .filter_map(|hash| hash_to_color.get(hash))
                .flat_map(|color| {
                    // TODO: remove this clone
                    query_colors.get(color).unwrap().clone().into_iter()
                })
                .for_each(|dataset| {
                    // TODO: collect the flat_map into a Counter, and remove more
                    //       than one at a time...
                    counter.entry(dataset).and_modify(|e| {
                        if *e > 0 {
                            *e -= 1
                        }
                    });
                });

            counter.remove(&dataset_id);
        }
        Ok(matches)
    }

    fn update(mut self, collection: CollectionSet) -> Result<module::RevIndex> {
        // TODO: verify new collection manifest is a superset of current one,
        //       and the initial chunk is the same
        let to_skip = self.collection.check_superset(&collection)?;

        // process the remainder
        let processed_sigs = AtomicUsize::new(0);

        self.collection = Arc::new(collection);

        self.collection
            .par_iter()
            .skip(to_skip)
            .for_each(|(dataset_id, _)| {
                let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                if i % 1000 == 0 {
                    info!("Processed {} reference sigs", i);
                }

                self.map_hashes_colors(dataset_id as Idx);
            });

        self.save_collection().expect("Error saving collection");

        info!("Compact SSTs");
        self.compact();

        info!(
            "Processed additional {} reference sigs",
            processed_sigs.into_inner()
        );

        Ok(module::RevIndex::Plain(self))
    }

    fn check(&self, quick: bool) -> DbStats {
        stats_for_cf(self.db.clone(), HASHES, true, quick)
    }

    fn compact(&self) {
        for cf_name in [HASHES, METADATA] {
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
}

fn cf_descriptors() -> Vec<ColumnFamilyDescriptor> {
    let mut cfopts = Options::default();
    cfopts.set_max_write_buffer_number(16);
    cfopts.set_merge_operator_associative("datasets operator", merge_datasets);
    cfopts.set_min_write_buffer_number_to_merge(10);

    // Updated default from
    // https://github.com/facebook/rocksdb/wiki/Setup-Options-and-Basic-Tuning#other-general-options
    cfopts.set_level_compaction_dynamic_level_bytes(true);

    let cf_hashes = ColumnFamilyDescriptor::new(HASHES, cfopts);

    let mut cfopts = Options::default();
    cfopts.set_max_write_buffer_number(16);
    // Updated default
    cfopts.set_level_compaction_dynamic_level_bytes(true);
    //cfopts.set_merge_operator_associative("colors operator", merge_colors);

    let cf_metadata = ColumnFamilyDescriptor::new(METADATA, cfopts);

    let mut cfopts = Options::default();
    cfopts.set_max_write_buffer_number(16);
    // Updated default
    cfopts.set_level_compaction_dynamic_level_bytes(true);
    //cfopts.set_merge_operator_associative("colors operator", merge_colors);

    vec![cf_hashes, cf_metadata]
}
