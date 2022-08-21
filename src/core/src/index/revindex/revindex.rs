use std::hash::{BuildHasher, BuildHasherDefault, Hash, Hasher};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use byteorder::{LittleEndian, WriteBytesExt};
use log::{info, trace};
use rayon::prelude::*;
use rocksdb::{ColumnFamilyDescriptor, MergeOperands, Options};

use crate::index::revindex::mem_revindex::GatherResult;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;

use crate::index::revindex::prepare_query;
use crate::index::revindex::{
    self as module, sig_save_to_db, stats_for_cf, Color, DatasetID, Datasets, HashToColor,
    QueryColors, SigCounter, SignatureData, DB, HASHES, SIGS,
};

fn compute_color(idxs: &Datasets) -> Color {
    let s = BuildHasherDefault::<twox_hash::Xxh3Hash128>::default();
    let mut hasher = s.build_hasher();
    /*
    // TODO: remove this...
    let mut sorted: Vec<_> = idxs.iter().collect();
    sorted.sort();
    */
    idxs.hash(&mut hasher);
    hasher.finish()
}

#[derive(Debug, Clone)]
pub struct RevIndex {
    db: Arc<DB>,
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
    pub fn create(path: &Path) -> module::RevIndex {
        let mut opts = module::RevIndex::db_options();
        opts.create_if_missing(true);
        opts.create_missing_column_families(true);

        // prepare column family descriptors
        let cfs = cf_descriptors();

        let db = Arc::new(DB::open_cf_descriptors(&opts, path, cfs).unwrap());

        module::RevIndex::Plain(Self { db })
    }

    pub fn open(path: &Path, read_only: bool) -> module::RevIndex {
        let opts = module::RevIndex::db_options();

        // prepare column family descriptors
        let cfs = cf_descriptors();

        let db = if read_only {
            Arc::new(DB::open_cf_descriptors_read_only(&opts, path, cfs, false).unwrap())
        } else {
            Arc::new(DB::open_cf_descriptors(&opts, path, cfs).unwrap())
        };

        module::RevIndex::Plain(Self { db })
    }

    fn map_hashes_colors(
        &self,
        dataset_id: DatasetID,
        filename: &PathBuf,
        threshold: f64,
        template: &Sketch,
        save_paths: bool,
    ) {
        let search_sig = Signature::from_path(&filename)
            .unwrap_or_else(|_| panic!("Error processing {:?}", filename))
            .swap_remove(0);

        let search_mh =
            prepare_query(&search_sig, template).expect("Couldn't find a compatible MinHash");

        let colors = Datasets::new(&[dataset_id]).as_bytes().unwrap();

        let cf_hashes = self.db.cf_handle(HASHES).unwrap();

        let matched = search_mh.mins();
        let size = matched.len() as u64;
        if !matched.is_empty() || size > threshold as u64 {
            // FIXME threshold is f64
            let mut hash_bytes = [0u8; 8];
            for hash in matched {
                (&mut hash_bytes[..])
                    .write_u64::<LittleEndian>(hash)
                    .expect("error writing bytes");
                self.db
                    .merge_cf(&cf_hashes, &hash_bytes[..], colors.as_slice())
                    .expect("error merging");
            }
        }

        sig_save_to_db(
            self.db.clone(),
            search_sig,
            search_mh,
            size,
            threshold,
            save_paths,
            filename,
            dataset_id,
        );
    }

    pub fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
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

    pub fn prepare_gather_counters(
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
            .zip(self.db.multi_get_cf(hashes_iter).into_iter())
            .filter_map(|(k, r)| {
                let raw = r.ok().unwrap_or(None);
                raw.map(|raw| {
                    let new_vals = Datasets::from_slice(&raw).unwrap();
                    let color = compute_color(&new_vals);
                    query_colors
                        .entry(color)
                        .or_insert_with(|| new_vals.clone());
                    counter.update(new_vals.into_iter());
                    (*k, color)
                })
            })
            .collect();

        (counter, query_colors, hash_to_colors)
    }

    pub fn matches_from_counter(
        &self,
        counter: SigCounter,
        threshold: usize,
    ) -> Vec<(String, usize)> {
        let cf_sigs = self.db.cf_handle(SIGS).unwrap();

        let matches_iter = counter
            .most_common()
            .into_iter()
            .filter_map(|(dataset_id, size)| {
                if size >= threshold {
                    let mut v = vec![0_u8; 8];
                    (&mut v[..])
                        .write_u64::<LittleEndian>(dataset_id)
                        .expect("error writing bytes");
                    Some((&cf_sigs, v, size))
                } else {
                    None
                }
            });

        let matches_sizes = matches_iter.clone().map(|(_, _, v)| v);

        info!("Multi get matches");
        self.db
            .multi_get_cf(matches_iter.map(|(k, v, _)| (k, v)))
            .into_iter()
            .zip(matches_sizes)
            .filter_map(|(r, size)| r.ok().unwrap_or(None).map(|v| (v, size)))
            .filter_map(
                |(sigdata, size)| match SignatureData::from_slice(&sigdata).unwrap() {
                    SignatureData::Empty => None,
                    SignatureData::External(p) => Some((p, size)),
                    SignatureData::Internal(sig) => Some((sig.name(), size)),
                },
            )
            .collect()
    }

    pub fn gather(
        &self,
        mut counter: SigCounter,
        query_colors: QueryColors,
        hash_to_color: HashToColor,
        threshold: usize,
        orig_query: &KmerMinHash,
        template: &Sketch,
    ) -> Result<Vec<GatherResult>, Box<dyn std::error::Error>> {
        let mut match_size = usize::max_value();
        let mut matches = vec![];
        let mut key_bytes = [0u8; 8];
        //let mut query: KmerMinHashBTree = orig_query.clone().into();

        let cf_sigs = self.db.cf_handle(SIGS).unwrap();

        while match_size > threshold && !counter.is_empty() {
            trace!("counter len: {}", counter.len());
            trace!("match size: {}", match_size);

            let (dataset_id, size) = counter.k_most_common_ordered(1)[0];
            match_size = if size >= threshold { size } else { break };

            (&mut key_bytes[..])
                .write_u64::<LittleEndian>(dataset_id)
                .expect("error writing bytes");

            let match_sig = self
                .db
                .get_cf(&cf_sigs, &key_bytes[..])
                .ok()
                .map(
                    |sigdata| match SignatureData::from_slice(&(sigdata.unwrap())).unwrap() {
                        SignatureData::Empty => todo!("throw error, empty sig"),
                        SignatureData::External(_p) => todo!("Load from external"),
                        SignatureData::Internal(sig) => sig,
                    },
                )
                .unwrap_or_else(|| panic!("Unknown dataset {}", dataset_id));

            let match_mh =
                prepare_query(&match_sig, template).expect("Couldn't find a compatible MinHash");

            // Calculate stats
            let f_orig_query = match_size as f64 / orig_query.size() as f64;
            let f_match = match_size as f64 / match_mh.size() as f64;
            let name = match_sig.name();
            let unique_intersect_bp = match_mh.scaled() as usize * match_size;
            let gather_result_rank = matches.len();

            let (intersect_orig, _) = match_mh.intersection_size(orig_query)?;
            let intersect_bp = (match_mh.scaled() as u64 * intersect_orig) as usize;

            let f_unique_to_query = intersect_orig as f64 / orig_query.size() as f64;
            let match_ = match_sig.clone();
            let md5 = match_sig.md5sum();

            // TODO: all of these
            let filename = "".into();
            let f_unique_weighted = 0.;
            let average_abund = 0;
            let median_abund = 0;
            let std_abund = 0;
            let f_match_orig = 0.;
            let remaining_bp = 0;

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
                .match_(match_)
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

    pub fn index(
        &self,
        index_sigs: Vec<PathBuf>,
        template: &Sketch,
        threshold: f64,
        save_paths: bool,
    ) {
        let processed_sigs = AtomicUsize::new(0);

        index_sigs
            .par_iter()
            .enumerate()
            .for_each(|(dataset_id, filename)| {
                let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                if i % 1000 == 0 {
                    info!("Processed {} reference sigs", i);
                }

                self.map_hashes_colors(
                    dataset_id as DatasetID,
                    filename,
                    threshold,
                    template,
                    save_paths,
                );
            });
        info!("Processed {} reference sigs", processed_sigs.into_inner());
    }

    pub fn check(&self, quick: bool) {
        stats_for_cf(self.db.clone(), HASHES, true, quick);
        info!("");
        stats_for_cf(self.db.clone(), SIGS, false, quick);
    }

    pub fn compact(&self) {
        for cf_name in [HASHES, SIGS] {
            let cf = self.db.cf_handle(cf_name).unwrap();
            self.db.compact_range_cf(&cf, None::<&[u8]>, None::<&[u8]>)
        }
    }

    pub fn flush(&self) -> Result<(), Box<dyn std::error::Error>> {
        self.db.flush_wal(true)?;

        for cf_name in [HASHES, SIGS] {
            let cf = self.db.cf_handle(cf_name).unwrap();
            self.db.flush_cf(&cf)?;
        }

        Ok(())
    }

    pub fn convert(&self, output_db: module::RevIndex) -> Result<(), Box<dyn std::error::Error>> {
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

    let cf_sigs = ColumnFamilyDescriptor::new(SIGS, cfopts);

    vec![cf_hashes, cf_sigs]
}
