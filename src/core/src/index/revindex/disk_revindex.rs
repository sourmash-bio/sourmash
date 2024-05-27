use std::hash::{BuildHasher, BuildHasherDefault, Hash, Hasher};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;

use byteorder::{LittleEndian, WriteBytesExt};
use log::{info, trace};
use rayon::prelude::*;
use rocksdb::{ColumnFamilyDescriptor, MergeOperands, WriteBatchWithTransaction};

use crate::collection::{Collection, CollectionSet};
use crate::encodings::{Color, Idx};
use crate::index::revindex::{
    self as module, stats_for_cf, Datasets, DbStats, HashToColor, QueryColors, RevIndexOps, DB,
    HASHES, MANIFEST, METADATA, STORAGE_SPEC, VERSION,
};
use crate::index::{calculate_gather_stats, GatherResult, SigCounter};
use crate::manifest::Manifest;
use crate::prelude::*;
use crate::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use crate::storage::{InnerStorage, Storage};
use crate::Result;

const DB_VERSION: u8 = 1;
static SOURMASH_MEM_CACHE: &str = "SOURMASH_MEM_CACHE";

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
        //opts.prepare_for_bulk_load();

        let cache_size: usize = std::env::var(SOURMASH_MEM_CACHE)
            .unwrap_or_else(|_| "1".into())
            .parse()
            .unwrap();
        // in bytes, (1024 << 20 == 1GiB)
        let cache = rocksdb::Cache::new_lru_cache(cache_size * (1024 << 20));

        // prepare column family descriptors
        let cfs = cf_descriptors(cache.clone());

        let db = Arc::new(DB::open_cf_descriptors(&opts, path, cfs).unwrap());

        let processed_sigs = AtomicUsize::new(0);

        let index = Self {
            db,
            collection: Arc::new(collection),
        };

        index.save_collection().expect("Error saving collection");

        info!("Starting indexing");
        index
            .collection
            .par_iter()
            .map(|(dataset_id, _)| {
                let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                if i % 1000 == 0 && i > 0 {
                    info!("Processed {} reference sigs", i);
                }

                if i % 10000 == 0 && i > 0 {
                    info!("Triggering manual compaction");
                    index.compact();
                    info!("Finished manual compaction");
                }

                (
                    index.map_hashes_colors(dataset_id as Idx),
                    dataset_id as Idx,
                )
            })
            .chunks(10)
            .for_each(|chunk| {
                //let mut batch = WriteBatchWithTransaction::<false>::default();
                let cf_hashes = index.db.cf_handle(HASHES).unwrap();

                let mut hash_bytes = [0u8; 8];
                for (hashes, dataset_id) in chunk {
                    let colors = Datasets::new(&[dataset_id]).as_bytes().unwrap();

                    for hash in hashes {
                        (&mut hash_bytes[..])
                            .write_u64::<LittleEndian>(hash)
                            .expect("error writing bytes");
                        index
                            .db
                            .merge_cf(&cf_hashes, &hash_bytes[..], colors.as_slice())
                            .expect("error merging");
                        //batch.merge_cf(&cf_hashes, &hash_bytes[..], colors.as_slice());
                    }
                }
                //index.db.write(batch).expect("error merging batch"); // Atomically commits the batch
            });

        info!("Compact SSTs");
        index.compact();
        info!("Processed {} reference sigs", processed_sigs.into_inner());

        Ok(module::RevIndex::Plain(index))
    }

    pub fn open<P: AsRef<Path>>(
        path: P,
        read_only: bool,
        storage_spec: Option<&str>,
    ) -> Result<module::RevIndex> {
        let mut opts = module::RevIndex::db_options();
        if !read_only {
            //opts.prepare_for_bulk_load();
        }

        let cache_size: usize = std::env::var(SOURMASH_MEM_CACHE)
            .unwrap_or_else(|_| "1".into())
            .parse()
            .unwrap();
        // in bytes, (1024 << 20 == 1GiB)
        let cache = rocksdb::Cache::new_lru_cache(cache_size * (1024 << 20));

        // prepare column family descriptors
        let cfs = cf_descriptors(cache);

        let db = if read_only {
            Arc::new(DB::open_cf_descriptors(&opts, path.as_ref(), cfs)?)
        } else {
            Arc::new(DB::open_cf_descriptors(&opts, path.as_ref(), cfs)?)
        };

        let collection = Arc::new(Self::load_collection_from_rocksdb(
            db.clone(),
            storage_spec,
        )?);

        Ok(module::RevIndex::Plain(Self { db, collection }))
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

    fn map_hashes_colors(&self, dataset_id: Idx) -> impl Iterator<Item = crate::HashIntoType> {
        let search_sig = self
            .collection
            .sig_for_dataset(dataset_id)
            .expect("Couldn't find a compatible Signature");
        let mh = search_sig.minhash().expect("Error extracting a minhash");

        mh.mins().into_iter()
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
        let mut query = KmerMinHashBTree::from(orig_query.clone());
        let mut sum_weighted_found = 0;
        let _selection = selection.unwrap_or_else(|| self.collection.selection());
        let mut orig_query_ds = orig_query.clone();
        let total_weighted_hashes = orig_query.sum_abunds();

        // or set this with user --track-abundance?
        let calc_abund_stats = orig_query.track_abundance();

        // todo: let user pass these options in
        let calc_ani_ci = false;
        let ani_confidence_interval_fraction = None;

        while match_size > threshold && !counter.is_empty() {
            trace!("counter len: {}", counter.len());
            trace!("match size: {}", match_size);

            let (dataset_id, size) = counter.k_most_common_ordered(1)[0];
            match_size = if size >= threshold { size } else { break };
            // handle special case where threshold was set to 0
            if match_size == 0 {
                break;
            }

            // this should downsample mh for us
            let match_sig = self.collection.sig_for_dataset(dataset_id)?;

            // get downsampled minhashes for comparison.
            let match_mh = match_sig.minhash().unwrap().clone();
            query = query.downsample_scaled(match_mh.scaled())?;
            orig_query_ds = orig_query_ds.downsample_scaled(match_mh.scaled())?;

            // just calculate essentials here
            let gather_result_rank = matches.len();

            let query_mh = KmerMinHash::from(query.clone());

            // grab the specific intersection:
            let isect = match_mh.intersection(&query_mh)?;
            let mut isect_mh = match_mh.clone();
            isect_mh.clear();
            isect_mh.add_many(&isect.0)?;

            // Calculate stats
            let gather_result = calculate_gather_stats(
                &orig_query_ds,
                KmerMinHash::from(query.clone()),
                match_sig,
                match_size,
                gather_result_rank,
                sum_weighted_found,
                total_weighted_hashes.try_into().unwrap(),
                calc_abund_stats,
                calc_ani_ci,
                ani_confidence_interval_fraction,
            )?;
            // keep track of the sum weighted found
            sum_weighted_found = gather_result.sum_weighted_found();
            matches.push(gather_result);

            trace!("Preparing counter for next round");
            // Prepare counter for finding the next match by decrementing
            // all hashes found in the current match in other datasets
            // TODO: not used at the moment, so just skip.
            query.remove_many(match_mh.iter_mins().copied())?; // is there a better way?

            // TODO: Use HashesToColors here instead. If not initialized,
            //       build it.
            isect
                .0
                .iter()
                .filter_map(|hash| hash_to_color.get(hash))
                .flat_map(|color| {
                    // TODO: remove this clone
                    query_colors.get(color).unwrap().clone().into_iter()
                })
                .for_each(|dataset| {
                    // TODO: collect the flat_map into a Counter, and remove more
                    //       than one at a time...
                    counter.entry(dataset).and_modify(|e| *e -= 1);
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

        self.save_collection().expect("Error saving collection");

        self.collection
            .par_iter()
            .skip(to_skip)
            .map(|(dataset_id, _)| {
                let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                if i % 1000 == 0 {
                    info!("Processed {} reference sigs", i);
                }

                (self.map_hashes_colors(dataset_id as Idx), dataset_id as Idx)
            })
            .chunks(10)
            .for_each(|chunk| {
                let mut batch = WriteBatchWithTransaction::<true>::default();
                let cf_hashes = self.db.cf_handle(HASHES).unwrap();

                let mut hash_bytes = [0u8; 8];
                for (hashes, dataset_id) in chunk {
                    let colors = Datasets::new(&[dataset_id]).as_bytes().unwrap();

                    for hash in hashes {
                        (&mut hash_bytes[..])
                            .write_u64::<LittleEndian>(hash)
                            .expect("error writing bytes");
                        batch.merge_cf(&cf_hashes, &hash_bytes[..], colors.as_slice());
                    }
                }
                self.db.write(batch).expect("error merging batch"); // Atomically commits the batch
            });

        info!("Compact SSTs");
        self.compact();

        info!(
            "Processed additional {} reference sigs",
            processed_sigs.into_inner()
        );

        Ok(module::RevIndex::Plain(self))
    }

    fn collection(&self) -> &CollectionSet {
        &self.collection
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

fn cf_descriptors(cache: rocksdb::Cache) -> Vec<ColumnFamilyDescriptor> {
    let mut cfopts = module::RevIndex::db_options();

    // following https://rocksdb.org/blog/2021/05/26/integrated-blob-db.html
    cfopts.set_enable_blob_files(true);
    // If empty or one dataset, avoid saving to blob store
    cfopts.set_min_blob_size(8);
    // TODO: set blob file size to write_buffer_size
    //cfopts.set_blob_file_size(cfopts.write_bufffer_size());
    cfopts.set_blob_file_size(0x4000000); // 64 MiB
    cfopts.set_enable_blob_gc(true);
    //
    cfopts.set_blob_compression_type(rocksdb::DBCompressionType::Zstd);

    cfopts.set_max_write_buffer_number(16);
    cfopts.set_merge_operator_associative("datasets operator", merge_datasets);
    cfopts.set_min_write_buffer_number_to_merge(10);

    // Updated default from
    // https://github.com/facebook/rocksdb/wiki/Setup-Options-and-Basic-Tuning#other-general-options
    cfopts.set_level_compaction_dynamic_level_bytes(true);
    //cfopts.prepare_for_bulk_load();

    let mut tfopts = rocksdb::BlockBasedOptions::default();
    //tfopts.set_index_type(rocksdb::BlockBasedIndexType::TwoLevelIndexSearch);
    tfopts.set_block_cache(&cache);
    tfopts.set_optimize_filters_for_memory(true);
    //tfopts.set_data_block_index_type(rocksdb::DataBlockIndexType::BinaryAndHash);
    // Keys for HASHES are HashIntoType, a u64
    //tfopts.set_hybrid_ribbon_filter(64.0, 2);

    // these are from db_options, not sure if overwritten if not here
    tfopts.set_block_size(0x4000000);
    tfopts.set_cache_index_and_filter_blocks(true);
    tfopts.set_pin_l0_filter_and_index_blocks_in_cache(true);
    tfopts.set_format_version(6);

    cfopts.set_block_based_table_factory(&tfopts);
    // Keys for HASHES are HashIntoType, a u64, and so 8 bytes
    //cfopts.set_prefix_extractor(rocksdb::SliceTransform::create_fixed_prefix(8));

    // 10GB for memory budget
    cfopts.optimize_level_style_compaction(10 * 1024 * 1024 * 1024);

    let cf_hashes = ColumnFamilyDescriptor::new(HASHES, cfopts);

    let mut cfopts = module::RevIndex::db_options();
    cfopts.set_max_write_buffer_number(16);
    // Updated default
    cfopts.set_level_compaction_dynamic_level_bytes(true);
    //cfopts.set_merge_operator_associative("colors operator", merge_colors);

    let cf_metadata = ColumnFamilyDescriptor::new(METADATA, cfopts);

    vec![cf_hashes, cf_metadata]
}
