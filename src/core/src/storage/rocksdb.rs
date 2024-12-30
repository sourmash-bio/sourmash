use std::sync::Arc;

use rocksdb::ColumnFamilyDescriptor;

use crate::storage::{Storage, StorageArgs, StorageError};
use crate::Result;

// Column families
pub(crate) const HASHES: &str = "hashes";
pub(crate) const COLORS: &str = "colors";
pub(crate) const METADATA: &str = "metadata";

// Column family for using rocksdb as a Storage
pub(crate) const STORAGE: &str = "storage";

pub(crate) const ALL_CFS: [&str; 3] = [HASHES, METADATA, STORAGE];

// Env var for controlling cache size
pub(crate) const SOURMASH_MEM_CACHE: &str = "SOURMASH_MEM_CACHE";

pub type DB = rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>;
//pub type DB = rocksdb::OptimisticTransactionDB<rocksdb::MultiThreaded>;

/// Store data in RocksDB
#[derive(Debug, Clone)]
pub struct RocksDBStorage {
    db: Arc<DB>,
}

impl RocksDBStorage {
    pub fn from_path(path: &str) -> Self {
        let mut opts = db_options();
        opts.create_if_missing(true);
        opts.create_missing_column_families(true);
        opts.prepare_for_bulk_load();

        let cache_size: usize = std::env::var(SOURMASH_MEM_CACHE)
            .unwrap_or_else(|_| "1".into())
            .parse()
            .unwrap();
        // in bytes, (1024 << 20 == 1GiB)
        let cache = rocksdb::Cache::new_lru_cache(cache_size * (1024 << 20));

        // prepare column family descriptors
        let cfs = cf_descriptors(cache.clone());

        let db = Arc::new(DB::open_cf_descriptors(&opts, path, cfs).unwrap());

        Self { db }
    }

    pub fn from_db(db: Arc<DB>) -> Self {
        Self { db: db.clone() }
    }
}

impl Storage for RocksDBStorage {
    fn save(&self, path: &str, content: &[u8]) -> Result<String> {
        let cf_storage = self.db.cf_handle(STORAGE).unwrap();
        // TODO(lirber): deal with conflict for path?
        self.db.put_cf(&cf_storage, path.as_bytes(), content)?;
        Ok(path.into())
    }

    fn load(&self, path: &str) -> Result<Vec<u8>> {
        let cf_storage = self.db.cf_handle(STORAGE).unwrap();
        let data = self.db.get_cf(&cf_storage, path.as_bytes())?;
        data.ok_or_else(|| StorageError::DataReadError(path.into()).into())
    }

    fn args(&self) -> StorageArgs {
        unimplemented!()
    }

    fn spec(&self) -> String {
        format!("rocksdb://{}", self.db.path().display())
    }
}

pub(crate) fn cf_descriptors(cache: rocksdb::Cache) -> Vec<ColumnFamilyDescriptor> {
    let mut cfopts = db_options();

    /*
    // following https://rocksdb.org/blog/2021/05/26/integrated-blob-db.html
    cfopts.set_enable_blob_files(true);
    // If empty or one dataset, avoid saving to blob store
    cfopts.set_min_blob_size(8);
    // TODO: set blob file size to write_buffer_size
    //cfopts.set_blob_file_size(cfopts.write_bufffer_size());
    cfopts.set_blob_file_size(0x4000000); // 64 MiB
    cfopts.set_enable_blob_gc(true);
    cfopts.set_blob_compression_type(rocksdb::DBCompressionType::Zstd);
    */

    cfopts.set_max_write_buffer_number(16);
    cfopts.set_merge_operator_associative(
        "datasets operator",
        crate::index::revindex::disk_revindex::merge_datasets,
    );
    cfopts.set_min_write_buffer_number_to_merge(10);

    // Updated default from
    // https://github.com/facebook/rocksdb/wiki/Setup-Options-and-Basic-Tuning#other-general-options
    cfopts.set_level_compaction_dynamic_level_bytes(true);

    let mut tfopts = rocksdb::BlockBasedOptions::default();
    //tfopts.set_index_type(rocksdb::BlockBasedIndexType::TwoLevelIndexSearch);
    tfopts.set_block_cache(&cache);
    tfopts.set_optimize_filters_for_memory(true);
    //tfopts.set_data_block_index_type(rocksdb::DataBlockIndexType::BinaryAndHash);
    // Keys for HASHES are HashIntoType, a u64
    //tfopts.set_hybrid_ribbon_filter(64.0, 2);

    // these are from db_options, not sure if overwritten if not here
    //tfopts.set_block_size(0x4000000); // 64 MiB
    tfopts.set_block_size(16 * 1024);
    tfopts.set_cache_index_and_filter_blocks(true);
    tfopts.set_pin_l0_filter_and_index_blocks_in_cache(true);
    tfopts.set_format_version(6);

    cfopts.set_block_based_table_factory(&tfopts);
    // Keys for HASHES are HashIntoType, a u64, and so 8 bytes
    //cfopts.set_prefix_extractor(rocksdb::SliceTransform::create_fixed_prefix(8));

    // 10GB for memory budget
    //cfopts.optimize_level_style_compaction(10 * 1024 * 1024 * 1024);

    let cf_hashes = ColumnFamilyDescriptor::new(HASHES, cfopts);

    let mut cfopts = db_options();
    cfopts.set_max_write_buffer_number(16);
    cfopts.set_merge_operator_associative(
        "datasets operator",
        crate::index::revindex::disk_revindex::merge_datasets,
    );
    // Updated default
    cfopts.set_level_compaction_dynamic_level_bytes(true);

    let cf_metadata = ColumnFamilyDescriptor::new(METADATA, cfopts);

    let mut cfopts = db_options();
    cfopts.set_max_write_buffer_number(16);
    // Updated default
    cfopts.set_level_compaction_dynamic_level_bytes(true);

    let cf_storage = ColumnFamilyDescriptor::new(STORAGE, cfopts);

    let mut cfopts = db_options();
    cfopts.set_max_write_buffer_number(16);
    // Updated default
    cfopts.set_level_compaction_dynamic_level_bytes(true);

    vec![cf_hashes, cf_metadata, cf_storage]
}

pub(crate) fn db_options() -> rocksdb::Options {
    let mut opts = rocksdb::Options::default();
    opts.set_max_open_files(500);

    // Updated defaults from
    // https://github.com/facebook/rocksdb/wiki/Setup-Options-and-Basic-Tuning#other-general-options
    opts.set_bytes_per_sync(1048576);
    let mut block_opts = rocksdb::BlockBasedOptions::default();
    block_opts.set_block_size(16 * 1024);
    //block_opts.set_block_size(0x4000000); // 64 MiB
    block_opts.set_cache_index_and_filter_blocks(true);
    block_opts.set_pin_l0_filter_and_index_blocks_in_cache(true);
    block_opts.set_format_version(6);
    opts.set_block_based_table_factory(&block_opts);
    // End of updated defaults

    opts.increase_parallelism(rayon::current_num_threads() as i32);
    opts.set_max_background_jobs(rayon::current_num_threads() as i32);
    opts.set_stats_dump_period_sec(300);
    // opts.optimize_level_style_compaction();
    // opts.optimize_universal_style_compaction();

    opts.set_bottommost_compression_type(rocksdb::DBCompressionType::Zstd);
    opts.set_bottommost_zstd_max_train_bytes(1024 << 10, true); // 1MiB

    opts
}
