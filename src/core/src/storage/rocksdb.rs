use std::sync::Arc;

use rocksdb::{ColumnFamilyDescriptor, Options};

use crate::storage::{Storage, StorageArgs, StorageError};
use crate::Result;

// Column families
pub(crate) const HASHES: &str = "hashes";
pub(crate) const COLORS: &str = "colors";
pub(crate) const METADATA: &str = "metadata";

// Column family for using rocksdb as a Storage
pub(crate) const STORAGE: &str = "storage";

pub(crate) const ALL_CFS: [&str; 3] = [HASHES, METADATA, STORAGE];

pub type DB = rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>;

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

        // prepare column family descriptors
        let cfs = cf_descriptors();

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

pub(crate) fn cf_descriptors() -> Vec<ColumnFamilyDescriptor> {
    let mut cfopts = Options::default();
    cfopts.set_max_write_buffer_number(16);
    cfopts.set_merge_operator_associative(
        "datasets operator",
        crate::index::revindex::disk_revindex::merge_datasets,
    );
    cfopts.set_min_write_buffer_number_to_merge(10);

    // Updated default from
    // https://github.com/facebook/rocksdb/wiki/Setup-Options-and-Basic-Tuning#other-general-options
    cfopts.set_level_compaction_dynamic_level_bytes(true);

    let cf_hashes = ColumnFamilyDescriptor::new(HASHES, cfopts);

    let mut cfopts = Options::default();
    cfopts.set_max_write_buffer_number(16);
    cfopts.set_merge_operator_associative(
        "datasets operator",
        crate::index::revindex::disk_revindex::merge_datasets,
    );
    // Updated default
    cfopts.set_level_compaction_dynamic_level_bytes(true);

    let cf_metadata = ColumnFamilyDescriptor::new(METADATA, cfopts);

    let mut cfopts = Options::default();
    cfopts.set_max_write_buffer_number(16);
    // Updated default
    cfopts.set_level_compaction_dynamic_level_bytes(true);

    let cf_storage = ColumnFamilyDescriptor::new(STORAGE, cfopts);

    let mut cfopts = Options::default();
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
    block_opts.set_cache_index_and_filter_blocks(true);
    block_opts.set_pin_l0_filter_and_index_blocks_in_cache(true);
    block_opts.set_format_version(6);
    opts.set_block_based_table_factory(&block_opts);
    // End of updated defaults

    opts.increase_parallelism(rayon::current_num_threads() as i32);
    //opts.max_background_jobs = 6;
    // opts.optimize_level_style_compaction();
    // opts.optimize_universal_style_compaction();

    opts
}
