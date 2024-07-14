use std::sync::Arc;

use crate::storage::{Storage, StorageArgs, StorageError};
use crate::{Error, Result};

// Column family for using rocksdb as a Storage
pub(crate) const STORAGE: &str = "storage";

/// Store data in RocksDB
#[derive(Debug, Clone)]
pub struct RocksDBStorage {
    db: Arc<crate::index::revindex::DB>,
}

impl RocksDBStorage {
    pub fn from_path(path: &str) -> Self {
        let mut opts = crate::index::revindex::RevIndex::db_options();
        opts.create_if_missing(true);
        opts.create_missing_column_families(true);
        opts.prepare_for_bulk_load();

        // prepare column family descriptors
        let cfs = crate::index::revindex::disk_revindex::cf_descriptors();

        let db =
            Arc::new(crate::index::revindex::DB::open_cf_descriptors(&opts, path, cfs).unwrap());

        Self { db }
    }

    pub fn from_db(db: Arc<crate::index::revindex::DB>) -> Self {
        Self { db: db.clone() }
    }
}

impl Storage for RocksDBStorage {
    fn save(&self, path: &str, content: &[u8]) -> Result<String> {
        let cf_storage = self.db.cf_handle(STORAGE).unwrap();
        // TODO(lirber): deal with conflict for path?
        self.db.put_cf(&cf_storage, path.as_bytes(), &content[..])?;
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
        "rocksdb://".into()
    }
}
