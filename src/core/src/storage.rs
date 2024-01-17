use std::collections::{BTreeMap, HashMap};
use std::ffi::OsStr;
use std::fs::{DirBuilder, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::ops::Deref;
use std::sync::{Arc, RwLock};

use camino::Utf8Path as Path;
use camino::Utf8PathBuf as PathBuf;
use once_cell::sync::OnceCell;
use serde::{Deserialize, Serialize};
use thiserror::Error;
use typed_builder::TypedBuilder;

use crate::errors::ReadDataError;
use crate::prelude::*;
use crate::signature::SigsTrait;
use crate::sketch::Sketch;
use crate::{Error, Result};

/// An abstraction for any place where we can store data.
pub trait Storage {
    /// Save bytes into path
    fn save(&self, path: &str, content: &[u8]) -> Result<String>;

    /// Load bytes from path
    fn load(&self, path: &str) -> Result<Vec<u8>>;

    /// Args for initializing a new Storage
    fn args(&self) -> StorageArgs;

    /// Load signature from internal path
    fn load_sig(&self, path: &str) -> Result<SigStore>;

    /// Return a spec for creating/opening a storage
    fn spec(&self) -> String;

    /// Save signature to internal path
    fn save_sig(&self, path: &str, sig: Signature) -> Result<String> {
        let mut buffer = vec![];
        {
            sig.to_writer(&mut buffer).unwrap();
        }
        self.save(path, &buffer)
    }
}

#[derive(Debug, Error)]
pub enum StorageError {
    #[error("Path can't be empty")]
    EmptyPathError,

    #[error("Path not found: {0}")]
    PathNotFoundError(String),

    #[error("Error reading data from {0}")]
    DataReadError(String),
}

#[derive(Clone)]
pub struct InnerStorage(Arc<RwLock<dyn Storage + Send + Sync + 'static>>);

#[derive(TypedBuilder, Default, Clone)]
pub struct SigStore {
    #[builder(setter(into))]
    filename: String,

    #[builder(setter(into))]
    name: String,

    #[builder(setter(into))]
    metadata: String,

    storage: Option<InnerStorage>,

    #[builder(setter(into), default)]
    data: OnceCell<Signature>,
}

#[derive(Serialize, Deserialize)]
pub(crate) struct StorageInfo {
    pub backend: String,
    pub args: StorageArgs,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(untagged)]
pub enum StorageArgs {
    FSStorage { path: String },
}

/// Store files locally into a directory
#[derive(TypedBuilder, Debug, Clone, Default)]
pub struct FSStorage {
    /// absolute path for the directory where data is saved.
    fullpath: PathBuf,
    subdir: String,
}

#[ouroboros::self_referencing]
pub struct ZipStorage {
    mapping: Option<memmap2::Mmap>,

    #[borrows(mapping)]
    #[covariant]
    archive: piz::ZipArchive<'this>,

    subdir: Option<String>,
    path: Option<PathBuf>,

    #[borrows(archive)]
    #[covariant]
    metadata: Metadata<'this>,
}

/// Store data in memory (no permanent storage)
#[derive(TypedBuilder, Debug, Clone, Default)]
pub struct MemStorage {
    //store: HashMap<String, Vec<u8>>,
    sigs: Arc<RwLock<HashMap<String, SigStore>>>,
}

pub type Metadata<'a> = BTreeMap<&'a OsStr, &'a piz::read::FileMetadata<'a>>;

// =========================================

impl InnerStorage {
    pub fn new(inner: impl Storage + Send + Sync + 'static) -> InnerStorage {
        InnerStorage(Arc::new(RwLock::new(inner)))
    }

    pub fn from_spec(spec: String) -> Result<Self> {
        Ok(match spec {
            x if x.starts_with("fs") => {
                let path = x.split("://").last().expect("not a valid path");
                InnerStorage::new(FSStorage::new("", path))
            }
            x if x.starts_with("memory") => InnerStorage::new(MemStorage::new()),
            x if x.starts_with("zip") => {
                let path = x.split("://").last().expect("not a valid path");
                InnerStorage::new(ZipStorage::from_file(path)?)
            }
            _ => todo!("storage not supported, throw error"),
        })
    }
}

impl Storage for InnerStorage {
    fn save(&self, path: &str, content: &[u8]) -> Result<String> {
        self.0.save(path, content)
    }

    fn load(&self, path: &str) -> Result<Vec<u8>> {
        self.0.load(path)
    }

    fn args(&self) -> StorageArgs {
        self.0.args()
    }

    fn load_sig(&self, path: &str) -> Result<SigStore> {
        let mut store = self.0.load_sig(path)?;
        store.storage = Some(self.clone());
        Ok(store)
    }

    fn spec(&self) -> String {
        self.0.spec()
    }
}

impl From<&StorageArgs> for FSStorage {
    fn from(other: &StorageArgs) -> FSStorage {
        match other {
            StorageArgs::FSStorage { path } => {
                let mut fullpath = PathBuf::new();
                fullpath.push(".");
                fullpath.push(path);

                FSStorage {
                    fullpath,
                    subdir: path.clone(),
                }
            }
        }
    }
}

impl<L> Storage for RwLock<L>
where
    L: ?Sized + Storage,
{
    fn save(&self, path: &str, content: &[u8]) -> Result<String> {
        self.read().unwrap().save(path, content)
    }

    fn load(&self, path: &str) -> Result<Vec<u8>> {
        self.read().unwrap().load(path)
    }

    fn args(&self) -> StorageArgs {
        self.read().unwrap().args()
    }

    fn load_sig(&self, path: &str) -> Result<SigStore> {
        self.read().unwrap().load_sig(path)
    }

    fn spec(&self) -> String {
        self.read().unwrap().spec()
    }
}

impl FSStorage {
    pub fn new(location: &str, subdir: &str) -> FSStorage {
        let mut fullpath = PathBuf::new();
        fullpath.push(location);
        fullpath.push(subdir);

        FSStorage {
            fullpath,
            subdir: subdir.into(),
        }
    }

    pub fn set_base(&mut self, location: &str) {
        let mut fullpath = PathBuf::new();
        fullpath.push(location);
        fullpath.push(&self.subdir);
        self.fullpath = fullpath;
    }
}

impl Storage for FSStorage {
    fn save(&self, path: &str, content: &[u8]) -> Result<String> {
        if path.is_empty() {
            return Err(StorageError::EmptyPathError.into());
        }

        let fpath = self.fullpath.join(path);
        DirBuilder::new()
            .recursive(true)
            .create(fpath.parent().unwrap())?;

        let file = File::create(&fpath)?;
        let mut buf_writer = BufWriter::new(file);
        buf_writer.write_all(content)?;
        Ok(path.into())
    }

    fn load(&self, path: &str) -> Result<Vec<u8>> {
        let path = self.fullpath.join(path);
        let file = File::open(path)?;
        let mut buf_reader = BufReader::new(file);
        let mut contents = Vec::new();
        buf_reader.read_to_end(&mut contents)?;
        Ok(contents)
    }

    fn args(&self) -> StorageArgs {
        StorageArgs::FSStorage {
            path: self.subdir.clone(),
        }
    }

    fn load_sig(&self, path: &str) -> Result<SigStore> {
        let raw = self.load(path)?;
        let sig = Signature::from_reader(&mut &raw[..])?
            // TODO: select the right sig?
            .swap_remove(0);

        Ok(sig.into())
    }

    fn spec(&self) -> String {
        format!("fs://{}", self.subdir)
    }
}

fn lookup<'a, P: AsRef<Path>>(
    metadata: &'a Metadata,
    path: P,
) -> Result<&'a piz::read::FileMetadata<'a>> {
    let path = path.as_ref();
    metadata
        .get(&path.as_os_str())
        .ok_or_else(|| StorageError::PathNotFoundError(path.to_string()).into())
        .map(|entry| *entry)
}

fn find_subdirs<'a>(archive: &'a piz::ZipArchive<'a>) -> Result<Option<String>> {
    let subdirs: Vec<_> = archive
        .entries()
        .iter()
        .filter(|entry| entry.is_dir())
        .collect();
    if subdirs.len() == 1 {
        Ok(Some(subdirs[0].path.as_str().into()))
    } else {
        Ok(None)
    }
}

impl Storage for ZipStorage {
    fn save(&self, _path: &str, _content: &[u8]) -> Result<String> {
        unimplemented!();
    }

    fn load(&self, path: &str) -> Result<Vec<u8>> {
        let metadata = self.borrow_metadata();

        let entry = lookup(metadata, path).or_else(|_| {
            if let Some(subdir) = self.borrow_subdir() {
                lookup(metadata, subdir.to_owned() + path)
                    .map_err(|_| StorageError::PathNotFoundError(path.into()))
            } else {
                Err(StorageError::PathNotFoundError(path.into()))
            }
        })?;

        let mut reader = BufReader::new(
            self.borrow_archive()
                .read(entry)
                .map_err(|_| StorageError::DataReadError(path.into()))?,
        );
        let mut contents = Vec::new();
        reader.read_to_end(&mut contents)?;

        Ok(contents)
    }

    fn args(&self) -> StorageArgs {
        unimplemented!();
    }

    fn load_sig(&self, path: &str) -> Result<SigStore> {
        let raw = self.load(path)?;
        let sig = Signature::from_reader(&mut &raw[..])?
            // TODO: select the right sig?
            .swap_remove(0);

        Ok(sig.into())
    }

    fn spec(&self) -> String {
        format!("zip://{}", self.path().unwrap_or("".into()))
    }
}

impl ZipStorage {
    pub fn from_file<P: AsRef<Path>>(location: P) -> Result<Self> {
        let zip_file = File::open(location.as_ref())?;
        let mapping = unsafe { memmap2::Mmap::map(&zip_file)? };

        let mut storage = ZipStorageBuilder {
            mapping: Some(mapping),
            archive_builder: |mapping: &Option<memmap2::Mmap>| {
                piz::ZipArchive::new(mapping.as_ref().unwrap()).unwrap()
            },
            metadata_builder: |archive: &piz::ZipArchive| {
                archive
                    .entries()
                    .iter()
                    .map(|entry| (entry.path.as_os_str(), entry))
                    .collect()
            },
            subdir: None,
            path: Some(location.as_ref().into()),
        }
        .build();

        let subdir = find_subdirs(storage.borrow_archive())?;
        storage.with_mut(|fields| *fields.subdir = subdir);

        Ok(storage)
    }

    pub fn path(&self) -> Option<PathBuf> {
        self.borrow_path().clone()
    }

    pub fn subdir(&self) -> Option<String> {
        self.borrow_subdir().clone()
    }

    pub fn set_subdir(&mut self, path: String) {
        self.with_mut(|fields| *fields.subdir = Some(path))
    }

    pub fn list_sbts(&self) -> Result<Vec<String>> {
        Ok(self
            .borrow_archive()
            .entries()
            .iter()
            .filter_map(|entry| {
                let path = entry.path.as_str();
                if path.ends_with(".sbt.json") {
                    Some(path.into())
                } else {
                    None
                }
            })
            .collect())
    }

    pub fn filenames(&self) -> Result<Vec<String>> {
        Ok(self
            .borrow_archive()
            .entries()
            .iter()
            .map(|entry| entry.path.as_str().into())
            .collect())
    }
}

impl SigStore {
    pub fn new_with_storage(sig: Signature, storage: InnerStorage) -> Self {
        let name = sig.name();
        let filename = sig.filename();

        SigStore::builder()
            .name(name)
            .filename(filename)
            .data(sig)
            .metadata("")
            .storage(Some(storage))
            .build()
    }

    pub fn name(&self) -> String {
        self.name.clone()
    }
}

impl Select for SigStore {
    fn select(mut self, selection: &Selection) -> Result<Self> {
        // TODO: find better error (perhaps ReadDataError or similar?)
        let sig = self.data.take().ok_or(Error::MismatchKSizes)?;
        self.data = OnceCell::with_value(sig.select(selection)?);

        Ok(self)
    }
}

impl std::fmt::Debug for SigStore {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SigStore [filename: {}, name: {}, metadata: {}]",
            self.filename, self.name, self.metadata
        )
    }
}

impl ReadData<Signature> for SigStore {
    fn data(&self) -> Result<&Signature> {
        if let Some(sig) = self.data.get() {
            Ok(sig)
        } else if let Some(storage) = &self.storage {
            let sig = self.data.get_or_init(|| {
                let raw = storage.load(&self.filename).unwrap();
                Signature::from_reader(&mut &raw[..])
                    .unwrap()
                    // TODO: select the right sig?
                    .swap_remove(0)
            });

            Ok(sig)
        } else {
            Err(ReadDataError::LoadError.into())
        }
    }
}

impl SigStore {
    pub fn save(&self, path: &str) -> Result<String> {
        if let Some(storage) = &self.storage {
            if let Some(data) = self.data.get() {
                let mut buffer = Vec::new();
                data.to_writer(&mut buffer)?;

                Ok(storage.save(path, &buffer)?)
            } else {
                unimplemented!()
            }
        } else {
            unimplemented!()
        }
    }
}

impl From<SigStore> for Signature {
    fn from(other: SigStore) -> Signature {
        other.data.get().unwrap().to_owned()
    }
}

impl Deref for SigStore {
    type Target = Signature;

    fn deref(&self) -> &Signature {
        self.data.get().unwrap()
    }
}

impl From<Signature> for SigStore {
    fn from(other: Signature) -> SigStore {
        let name = other.name();
        let filename = other.filename();

        SigStore::builder()
            .name(name)
            .filename(filename)
            .data(other)
            .metadata("")
            .storage(None)
            .build()
    }
}

impl Comparable<SigStore> for SigStore {
    fn similarity(&self, other: &SigStore) -> f64 {
        let ng: &Signature = self.data().unwrap();
        let ong: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &ng.signatures[0] {
            if let Sketch::MinHash(omh) = &ong.signatures[0] {
                return mh.similarity(omh, true, false).unwrap();
            }
        }

        unimplemented!()
    }

    fn containment(&self, other: &SigStore) -> f64 {
        let ng: &Signature = self.data().unwrap();
        let ong: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &ng.signatures[0] {
            if let Sketch::MinHash(omh) = &ong.signatures[0] {
                let common = mh.count_common(omh, false).unwrap();
                let size = mh.size();
                return common as f64 / size as f64;
            }
        }
        unimplemented!()
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct DatasetInfo {
    pub filename: String,
    pub name: String,
    pub metadata: String,
}
impl From<DatasetInfo> for SigStore {
    fn from(other: DatasetInfo) -> SigStore {
        SigStore {
            filename: other.filename,
            name: other.name,
            metadata: other.metadata,
            storage: None,
            data: OnceCell::new(),
        }
    }
}

impl Comparable<Signature> for Signature {
    fn similarity(&self, other: &Signature) -> f64 {
        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &self.signatures[0] {
            if let Sketch::MinHash(omh) = &other.signatures[0] {
                return mh.similarity(omh, true, false).unwrap();
            }
        }
        unimplemented!()
    }

    fn containment(&self, other: &Signature) -> f64 {
        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &self.signatures[0] {
            if let Sketch::MinHash(omh) = &other.signatures[0] {
                let common = mh.count_common(omh, false).unwrap();
                let size = mh.size();
                return common as f64 / size as f64;
            }
        }
        unimplemented!()
    }
}

impl MemStorage {
    pub fn new() -> Self {
        Self {
            sigs: Arc::new(RwLock::new(HashMap::default())),
        }
    }
}

impl Storage for MemStorage {
    fn save(&self, _path: &str, _content: &[u8]) -> Result<String> {
        unimplemented!()
    }

    fn load(&self, _path: &str) -> Result<Vec<u8>> {
        unimplemented!()
    }

    fn args(&self) -> StorageArgs {
        unimplemented!()
    }

    fn load_sig(&self, path: &str) -> Result<SigStore> {
        Ok(self.sigs.read().unwrap().get(path).unwrap().clone())
    }

    fn save_sig(&self, path: &str, sig: Signature) -> Result<String> {
        // side-step saving to store
        let sig_store: SigStore = sig.into();
        self.sigs.write().unwrap().insert(path.into(), sig_store);
        Ok(path.into())
    }

    fn spec(&self) -> String {
        "memory://".into()
    }
}
