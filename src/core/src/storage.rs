use std::collections::BTreeMap;
use std::ffi::OsStr;
use std::fs::{DirBuilder, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

use serde::{Deserialize, Serialize};
use thiserror::Error;
use typed_builder::TypedBuilder;

use crate::Error;

/// An abstraction for any place where we can store data.
pub trait Storage {
    /// Save bytes into path
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error>;

    /// Load bytes from path
    fn load(&self, path: &str) -> Result<Vec<u8>, Error>;

    /// Args for initializing a new Storage
    fn args(&self) -> StorageArgs;
}

#[derive(Clone)]
pub struct InnerStorage(Arc<Mutex<dyn Storage>>);

impl InnerStorage {
    pub fn new(inner: impl Storage + 'static) -> InnerStorage {
        InnerStorage(Arc::new(Mutex::new(inner)))
    }
}

impl Storage for InnerStorage {
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error> {
        self.0.save(path, content)
    }
    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
        self.0.load(path)
    }
    fn args(&self) -> StorageArgs {
        self.0.args()
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

impl<L> Storage for Mutex<L>
where
    L: ?Sized + Storage,
{
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error> {
        self.lock().unwrap().save(path, content)
    }

    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
        self.lock().unwrap().load(path)
    }

    fn args(&self) -> StorageArgs {
        self.lock().unwrap().args()
    }
}

/// Store files locally into a directory
#[derive(TypedBuilder, Debug, Clone, Default)]
pub struct FSStorage {
    /// absolute path for the directory where data is saved.
    fullpath: PathBuf,
    subdir: String,
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
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error> {
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

    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
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
}

#[ouroboros::self_referencing]
pub struct ZipStorage {
    mapping: Option<memmap2::Mmap>,

    #[borrows(mapping)]
    #[covariant]
    archive: piz::ZipArchive<'this>,

    subdir: Option<String>,
    path: Option<String>,

    #[borrows(archive)]
    #[covariant]
    metadata: Metadata<'this>,
}

pub type Metadata<'a> = BTreeMap<&'a OsStr, &'a piz::read::FileMetadata<'a>>;

fn lookup<'a, P: AsRef<Path>>(
    metadata: &'a Metadata,
    path: P,
) -> Result<&'a piz::read::FileMetadata<'a>, Error> {
    let path = path.as_ref();
    metadata
        .get(&path.as_os_str())
        .ok_or_else(|| StorageError::PathNotFoundError(path.to_str().unwrap().into()).into())
        .map(|entry| *entry)
}

fn find_subdirs<'a>(archive: &'a piz::ZipArchive<'a>) -> Result<Option<String>, Error> {
    let subdirs: Vec<_> = archive
        .entries()
        .iter()
        .filter(|entry| entry.is_dir())
        .collect();
    if subdirs.len() == 1 {
        Ok(Some(
            subdirs[0]
                .path
                .to_str()
                .expect("Error converting path")
                .into(),
        ))
    } else {
        Ok(None)
    }
}

impl Storage for ZipStorage {
    fn save(&self, _path: &str, _content: &[u8]) -> Result<String, Error> {
        unimplemented!();
    }

    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
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
}

impl ZipStorage {
    pub fn from_file(location: &str) -> Result<Self, Error> {
        let zip_file = File::open(location)?;
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
            path: Some(location.to_owned()),
        }
        .build();

        let subdir = find_subdirs(storage.borrow_archive())?;
        storage.with_mut(|fields| *fields.subdir = subdir);

        Ok(storage)
    }

    pub fn path(&self) -> Option<String> {
        self.borrow_path().clone()
    }

    pub fn subdir(&self) -> Option<String> {
        self.borrow_subdir().clone()
    }

    pub fn set_subdir(&mut self, path: String) {
        self.with_mut(|fields| *fields.subdir = Some(path))
    }

    pub fn list_sbts(&self) -> Result<Vec<String>, Error> {
        Ok(self
            .borrow_archive()
            .entries()
            .iter()
            .filter_map(|entry| {
                let path = entry.path.to_str().expect("Error converting path");
                if path.ends_with(".sbt.json") {
                    Some(path.into())
                } else {
                    None
                }
            })
            .collect())
    }

    pub fn filenames(&self) -> Result<Vec<String>, Error> {
        Ok(self
            .borrow_archive()
            .entries()
            .iter()
            .map(|entry| entry.path.to_str().expect("Error converting path").into())
            .collect())
    }
}
