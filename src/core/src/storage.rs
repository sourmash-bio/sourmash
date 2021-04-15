use std::fs::{DirBuilder, File};
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::PathBuf;
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

pub struct ZipStorage<'a> {
    mapping: Option<memmap2::Mmap>,
    archive: Option<piz::ZipArchive<'a>>,
    //metadata: piz::read::DirectoryContents<'a>,
}

fn load_from_archive<'a>(archive: &'a piz::ZipArchive<'a>, path: &str) -> Result<Vec<u8>, Error> {
    use piz::read::FileTree;

    // FIXME error
    let tree = piz::read::as_tree(archive.entries()).map_err(|_| StorageError::EmptyPathError)?;
    // FIXME error
    let entry = tree
        .lookup(path)
        .map_err(|_| StorageError::EmptyPathError)?;

    // FIXME error
    let mut reader = BufReader::new(
        archive
            .read(entry)
            .map_err(|_| StorageError::EmptyPathError)?,
    );
    let mut contents = Vec::new();
    reader.read_to_end(&mut contents)?;

    Ok(contents)
}

impl<'a> Storage for ZipStorage<'a> {
    fn save(&self, _path: &str, _content: &[u8]) -> Result<String, Error> {
        unimplemented!();
    }

    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
        if let Some(archive) = &self.archive {
            load_from_archive(archive, path)
        } else {
            //FIXME
            let archive = piz::ZipArchive::new((&self.mapping.as_ref()).unwrap())
                .map_err(|_| StorageError::EmptyPathError)?;
            load_from_archive(&archive, path)
        }
    }

    fn args(&self) -> StorageArgs {
        unimplemented!();
    }
}

impl<'a> ZipStorage<'a> {
    pub fn new(location: &str) -> Result<Self, Error> {
        let zip_file = File::open(location)?;
        let mapping = unsafe { memmap2::Mmap::map(&zip_file)? };

        //FIXME
        //let archive = piz::ZipArchive::new(&mapping).map_err(|_| StorageError::EmptyPathError)?;

        //FIXME
        //  let tree =
        //      piz::read::as_tree(archive.entries()).map_err(|_| StorageError::EmptyPathError)?;

        Ok(Self {
            mapping: Some(mapping),
            archive: None,
            //metadata: tree,
        })
    }

    pub fn from_slice(mapping: &'a [u8]) -> Result<Self, Error> {
        //FIXME
        let archive = piz::ZipArchive::new(mapping).map_err(|_| StorageError::EmptyPathError)?;

        //FIXME
        //let entries: Vec<_> = archive.entries().iter().map(|x| x.to_owned()).collect();
        //let tree =
        //    piz::read::as_tree(entries.as_slice()).map_err(|_| StorageError::EmptyPathError)?;

        Ok(Self {
            archive: Some(archive),
            mapping: None,
            /*            metadata: archive
            .as_tree()
            .map_err(|_| StorageError::EmptyPathError)?, */
        })
    }
}

#[derive(Default)]
pub struct MemStorage {
    storage: Mutex<HashMap<String, Vec<u8>>>,
}

impl MemStorage {}

impl Storage for MemStorage {
    fn save(&self, path: &str, content: &[u8]) -> Result<String, Error> {
        let mut lock = self.storage.lock().unwrap();
        lock.insert(path.into(), content.into());
        Ok(path.into())
    }

    fn load(&self, path: &str) -> Result<Vec<u8>, Error> {
        let lock = self.storage.lock().unwrap();
        let v = lock.get(path).ok_or(ReadDataError::LoadError)?;
        Ok(v.clone())
    }

    fn args(&self) -> StorageArgs {
        unimplemented!()
    }
}
