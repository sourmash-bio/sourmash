use std::fs::File;

use camino::Utf8Path as Path;
use camino::Utf8PathBuf as PathBuf;
use rc_zip_sync::{ArchiveHandle, ReadZip};

use crate::prelude::*;
use crate::storage::{StorageError,SigStore, StorageArgs};
use crate::Result;

/// Store files in a zip file.
#[ouroboros::self_referencing]
pub struct ZipStorage {
    file: std::fs::File,

    #[borrows(file)]
    #[covariant]
    archive: ArchiveHandle<'this, std::fs::File>,

    subdir: Option<String>,
    path: Option<PathBuf>,
}

impl Storage for ZipStorage {
    fn save(&self, _path: &str, _content: &[u8]) -> Result<String> {
        unimplemented!();
    }

    fn load(&self, path: &str) -> Result<Vec<u8>> {
        let archive = self.borrow_archive();
        if let Some(entry) = archive.by_name(path) {
            return Ok(entry.bytes()?);
        }

        if let Some(subdir) = &self.borrow_subdir() {
            if let Some(entry) = archive.by_name(subdir.to_owned() + path) {
                return Ok(entry.bytes()?);
            }
        }

        Err(StorageError::PathNotFoundError(path.into()).into())
    }

    fn args(&self) -> StorageArgs {
        unimplemented!();
    }

    fn load_sig(&self, path: &str) -> Result<SigStore> {
        let raw = self.load(path)?;
        let mut vs = Signature::from_reader(&mut &raw[..])?;
        if vs.len() > 1 {
            unimplemented!("only one Signature currently allowed");
        }
        let sig = vs.swap_remove(0);

        Ok(sig.into())
    }

    fn spec(&self) -> String {
        format!("zip://{}", self.borrow_path().clone().unwrap_or("".into()))
    }
}

impl ZipStorage {
    pub fn from_file<P: AsRef<Path>>(location: P) -> Result<Self> {
        let file = File::open(location.as_ref())?;

        let mut storage = ZipStorageBuilder {
            file,
            archive_builder: |file: &std::fs::File| file.read_zip().expect("Error loading zipfile"),
            subdir: None,
            path: Some(location.as_ref().into()),
        }
        .build();

        let subdir = {
            let subdirs: Vec<_> = storage
                .borrow_archive()
                .entries()
                .filter(|entry| matches!(entry.kind(), rc_zip::parse::EntryKind::Directory))
                .collect();
            if subdirs.len() == 1 {
                Some(
                    subdirs[0]
                        .sanitized_name()
                        .expect("TODO throw right error")
                        .into(),
                )
            } else {
                None
            }
        };

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
            .filter_map(|entry| {
                let path = entry.sanitized_name().expect("TODO throw right error");
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
            .map(|entry| {
                entry
                    .sanitized_name()
                    .expect("TODO throw right error")
                    .into()
            })
            .collect())
    }
}
