use std::ops::{Deref, DerefMut};

use camino::Utf8Path as Path;
use camino::Utf8PathBuf as PathBuf;

use crate::encodings::Idx;
use crate::manifest::{Manifest, Record};
use crate::prelude::*;
use crate::signature::Signature;
use crate::storage::{FSStorage, InnerStorage, MemStorage, SigStore, Storage, ZipStorage};
use crate::{Error, Result};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub struct Collection {
    manifest: Manifest,
    storage: InnerStorage,
}

pub struct CollectionSet {
    collection: Collection,
}

impl Deref for CollectionSet {
    type Target = Collection;

    fn deref(&self) -> &Self::Target {
        &self.collection
    }
}

impl DerefMut for CollectionSet {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.collection
    }
}

impl TryFrom<Collection> for CollectionSet {
    type Error = crate::Error;

    fn try_from(collection: Collection) -> Result<Self> {
        let first = if let Some(first) = collection.manifest.first() {
            first
        } else {
            // empty collection is consistent ¯\_(ツ)_/¯
            return Ok(Self { collection });
        };

        collection
            .manifest
            .iter()
            .skip(1)
            .try_for_each(|c| first.check_compatible(c))?;

        Ok(Self { collection })
    }
}

impl CollectionSet {
    pub fn into_inner(self) -> Collection {
        self.collection
    }

    pub fn selection(&self) -> Selection {
        todo!("Extract selection from first sig")
    }
}

impl Collection {
    pub fn new(manifest: Manifest, storage: InnerStorage) -> Self {
        Self { manifest, storage }
    }

    pub fn iter(&self) -> impl Iterator<Item = (Idx, &Record)> {
        self.manifest.iter().enumerate().map(|(i, r)| (i as Idx, r))
    }

    #[cfg(feature = "parallel")]
    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = (Idx, &Record)> {
        self.manifest
            .par_iter()
            .enumerate()
            .map(|(i, r)| (i as Idx, r))
    }

    pub fn len(&self) -> usize {
        self.manifest.len()
    }

    pub fn is_empty(&self) -> bool {
        self.manifest.len() == 0
    }

    pub fn manifest(&self) -> &Manifest {
        &self.manifest
    }

    pub fn storage(&self) -> &InnerStorage {
        &self.storage
    }

    pub fn check_superset(&self, other: &Collection) -> Result<usize> {
        self.iter()
            .zip(other.iter())
            .all(|((id1, rec1), (id2, rec2))| id1 == id2 && rec1 == rec2)
            .then(|| self.len())
            // TODO: right error here
            .ok_or(Error::MismatchKSizes)
    }

    pub fn from_zipfile<P: AsRef<Path>>(zipfile: P) -> Result<Self> {
        let storage = ZipStorage::from_file(zipfile)?;
        // Load manifest from standard location in zipstorage
        let manifest = Manifest::from_reader(storage.load("SOURMASH-MANIFEST.csv")?.as_slice())?;
        Ok(Self {
            manifest,
            storage: InnerStorage::new(storage),
        })
    }

    pub fn from_sigs(sigs: Vec<Signature>) -> Result<Self> {
        let storage = MemStorage::new();

        #[cfg(feature = "parallel")]
        let iter = sigs.into_par_iter();

        #[cfg(not(feature = "parallel"))]
        let iter = sigs.into_iter();

        let records: Vec<_> = iter
            .enumerate()
            .flat_map(|(i, sig)| {
                let path = format!("{}", i);
                let mut record = Record::from_sig(&sig, &path);
                let path = storage.save_sig(&path, sig).expect("Error saving sig");
                record.iter_mut().for_each(|rec| {
                    rec.set_internal_location(path.clone().into());
                });
                record
            })
            .collect();

        Ok(Self {
            manifest: records.into(),
            storage: InnerStorage::new(storage),
        })
    }

    pub fn from_paths(paths: &[PathBuf]) -> Result<Self> {
        // TODO:
        // - figure out if there is a common path between sigs for FSStorage?

        Ok(Self {
            manifest: paths.into(),
            storage: InnerStorage::new(
                FSStorage::builder()
                    .fullpath("".into())
                    .subdir("".into())
                    .build(),
            ),
        })
    }

    pub fn record_for_dataset(&self, dataset_id: Idx) -> Result<&Record> {
        Ok(&self.manifest[dataset_id as usize])
    }

    pub fn sig_for_dataset(&self, dataset_id: Idx) -> Result<SigStore> {
        let match_path = if self.manifest.is_empty() {
            ""
        } else {
            self.manifest[dataset_id as usize]
                .internal_location()
                .as_str()
        };

        let selection = Selection::from_record(&self.manifest[dataset_id as usize])?;
        let sig = self.storage.load_sig(match_path)?.select(&selection)?;
        assert_eq!(sig.signatures.len(), 1);
        Ok(sig)
    }
}

impl Select for Collection {
    fn select(mut self, selection: &Selection) -> Result<Self> {
        self.manifest = self.manifest.select(selection)?;
        Ok(self)
    }
}
