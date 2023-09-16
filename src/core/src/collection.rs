use std::ops::{Deref, DerefMut};

use camino::Utf8Path as Path;

use crate::encodings::Idx;
use crate::manifest::{Manifest, Record};
use crate::prelude::*;
use crate::signature::Signature;
use crate::storage::{FSStorage, InnerStorage, MemStorage, SigStore, Storage, ZipStorage};
use crate::Result;

pub struct Collection {
    pub(crate) manifest: Manifest,
    pub(crate) storage: InnerStorage,
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
}

impl Collection {
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

        let mut records = vec![];
        for (i, sig) in sigs.into_iter().enumerate() {
            let path = format!("{}", i);
            let mut record = Record::from_sig(&sig, &path);
            let path = storage.save_sig(&path, sig)?;
            record.iter_mut().for_each(|rec| {
                rec.set_internal_location(path.clone().into());
            });
            records.extend(record);
        }

        Ok(Self {
            manifest: records.into(),
            storage: InnerStorage::new(storage),
        })
    }

    pub fn from_paths<P: AsRef<Path>>(paths: &[P]) -> Result<Self> {
        // TODO:
        // - Build manifest from paths
        //   - Might need to load the data?
        // - Use FSStorage (figure out if there is a common path between sigs?)
        let records: Vec<Record> = paths
            .iter()
            .flat_map(|p| {
                let recs: Vec<Record> = Signature::from_path(p.as_ref())
                    .unwrap_or_else(|_| panic!("Error processing {:?}", p.as_ref()))
                    .into_iter()
                    .flat_map(|v| Record::from_sig(&v, p.as_ref().as_str()))
                    .collect();
                recs
            })
            //.map(|p| self.collection().storage.load_sig(p.as_str())?.into())
            .collect();

        Ok(Self {
            manifest: records.into(),
            storage: InnerStorage::new(
                FSStorage::builder()
                    .fullpath("".into())
                    .subdir("".into())
                    .build(),
            ),
        })
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
