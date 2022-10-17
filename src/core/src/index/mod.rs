//! # Indexing structures for fast similarity search
//!
//! An index organizes signatures to allow for fast similarity search.
//! Some indices also support containment searches.

pub mod linear;
pub mod revindex;

pub mod search;

use std::ops::Deref;
use std::path::Path;

use once_cell::sync::OnceCell;
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::encodings::HashFunctions;
use crate::errors::ReadDataError;
use crate::index::search::{search_minhashes, search_minhashes_containment};
use crate::prelude::*;
use crate::signature::SigsTrait;
use crate::sketch::Sketch;
use crate::storage::{InnerStorage, Storage};
use crate::Error;

#[derive(Default)]
pub struct Selection {
    ksize: Option<u32>,
    abund: Option<bool>,
    num: Option<u32>,
    scaled: Option<u32>,
    containment: Option<bool>,
    moltype: Option<HashFunctions>,
}

impl Selection {
    pub fn ksize(&self) -> Option<u32> {
        self.ksize
    }

    pub fn set_ksize(&mut self, ksize: u32) {
        self.ksize = Some(ksize);
    }

    pub fn abund(&self) -> Option<bool> {
        self.abund
    }

    pub fn set_abund(&mut self, value: bool) {
        self.abund = Some(value);
    }

    pub fn num(&self) -> Option<u32> {
        self.num
    }

    pub fn set_num(&mut self, num: u32) {
        self.num = Some(num);
    }

    pub fn scaled(&self) -> Option<u32> {
        self.scaled
    }

    pub fn set_scaled(&mut self, scaled: u32) {
        self.scaled = Some(scaled);
    }

    pub fn containment(&self) -> Option<bool> {
        self.containment
    }

    pub fn set_containment(&mut self, containment: bool) {
        self.containment = Some(containment);
    }

    pub fn moltype(&self) -> Option<HashFunctions> {
        self.moltype
    }

    pub fn set_moltype(&mut self, value: HashFunctions) {
        self.moltype = Some(value);
    }
}

pub trait Index<'a> {
    type Item: Comparable<Self::Item>;
    //type SignatureIterator: Iterator<Item = Self::Item>;

    fn find<F>(
        &self,
        search_fn: F,
        sig: &Self::Item,
        threshold: f64,
    ) -> Result<Vec<&Self::Item>, Error>
    where
        F: Fn(&dyn Comparable<Self::Item>, &Self::Item, f64) -> bool,
    {
        Ok(self
            .signature_refs()
            .into_iter()
            .flat_map(|node| {
                if search_fn(&node, sig, threshold) {
                    Some(node)
                } else {
                    None
                }
            })
            .collect())
    }

    fn search(
        &self,
        sig: &Self::Item,
        threshold: f64,
        containment: bool,
    ) -> Result<Vec<&Self::Item>, Error> {
        if containment {
            self.find(search_minhashes_containment, sig, threshold)
        } else {
            self.find(search_minhashes, sig, threshold)
        }
    }

    //fn gather(&self, sig: &Self::Item, threshold: f64) -> Result<Vec<&Self::Item>, Error>;

    fn insert(&mut self, node: Self::Item) -> Result<(), Error>;

    fn batch_insert(&mut self, nodes: Vec<Self::Item>) -> Result<(), Error> {
        for node in nodes {
            self.insert(node)?;
        }

        Ok(())
    }

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error>;

    fn load<P: AsRef<Path>>(path: P) -> Result<(), Error>;

    fn signatures(&self) -> Vec<Self::Item>;

    fn signature_refs(&self) -> Vec<&Self::Item>;

    fn len(&self) -> usize {
        self.signature_refs().len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /*
    fn iter_signatures(&self) -> Self::SignatureIterator;
    */
}

impl<'a, N, L> Comparable<L> for &'a N
where
    N: Comparable<L>,
{
    fn similarity(&self, other: &L) -> f64 {
        (*self).similarity(other)
    }

    fn containment(&self, other: &L) -> f64 {
        (*self).containment(other)
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct DatasetInfo {
    pub filename: String,
    pub name: String,
    pub metadata: String,
}

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

impl SigStore {
    pub fn name(&self) -> String {
        self.name.clone()
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
    fn data(&self) -> Result<&Signature, Error> {
        if let Some(sig) = self.data.get() {
            Ok(sig)
        } else if let Some(storage) = &self.storage {
            let sig = self.data.get_or_init(|| {
                let raw = storage.load(&self.filename).unwrap();
                let sigs: Result<Vec<Signature>, _> = serde_json::from_reader(&mut &raw[..]);
                if let Ok(sigs) = sigs {
                    // TODO: select the right sig?
                    sigs[0].to_owned()
                } else {
                    let sig: Signature = serde_json::from_reader(&mut &raw[..]).unwrap();
                    sig
                }
            });

            Ok(sig)
        } else {
            Err(ReadDataError::LoadError.into())
        }
    }
}

impl SigStore {
    pub fn save(&self, path: &str) -> Result<String, Error> {
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

impl SigStore {
    pub fn count_common(&self, other: &SigStore) -> u64 {
        let ng: &Signature = self.data().unwrap();
        let ong: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &ng.signatures[0] {
            if let Sketch::MinHash(omh) = &ong.signatures[0] {
                return mh.count_common(omh, false).unwrap() as u64;
            }
        }
        unimplemented!();
    }

    pub fn mins(&self) -> Vec<u64> {
        let ng: &Signature = self.data().unwrap();

        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &ng.signatures[0] {
            mh.mins()
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

        /* FIXME: bring back after boomphf changes
        if let Sketch::UKHS(mh) = &ng.signatures[0] {
            if let Sketch::UKHS(omh) = &ong.signatures[0] {
                return 1. - mh.distance(&omh);
            }
        }
        */

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

impl Comparable<Signature> for Signature {
    fn similarity(&self, other: &Signature) -> f64 {
        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &self.signatures[0] {
            if let Sketch::MinHash(omh) = &other.signatures[0] {
                return mh.similarity(omh, true, false).unwrap();
            }
        }

        /* FIXME: bring back after boomphf changes
        if let Sketch::UKHS(mh) = &self.signatures[0] {
            if let Sketch::UKHS(omh) = &other.signatures[0] {
                return 1. - mh.distance(&omh);
            }
        }
        */

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
