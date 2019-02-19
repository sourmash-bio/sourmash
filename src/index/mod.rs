//! # Indexing structures for fast similarity search
//!
//! An index organizes signatures to allow for fast similarity search.
//! Some indices also support containment searches.

pub mod bigsi;
pub mod linear;
pub mod sbt;

pub mod storage;

pub mod search;

use std::path::Path;
use std::rc::Rc;

use serde_derive::{Deserialize, Serialize};

use failure::Error;
use lazy_init::Lazy;
use typed_builder::TypedBuilder;

use crate::index::sbt::{Node, SBT};
use crate::index::search::{search_minhashes, search_minhashes_containment};
use crate::index::storage::{ReadData, ReadDataError, Storage};
use crate::signature::Signature;
use crate::sketch::nodegraph::Nodegraph;
use crate::sketch::ukhs::{FlatUKHS, UKHSTrait};
use crate::sketch::Sketch;

pub type MHBT = SBT<Node<Nodegraph>, Dataset<Signature>>;
pub type UKHSTree = SBT<Node<FlatUKHS>, Dataset<Signature>>;

pub trait Index {
    type Item;

    fn find<F>(
        &self,
        search_fn: F,
        sig: &Self::Item,
        threshold: f64,
    ) -> Result<Vec<&Self::Item>, Error>
    where
        F: Fn(&dyn Comparable<Self::Item>, &Self::Item, f64) -> bool;

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

    fn insert(&mut self, node: &Self::Item) -> Result<(), Error>;

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error>;

    fn load<P: AsRef<Path>>(path: P) -> Result<(), Error>;

    fn datasets(&self) -> Vec<Self::Item>;
}

// TODO: split into two traits, Similarity and Containment?
pub trait Comparable<O> {
    fn similarity(&self, other: &O) -> f64;
    fn containment(&self, other: &O) -> f64;
}

impl<'a, N, L> Comparable<L> for &'a N
where
    N: Comparable<L>,
{
    fn similarity(&self, other: &L) -> f64 {
        (*self).similarity(&other)
    }

    fn containment(&self, other: &L) -> f64 {
        (*self).containment(&other)
    }
}

#[derive(Serialize, Deserialize, Debug)]
pub struct DatasetInfo {
    pub filename: String,
    pub name: String,
    pub metadata: String,
}

#[derive(TypedBuilder, Default, Clone)]
pub struct Dataset<T>
where
    T: std::marker::Sync,
{
    pub(crate) filename: String,
    pub(crate) name: String,
    pub(crate) metadata: String,

    pub(crate) storage: Option<Rc<dyn Storage>>,

    pub(crate) data: Rc<Lazy<T>>,
}

impl<T> Dataset<T>
where
    T: std::marker::Sync + Default,
{
    pub fn name(&self) -> String {
        self.name.clone()
    }
}

impl<T> std::fmt::Debug for Dataset<T>
where
    T: std::marker::Sync,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Dataset [filename: {}, name: {}, metadata: {}]",
            self.filename, self.name, self.metadata
        )
    }
}

impl ReadData<Signature> for Dataset<Signature> {
    fn data(&self) -> Result<&Signature, Error> {
        if let Some(sig) = self.data.get() {
            Ok(sig)
        } else if let Some(storage) = &self.storage {
            let sig = self.data.get_or_create(|| {
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

impl Dataset<Signature> {
    pub fn count_common(&self, other: &Dataset<Signature>) -> u64 {
        let ng: &Signature = self.data().unwrap();
        let ong: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &ng.signatures[0] {
            if let Sketch::MinHash(omh) = &ong.signatures[0] {
                return mh.count_common(&omh).unwrap() as u64;
            }
        }
        unimplemented!();
    }

    pub fn mins(&self) -> Vec<u64> {
        let ng: &Signature = self.data().unwrap();

        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &ng.signatures[0] {
            mh.mins.to_vec()
        } else {
            unimplemented!()
        }
    }
}

impl From<Dataset<Signature>> for Signature {
    fn from(other: Dataset<Signature>) -> Signature {
        other.data.get().unwrap().to_owned()
    }
}

impl From<Signature> for Dataset<Signature> {
    fn from(other: Signature) -> Dataset<Signature> {
        let name = other.name();
        let filename = other.filename();

        let data = Lazy::new();
        data.get_or_create(|| other);

        Dataset::builder()
            .name(name)
            .filename(filename)
            .data(data)
            .metadata("")
            .storage(None)
            .build()
    }
}

impl Comparable<Dataset<Signature>> for Dataset<Signature> {
    fn similarity(&self, other: &Dataset<Signature>) -> f64 {
        let ng: &Signature = self.data().unwrap();
        let ong: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &ng.signatures[0] {
            if let Sketch::MinHash(omh) = &ong.signatures[0] {
                return mh.compare(&omh).unwrap();
            }
        }

        if let Sketch::UKHS(mh) = &ng.signatures[0] {
            if let Sketch::UKHS(omh) = &ong.signatures[0] {
                return 1. - mh.distance(&omh);
            }
        }

        unimplemented!()
    }

    fn containment(&self, other: &Dataset<Signature>) -> f64 {
        let ng: &Signature = self.data().unwrap();
        let ong: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        // TODO: better matching here, what if it is not a mh?
        if let Sketch::MinHash(mh) = &ng.signatures[0] {
            if let Sketch::MinHash(omh) = &ong.signatures[0] {
                let common = mh.count_common(&omh).unwrap();
                let size = mh.mins.len();
                return common as f64 / size as f64;
            }
        }
        unimplemented!()
    }
}
