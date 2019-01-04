pub mod sbt;

pub mod storage;

pub mod nodegraph;

pub mod linear;

pub mod search;

use std::path::Path;
use std::rc::Rc;

use serde_derive::Deserialize;

use derive_builder::Builder;
use failure::Error;
use lazy_init::Lazy;

use crate::index::storage::{ReadData, Storage};
use crate::Signature;

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

    fn insert(&mut self, node: &Self::Item);

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error>;

    fn load<P: AsRef<Path>>(path: P) -> Result<(), Error>;
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

#[derive(Deserialize)]
pub struct LeafInfo {
    pub filename: String,
    pub name: String,
    pub metadata: String,
}

#[derive(Builder, Default, Clone)]
pub struct Leaf<T>
where
    T: std::marker::Sync,
{
    pub(crate) filename: String,
    pub(crate) name: String,
    pub(crate) metadata: String,

    pub(crate) storage: Option<Rc<dyn Storage>>,

    pub(crate) data: Rc<Lazy<T>>,
}

impl<T> std::fmt::Debug for Leaf<T>
where
    T: std::marker::Sync,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Leaf [filename: {}, name: {}, metadata: {}]",
            self.filename, self.name, self.metadata
        )
    }
}

impl<S: Storage + ?Sized> ReadData<Signature, S> for Leaf<Signature> {
    fn data(&self, storage: &S) -> Result<&Signature, Error> {
        let sig = self.data.get_or_create(|| {
            let raw = storage.load(&self.filename).unwrap();
            let sigs: Vec<Signature> = serde_json::from_reader(&mut &raw[..]).unwrap();
            // TODO: select the right sig?
            sigs[0].to_owned()
        });

        Ok(sig)
    }
}

impl Leaf<Signature> {
    pub fn count_common(&self, other: &Leaf<Signature>) -> u64 {
        if let Some(storage) = &self.storage {
            let ng: &Signature = self.data(&**storage).unwrap();
            let ong: &Signature = other.data(&**storage).unwrap();

            // TODO: select the right signatures...
            ng.signatures[0].count_common(&ong.signatures[0]).unwrap() as u64
        } else {
            0
        }
    }

    pub fn mins(&self) -> Vec<u64> {
        if let Some(storage) = &self.storage {
            let ng: &Signature = self.data(&**storage).unwrap();
            ng.signatures[0].mins.to_vec()
        } else {
            Vec::new()
        }
    }
}

impl Comparable<Leaf<Signature>> for Leaf<Signature> {
    fn similarity(&self, other: &Leaf<Signature>) -> f64 {
        if let Some(storage) = &self.storage {
            let ng: &Signature = self.data(&**storage).unwrap();
            let ong: &Signature = other.data(&**storage).unwrap();

            // TODO: select the right signatures...
            ng.signatures[0].compare(&ong.signatures[0]).unwrap()
        } else {
            // TODO: in this case storage is not set up,
            // so we should throw an error?
            0.0
        }
    }

    fn containment(&self, other: &Leaf<Signature>) -> f64 {
        if let Some(storage) = &self.storage {
            let ng: &Signature = self.data(&**storage).unwrap();
            let ong: &Signature = other.data(&**storage).unwrap();

            // TODO: select the right signatures...
            let common = ng.signatures[0].count_common(&ong.signatures[0]).unwrap();
            let size = ng.signatures[0].mins.len();
            common as f64 / size as f64
        } else {
            // TODO: in this case storage is not set up,
            // so we should throw an error?
            0.0
        }
    }
}
