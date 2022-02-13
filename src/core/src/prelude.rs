use std::io::Write;

use crate::Error;

pub use crate::signature::Signature;
pub use crate::storage::Storage;

pub trait ToWriter {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write;
}

pub trait Update<O> {
    fn update(&self, other: &mut O) -> Result<(), Error>;
}

pub trait FromFactory<N> {
    fn factory(&self, name: &str) -> Result<N, Error>;
}

/// Implemented by anything that wants to read specific data from a storage.
pub trait ReadData<D> {
    fn data(&self) -> Result<&D, Error>;
}

// TODO: split into two traits, Similarity and Containment?
pub trait Comparable<O> {
    fn similarity(&self, other: &O) -> f64;
    fn containment(&self, other: &O) -> f64;
}
