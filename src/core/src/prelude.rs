use std::io::Write;

use crate::Result;

pub use crate::selection::{Select, Selection};
pub use crate::signature::Signature;
pub use crate::storage::Storage;

pub trait ToWriter {
    fn to_writer<W>(&self, writer: &mut W) -> Result<()>
    where
        W: Write;
}

pub trait Update<O> {
    fn update(&self, other: &mut O) -> Result<()>;
}

pub trait FromFactory<N> {
    fn factory(&self, name: &str) -> Result<N>;
}

/// Implemented by anything that wants to read specific data from a storage.
pub trait ReadData<D> {
    fn data(&self) -> Result<&D>;
}

// TODO: split into two traits, Similarity and Containment?
pub trait Comparable<O> {
    fn similarity(&self, other: &O) -> f64;
    fn containment(&self, other: &O) -> f64;
}
