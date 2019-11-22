use std::io::Write;

use failure::Error;

use crate::index::sbt::{FromFactory, Node, Update, SBT};
use crate::index::storage::{ReadData, ReadDataError, ToWriter};
use crate::index::{Comparable, Dataset};
use crate::signature::{Signature, SigsTrait};
use crate::sketch::nodegraph::Nodegraph;
use crate::sketch::Sketch;

impl ToWriter for Nodegraph {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write,
    {
        self.save_to_writer(writer)
    }
}

impl<L> FromFactory<Node<Nodegraph>> for SBT<Node<Nodegraph>, L> {
    fn factory(&self, _name: &str) -> Result<Node<Nodegraph>, Error> {
        unimplemented!()
    }
}

impl Update<Node<Nodegraph>> for Node<Nodegraph> {
    fn update(&self, _other: &mut Node<Nodegraph>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Update<Node<Nodegraph>> for Dataset<Signature> {
    fn update(&self, _other: &mut Node<Nodegraph>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Comparable<Node<Nodegraph>> for Node<Nodegraph> {
    fn similarity(&self, other: &Node<Nodegraph>) -> f64 {
        let ng: &Nodegraph = self.data().unwrap();
        let ong: &Nodegraph = other.data().unwrap();
        ng.similarity(&ong)
    }

    fn containment(&self, other: &Node<Nodegraph>) -> f64 {
        let ng: &Nodegraph = self.data().unwrap();
        let ong: &Nodegraph = other.data().unwrap();
        ng.containment(&ong)
    }
}

impl Comparable<Dataset<Signature>> for Node<Nodegraph> {
    fn similarity(&self, other: &Dataset<Signature>) -> f64 {
        let ng: &Nodegraph = self.data().unwrap();
        let oth: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &oth.signatures[0] {
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig.mins.iter().map(|h| ng.get(*h)).sum();

            let min_n_below = self.metadata["min_n_below"] as f64;

            // This overestimates the similarity, but better than truncating too
            // soon and losing matches
            matches as f64 / min_n_below
        } else {
            //TODO what if it is not a minhash?
            unimplemented!()
        }
    }

    fn containment(&self, other: &Dataset<Signature>) -> f64 {
        let ng: &Nodegraph = self.data().unwrap();
        let oth: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &oth.signatures[0] {
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig.mins.iter().map(|h| ng.get(*h)).sum();

            matches as f64 / sig.size() as f64
        } else {
            //TODO what if it is not a minhash?
            unimplemented!()
        }
    }
}

impl ReadData<Nodegraph> for Node<Nodegraph> {
    fn data(&self) -> Result<&Nodegraph, Error> {
        if let Some(storage) = &self.storage {
            Ok(self.data.get_or_create(|| {
                let raw = storage.load(&self.filename).unwrap();
                Nodegraph::from_reader(&mut &raw[..]).unwrap()
            }))
        } else if let Some(data) = self.data.get() {
            Ok(data)
        } else {
            Err(ReadDataError::LoadError.into())
        }
    }
}
