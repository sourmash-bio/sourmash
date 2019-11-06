use std::io::{Read, Write};

use failure::Error;
use mqf::MQF;

use crate::index::sbt::{FromFactory, Node, Update, SBT};
use crate::index::storage::{ReadData, ReadDataError, ToWriter};
use crate::index::{Comparable, Dataset};
use crate::signature::{Signature, SigsTrait};
use crate::sketch::Sketch;

impl ToWriter for MQF {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write,
    {
        // TODO: using tempfile for now, but ideally want to avoid that
        let mut tmpfile = tempfile::NamedTempFile::new()?;
        self.serialize(tmpfile.path()).unwrap(); // TODO: convert this to a proper error

        let mut buffer = Vec::new();
        tmpfile.read_to_end(&mut buffer)?;
        writer.write_all(&buffer)?;

        Ok(())
    }
}

impl ReadData<MQF> for Node<MQF> {
    fn data(&self) -> Result<&MQF, Error> {
        if let Some(storage) = &self.storage {
            Ok(self.data.get_or_create(|| {
                let raw = storage.load(&self.filename).unwrap();

                // TODO: using tempfile for now, but ideally want to avoid that
                let mut tmpfile = tempfile::NamedTempFile::new().unwrap();
                tmpfile.write_all(&mut &raw[..]).unwrap();

                MQF::deserialize(tmpfile.path()).unwrap()
            }))
        } else if let Some(data) = self.data.get() {
            Ok(data)
        } else {
            Err(ReadDataError::LoadError.into())
        }
    }
}

impl<L> FromFactory<Node<MQF>> for SBT<Node<MQF>, L> {
    fn factory(&self, _name: &str) -> Result<Node<MQF>, Error> {
        unimplemented!()
    }
}

impl Update<Node<MQF>> for Node<MQF> {
    fn update(&self, _other: &mut Node<MQF>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Update<Node<MQF>> for Dataset<Signature> {
    fn update(&self, _other: &mut Node<MQF>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Comparable<Node<MQF>> for Node<MQF> {
    fn similarity(&self, other: &Node<MQF>) -> f64 {
        let ng: &MQF = self.data().unwrap();
        let ong: &MQF = other.data().unwrap();
        unimplemented!();
        //ng.similarity(&ong)
    }

    fn containment(&self, other: &Node<MQF>) -> f64 {
        let ng: &MQF = self.data().unwrap();
        let ong: &MQF = other.data().unwrap();
        unimplemented!();
        //ng.containment(&ong)
    }
}

impl Comparable<Dataset<Signature>> for Node<MQF> {
    fn similarity(&self, other: &Dataset<Signature>) -> f64 {
        let ng: &MQF = self.data().unwrap();
        let oth: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &oth.signatures[0] {
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig
                .mins
                .iter()
                .filter(|h| dbg!(ng.count_key(**h % u64::pow(2, 26))) > 0)
                //.filter(|h| dbg!(ng.count_key(**h)) > 0)
                .count();

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
        let ng: &MQF = self.data().unwrap();
        let oth: &Signature = other.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &oth.signatures[0] {
            if sig.size() == 0 {
                return 0.0;
            }

            let matches: usize = sig
                .mins
                .iter()
                .filter(|h| ng.count_key(**h % u64::pow(2, 26)) > 0)
                //.filter(|h| ng.count_key(**h) > 0)
                .count();

            matches as f64 / sig.size() as f64
        } else {
            //TODO what if it is not a minhash?
            unimplemented!()
        }
    }
}
