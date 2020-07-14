use std::io::{Read, Write};

use mqf::MQF;

use crate::Error;
use crate::index::sbt::{FromFactory, Node, Update, SBT};
use crate::index::storage::{ReadData, ReadDataError, ToWriter};
use crate::index::Comparable;
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
                tmpfile.write_all(&raw[..]).unwrap();

                MQF::deserialize(tmpfile.path()).unwrap()
            }))
        } else if let Some(data) = self.data.get() {
            Ok(data)
        } else {
            Err(ReadDataError::LoadError.into())
        }
    }
}

impl<L: Sync> FromFactory<Node<MQF>> for SBT<Node<MQF>, L> {
    fn factory(&self, _name: &str) -> Result<Node<MQF>, Error> {
        unimplemented!()
    }
}

impl Update<Node<MQF>> for Node<MQF> {
    fn update(&self, _other: &mut Node<MQF>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Update<Node<MQF>> for Signature {
    fn update(&self, _other: &mut Node<MQF>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Comparable<Node<MQF>> for Node<MQF> {
    fn similarity(&self, other: &Node<MQF>) -> f64 {
        let _ng: &MQF = self.data().unwrap();
        let _ong: &MQF = other.data().unwrap();
        unimplemented!();
        //ng.similarity(&ong)
    }

    fn containment(&self, other: &Node<MQF>) -> f64 {
        let _ng: &MQF = self.data().unwrap();
        let _ong: &MQF = other.data().unwrap();
        unimplemented!();
        //ng.containment(&ong)
    }
}

impl Comparable<Signature> for Node<MQF> {
    fn similarity(&self, other: &Signature) -> f64 {
        let ng: &MQF = self.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &other.signatures[0] {
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

    fn containment(&self, other: &Signature) -> f64 {
        let ng: &MQF = self.data().unwrap();

        // TODO: select the right signatures...
        if let Sketch::MinHash(sig) = &other.signatures[0] {
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

/* FIXME: bring back after MQF works on macOS and Windows
#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::{BufReader, Seek, SeekFrom};
    use std::path::PathBuf;
    use std::rc::Rc;
    use tempfile;

    use assert_matches::assert_matches;
    use lazy_init::Lazy;

    use super::{scaffold, Factory};

    use crate::index::linear::LinearIndex;
    use crate::index::search::{search_minhashes, search_minhashes_containment};
    use crate::index::storage::ReadData;
    use crate::index::{Index, SigStore, MHBT};
    use crate::signature::Signature;

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn load_mhmt() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/v5_mhmt.sbt.json");

        let mut sbt = crate::index::MHMT::from_path(filename).expect("Loading error");

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/.sbt.v3/60f7e23c24a8d94791cc7a8680c493f9");

        let mut reader = BufReader::new(File::open(filename).unwrap());
        let sigs = Signature::load_signatures(&mut reader, 31, Some("DNA".into()), None).unwrap();
        let sig_data = sigs[0].clone();

        let data = Lazy::new();
        data.get_or_create(|| sig_data);

        let leaf = SigStore::builder()
            .data(Rc::new(data))
            .filename("")
            .name("")
            .metadata("")
            .storage(None)
            .build();

        let results = sbt.find(search_minhashes, &leaf, 0.5).unwrap();
        //assert_eq!(results.len(), 1);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = sbt.find(search_minhashes, &leaf, 0.1).unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let mut linear = LinearIndex::builder().storage(sbt.storage()).build();
        for l in &sbt.leaves {
            linear.insert(l.1.data().unwrap().clone()).unwrap();
        }

        println!(
            "linear leaves {:?} {:?}",
            linear.datasets.len(),
            linear.datasets
        );

        let results = linear.find(search_minhashes, &leaf, 0.5).unwrap();
        assert_eq!(results.len(), 1);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear.find(search_minhashes, &leaf, 0.1).unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear
            .find(search_minhashes_containment, &leaf, 0.5)
            .unwrap();
        assert_eq!(results.len(), 2);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);

        let results = linear
            .find(search_minhashes_containment, &leaf, 0.1)
            .unwrap();
        assert_eq!(results.len(), 4);
        println!("results: {:?}", results);
        println!("leaf: {:?}", leaf);
    }
    */
}
