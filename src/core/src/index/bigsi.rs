use std::collections::HashMap;
use std::path::Path;

use fixedbitset::FixedBitSet;
use thiserror::Error;
use typed_builder::TypedBuilder;

use crate::index::Index;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::nodegraph::Nodegraph;
use crate::sketch::Sketch;
use crate::Error;
use crate::HashIntoType;

#[derive(Clone, TypedBuilder)]
pub struct BIGSI<L> {
    matrix: Vec<FixedBitSet>,
    ksize: usize,
    datasets: Vec<L>,
    //#[builder(setter(skip))]
    //storage: Rc<dyn Storage>,
}

#[derive(Debug, Error)]
pub enum BIGSIError {
    #[error("BIGSI doesn't support this method")]
    MethodDisabled,
}

impl<L> BIGSI<L> {
    pub fn new(bf_size: usize, ksize: usize) -> BIGSI<L> {
        let mut matrix = Vec::with_capacity(bf_size);
        for _ in 0..bf_size {
            // TODO: figure initial capacity for each row
            matrix.push(FixedBitSet::with_capacity(100));
        }

        BIGSI {
            matrix,
            ksize,
            datasets: Vec::new(),
        }
    }
}

impl BIGSI<Signature> {
    pub fn add(&mut self, dataset: Signature) {
        let mut ng = Nodegraph::new(&[self.matrix.len()], self.ksize);

        // TODO: select correct minhash
        if let Sketch::MinHash(mh) = &dataset.signatures[0] {
            for h in mh.mins() {
                ng.count(h);
            }
        } else {
            // TODO: what if it is not a mh?
            unimplemented!()
        }

        self.datasets.push(dataset);
        let col = self.datasets.len() - 1;

        let bs = ng.into_bitsets();
        for pos in bs[0].ones() {
            let bs = &mut self.matrix[pos];
            if bs.len() == col {
                bs.grow(col + col / 2);
            }
            bs.insert(col);
        }
    }

    pub fn query(&self, hash: HashIntoType) -> impl Iterator<Item = usize> + '_ {
        let pos = hash as usize % self.matrix.len();
        let bs = &self.matrix[pos];
        bs.ones()
    }

    pub fn query_datasets(&self, hash: HashIntoType) -> impl Iterator<Item = Signature> + '_ {
        self.query(hash).map(move |pos| self.datasets[pos].clone())
    }
}

impl<'a> Index<'a> for BIGSI<Signature> {
    type Item = Signature;
    //type SignatureIterator = std::slice::Iter<'a, Self::Item>;

    fn search(
        &self,
        sig: &Self::Item,
        threshold: f64,
        containment: bool,
    ) -> Result<Vec<&Self::Item>, Error> {
        let mut results = Vec::new();

        //TODO: still assuming one mh in the signature!
        if let Sketch::MinHash(hashes) = &sig.signatures[0] {
            let mut counter: HashMap<usize, usize> = HashMap::with_capacity(hashes.size());

            for hash in hashes.mins() {
                self.query(hash).for_each(|dataset_idx| {
                    let idx = counter.entry(dataset_idx).or_insert(0);
                    *idx += 1;
                });
            }

            for (idx, count) in counter {
                let match_sig = &self.datasets[idx];
                //TODO: still assuming one mh in the signature!
                let match_mh = match_sig.signatures[0].size();

                let score = if containment {
                    count as f64 / hashes.size() as f64
                } else {
                    count as f64 / (hashes.size() + match_mh - count) as f64
                };

                if score >= threshold {
                    results.push(match_sig)
                }
            }

            Ok(results)
        } else {
            // TODO: what if it is not a minhash?
            unimplemented!()
        }
    }

    fn insert(&mut self, node: Self::Item) -> Result<(), Error> {
        self.add(node);
        Ok(())
    }

    fn save<P: AsRef<Path>>(&self, _path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!()
    }

    fn signatures(&self) -> Vec<Self::Item> {
        unimplemented!()
    }

    fn signature_refs(&self) -> Vec<&Self::Item> {
        unimplemented!()
    }

    /*
    fn iter_signatures(&'a self) -> Self::SignatureIterator {
        self.datasets.iter()
    }
    */
}

#[cfg(test)]
mod test {
    use std::convert::TryInto;
    use std::fs::File;
    use std::io::BufReader;
    use std::path::PathBuf;

    use super::BIGSI;

    use crate::index::SigStore;
    use crate::index::{Index, MHBT};
    use crate::signature::Signature;

    #[test]
    fn bigsi_sbt_oracle() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/v5.sbt.json");

        let sbt = MHBT::from_path(filename).expect("Loading error");

        let mut bigsi = BIGSI::new(10000, 10);
        let datasets = sbt.signatures();

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/60f7e23c24a8d94791cc7a8680c493f9");

        let mut reader = BufReader::new(File::open(filename).unwrap());
        let sigs = Signature::load_signatures(
            &mut reader,
            Some(31),
            Some("DNA".try_into().unwrap()),
            None,
        )
        .unwrap();
        let sig_data = sigs[0].clone();

        let leaf: SigStore<_> = sig_data.into();

        for l in datasets {
            bigsi.insert(l).expect("insertion error!");
        }

        let results_sbt = sbt.search(&leaf, 0.5, false).unwrap();
        assert_eq!(results_sbt.len(), 1);

        let data = leaf.data.get().unwrap();
        let results_bigsi = bigsi.search(data, 0.5, false).unwrap();
        assert_eq!(results_bigsi.len(), 1);

        assert_eq!(results_sbt.len(), results_bigsi.len());

        let results_sbt = sbt.search(&leaf, 0.1, false).unwrap();
        assert_eq!(results_sbt.len(), 2);

        let data = leaf.data.get().unwrap();
        let results_bigsi = bigsi.search(data, 0.1, false).unwrap();
        assert_eq!(results_bigsi.len(), 2);

        assert_eq!(results_sbt.len(), results_bigsi.len());
    }
}
