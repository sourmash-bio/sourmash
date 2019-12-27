use std::collections::HashMap;
use std::mem;
use std::rc::Rc;

use failure::Error;
use lazy_init::Lazy;

use crate::index::sbt::{FromFactory, Node, Update, SBT};
use crate::index::storage::{ReadData, ReadDataError};
use crate::index::Comparable;
use crate::signature::Signature;
use crate::sketch::ukhs::{FlatUKHS, UKHSTrait};
use crate::sketch::Sketch;

impl<L: Sync> FromFactory<Node<FlatUKHS>> for SBT<Node<FlatUKHS>, L> {
    fn factory(&self, name: &str) -> Result<Node<FlatUKHS>, Error> {
        let data = Lazy::new();
        // TODO: don't hardcode this!
        data.get_or_create(|| FlatUKHS::new(9, 31).unwrap());

        Ok(Node {
            name: name.into(),
            filename: name.into(),
            metadata: HashMap::default(),
            storage: self.storage.clone(),
            data: Rc::new(data),
        })
    }
}

impl Update<Node<FlatUKHS>> for Node<FlatUKHS> {
    fn update(&self, _other: &mut Node<FlatUKHS>) -> Result<(), Error> {
        unimplemented!();
    }
}

impl Update<Node<FlatUKHS>> for Signature {
    fn update(&self, other: &mut Node<FlatUKHS>) -> Result<(), Error> {
        let sigs = if self.signatures.len() > 1 {
            // TODO: select the right signatures...
            unimplemented!()
        } else {
            &self.signatures[0]
        };

        if let Sketch::UKHS(sig) = sigs {
            let mut data: FlatUKHS = other.data()?.clone();
            data.merge(sig);

            // TODO update metadata?

            let new_data = Lazy::new();
            new_data.get_or_create(|| data);

            mem::replace(&mut other.data, Rc::new(new_data));
            return Ok(());
        }
        unimplemented!()
    }
}

impl Comparable<Node<FlatUKHS>> for Node<FlatUKHS> {
    fn similarity(&self, other: &Node<FlatUKHS>) -> f64 {
        let o_sig: &FlatUKHS = other.data().unwrap();
        let me_sig: &FlatUKHS = self.data().unwrap();
        1.0 - me_sig.distance(o_sig)
    }

    fn containment(&self, _other: &Node<FlatUKHS>) -> f64 {
        unimplemented!();
    }
}

impl Comparable<Signature> for Node<FlatUKHS> {
    fn similarity(&self, other: &Signature) -> f64 {
        if other.signatures.len() > 1 {
            // TODO: select the right signatures...
            unimplemented!()
        } else if let Sketch::UKHS(o_sig) = &other.signatures[0] {
            // This is doing a variation of Weighted Jaccard.
            // The internal nodes are built with max(l_i, r_i) for each
            // left and right children, so if we do a WJ similarity directly
            // we will underestimate it.
            //
            // Instead of doing sum(mins) / sum(max), we do instead
            // sum(mins) / sum(sig_buckets), which might overestimate
            // but will never loose any leaf.
            // TODO: is this weighted containment?
            //
            // TODO: Still need to test with larger trees! Does it save any
            // comparisons?
            let me_sig: &FlatUKHS = self.data().unwrap();
            let mins: u64 = me_sig
                .buckets()
                .zip(o_sig.buckets())
                .map(|(a, b)| u64::min(*a, *b))
                .sum();
            let maxs: u64 = o_sig.buckets().sum();

            mins as f64 / maxs as f64
        } else {
            // TODO: sig[0] was not a UKHS
            unimplemented!()
        }
    }

    fn containment(&self, _other: &Signature) -> f64 {
        unimplemented!();
    }
}

impl ReadData<FlatUKHS> for Node<FlatUKHS> {
    fn data(&self) -> Result<&FlatUKHS, Error> {
        if let Some(storage) = &self.storage {
            Ok(self.data.get_or_create(|| {
                let raw = storage.load(&self.filename).unwrap();
                FlatUKHS::from_reader(&mut &raw[..]).unwrap()
            }))
        } else if let Some(data) = self.data.get() {
            Ok(data)
        } else {
            Err(ReadDataError::LoadError.into())
        }
    }
}
