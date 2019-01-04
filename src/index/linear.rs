use std::path::Path;
use std::rc::Rc;

use derive_builder::Builder;
use failure::Error;

use crate::index::storage::Storage;
use crate::index::{Comparable, Index};

#[derive(Builder)]
pub struct LinearIndex<L> {
    //#[builder(setter(skip))]
    storage: Rc<dyn Storage>,

    #[builder(setter(skip))]
    pub(crate) leaves: Vec<L>,
}

impl<L> Index for LinearIndex<L>
where
    L: Clone + Comparable<L>,
{
    type Item = L;

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
            .leaves
            .iter()
            .flat_map(|node| {
                if search_fn(node, sig, threshold) {
                    Some(node)
                } else {
                    None
                }
            })
            .collect())
    }

    fn insert(&mut self, node: &L) {
        self.leaves.push(node.clone());
    }

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        Ok(())
    }

    fn load<P: AsRef<Path>>(path: P) -> Result<(), Error> {
        Ok(())
    }
}
