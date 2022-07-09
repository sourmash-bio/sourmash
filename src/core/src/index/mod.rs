//! # Indexing structures for fast similarity search
//!
//! An index organizes signatures to allow for fast similarity search.
//! Some indices also support containment searches.

pub mod linear;

#[cfg(not(target_arch = "wasm32"))]
#[cfg(feature = "branchwater")]
pub mod revindex;

pub mod search;

use std::path::Path;

use getset::{CopyGetters, Getters, Setters};

use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::encodings::{HashFunctions, Idx};
use crate::index::search::{search_minhashes, search_minhashes_containment};
use crate::manifest::Record;
use crate::picklist::Picklist;
use crate::prelude::*;
use crate::signature::SigsTrait;
use crate::sketch::Sketch;
use crate::Result;

#[derive(TypedBuilder, CopyGetters, Getters, Setters, Serialize, Deserialize, Debug, PartialEq)]
pub struct GatherResult {
    #[getset(get_copy = "pub")]
    intersect_bp: usize,

    #[getset(get_copy = "pub")]
    f_orig_query: f64,

    #[getset(get_copy = "pub")]
    f_match: f64,

    f_unique_to_query: f64,
    f_unique_weighted: f64,
    average_abund: usize,
    median_abund: usize,
    std_abund: usize,

    #[getset(get = "pub")]
    filename: String,

    #[getset(get = "pub")]
    name: String,

    #[getset(get = "pub")]
    md5: String,

    #[serde(skip)]
    match_: Signature,

    f_match_orig: f64,
    unique_intersect_bp: usize,
    gather_result_rank: usize,
    remaining_bp: usize,
}

impl GatherResult {
    pub fn get_match(&self) -> Signature {
        self.match_.clone()
    }
}

#[derive(Default, Debug)]
pub struct Selection {
    ksize: Option<u32>,
    abund: Option<bool>,
    num: Option<u32>,
    scaled: Option<u32>,
    containment: Option<bool>,
    moltype: Option<HashFunctions>,
    picklist: Option<Picklist>,
}

type SigCounter = counter::Counter<Idx>;

impl Selection {
    pub fn ksize(&self) -> Option<u32> {
        self.ksize
    }

    pub fn set_ksize(&mut self, ksize: u32) {
        self.ksize = Some(ksize);
    }

    pub fn abund(&self) -> Option<bool> {
        self.abund
    }

    pub fn set_abund(&mut self, value: bool) {
        self.abund = Some(value);
    }

    pub fn num(&self) -> Option<u32> {
        self.num
    }

    pub fn set_num(&mut self, num: u32) {
        self.num = Some(num);
    }

    pub fn scaled(&self) -> Option<u32> {
        self.scaled
    }

    pub fn set_scaled(&mut self, scaled: u32) {
        self.scaled = Some(scaled);
    }

    pub fn containment(&self) -> Option<bool> {
        self.containment
    }

    pub fn set_containment(&mut self, containment: bool) {
        self.containment = Some(containment);
    }

    pub fn moltype(&self) -> Option<HashFunctions> {
        self.moltype
    }

    pub fn set_moltype(&mut self, value: HashFunctions) {
        self.moltype = Some(value);
    }

    pub fn picklist(&self) -> Option<Picklist> {
        self.picklist.clone()
    }

    pub fn set_picklist(&mut self, value: Picklist) {
        self.picklist = Some(value);
    }

    pub fn from_template(template: &Sketch) -> Self {
        let (num, scaled) = match template {
            Sketch::MinHash(mh) => (Some(mh.num()), Some(mh.scaled() as u32)),
            Sketch::LargeMinHash(mh) => (Some(mh.num()), Some(mh.scaled() as u32)),
            _ => (None, None),
        };

        Selection {
            ksize: Some(template.ksize() as u32),
            abund: None,
            containment: None,
            //moltype: Some(template.hash_function()),
            moltype: None,
            num,
            picklist: None,
            scaled,
        }
    }

    pub fn from_record(row: &Record) -> Result<Self> {
        Ok(Self {
            ksize: Some(*row.ksize()),
            abund: Some(*row.with_abundance()),
            moltype: Some(row.moltype()),
            num: None,
            scaled: None,
            containment: None,
            picklist: None,
        })
    }
}

pub trait Index<'a> {
    type Item: Comparable<Self::Item>;
    //type SignatureIterator: Iterator<Item = Self::Item>;

    fn find<F>(&self, search_fn: F, sig: &Self::Item, threshold: f64) -> Result<Vec<&Self::Item>>
    where
        F: Fn(&dyn Comparable<Self::Item>, &Self::Item, f64) -> bool,
    {
        Ok(self
            .signature_refs()
            .into_iter()
            .flat_map(|node| {
                if search_fn(&node, sig, threshold) {
                    Some(node)
                } else {
                    None
                }
            })
            .collect())
    }

    fn search(
        &self,
        sig: &Self::Item,
        threshold: f64,
        containment: bool,
    ) -> Result<Vec<&Self::Item>> {
        if containment {
            self.find(search_minhashes_containment, sig, threshold)
        } else {
            self.find(search_minhashes, sig, threshold)
        }
    }

    //fn gather(&self, sig: &Self::Item, threshold: f64) -> Result<Vec<&Self::Item>>;

    fn insert(&mut self, node: Self::Item) -> Result<()>;

    fn batch_insert(&mut self, nodes: Vec<Self::Item>) -> Result<()> {
        for node in nodes {
            self.insert(node)?;
        }

        Ok(())
    }

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<()>;

    fn load<P: AsRef<Path>>(path: P) -> Result<()>;

    fn signatures(&self) -> Vec<Self::Item>;

    fn signature_refs(&self) -> Vec<&Self::Item>;

    fn len(&self) -> usize {
        self.signature_refs().len()
    }

    fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /*
    fn iter_signatures(&self) -> Self::SignatureIterator;
    */
}

impl<'a, N, L> Comparable<L> for &'a N
where
    N: Comparable<L>,
{
    fn similarity(&self, other: &L) -> f64 {
        (*self).similarity(other)
    }

    fn containment(&self, other: &L) -> f64 {
        (*self).containment(other)
    }
}
