use std::fs::File;
use std::io;
use std::iter::Iterator;
use std::path::Path;
use std::str;
use std::ffi::CString;

use cfg_if::cfg_if;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use once_cell::sync::Lazy;
use serde::de::Deserializer;
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

use crate::index::storage::ToWriter;
use crate::sketch::minhash::{KmerMinHash, HashFunctions, max_hash_for_scaled};
use crate::sketch::Sketch;
use crate::Error;

impl Serialize for LcaDB {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = 12;

        let mut partial = serializer.serialize_struct("LcaDB", n_fields)?;
        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("scaled", &self.scaled)?;
        partial.serialize_field("filename", &self.filename.to_string())?;
        partial.serialize_field("moltype", &self.moltype.to_string())?;
        // partial.serialize_field("_next_index", &self._next_index)?;
        // partial.serialize_field("_next_lid", &self._next_lid)?;
        // partial.serialize_field("ident_to_name", &self.ident_to_name)?;
        // partial.serialize_field("ident_to_idx", &self.ident_to_idx)?;
        // partial.serialize_field("idx_to_lid", &self.idx_to_lid)?;
        // partial.serialize_field("lineage_to_lid", &self._next_lid)?;
        // partial.serialize_field("lid_to_lineage", &self.lid_to_lineage)?;
        // partial.serialize_field("hashval_to_idx", &self.hashval_to_idx)?;

        partial.end()
    }
}

impl<'de> Deserialize<'de> for LcaDB {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct TempLcaDB {
            ksize: u32,
            scaled: u64,
            filename: String,
            moltype: String,
            _next_index: u32,
            _next_lid: u32,
            ident_to_name: Option<Vec<u32>>,
            ident_to_idx: Option<Vec<u32>>,
            idx_to_lid: Option<Vec<u32>>,
            lineage_to_lid: Option<Vec<u32>>,
            lid_to_lineage: Option<Vec<u32>>,
            hashval_to_idx: Option<Vec<u32>>,
        }

        let templcadb = TempLcaDB::deserialize(deserializer)?;

        Ok(LcaDB {
            ksize: templcadb.ksize,
            scaled: templcadb.scaled,
            filename: templcadb.filename,
            moltype: templcadb.moltype,
            _next_index: templcadb._next_index,
            _next_lid: templcadb._next_lid,
            ident_to_name: templcadb.ident_to_name,
            ident_to_idx: templcadb.ident_to_idx,
            idx_to_lid: templcadb.idx_to_lid,
            lineage_to_lid: templcadb.lineage_to_lid,
            lid_to_lineage: templcadb.lid_to_lineage,
            hashval_to_idx: templcadb.hashval_to_idx,
        })
    }
}

#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(Debug, TypedBuilder)]
pub struct LcaDB {
    ksize: u32,
    scaled: u64,
    moltype: String,

    #[builder(default = "filename".to_string())]
    filename: String,

    #[builder(default = 0u32)]
    _next_index: u32,

    #[builder(default = 0u32)]
    _next_lid: u32,

    #[builder(default)]
    ident_to_name: Option<Vec<u32>>,

    #[builder(default)]
    ident_to_idx: Option<Vec<u32>>,

    #[builder(default)]
    idx_to_lid: Option<Vec<u32>>,

    #[builder(default)]
    lineage_to_lid: Option<Vec<u32>>,

    #[builder(default)]
    lid_to_lineage: Option<Vec<u32>>,

    #[builder(default)]
    hashval_to_idx: Option<Vec<u32>>, //HashSet as keys
}

impl Default for LcaDB {
    fn default() -> LcaDB {
        LcaDB {
            ksize: 21,
            scaled: 0,
            filename: "filename".to_string(),
            moltype: "DNA".to_string(),

            _next_index: 0,
            _next_lid: 0,
            ident_to_name: None,
            ident_to_idx: None,
            idx_to_lid: None,
            lineage_to_lid: None,
            lid_to_lineage: None,
            hashval_to_idx: None,
        }
    }
}

impl LcaDB {
    pub fn new(
        ksize: u32,
        scaled: u64,
        filename: &[u8],
        moltype: &[u8],
    ) -> LcaDB {

        let filename = str::from_utf8(filename).unwrap().to_string();
        let moltype = str::from_utf8(moltype).unwrap().to_string();
        
        let _next_index = 0;
        let _next_lid = 0;
        let ident_to_name = None;
        let ident_to_idx = None;
        let idx_to_lid = None;
        let lineage_to_lid = None;
        let lid_to_lineage = None;
        let hashval_to_idx = None;

        LcaDB {
            ksize,
            scaled,
            filename,
            moltype,
            _next_index,
            _next_lid,
            ident_to_name,
            ident_to_idx,
            idx_to_lid,
            lineage_to_lid,
            lid_to_lineage,
            hashval_to_idx,
        }
    }

    pub fn ksize(&self) -> u32 {
        self.ksize
    }

    pub fn scaled(&self) -> u64 {
        self.scaled
    }

    pub fn filename(&self) -> String {
        self.filename.clone()
    }

    pub fn moltype(&self) -> String {
        self.moltype.clone()
    }

    pub fn _next_index(&self) -> u32 {
        self._next_index
    }

    pub fn _next_lid(&self) -> u32 {
        self._next_lid
    }

    pub fn ident_to_name(&self) -> Option<Vec<u32>> {
        self.ident_to_name.clone()
    }

    pub fn ident_to_idx(&self) -> Option<Vec<u32>> {
        self.ident_to_idx.clone()
    }

    pub fn idx_to_lid(&self) -> Option<Vec<u32>> {
        self.idx_to_lid.clone()
    }

    pub fn lineage_to_lid(&self) -> Option<Vec<u32>> {
        self.lineage_to_lid.clone()
    }

    pub fn lid_to_lineage(&self) -> Option<Vec<u32>> {
        self.lid_to_lineage.clone()
    }

    pub fn hashval_to_idx(&self) -> Option<Vec<u32>> {
        self.hashval_to_idx.clone()
    }

    pub fn signatures(&self) {
        let protein = String::from("protein");
        let hp = String::from("hp");
        let dayhoff = String::from("dayhoff");
        let moltype = self.moltype();
        let hash_function = match moltype {
            protein => HashFunctions::murmur64_protein,
            hp => HashFunctions::murmur64_hp,
            dayhoff => HashFunctions::murmur64_dayhoff,
        };
        // defaults for the rest of the parameters
        let seed = 42;
        let track_abundance = false;
        let max_hash = max_hash_for_scaled(self.scaled()).unwrap();
        // max_hash_for_scaled(scaled: u64) -> Option<u64>
        let mh = KmerMinHash::new(0, self.ksize, hash_function, 42, max_hash, track_abundance);
    }
}