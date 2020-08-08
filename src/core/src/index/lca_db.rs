use std::str;
use std::ffi::CStr;
use std::collections::{HashMap, BTreeMap};
use std::os::raw::c_char;

#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::de::Deserializer;
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

use crate::sketch::Sketch;
use crate::index::Index;
use crate::signature::Signature;
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
        partial.serialize_field("_next_index", &self._next_index)?;
        partial.serialize_field("_next_lid", &self._next_lid)?;
        partial.serialize_field("ident_to_name", &self.ident_to_name)?;
        partial.serialize_field("ident_to_idx", &self.ident_to_idx)?;
        partial.serialize_field("idx_to_lid", &self.idx_to_lid)?;
        partial.serialize_field("lineage_to_lid", &self._next_lid)?;
        partial.serialize_field("lid_to_lineage", &self.lid_to_lineage)?;
        partial.serialize_field("hashval_to_idx", &self.hashval_to_idx)?;

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
            ident_to_name: HashMap<String, String>,
            ident_to_idx: HashMap<String, u32>,
            idx_to_lid: HashMap<u32, u32>,
            lineage_to_lid: HashMap<LineagePairs, u32>,
            lid_to_lineage: HashMap<u32, LineagePairs>,
            hashval_to_idx: HashMap<u64, Vec<u32>>,
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

pub type LineagePairs = BTreeMap<String, String>; 

// impl<'de> Deserialize<'de> for LineagePairs {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//     where
//         D: Deserializer<'de>,
//     {
//         #[derive(Deserialize)]
//         type TempLineagePairs = BTreeMap<String, String>

//         let templineagepairs = TempLineagePairs::deserialize(deserializer)?;

//         Ok(LcaDB {
//             ksize: templcadb.ksize,
//             scaled: templcadb.scaled,
//             filename: templcadb.filename,
//             moltype: templcadb.moltype,
//             _next_index: templcadb._next_index,
//             _next_lid: templcadb._next_lid,
//             ident_to_name: templcadb.ident_to_name,
//             ident_to_idx: templcadb.ident_to_idx,
//             idx_to_lid: templcadb.idx_to_lid,
//             lineage_to_lid: templcadb.lineage_to_lid,
//             lid_to_lineage: templcadb.lid_to_lineage,
//             hashval_to_idx: templcadb.hashval_to_idx,
//         })
//     }
// }

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

    #[builder(default = HashMap::new())]
    ident_to_name: HashMap<String, String>,

    #[builder(default = HashMap::new())]
    ident_to_idx: HashMap<String, u32>,

    #[builder(default = HashMap::new())]
    idx_to_lid: HashMap<u32, u32>,

    #[builder(default = HashMap::new())]
    lineage_to_lid: HashMap<LineagePairs, u32>,

    #[builder(default = HashMap::new())]
    lid_to_lineage: HashMap<u32, LineagePairs>,

    #[builder(default = HashMap::new())]
    hashval_to_idx: HashMap<u64, Vec<u32>>, //HashSet as keys
}

impl Default for LcaDB {
    fn default() -> LcaDB {
        LcaDB {
            ksize: 32,
            scaled: 1,
            filename: "filename".to_string(),
            moltype: "DNA".to_string(),

            _next_index: 0,
            _next_lid: 0,
            ident_to_name: HashMap::new(),
            ident_to_idx: HashMap::new(),
            idx_to_lid: HashMap::new(),
            lineage_to_lid: HashMap::new(),
            lid_to_lineage: HashMap::new(),
            hashval_to_idx: HashMap::new(),
        }
    }
}

impl LcaDB {
    pub fn new() -> LcaDB {

        let ksize = 32;
        let scaled = 1;
        let filename = "filename".to_string();
        let moltype = "DNA".to_string();
        
        let _next_index = 0;
        let _next_lid = 0;
        let ident_to_name = HashMap::new();
        let ident_to_idx = HashMap::new();
        let idx_to_lid = HashMap::new();
        let lineage_to_lid = HashMap::new();
        let lid_to_lineage = HashMap::new();
        let hashval_to_idx = HashMap::new();

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

    pub fn ident_to_name(&self) -> HashMap<String, String> {
        self.ident_to_name.clone()
    }

    pub fn ident_to_idx(&self) -> HashMap<String, u32> {
        self.ident_to_idx.clone()
    }

    pub fn idx_to_lid(&self) -> HashMap<u32, u32> {
        self.idx_to_lid.clone()
    }

    pub fn lineage_to_lid(&self) -> HashMap<LineagePairs, u32> {
        self.lineage_to_lid.clone()
    }

    pub fn lid_to_lineage(&self) -> HashMap<u32, LineagePairs> {
        self.lid_to_lineage.clone()
    }

    pub fn hashval_to_idx(&self) -> HashMap<u64, Vec<u32>> {
        self.hashval_to_idx.clone()
    }


    // insert functions (possibly not neccessary but idk yet)
    pub fn insert_ident_to_name(&mut self, ident: String, name: String) {
        self.ident_to_name.insert(ident, name);
    }

    pub fn insert_ident_to_idx(&mut self, ident: String, idx: u32) {
        self.ident_to_idx.insert(ident, idx);
    }

    pub fn insert_idx_to_lid(&mut self, idx: u32, lid: u32) {
        self.idx_to_lid.insert(idx, lid);
    }

    pub fn insert_lineage_to_lid(&mut self, lineage: LineagePairs, lid: u32) {
        self.lineage_to_lid.insert(lineage, lid);
    }

    pub fn insert_lid_to_lineage(&mut self, lid: u32, lineage: LineagePairs) {
        self.lid_to_lineage.insert(lid, lineage);
    }

    pub fn insert_hashval_to_idx(&mut self, hashval: u64, idx: Vec<u32>) {
        self.hashval_to_idx.insert(hashval, idx);
    }

    pub fn insert_hashval_to_idx_vec(&mut self, hashval: u64, idx: u32) {
        if let Some(idx_vec) = self.hashval_to_idx.get_mut(&hashval) {
            idx_vec.push(idx);
        }
    }

    pub fn _get_ident_index(&mut self, ident: &String, fail_on_duplicate: bool) -> u32 {
        // Get (create if nec) a unique int id, idx, for each identifier.
        if fail_on_duplicate {
            assert!(!self.ident_to_idx.contains_key(ident));    // should be no duplicate identities
        }

        let idx = match self.ident_to_idx.get(ident) {
            Some(i) => *i,
            None => {
                let i = self._next_index;
                self._next_index += 1;
                self.ident_to_idx.insert(ident.to_string(), i);
                i
            },
        };
        idx
    }

    pub fn _get_lineage_id(&mut self, lineage: &LineagePairs) -> u32 {
        // "Get (create if nec) a unique lineage ID for each LineagePair tuples."
        
        // does one exist already?
        // nope - create one. Increment next_lid.
        if !self.lineage_to_lid.contains_key(&lineage) {
            let lid = self._next_lid;
            self._next_lid += 1;

            // build mappings
            self.lineage_to_lid.insert(lineage.clone(), lid);
            self.lid_to_lineage.insert(lid, lineage.clone());
            lid
        } else {
            let lid = self.lineage_to_lid.get(&lineage);
            *lid.unwrap()
        }
    }

    pub fn insert(&mut self, sig: &Signature, ident_opt: &[u8], lineage: &LineagePairs) -> Result<u32, Error> {
        // set ident
        let ident = if ident_opt != b"" {
            String::from_utf8(ident_opt.to_vec()).unwrap()
        } else {
            sig.name()
        };

        if let Sketch::MinHash(minhash) = &sig.signatures[0] {
    
            // check for errors
                
            assert!(!self.ident_to_name().contains_key(&ident));
                //error: raise ValueError("signature {} is already in this LCA db.".format(ident))
            
            // implement self.cache property and _invalidate_cache() method
            // self._invalidate_cache()
                
            // downsample to specified scaled; this has the side effect of
            // making sure they're all at the same scaled value!
            let minhash = minhash.downsample_scaled(self.scaled()).unwrap();
                // error if it is a scaled signature
    
            // store name
            self.ident_to_name.insert(ident.clone(), sig.name());
    
            // identifier -> integer index (idx)
            let idx = self._get_ident_index(&ident, true);
    
            if lineage.len() > 0 {
                // (LineagePairs*) -> integer lineage ids (lids)
                let lid = self._get_lineage_id(lineage);
    
                // map idx to lid as well.
                if !self.idx_to_lid().contains_key(&idx) {
                    self.idx_to_lid.insert(idx, lid);
                }
            }
    
            // append idx to each hashval's idx vector
            for hashval in minhash.mins() {
                if self.hashval_to_idx().contains_key(&hashval) {
                    self.hashval_to_idx.get_mut(&hashval).unwrap().push(idx);
                } else {
                    let mut idx_vec: Vec<u32> = Vec::with_capacity(1);
                    idx_vec.push(idx);
                    self.hashval_to_idx.insert(hashval, idx_vec);
                }
            }
        
            Ok(minhash.num())
        }
        else {
            unimplemented!()
        }
    }

    pub fn save(&self, filename: &[u8]) {
        unimplemented!()
    }

    pub fn load(filename: &[u8]) {
        unimplemented!()
    }

    // pub fn _signatures(&self) {
    //     let protein = String::from("protein");
    //     let hp = String::from("hp");
    //     let dayhoff = String::from("dayhoff");
    //     let moltype = self.moltype();

    //     let hash_function = match moltype {
    //         protein => HashFunctions::murmur64_protein,
    //         hp => HashFunctions::murmur64_hp,
    //         dayhoff => HashFunctions::murmur64_dayhoff,
    //     };
    //     // defaults for the rest of the parameters
    //     let seed = 42;
    //     let track_abundance = false;
    //     let max_hash = max_hash_for_scaled(self.scaled()).unwrap();

    //     let mh = KmerMinHash::new(0, self.ksize, hash_function, 42, max_hash, track_abundance);
    // }
}

// impl<'a> Index<'a> for LcaDB {
//     type Item = Signature;

//     fn insert(
//         &mut self, 
//         node: Self::Item, 
//     ) -> Result<(), Error> {
//         unimplemented!()
//     }

//     fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
//         unimplemented!()
//     }

//     fn load<P: AsRef<Path>>(path: P) -> Result<(), Error> {
//         unimplemented!()
//     }

//     fn signatures(&self) -> Vec<Self::Item> {
//         unimplemented!()
//     }

//     fn signature_refs(&self) -> Vec<&Self::Item> {
//         unimplemented!()
//     }
// }



// TESTING FUNCTIONS
#[cfg(test)]                                                                                   
mod test { 
    use crate::cmd::ComputeParameters;
    use std::fs::File;
    use std::io::{BufReader, Read, Seek, SeekFrom};
    use std::path::PathBuf;
    use std::collections::{HashMap, BTreeMap};

    use crate::signature::Signature;
    use crate::index::lca_db::{LineagePairs, LcaDB};

    
        #[test]
        fn lca_roundtrip() {
            let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
            filename.push("../../tests/test-data/lca/delmont-1.lca.json");
    
            // let lcadb = LcaDB::load(filename).unwrap();
    
            // let mut tmpfile = tempfile::NamedTempFile::new().unwrap();
            // lcadb.save(tmpfile.path()).unwrap();
    
            // tmpfile.seek(SeekFrom::Start(0)).unwrap();
    
            // let lcadb_2 = LcaDB::load(tmpfile.path()).unwrap();
    
            // assert_eq!(lcadb.ksize, lcadb_2.ksize);
            // assert_eq!(lcadb.scaled, lcadb_2.scaled);
            // assert_eq!(lcadb.moltype, lcadb_2.moltype);
            // assert_eq!(lcadb.idx_to_lid, lcadb_2.idx_to_lid);
            // assert_eq!(lcadb.hashval_to_idx, lcadb_2.hashval_to_idx);
        }
    #[test]
    fn build_default_struct() {
        let lca_db = LcaDB::new();

        assert!(lca_db.ksize() == 32);
        assert!(lca_db.scaled() == 1);
        assert!(lca_db.filename() == "filename".to_string());
        assert!(lca_db.moltype() == "DNA".to_string());
        assert!(lca_db._next_index() as u32 == 0);
        assert!(lca_db._next_lid() as u32 == 0);
        assert!(lca_db.ident_to_name() == HashMap::new());
        assert!(lca_db.ident_to_idx() == HashMap::new());
        assert!(lca_db.idx_to_lid() == HashMap::new());
        assert!(lca_db.lineage_to_lid() == HashMap::new());
        assert!(lca_db.lid_to_lineage() == HashMap::new());
        assert!(lca_db.hashval_to_idx() == HashMap::new());
    }

    #[test]
    fn test_insert() {
        let mut lca_db = LcaDB::new();

        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10+s11.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        let mut lineage: LineagePairs = BTreeMap::new();
        lineage.insert("name1".to_string(), "rank1".to_string()); 
        lineage.insert("name2".to_string(), "rank2".to_string());
        
        lca_db.insert(&sigs[0], b"erik", &lineage);

        println!("{:?}", lca_db);
                                            
        assert!(lca_db.ksize() == 32);
        assert!(lca_db.scaled() == 1);
        assert!(lca_db.filename() == "filename".to_string());
        assert!(lca_db.moltype() == "DNA".to_string());
        assert!(lca_db._next_index() as u32 == 1);
        assert!(lca_db._next_lid() as u32 == 1);

        let mut ident_to_name2 = HashMap::with_capacity(1);
        ident_to_name2.insert("erik".to_string(), "s10+s11".to_string());
        assert!(lca_db.ident_to_name() == ident_to_name2);

        let mut lid_to_lineage2 = HashMap::with_capacity(1);
        lid_to_lineage2.insert(0, lineage);
        assert!(lca_db.lid_to_lineage() == lid_to_lineage2);
    }

    // #[test]
    // fn test_save() {
    //     let mut lca_db2 = LcaDB::new(32, 1, b"", b"DNA");
    //     let mh = KmerMinHash::new(0, 32, HashFunctions::murmur64_protein, 42, 10000, false);
    //     let sig = Signature::default();
    //     let lineage: Vec<LineagePair> = vec![LineagePair::new("name1".to_string(), "rank1".to_string()), 
    //                                         LineagePair::new("name2".to_string(), "rank2".to_string())];
                                            
    //     let lca_db = lca_db2.save(&mh, b"eriksjson");

                                            
    //     assert!(false);
    //     assert!(lca_db.ksize() == 32);
    //     assert!(lca_db.scaled() == 1);
    //     assert!(lca_db.filename() == "");
    //     assert!(lca_db.moltype() == "DNA");
    //     assert!(lca_db._next_index() as u32 == 0);
    //     assert!(lca_db._next_lid() as u32 == 0);
    //     assert!(lca_db.ident_to_name() == HashMap::new());
    //     assert!(lca_db.ident_to_idx() == HashMap::new());
    //     assert!(lca_db.idx_to_lid() == HashMap::new());
    //     assert!(lca_db.lineage_to_lid() == HashMap::new());
    //     assert!(lca_db.lid_to_lineage() == HashMap::new());
    //     assert!(lca_db.hashval_to_idx() == HashMap::new());
    // }
}