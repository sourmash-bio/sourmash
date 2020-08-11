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

use std::io;
use std::fs::File;
use std::path::Path;
use std::io::BufReader;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

use crate::Error;
use crate::sketch::Sketch;
use crate::index::Index;
use crate::signature::Signature;
use crate::sketch::minhash::{KmerMinHash, HashFunctions, max_hash_for_scaled};
use crate::cmd::ComputeParameters;
use std::convert::TryFrom;


pub type LineagePairs = BTreeMap<String, String>; 

impl Serialize for LcaDB {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = 12;

        // TODO: MAKE THIS CONVERSION PRETTIER/FASTER
        // convert lineage in lid_to_lineage
        let mut lid_to_lineage_vec: HashMap<u32, Vec<Vec<String>>> = HashMap::new();
        for pair in self.lid_to_lineage.clone() {

            let mut lineage_vec_vec: Vec<Vec<String>> = [].to_vec();

            for pair2 in pair.1 {
                // pair.1 = BTreeMap<String, String>
                // pair2 = (String, String)
                let lineage_vec: Vec<String> = vec![pair2.0, pair2.1];
                lineage_vec_vec.push(lineage_vec);
            }

            lid_to_lineage_vec.insert(pair.0, lineage_vec_vec);
        }

        let mut partial = serializer.serialize_struct("LcaDB", n_fields)?;
        partial.serialize_field("version", "2.1")?;
        partial.serialize_field("type", "sourmash_lca")?;
        partial.serialize_field("license", "CC0")?;

        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("scaled", &self.scaled)?;
        // partial.serialize_field("moltype", &self.moltype.to_string())?;

        partial.serialize_field("lid_to_lineage", &lid_to_lineage_vec)?;
        partial.serialize_field("hashval_to_idx", &self.hashval_to_idx)?;
        partial.serialize_field("ident_to_name", &self.ident_to_name)?;
        partial.serialize_field("ident_to_idx", &self.ident_to_idx)?;
        partial.serialize_field("idx_to_lid", &self.idx_to_lid)?;

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
            version: String,
            r#type: String,
            license: String,

            ksize: u32,
            scaled: u64,
            // moltype: String,

            lid_to_lineage: HashMap<u32, Vec<Vec<String>>>,
            hashval_to_idx: HashMap<u64, Vec<u32>>,
            ident_to_name: HashMap<String, String>,
            ident_to_idx: HashMap<String, u32>,
            idx_to_lid: HashMap<u32, u32>,
        }

        let templcadb = TempLcaDB::deserialize(deserializer)?;

        // TODO: MAKE THIS CONVERSION PRETTIER/FASTER
        // convert lineage in lid_to_lineage
        let mut lid_to_lineage_map: HashMap<u32, LineagePairs> = HashMap::new();
        for pair in templcadb.lid_to_lineage {
            let lid = pair.0;
            let lineage_vec = pair.1;

            let mut lineage_map: LineagePairs = BTreeMap::new();
            for vec_pair in lineage_vec {
                // lineage_vec = Vec<Vec<String>>
                lineage_map.insert(vec_pair[0].clone(), vec_pair[1].clone());
            }
            lid_to_lineage_map.insert(lid, lineage_map);
        }

        // TODO: MAKE THIS CONVERSION PRETTIER/FASTER
        // invert lid_to_lineage
        let mut lineage_to_lid: HashMap<LineagePairs, u32> = HashMap::new();
        if lid_to_lineage_map.len() > 0 {
            for pair in &lid_to_lineage_map {
                lineage_to_lid.insert(pair.1.clone(), pair.0.clone());
            }
        }

        // let moltype = templcadb.moltype;
        let moltype = "DNA".to_string();

        let _next_index = templcadb.ident_to_idx.len() as u32;
        let _next_lid = templcadb.idx_to_lid.len() as u32;

        Ok(LcaDB {
            ksize: templcadb.ksize,
            scaled: templcadb.scaled,
            filename: "".to_string(),
            moltype: moltype,
            _next_index: _next_index,
            _next_lid: _next_lid,
            ident_to_name: templcadb.ident_to_name,
            ident_to_idx: templcadb.ident_to_idx,
            idx_to_lid: templcadb.idx_to_lid,
            lineage_to_lid: lineage_to_lid,
            lid_to_lineage: lid_to_lineage_map,
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

    pub fn idx_to_ident(&self) -> HashMap<u32, String> {
        let mut idx_to_ident: HashMap<u32, String> = HashMap::new();
        for (ident, idx) in &self.ident_to_idx {
            if !idx_to_ident.contains_key(&idx) {
                idx_to_ident.insert(*idx, ident.to_string());
            }
        }
        return idx_to_ident;
    }

    pub fn select(&self, ksize: Option<u32>, moltype: Option<&str>) -> Result<(), Error> {
        if ksize != None {
            if self.ksize != ksize.unwrap() {
                return Err(Error::MismatchKSizes);
            }
        } 

        if moltype != None {
            if self.moltype != moltype.unwrap().to_string() {
                // TODO: ask if this is the right error code
                return Err(Error::MismatchDNAProt);
            }
        } 

        Ok(())
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

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        // let filename = str::from_utf8(filename);

        let file = File::create(path)?;
        serde_json::to_writer(file, self)?;
        Ok(())
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<LcaDB, Error> {
        let mut reader = io::BufReader::new(File::open(path)?);
        let lca_db = serde_json::from_reader(reader)?;
        // TODO: set self.filename from end of path
        Ok(lca_db)
    }

    // Return all of the signatures in this LCA database.
    pub fn signatures(&self) -> Result<Vec<Signature>, Error>{
        let mut sigs: Vec<Signature> = Vec::new();
        for v in self._signatures().values() {
            sigs.push(v.clone());
        }
        Ok(sigs)
    }

    pub fn _signatures(&self) -> HashMap<u32, Signature> {

        // let hash_function: hash_function = match self.moltype.as_str() {
        //     "protein" => HashFunctions::murmur64_protein,
        //     "hp" => HashFunctions::murmur64_hp,
        //     "dayhoff" => HashFunctions::murmur64_dayhoff,
        // };
        let hash_function = HashFunctions::try_from(self.moltype.as_str()).unwrap();

        // defaults for the rest of the parameters
        let seed = 42;
        let track_abundance = false;
        let max_hash = max_hash_for_scaled(self.scaled()).unwrap();

        let mh = KmerMinHash::new(0, self.ksize, hash_function, 42, max_hash, track_abundance);

        println!("creating signatures for LCA DB...");
        let mut mhd: HashMap<u32, KmerMinHash> = HashMap::new();
        let mut temp_vals: HashMap<u32, Vec<u64>> = HashMap::new();

        // invert hashval_to_idx to temp_vals (idx to hashval)
        for (hashval, idlist) in &self.hashval_to_idx {
            for idx in idlist {
                // set temp_vals and get temp_hashes
                let mut temp_hashes: &mut Vec<u64> = match temp_vals.get_mut(idx) {
                    Some(s) => {
                        s.push(*hashval);
                        s
                    },
                    None => {
                        temp_vals.insert(*idx, vec![*hashval]);
                        temp_vals.get_mut(idx).unwrap()
                    },
                };

                // set minhash and/or add hashes
                if temp_hashes.len() > 50 {
                    match mhd.get_mut(idx) {
                        Some(s) => s.add_many(temp_hashes).unwrap(),
                        None => {
                            let mut mh_temp = mh.copy_and_clear().unwrap();
                            mh_temp.add_many(temp_hashes).unwrap();
                            mhd.insert(*idx, mh_temp);
                        },
                    }

                    temp_hashes = &mut vec![];
                }
            }
        }
        // loop temp_vals again to add remaining hashes
        // (each list of hashes is smaller than 50 items)
        for (idx, hashes) in temp_vals {
            match mhd.get_mut(&idx) {
                Some(s) => s.add_many(&hashes).unwrap(),
                None => {
                    let mut mh_temp = mh.copy_and_clear().unwrap();
                    mh_temp.add_many(&hashes).unwrap();
                    mhd.insert(idx, mh_temp);
                },
            }
        }

        let mut sigd: HashMap<u32, Signature> = HashMap::new();
        for (idx, mh) in mhd {
            // get name from idx
            let ident = &self.idx_to_ident()[&idx];
            let name = &self.ident_to_name[ident];

            // make sig and set name
            let mut sig = Signature::default();
            sig.set_name(name.as_str());
            sig.push(Sketch::MinHash(mh.clone()));

            // set sig in HashMap
            if sigd.contains_key(&idx) {
                *sigd.get_mut(&idx).unwrap() = sig;
            } else {
                sigd.insert(idx, sig);
            }
        }

        println!("=> {} signatures!", sigd.len());
        sigd
    }

    // Get a list of lineages for this hashval.
    pub fn get_lineage_assignments(&self, hashval: u64) -> Result<Vec<LineagePairs>, Error> {
        let mut x: Vec<LineagePairs> = Vec::new();

        let idx_list = match self.hashval_to_idx.get(&hashval) {
            Some(s) => s.to_vec(),
            None => Vec::new(),
        };
        for idx in idx_list {
            let lid = self.idx_to_lid.get(&idx);

            if lid != None {
                let lineage = &self.lid_to_lineage[&lid.unwrap()];
                x.push(lineage.clone());
            }
        }
        Ok(x)
    }
}



// TESTINGGGGGG!
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
    fn lca_signatures() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/lca/delmont-1.lca.json");

        let mut lca_db = LcaDB::load(filename.as_path()).unwrap();

        // get signature to add
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10+s11.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let signatures: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        // construct lineage
        let mut lineage: LineagePairs = BTreeMap::new();
        lineage.insert("name1".to_string(), "rank1".to_string()); 
        lineage.insert("name2".to_string(), "rank2".to_string());
        
        // add new sigs and lineage
        lca_db.insert(&signatures[0], b"erik", &lineage).unwrap();

        let sigs = lca_db.signatures().unwrap();
        dbg!(&sigs);

        // should be 3. 2 from the json file and 1 added using .insert(...)
        assert!(sigs.len() == 3);
    }

    #[test]
    fn lca_get_lineage_assignments() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/lca/both.lca.json");

        let lca_db = LcaDB::load(filename.as_path()).unwrap();

        let hashval = 270342483362557;

        let lineages = lca_db.get_lineage_assignments(hashval).unwrap();

        dbg!(&lineages);
        assert!(lineages.len() == 1);
    }

    #[test]
    fn lca_select() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/lca/delmont-1.lca.json");

        let lca_db = LcaDB::load(filename.as_path()).unwrap();

        // actual values: ksize = 31 and moltype = "DNA"

        lca_db.select(Some(31), Some("DNA")).unwrap();
        lca_db.select(Some(31), Some("protein")).unwrap_err();
        lca_db.select(Some(32), Some("DNA")).unwrap_err();
        lca_db.select(Some(32), Some("protein")).unwrap_err();
    }

    #[test]
    fn lca_roundtrip() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/lca/delmont-1.lca.json");

        let lcadb = LcaDB::load(filename.as_path()).unwrap();

        let mut tmpfile = tempfile::NamedTempFile::new().unwrap();
        lcadb.save(tmpfile.path()).unwrap();

        tmpfile.seek(SeekFrom::Start(0)).unwrap();

        let lcadb_2 = LcaDB::load(tmpfile.path()).unwrap();

        assert_eq!(lcadb.ksize, lcadb_2.ksize);
        assert_eq!(lcadb.scaled, lcadb_2.scaled);
        assert_eq!(lcadb.moltype, lcadb_2.moltype);

        assert_eq!(lcadb._next_index, lcadb_2._next_index);
        assert_eq!(lcadb._next_lid, lcadb_2._next_lid);
        assert_eq!(lcadb.ident_to_name, lcadb_2.ident_to_name);
        assert_eq!(lcadb.ident_to_idx, lcadb_2.ident_to_idx);
        assert_eq!(lcadb.idx_to_lid, lcadb_2.idx_to_lid);
        assert_eq!(lcadb.lineage_to_lid, lcadb_2.lineage_to_lid);
        assert_eq!(lcadb.lid_to_lineage, lcadb_2.lid_to_lineage);
        assert_eq!(lcadb.hashval_to_idx, lcadb_2.hashval_to_idx);
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

        // construct lineage
        let mut lineage: LineagePairs = BTreeMap::new();
        lineage.insert("name1".to_string(), "rank1".to_string()); 
        lineage.insert("name2".to_string(), "rank2".to_string());
        
        lca_db.insert(&sigs[0], b"erik", &lineage).unwrap();

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
}