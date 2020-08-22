use std::str;
use std::collections::{HashMap, BTreeMap};

#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::de::Deserializer;
use serde::ser::{SerializeStruct, Serializer};
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use std::io;
use std::fs::File;
use std::path::Path;
use std::ffi::OsStr;
use log::debug;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

use crate::Error;
use crate::sketch::Sketch;
use crate::signature::Signature;
use crate::sketch::minhash::{KmerMinHash, HashFunctions, max_hash_for_scaled, scaled_for_max_hash};
use std::convert::TryFrom;
use std::cmp::Ordering::Equal;


pub type Lineage = BTreeMap<String, String>; 

pub fn lineage_to_vec(lineage: Lineage) -> Vec<Vec<String>> {
    let mut lineage_vec: Vec<Vec<String>> = Vec::new();
        if lineage.contains_key(&"superkingdom".to_string()) {
            lineage_vec.push(vec!["superkingdom".to_string(), lineage.get("superkingdom").unwrap().to_string()]);
        }
        if lineage.contains_key(&"phylum".to_string()) {
            lineage_vec.push(vec!["phylum".to_string(), lineage.get("phylum").unwrap().to_string()]);
        }
        if lineage.contains_key(&"class".to_string()) {
            lineage_vec.push(vec!["class".to_string(), lineage.get("class").unwrap().to_string()]);
        }
        if lineage.contains_key(&"order".to_string()) {
            lineage_vec.push(vec!["order".to_string(), lineage.get("order").unwrap().to_string()]);
        }
        if lineage.contains_key(&"family".to_string()) {
            lineage_vec.push(vec!["family".to_string(), lineage.get("family").unwrap().to_string()]);
        }
        if lineage.contains_key(&"genus".to_string()) {
            lineage_vec.push(vec!["genus".to_string(), lineage.get("genus").unwrap().to_string()]);
        }
        if lineage.contains_key(&"species".to_string()) {
            lineage_vec.push(vec!["species".to_string(), lineage.get("species").unwrap().to_string()]);
        }
        if lineage.contains_key(&"strain".to_string()) {
            lineage_vec.push(vec!["strain".to_string(), lineage.get("strain").unwrap().to_string()]);
        }

        // for tests that dont have normal lineages
        if lineage_vec.len() == 0 && lineage.len() != 0 {
            for (rank, name) in lineage {
                lineage_vec.push(vec![rank, name]);
            }
        }
        lineage_vec

}

impl Serialize for LcaDB {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = 11;

        // TODO: MAKE THIS CONVERSION PRETTIER/FASTER
        // convert lineage in lid_to_lineage
        let mut lid_to_lineage_vec: HashMap<u32, Vec<Vec<String>>> = HashMap::new();
        for (lid, lineage) in self.lid_to_lineage.clone() {
            let lineage_vec = lineage_to_vec(lineage);
            lid_to_lineage_vec.insert(lid, lineage_vec);
        }

        let mut partial = serializer.serialize_struct("LcaDB", n_fields)?;
        partial.serialize_field("version", "2.1")?;
        partial.serialize_field("type", "sourmash_lca")?;
        partial.serialize_field("license", "CC0")?;

        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("scaled", &self.scaled)?;
        partial.serialize_field("moltype", &self.moltype.to_string())?;

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
        #[derive(TypedBuilder, Deserialize)]
        struct TempLcaDB {
            #[serde(default)]
            version: String,

            #[serde(default)]
            #[serde(rename(deserialize = "type"))]
            class: String,

            #[serde(default)]
            license: String,

            ksize: u32,
            scaled: u64,

            #[serde(default)]
            moltype: String,

            #[serde(default = "missing_field")]
            lid_to_lineage: HashMap<u32, Vec<Vec<String>>>,

            hashval_to_idx: HashMap<u64, Vec<u32>>,
            ident_to_name: HashMap<String, String>,
            idx_to_lid: HashMap<u32, u32>,
            ident_to_idx: HashMap<String, u32>,
        }

        let templcadb = TempLcaDB::deserialize(deserializer)?;

        // error handling
        if templcadb.class != "sourmash_lca" {
            panic!("database file is not an LCA db.")
        }

        // TODO: should raise error also if lid_to_lineage isnt parsed
        let version = templcadb.version.parse::<f32>().unwrap();
        if version < 2.0 {
            panic!("Error! This is an old-style LCA DB. You'll need to rebuild or download a newer one.")
        }

        // TODO: MAKE THIS CONVERSION PRETTIER/FASTER
        // convert lineage in lid_to_lineage
        let mut lid_to_lineage_map: HashMap<u32, Lineage> = HashMap::new();
        for pair in templcadb.lid_to_lineage {
            let lid = pair.0;
            let lineage_vec = pair.1;

            let mut lineage_map: Lineage = BTreeMap::new();
            for vec_pair in lineage_vec.into_iter().rev() {
                // lineage_vec = Vec<Vec<String>>
                lineage_map.insert(vec_pair[0].clone(), vec_pair[1].clone());
            }
            lid_to_lineage_map.insert(lid, lineage_map);
        }

        // TODO: MAKE THIS CONVERSION PRETTIER/FASTER
        // invert lid_to_lineage
        let mut lineage_to_lid: HashMap<Lineage, u32> = HashMap::new();
        if lid_to_lineage_map.len() > 0 {
            for pair in &lid_to_lineage_map {
                lineage_to_lid.insert(pair.1.clone(), pair.0.clone());
            }
        }

        // set moltype or default to "DNA"
        let mut moltype = "DNA".to_string();
        if templcadb.moltype != "" {
            moltype = templcadb.moltype.to_string();
        }

        let _next_index = templcadb.ident_to_idx.len() as u32;
        let _next_lid = templcadb.idx_to_lid.len() as u32;

        Ok(LcaDB {
            version: templcadb.version.to_string(),
            class: templcadb.class.to_string(),
            license: templcadb.license.to_string(),
            ksize: templcadb.ksize,
            scaled: templcadb.scaled,
            filename: "".to_string(),
            moltype: moltype.to_string(),
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

pub fn missing_field() -> HashMap<u32, Vec<Vec<String>>> {
    panic!("Error! This is an old-style LCA DB. You'll need to rebuild or download a newer one.");
}

#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(Debug, TypedBuilder)]
pub struct LcaDB {
    version: String,
    class: String,
    license: String,

    ksize: u32,
    scaled: u64,
    moltype: String,
    filename: String,

    _next_index: u32,
    _next_lid: u32,
    ident_to_name: HashMap<String, String>,
    ident_to_idx: HashMap<String, u32>,
    idx_to_lid: HashMap<u32, u32>,
    lineage_to_lid: HashMap<Lineage, u32>,
    lid_to_lineage: HashMap<u32, Lineage>,
    hashval_to_idx: HashMap<u64, Vec<u32>>,
}

impl Default for LcaDB {
    fn default() -> LcaDB {
        LcaDB {
            version: "2.1".to_string(),
            class: "sourmash_lca".to_string(),
            license: "CC0".to_string(),

            ksize: 32,
            scaled: 1,
            filename: "".to_string(),
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
        let version = "2.1".to_string();
        let class = "sourmash_lca".to_string();
        let license = "CC0".to_string();

        let ksize = 32;
        let scaled = 1;
        let filename = "".to_string();
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
            version,
            class,
            license,
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

    pub fn new_with_params(ksize: u32, scaled: u64, filename: &str, moltype: &str) -> LcaDB {
        let version = "2.1".to_string();
        let class = "sourmash_lca".to_string();
        let license = "CC0".to_string();

        let filename = filename.to_string();
        let moltype = moltype.to_string();

        let _next_index = 0;
        let _next_lid = 0;
        let ident_to_name = HashMap::new();
        let ident_to_idx = HashMap::new();
        let idx_to_lid = HashMap::new();
        let lineage_to_lid = HashMap::new();
        let lid_to_lineage = HashMap::new();
        let hashval_to_idx = HashMap::new();

        LcaDB {
            version,
            class,
            license,
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

    pub fn lineage_to_lid(&self) -> HashMap<Lineage, u32> {
        self.lineage_to_lid.clone()
    }

    pub fn lid_to_lineage(&self) -> HashMap<u32, Lineage> {
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

    pub fn _get_lineage_id(&mut self, lineage: &Lineage) -> u32 {
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

    pub fn insert(&mut self, sig: &Signature, ident_opt: Option<&str>, lineage_opt: Option<Lineage>) -> Result<u32, Error> {
        // set ident
        let ident = match ident_opt {
            Some(s) => s.to_string(),
            None => sig.name(),
        };

        if let Sketch::MinHash(minhash) = &sig.signatures[0] {
    
            // check for errors
                
            if self.ident_to_name.contains_key(&ident) {
                return Err(Error::DuplicateSignature { ident });
            }
            
            // implement self.cache property and _invalidate_cache() method
            // self._invalidate_cache()
                
            // downsample to specified scaled; this has the side effect of
            // making sure they're all at the same scaled value!
            let minhash = minhash.downsample_scaled(self.scaled).unwrap();
    
            // store name
            self.ident_to_name.insert(ident.clone(), sig.name());
    
            // identifier -> integer index (idx)
            let idx = self._get_ident_index(&ident, true);
    
            if let Some(lineage) = lineage_opt {
                // (Lineage*) -> integer lineage ids (lids)
                let lid = self._get_lineage_id(&lineage);
    
                // map idx to lid as well.
                self.idx_to_lid.insert(idx, lid);   
            }
    
            // append idx to each hashval's idx vector
            for hashval in minhash.mins() {
                if self.hashval_to_idx.contains_key(&hashval) {
                    self.hashval_to_idx.get_mut(&hashval).unwrap().push(idx);
                } else {
                    self.hashval_to_idx.insert(hashval, vec![idx]);
                }
            }
        
            Ok(minhash.mins().len() as u32)
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

    pub fn load<P: AsRef<Path>>(path: P, filename: &str) -> Result<LcaDB, Error> {
        let reader = io::BufReader::new(File::open(&path)?);
        let mut lca_db: LcaDB = LcaDB::from_reader(reader)?;

        lca_db.filename = filename.to_string();

        Ok(lca_db)
    }

    pub fn load_db<R>(
        buf: R,
        filename: &str
    ) -> Result<LcaDB, Error>
    where
        R: io::Read,
    { 
        let mut lca_db = LcaDB::from_reader(buf).unwrap();

        lca_db.filename = filename.to_string();

        Ok(lca_db)
    }

    pub fn from_reader<R>(rdr: R) -> Result<LcaDB, Error>
    where
        R: io::Read,
    {
        let (rdr, _format) = niffler::get_reader(Box::new(rdr))?;

        let dbs: LcaDB = serde_json::from_reader(rdr)?;
        Ok(dbs)
    }

    // Return all of the signatures in this LCA database.
    pub fn signatures(&self) -> Result<Vec<Signature>, Error> {
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
        let max_hash = max_hash_for_scaled(self.scaled).unwrap();

        let mh = KmerMinHash::new(0, self.ksize, hash_function, seed, max_hash, track_abundance);

        debug!("creating signatures for LCA DB...");
        let mut mhd: HashMap<u32, KmerMinHash> = HashMap::new();
        let mut temp_vals: HashMap<u32, Vec<u64>> = HashMap::new();

        // invert hashval_to_idx to temp_vals (idx to hashval)
        for (hashval, idlist) in &self.hashval_to_idx {
            for idx in idlist {
                // set temp_vals and get temp_hashes
                let temp_hashes: &mut Vec<u64> = match temp_vals.get_mut(idx) {
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

                    *temp_vals.get_mut(idx).unwrap() = vec![];
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
            sigd.insert(idx, sig);
        }

        debug!("=> {} signatures!", sigd.len());
        sigd
    }

    // Get a list of lineages for this hashval.
    pub fn get_lineage_assignments(&self, hashval: u64) -> Result<Vec<Lineage>, Error> {
        let mut x: Vec<Lineage> = Vec::new();

        let idx_list = match self.hashval_to_idx.get(&hashval) {
            Some(s) => s.to_vec(),
            None => Vec::new(),
        };
        
        for idx in idx_list {
            // dbg!(&self.idx_to_lid);
            let lid = self.idx_to_lid.get(&idx);

            if lid != None {
                let lineage = &self.lid_to_lineage[&lid.unwrap()];
                x.push(lineage.clone());
            }
        }
        Ok(x)
    }

    // default: containment = false, ignore_scaled = false, return_one = false (used for gather(...))
    pub fn _find_signatures(
        &self, 
        mh: &KmerMinHash, 
        threshold: f32, 
        containment: bool, 
        ignore_scaled: bool,
        return_one: bool,
    ) -> Result<Vec<(f32, Signature, String)>, Error> {
        // make sure we're looking at the same scaled value as database
        let mh_scaled = scaled_for_max_hash(mh.max_hash());

        if self.scaled > mh_scaled {
            let mh = mh.downsample_scaled(self.scaled).unwrap();

        } else if self.scaled < mh_scaled && !ignore_scaled {
            // note that containment can be calculated w/o matching scaled.
            return Err(Error::MismatchScaled);
            // panic!("lca db scaled is {} vs query {}; must downsample", self.scaled, mh_scaled);
        }

        // collect matching hashes for the query
        let mut counter = HashMap::new();
        let query_mins = mh.mins();


        for hashval in &query_mins {
            let idlist = match self.hashval_to_idx.get(hashval) {
                Some(s) => s.to_vec(),
                None => Vec::new(),
            };
            for idx in idlist {
                *counter.entry(idx).or_insert(0) += 1;
            }
        }

        debug!("number of matching signatures for hashes: {}", &counter.len());
        let mut most_common: Vec<(&u32, &u32)> = Vec::new();
        for (idx, count) in &counter {
            most_common.push((idx, count));
        }

        // sort by most common
        most_common.sort_by(|a, b| b.1.cmp(a.1));

        let mut result: Vec<(f32, Signature, String)> = Vec::new();
        let sigmap = self._signatures();
        let mut score: f32;

        for (idx, count) in most_common {
            let match_sig = sigmap.get(&idx).unwrap().clone();
            if let Sketch::MinHash(match_mh) = &match_sig.signatures[0] {
                let match_size = match_mh.mins().len();

                // calculate the containment or similarity
                if containment {
                    score = *count as f32 / query_mins.len() as f32;
                } else {
                    // query_mins is size of query signature
                    // match_size is size of match signature
                    // count is overlap
                    score = *count as f32 / (query_mins.len() as f32 + match_size as f32 - *count as f32);
                }

                if score >= threshold {
                    result.push((score, match_sig, self.filename.clone()));
                }

                // replace this by using a generator
                if return_one {
                    break;
                }
            } else {
                unimplemented!()
            }
        }
        Ok(result)
    }

    pub fn downsample_scaled(&mut self, scaled: u64) -> Result<(), Error> {
        if scaled == self.scaled {
            return Ok(());
        } else if scaled < self.scaled {
            // return Err(Error::MismatchScaled);
            panic!("cannot decrease scaled from {} to {}", self.scaled, scaled);
        }

        // self._invalidate_cache()

        let max_hash = max_hash_for_scaled(scaled).unwrap();

        // filter out all hashes over max_hash in value.
        let mut new_hashval_to_inx: HashMap<u64, Vec<u32>> = HashMap::new();
        for (hashval, idx) in self.hashval_to_idx.clone() {
            if hashval < max_hash {
                new_hashval_to_inx.insert(hashval, idx);
            }
        }
        self.hashval_to_idx = new_hashval_to_inx;
        self.scaled = scaled;
        Ok(())
    }

    //default: threshold_bp = 0
    pub fn gather(
        &self, 
        query: Signature, 
        threshold_bp: f32,
    ) -> Vec<(f32, Signature, String)> {

        let mut results: Vec<(f32, Signature, String)> = Vec::new();

        if query.signatures.len() == 0 {
            return results;
        }

        if let Sketch::MinHash(mh) = &query.signatures[0] {

            let threshold = threshold_bp / (mh.mins().len() as f32 * self.scaled as f32);

            for (score, sig, filename) in self._find_signatures(mh, threshold, true, true, true).unwrap() {
                if score > 0.0 {
                    results.push((score, sig, filename));
                }
            }
            
            return results;
        } else {
            unimplemented!()
        }

    }

    // default: do_containment = false, best_only = false, ignore_abundance = false.
    // "Return the match with the best Jaccard containment in the database."
    pub fn search(
        &self, 
        query: Signature, 
        threshold: f32, 
        do_containment: bool, 
        ignore_abundance: bool,
    ) -> Result<Vec<(f32, Signature, String)>, Error> {

        let mut results: Vec<(f32, Signature, String)> = Vec::new();

        if query.signatures.len() == 0 {
            return Ok(results);
        }

        if let Sketch::MinHash(mh) = &query.signatures[0] {
            let mut mh = mh.clone();
            if ignore_abundance {
                mh.disable_abundance();
            }

            for (score, sig, filename) in self._find_signatures(&mh, threshold, do_containment, false, false)? {
                results.push((score, sig, filename));
            }

            // sort by score
            results.sort_by(|a, b| b.0.partial_cmp(&a.0).unwrap_or(Equal));
            Ok(results)
        } else {
            unimplemented!()
        }
    }
}



// TESTINGGGGGG!
#[cfg(test)]                                                                                   
mod test { 
    use std::fs::File;
    use std::io::{BufReader, Seek, SeekFrom};
    use std::path::PathBuf;
    use std::collections::{HashMap, BTreeMap};

    use crate::signature::Signature;
    use crate::index::lca_db::{Lineage, LcaDB};


    #[test]
    fn lca_search() {
        let mut lca_db = LcaDB::new();

        // get signature to add
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10.fa.gz.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let signatures: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        lca_db.insert(&signatures[0], Some("erik"), None).unwrap();

        // get second signature to add
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10+s11.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let signatures: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        lca_db.insert(&signatures[0], Some("erik2"), None).unwrap();
        
        let tup = lca_db.search(signatures[0].clone(), 0.0, false, false);

        dbg!(&tup[1].2);

        // should be 1 result from one of the .insert(...)
        assert!(tup.len() == 2);
        assert!(tup[0].0 == 1.0); // testing the score it got
        assert!(tup[1].0 == 0.31752306); // testing the score the second sig got
        assert!(tup[0].2 == "".to_string()); // testing the filename it returned
    }

    #[test]
    fn lca_gather() {
        //load lca_db
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/lca/delmont-1.lca.json");

        let mut lca_db = LcaDB::load(filename.as_path()).unwrap();

        // get signature to add
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10+s11.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let signatures: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        lca_db.insert(&signatures[0], Some("erik"), None).unwrap();
        lca_db.insert(&signatures[0], Some("erik2"), None).unwrap();
        
        let tup = lca_db.gather(signatures[0].clone(), 0.0);

        dbg!(&tup[0].2);

        // should be 2 identicle results from both .insert(...)
        assert!(tup.len() == 1);
        assert!(tup[0].0 == 0.39); // testing the score it got
        assert!(tup[0].2 == "delmont-1.lca.json".to_string()); // testing the filename it returned
    }

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
        let mut lineage: Lineage = BTreeMap::new();
        lineage.insert("name1".to_string(), "rank1".to_string()); 
        lineage.insert("name2".to_string(), "rank2".to_string());
        
        // add new sigs and lineage
        lca_db.insert(&signatures[0], None, Some(&lineage)).unwrap();

        let sigs = lca_db.signatures().unwrap();
        dbg!(&sigs.len());

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

        dbg!(&lca_db.moltype);

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

        assert_eq!(lcadb.class, lcadb_2.class);
        assert_eq!(lcadb.license, lcadb_2.license);

        assert_eq!(lcadb.ksize, lcadb_2.ksize);
        assert_eq!(lcadb.scaled, lcadb_2.scaled);
        assert_eq!(lcadb.moltype, lcadb_2.moltype);
        assert_eq!(lcadb.filename, "delmont-1.lca.json".to_string());

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
        assert!(lca_db.filename() == "".to_string());
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
        let mut lineage: Lineage = BTreeMap::new();
        lineage.insert("name1".to_string(), "rank1".to_string()); 
        lineage.insert("name2".to_string(), "rank2".to_string());
        
        lca_db.insert(&sigs[0], Some("erik"), Some(&lineage)).unwrap();

        // println!("{:?}", lca_db);
                                            
        assert!(lca_db.ksize() == 32);
        assert!(lca_db.scaled() == 1);
        assert!(lca_db.filename() == "".to_string());
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