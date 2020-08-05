use std::ffi::CStr;
use std::os::raw::c_char;
use std::str;
use std::vec::Vec;
use std::slice;
use std::convert::From;
use crate::Error;

use crate::signature::{Signature, SigsTrait};
use crate::sketch::Sketch;
use crate::index::Index;
use crate::ffi::utils::{ForeignObject};
use crate::index::lca_db::{LcaDB, LineagePair};
use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::signature::SourmashSignature;


pub struct SourmashLcaDatabase;

impl ForeignObject for SourmashLcaDatabase {
    type RustObject = LcaDB;
}

#[no_mangle]
pub unsafe extern "C" fn LcaDB_new(
    ksize: u32,
    scaled: u64,
    filename: *const c_char,
    moltype: *const c_char,
) -> *mut SourmashLcaDatabase {
    let filename = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };
    
    let moltype = {
        assert!(!moltype.is_null());

        CStr::from_ptr(moltype)
    };

    let lca_db = LcaDB::new(ksize, scaled, filename.to_bytes(), moltype.to_bytes());
    SourmashLcaDatabase::from_rust(lca_db)
}

#[derive(Clone)]
#[repr(C)]
pub struct AcceptedLineagePair {
    name: *const c_char,
    rank: *const c_char,
}

impl AcceptedLineagePair {
    pub fn name(&self) -> *const c_char {
        self.name
    }
    pub fn rank(&self) -> *const c_char {
        self.rank
    }
}

ffi_fn! {
unsafe fn LcaDB_insert(
    ptr: *mut SourmashLcaDatabase, 
    sig_ptr:  *const SourmashSignature, 
    ident_char: *const c_char,
    lineage_ptr: *const AcceptedLineagePair,
    insize: usize,
) -> Result<u32> {

    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    let sig = SourmashSignature::as_rust(sig_ptr);
    let ident = if ident_char.is_null() {
        sig.name()
    } else {
        LcaDB::c_char_to_string(ident_char)
    };

    if let Sketch::MinHash(minhash) = &sig.signatures[0] {

        // set lineage
        let lineage_str = {
            assert!(!lineage_ptr.is_null());
            slice::from_raw_parts(lineage_ptr, insize)
        };
        let mut lineage: Vec<LineagePair> = Vec::new();
        for l in lineage_str {
            lineage.push(l.into());
        }

        // check for errors
            
        assert!(lca_db.ident_to_name().contains_key(&ident));
            //error: raise ValueError("signature {} is already in this LCA db.".format(ident))
        
        // implement self.cache property and _invalidate_cache() method
        // self._invalidate_cache()
            
        // downsample to specified scaled; this has the side effect of
        // making sure they're all at the same scaled value!
        let minhash = minhash.downsample_scaled(lca_db.scaled()).unwrap();
            // error if it is a scaled signature

        // store name
        lca_db.insert_ident_to_name(ident.clone(), sig.name());

        // identifier -> integer index (idx)
        let idx = lca_db._get_ident_index(ident.clone(), true);

        if lineage.len() > 0 {
            // (LineagePairs*) -> integer lineage ids (lids)
            let lid = lca_db._get_lineage_id(lineage);

            // map idx to lid as well.
            if !lca_db.idx_to_lid().contains_key(&idx) {
                lca_db.insert_idx_to_lid(idx, lid);
            }
        }

        // append idx to each idx
        for hashval in minhash.mins() {
            lca_db.insert_hashval_to_idx_vec(hashval, idx);
        }

        println!("\n\n\n{:?}\n\n\n", minhash);

        Ok(minhash.num())
    }
    else {
        unimplemented!()
    }
}
}

#[no_mangle]
pub unsafe extern "C" fn LcaDB_select(
    ptr: *mut SourmashLcaDatabase, 
    ksize: u32,
    moltype_char: *const c_char,
) -> bool {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    
    let moltype = LcaDB::c_char_to_string(moltype_char);
    println!("ksize = {}, moltype_char = {}", ksize, moltype);
    println!("ksize = {}, moltype_char = {}", lca_db.ksize(), lca_db.moltype());

    if ksize != 0 && lca_db.ksize() != ksize {
        return false
    }
    if moltype != "".to_string() && moltype != lca_db.moltype() {
        return false
    }
    return true
}

ffi_fn! {
unsafe fn LcaDB_save(
    ptr: *mut SourmashLcaDatabase, 
    db_name_char: *const c_char,
) -> Result<()> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    let db_name = LcaDB::c_char_to_string(db_name_char);

    let st = serde_json::to_string(lca_db).unwrap();
    println!("\n\nRUST:\n{}\n\n", st);
    Ok(())
}
}

// #[no_mangle]
// pub unsafe extern "C" fn LcaDB_temp(ptr: *mut SourmashLcaDatabase) -> *const SourmashLcaDatabase {
//     let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

//     println!("\n");
//     let hashes = lca_db.hashval_to_idx();
//     let temp = hashes.unwrap_or_else("nothing");
//     for h in temp {
//         println!("{}", h);
//     }
//     println!("\n");
//     SourmashLcaDatabase::from_ref(lca_db)
// }

// #[no_mangle]
// pub unsafe extern "C" fn LcaDB_signatures(ptr: *mut SourmashLcaDatabase) {
//     let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
//     lca_db.signatures()
// }
