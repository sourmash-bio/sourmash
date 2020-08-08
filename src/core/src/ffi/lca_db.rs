use std::ffi::CStr;
use std::os::raw::c_char;

use crate::ffi::utils::{ForeignObject};
use crate::index::lca_db::{LcaDB, LineagePairs};
use crate::ffi::signature::SourmashSignature;
use std::collections::BTreeMap;
use serde_json::{Result, Value};


pub struct SourmashLcaDatabase;

impl ForeignObject for SourmashLcaDatabase {
    type RustObject = LcaDB;
}

#[no_mangle]
pub unsafe extern "C" fn LcaDB_new() -> *mut SourmashLcaDatabase {
    let lca_db = LcaDB::new();
    SourmashLcaDatabase::from_rust(lca_db)
}

ffi_fn! {
unsafe fn LcaDB_insert(
    ptr: *mut SourmashLcaDatabase,
    sig_ptr:  *const SourmashSignature,
    ident_char: *const c_char,
    lineage_char: *const c_char,
) -> Result<u32> {

    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    let sig = SourmashSignature::as_rust(sig_ptr);
    let ident = {
        assert!(!ident_char.is_null());
        CStr::from_ptr(ident_char)
    };
    let data = {
        assert!(!lineage_char.is_null());
        CStr::from_ptr(lineage_char)
    };

    //set lineage
    let lineage: LineagePairs = serde_json::from_str(data.to_str().unwrap())?;
    println!("\n\n\nRUST:\nlca_db = {:?}", lca_db);
    

    lca_db.insert(sig, ident.to_bytes(), &lineage)

    // if let Sketch::MinHash(minhash) = &sig.signatures[0] {

    //     // set lineage
    //     let lineage_str = {
    //         assert!(!lineage_ptr.is_null());
    //         slice::from_raw_parts(lineage_ptr, insize)
    //     };

    //     // check for errors
            
    //     assert!(lca_db.ident_to_name().contains_key(&ident));
    //         //error: raise ValueError("signature {} is already in this LCA db.".format(ident))
        
    //     // implement self.cache property and _invalidate_cache() method
    //     // self._invalidate_cache()
            
    //     // downsample to specified scaled; this has the side effect of
    //     // making sure they're all at the same scaled value!
    //     let minhash = minhash.downsample_scaled(lca_db.scaled()).unwrap();
    //         // error if it is a scaled signature

    //     // store name
    //     lca_db.insert_ident_to_name(ident.clone(), sig.name());

    //     // identifier -> integer index (idx)
    //     let idx = lca_db._get_ident_index(ident.clone(), true);

    //     if lineage.len() > 0 {
    //         // (LineagePairs*) -> integer lineage ids (lids)
    //         let lid = lca_db._get_lineage_id(lineage);

    //         // map idx to lid as well.
    //         if !lca_db.idx_to_lid().contains_key(&idx) {
    //             lca_db.insert_idx_to_lid(idx, lid);
    //         }
    //     }

    //     // append idx to each idx
    //     for hashval in minhash.mins() {
    //         lca_db.insert_hashval_to_idx_vec(hashval, idx);
    //     }

    //     println!("\n\n\n{:?}\n\n\n", minhash);

    //     Ok(minhash.num())
    // }
    // else {
    //     unimplemented!()
    // }
}
}

// #[no_mangle]
// pub unsafe extern "C" fn LcaDB_select(
//     ptr: *mut SourmashLcaDatabase, 
//     ksize: u32,
//     moltype_char: *const c_char,
// ) -> bool {
//     let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    
//     let moltype = LcaDB::c_char_to_string(moltype_char);
//     println!("ksize = {}, moltype_char = {}", ksize, moltype);
//     println!("ksize = {}, moltype_char = {}", lca_db.ksize(), lca_db.moltype());

//     if ksize != 0 && lca_db.ksize() != ksize {
//         return false
//     }
//     if moltype != "".to_string() && moltype != lca_db.moltype() {
//         return false
//     }
//     return true
// }

// ffi_fn! {
// unsafe fn LcaDB_save(
//     ptr: *mut SourmashLcaDatabase, 
//     db_name_char: *const c_char,
// ) -> Result<()> {
//     let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
//     let db_name = LcaDB::c_char_to_string(db_name_char);

//     let st = serde_json::to_string(lca_db).unwrap();
//     println!("\n\nRUST:\n{}\n\n", st);
//     Ok(())
// }
// }

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
