use std::ffi::CStr;
use std::os::raw::c_char;

use crate::ffi::utils::{ForeignObject};
use crate::index::lca_db::{LcaDB, LineagePairs};
use crate::ffi::signature::SourmashSignature;
use std::collections::BTreeMap;


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
    dbg!(&lca_db.lid_to_lineage());
    

    lca_db.insert(sig, ident.to_bytes(), &lineage)
}
}

// #[no_mangle]
// pub unsafe extern "C" fn LcaDB_select(
//     ptr: *mut SourmashLcaDatabase, 
//     ksize_opt: u32,
//     moltype_char: *const c_char,
// ) -> bool {
//     let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    
//     let moltype = if moltype_char.is_null() { 
//         None
//     } else {
//         Some(CStr::from_ptr(moltype_char).to_str().unwrap())
//     };

//     let ksize = match ksize_opt {
//         -1 => None,
//         >-1 => Some(ksize_opt)
//     };

//     lca_db.select(ksize, moltype).is_ok()
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
