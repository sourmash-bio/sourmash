use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::ffi::utils::{ForeignObject, SourmashStr};
use crate::sketch::lca_db::LcaDB;
use crate::ffi::minhash::SourmashKmerMinHash;

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

#[no_mangle]
pub unsafe extern "C" fn LcaDB_signatures(ptr: *mut SourmashLcaDatabase) {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    lca_db.signatures()
}
