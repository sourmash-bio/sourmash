use std::ffi::CStr;
use std::os::raw::c_char;
use std::path::Path;

use crate::ffi::utils::{ForeignObject};
use crate::index::lca_db::{LcaDB, Lineage, lineage_to_vec};
use crate::signature::Signature;
use crate::ffi::signature::SourmashSignature;
use std::collections::{HashMap, BTreeMap};


pub struct SourmashLcaDatabase;

impl ForeignObject for SourmashLcaDatabase {
    type RustObject = LcaDB;
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_new() -> *mut SourmashLcaDatabase {
    let lca_db = LcaDB::new();
    SourmashLcaDatabase::from_rust(lca_db)
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_new_with_params(
    ksize: u32,
    scaled: u64,
    filename_char: *const c_char,
    moltype_char: *const c_char,
) -> *mut SourmashLcaDatabase {
    let filename = {
        assert!(!filename_char.is_null());
        CStr::from_ptr(filename_char).to_str().unwrap()
    };
    let moltype = {
        assert!(!moltype_char.is_null());
        CStr::from_ptr(moltype_char).to_str().unwrap()
    };
    
    let lca_db = LcaDB::new_with_params(Some(ksize), Some(scaled), Some(filename), Some(moltype));

    SourmashLcaDatabase::from_rust(lca_db)
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_ksize(ptr: *const SourmashLcaDatabase) -> u32 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    lca_db.ksize()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_scaled(ptr: *const SourmashLcaDatabase) -> u64 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    lca_db.scaled()
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_moltype(ptr: *const SourmashLcaDatabase, size: *mut usize) -> *const u8 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let moltype = lca_db.moltype();

    let str_boxed = moltype.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let moltype = Box::into_raw(b) as *const u8;

    moltype
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_filename(ptr: *const SourmashLcaDatabase, size: *mut usize) -> *const u8 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let filename = lca_db.filename();

    let str_boxed = filename.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let filename = Box::into_raw(b) as *const u8;

    filename
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_ident_to_name(ptr: *const SourmashLcaDatabase, size: *mut usize) -> *const u8 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let buf = serde_json::to_string(&lca_db.ident_to_name()).unwrap();

    let str_boxed = buf.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let ident_to_name = Box::into_raw(b) as *const u8;

    ident_to_name
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_ident_to_idx(ptr: *const SourmashLcaDatabase, size: *mut usize) -> *const u8 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let buf = serde_json::to_string(&lca_db.ident_to_idx()).unwrap();

    let str_boxed = buf.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let ident_to_idx = Box::into_raw(b) as *const u8;

    ident_to_idx
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_idx_to_lid(ptr: *const SourmashLcaDatabase, size: *mut usize) -> *const u8 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);
    let buf = serde_json::to_string(&lca_db.idx_to_lid()).unwrap();

    let str_boxed = buf.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let idx_to_lid = Box::into_raw(b) as *const u8;

    idx_to_lid
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_lineage_to_lid(ptr: *const SourmashLcaDatabase, size: *mut usize) -> *const u8 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let mut lineage_to_lid_vec: HashMap<Vec<Vec<String>>, u32> = HashMap::new();
    for (lineage, lid) in lca_db.lineage_to_lid() {
        lineage_to_lid_vec.insert(lineage_to_vec(lineage), lid);
    }

    let buf = serde_json::to_string(&lineage_to_lid_vec).unwrap();

    let str_boxed = buf.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let lineage_to_lid = Box::into_raw(b) as *const u8;

    lineage_to_lid
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_lid_to_lineage(ptr: *const SourmashLcaDatabase, size: *mut usize) -> *const u8 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let mut lid_to_lineage_vec: HashMap<u32, Vec<Vec<String>>> = HashMap::new();
    for (lid, lineage) in lca_db.lid_to_lineage() {
        lid_to_lineage_vec.insert(lid, lineage_to_vec(lineage));
    }

    let buf = serde_json::to_string(&lid_to_lineage_vec).unwrap();

    let str_boxed = buf.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let lid_to_lineage = Box::into_raw(b) as *const u8;

    lid_to_lineage
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_hashval_to_idx(ptr: *const SourmashLcaDatabase, size: *mut usize) -> *const u8 {
    let lca_db = SourmashLcaDatabase::as_rust(ptr);

    let buf = serde_json::to_string(&lca_db.hashval_to_idx()).unwrap();

    let str_boxed = buf.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let hashval_to_idx = Box::into_raw(b) as *const u8;

    // dbg!(&hashval_to_idx);
    hashval_to_idx
}

#[no_mangle]
pub unsafe extern "C" fn lcadb_hashval_to_idx_len(ptr: *mut SourmashLcaDatabase) -> u32 {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    lca_db.hashval_to_idx().len() as u32
}

ffi_fn! {
unsafe fn lcadb_get_idx_from_hashval(
    ptr: *mut SourmashLcaDatabase, 
    hashval: u64, 
    size: *mut usize,
) -> Result<*const u32> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

    let buf: Vec<u32> = match lca_db.hashval_to_idx().get(&hashval) {
        Some(s) => s.to_vec(),
        None => vec![],
    };

    let mut x: Vec<u32> = Vec::new();
    for v in buf {
        x.push(v);
    }
    
    let b = x.into_boxed_slice();
    *size = b.len();
    let result = Box::into_raw(b) as *const u32;


    Ok(result)
}
}

ffi_fn! {
unsafe fn lcadb_get_lineage_from_idx(
    ptr: *mut SourmashLcaDatabase, 
    idx: u32, 
    size: *mut usize,
) -> Result<*const u8> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

    let lineage: Lineage = match lca_db.idx_to_lid().get(&idx) {
        Some(lid) => lca_db.lid_to_lineage()[&lid].clone(),
        None => BTreeMap::new(),
    };
  
    let buff = serde_json::to_string(&lineage_to_vec(lineage)).unwrap();

    let str_boxed = buff.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let result = Box::into_raw(b) as *const u8;


    Ok(result)
}
}

ffi_fn! {
unsafe fn lcadb_get_match_size(
    ptr: *mut SourmashLcaDatabase, 
    best_idx: u32,
) -> Result<u32> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

    let mut match_size = 0;
    for (hashval, idx) in lca_db.hashval_to_idx() {
        if idx.contains(&best_idx) {
            match_size += 1;
        }
    }

    Ok(match_size)
}
}

ffi_fn! {
unsafe fn lcadb_best_name(
    ptr: *mut SourmashLcaDatabase, 
    best_idx: u32,
    size: *mut usize,
) -> Result<*const u8> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

    let mut name = String::new();
    for (ident, idx) in lca_db.ident_to_idx() {
        if idx == best_idx {
            name = lca_db.ident_to_name()[&ident.clone()].clone();
        }
    }
    
    let str_boxed = name.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let result = Box::into_raw(b) as *const u8;
    Ok(result)
}
}

ffi_fn! {
unsafe fn make_assignments_helper(
    ptr: *mut SourmashLcaDatabase, 
    min_num: usize,
) -> Result<u32> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

    let mut match_size = 0;
    for (hashval, idxlist) in lca_db.hashval_to_idx() {
        if min_num > 0 && idxlist.len() < min_num {
            continue;
        }

        let idx_to_lid = lca_db.idx_to_lid();
        let mut assignments: HashMap<u64, Vec<Lineage>> = HashMap::new();
        for idx in idxlist {
            let lid = idx_to_lid.get(&idx);
            if lid != None {
                let lineage = lca_db.lid_to_lineage().get(&lid.unwrap()).unwrap().clone();
                let temp = assignments.get_mut(&hashval).unwrap();
                temp.push(lineage.clone());
            }
        }
    }

    Ok(match_size)
}
}

ffi_fn! {
unsafe fn lcadb_get_lineage_assignments(
    ptr: *mut SourmashLcaDatabase, 
    hashval: u64, 
    size: *mut usize,
) -> Result<*const u8> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

    let buf = lca_db.get_lineage_assignments(hashval).unwrap();

    let mut x: Vec<Vec<Vec<String>>> = Vec::new();
    for lineage in buf {
        x.push(lineage_to_vec(lineage));
    }
  
    let buff = serde_json::to_string(&x).unwrap();
    // dbg!(&buff);
    // println!("\n");

    let str_boxed = buff.into_boxed_str();
    let b = str_boxed.into_boxed_bytes();
    *size = b.len();
    let result = Box::into_raw(b) as *const u8;


    Ok(result)
}
}

ffi_fn! {
unsafe fn lcadb_downsample_scaled(
    ptr: *mut SourmashLcaDatabase, 
    new_scaled: u64,
) -> Result<()> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    
    lca_db.downsample_scaled(new_scaled)?;

    Ok(())
}
}

ffi_fn! {
unsafe fn lcadb_insert(
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

    //set lineage
    let data = CStr::from_ptr(lineage_char).to_str().unwrap();
    let lineage: Option<&Lineage> = None;
    let lineage = if data != "" {
        let lineage_temp: Lineage = serde_json::from_str(data)?;
        Some(lineage_temp)
    } else {
        None
    };
    
    lca_db.insert(sig, Some(ident.to_str().unwrap()), lineage)
}
}

// WORKS!
ffi_fn! {
unsafe fn lcadb_signatures(
    ptr: *mut SourmashLcaDatabase,
    size: *mut usize,
) -> Result<*mut *mut SourmashSignature> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);

    // signatures() -> Result<Vec<Signature>>
    let sigs = lca_db.signatures()?;

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*mut SourmashSignature> = sigs.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut SourmashSignature
    }).collect();

    
    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashSignature)
}
}

ffi_fn! {
unsafe fn lcadb_load(
    filepath_char: *const c_char,
) -> Result<*mut SourmashLcaDatabase> {

    let buf = {
        assert!(!filepath_char.is_null());
        CStr::from_ptr(filepath_char)
    };
    let filename = &buf.to_str()?;
    
    let path = Path::new(filename);
    dbg!(&path);
    let lca_db = LcaDB::load(path, filename).unwrap();

    Ok(SourmashLcaDatabase::from_rust(lca_db))
}
}

ffi_fn! {
unsafe fn lcadb_load_db(
    filepath_char: *const c_char,
) -> Result<*mut SourmashLcaDatabase> {

    let buf = {
        assert!(!filepath_char.is_null());
        CStr::from_ptr(filepath_char)
    };
    let filename = &buf.to_str()?;
    
    let (mut input, _) = niffler::from_path(filename)?;
    dbg!(&filename);
    let lca_db = LcaDB::load_db(input, filename).unwrap();

    Ok(SourmashLcaDatabase::from_rust(lca_db))
}
}

ffi_fn! {
unsafe fn lcadb_save(
    ptr: *mut SourmashLcaDatabase, 
    filepath_char: *const c_char,
) -> Result<()> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    let filepath = {
        assert!(!filepath_char.is_null());
        CStr::from_ptr(filepath_char)
    };
    let path = Path::new(filepath.to_str().unwrap());
    lca_db.save(path)
}
}

#[repr(C)]
pub struct FFISearchResults {
    score: f32,
    sig: *mut SourmashSignature,
    filename: *const u8,
    name_size: u32
}

impl FFISearchResults {
    pub fn new(score: f32, sig: Signature, filename: String) -> FFISearchResults {

        let str_boxed = filename.into_boxed_str();
        let b = str_boxed.into_boxed_bytes();
        let name_size = b.len() as u32;
        let filename = Box::into_raw(b) as *const u8;

        let sig = unsafe { SourmashSignature::from_rust(sig) };

        FFISearchResults {
            score,
            sig,
            filename,
            name_size,
        }
    }
}

ffi_fn! {
unsafe fn lcadb_search(
    ptr: *mut SourmashLcaDatabase,
    query: *const SourmashSignature,
    threshold: f32,
    do_containment: bool,
    ignore_abundance: bool,
    osize: *mut usize
) -> Result<*mut FFISearchResults> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    let sig = SourmashSignature::as_rust(query);

    let searchresults: Vec<(f32, Signature, String)> = lca_db.search(sig.clone(), threshold, do_containment, ignore_abundance)?;
    
    let mut result_vec: Vec<FFISearchResults> = Vec::new();
    for tup in searchresults {
        result_vec.push(FFISearchResults::new(tup.0, tup.1, tup.2));
    }
    
    let b = result_vec.into_boxed_slice();
    *osize = b.len();

    Ok(Box::into_raw(b) as *mut FFISearchResults)
}
}

ffi_fn! {
unsafe fn lcadb_gather(
    ptr: *mut SourmashLcaDatabase,
    query: *const SourmashSignature,
    threshold_bp: f32,
    osize: *mut usize
) -> Result<*mut FFISearchResults> {
    let lca_db = SourmashLcaDatabase::as_rust_mut(ptr);
    let sig = SourmashSignature::as_rust(query);

    let searchresults: Vec<(f32, Signature, String)> = lca_db.gather(sig.clone(), threshold_bp);

    let mut result_vec: Vec<FFISearchResults> = Vec::new();
    for tup in searchresults {
        result_vec.push(FFISearchResults::new(tup.0, tup.1, tup.2));
    }

    // let ptr_results: Vec<FFISearchResults> = searchresults.into_iter().map(|(a, b, c)| {
    //     Box::into_raw(Box::new( FFISearchResults::new(a, b, c) ))
    // }).collect();
    
    let b = result_vec.into_boxed_slice();
    *osize = b.len();

    Ok(Box::into_raw(b) as *mut FFISearchResults)
}
}
