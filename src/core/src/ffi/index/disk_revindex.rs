use std::ffi::CStr;
use std::os::raw::c_char;

// use crate::ffi::index::SourmashSearchResult;
use crate::ffi::minhash::SourmashKmerMinHash;
use crate::index::revindex::RevIndexOps;
use crate::ffi::signature::SourmashSignature;
use crate::ffi::utils::{ForeignObject};
use crate::index::revindex::RevIndex as BasicRevIndex;
use crate::index::revindex::disk_revindex::RevIndex as DDRevIndex;
// use crate::collection::Collection;
// use crate::index::Index;
// use crate::prelude::*;
use crate::signature::Signature;
use crate::sketch::minhash::KmerMinHash;
// use crate::sketch::Sketch;
// use crate::ScaledType;

pub struct SourmashDiskRevIndex;

impl ForeignObject for SourmashDiskRevIndex {
    type RustObject = BasicRevIndex;
}

ffi_fn! {
unsafe fn disk_revindex_new_from_rocksdb(
    path_ptr: *const c_char,
) -> Result<*mut SourmashDiskRevIndex> {
    // FIXME use buffer + len instead of cstr
    let rocksdb_path = {
        assert!(!path_ptr.is_null());
        CStr::from_ptr(path_ptr)
    }.to_str()?;

    let rocksdb: BasicRevIndex = DDRevIndex::open(
        rocksdb_path,
        true,
        None
    )?;

    Ok(SourmashDiskRevIndex::from_rust(rocksdb))
}
}

#[no_mangle]
pub unsafe extern "C" fn disk_revindex_free(ptr: *mut SourmashDiskRevIndex) {
    SourmashDiskRevIndex::drop(ptr);
}

/*
ffi_fn! {
unsafe fn revindex_search(
    ptr: *const SourmashRevIndex,
    sig_ptr: *const SourmashSignature,
    threshold: f64,
    do_containment: bool,
    _ignore_abundance: bool,
    size: *mut usize,
) -> Result<*const *const SourmashSearchResult> {
    let revindex = SourmashRevIndex::as_rust(ptr);
    let sig = SourmashSignature::as_rust(sig_ptr);

    if sig.signatures.is_empty() {
        *size = 0;
        return Ok(std::ptr::null::<*const SourmashSearchResult>());
    }

    let mh = if let Sketch::MinHash(mh) = &sig.signatures[0] {
        mh
    } else {
        // TODO: what if it is not a mh?
        unimplemented!()
    };

    let results: Vec<(f64, Signature, String)> = revindex
        .find_signatures(mh, threshold, do_containment, true)?
        .into_iter()
        .collect();

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*const SourmashSearchResult> = results
        .into_iter()
        .map(|x| Box::into_raw(Box::new(x)) as *const SourmashSearchResult)
        .collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *const *const SourmashSearchResult)
}
}

ffi_fn! {
unsafe fn revindex_gather(
    ptr: *const SourmashRevIndex,
    sig_ptr: *const SourmashSignature,
    threshold: f64,
    _do_containment: bool,
    _ignore_abundance: bool,
    size: *mut usize,
) -> Result<*const *const SourmashSearchResult> {
    let revindex = SourmashRevIndex::as_rust(ptr);
    let sig = SourmashSignature::as_rust(sig_ptr);

    if sig.signatures.is_empty() {
        *size = 0;
        return Ok(std::ptr::null::<*const SourmashSearchResult>());
    }

    let mh = if let Sketch::MinHash(mh) = &sig.signatures[0] {
        mh
    } else {
        // TODO: what if it is not a mh?
        unimplemented!()
    };

    // TODO: proper threshold calculation
    let threshold: usize = (threshold * (mh.size() as f64)) as _;

    let counter = revindex.counter_for_query(mh);
    dbg!(&counter);

    let results: Vec<(f64, Signature, String)> = revindex
        .gather(counter, threshold, mh)
        .unwrap() // TODO: proper error handling
        .into_iter()
        .map(|r| {
            let filename = r.filename().to_owned();
            let sig = r.get_match();
            (r.f_match(), sig, filename)
        })
        .collect();

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*const SourmashSearchResult> = results
        .into_iter()
        .map(|x| Box::into_raw(Box::new(x)) as *const SourmashSearchResult)
        .collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *const *const SourmashSearchResult)
}
}

#[no_mangle]
pub unsafe extern "C" fn revindex_scaled(ptr: *const SourmashRevIndex) -> ScaledType {
    let revindex = SourmashRevIndex::as_rust(ptr);
    if let Sketch::MinHash(mh) = revindex.template() {
        mh.scaled()
    } else {
        unimplemented!()
    }
}
*/ 

#[no_mangle]
pub unsafe extern "C" fn disk_revindex_len(ptr: *const SourmashDiskRevIndex) -> u64 {
    let revindex = SourmashDiskRevIndex::as_rust(ptr);
    revindex.collection().len() as u64
}

ffi_fn! {
unsafe fn disk_revindex_signatures(
    ptr: *const SourmashDiskRevIndex,
    size: *mut usize,
) -> Result<*mut *mut SourmashSignature> {
    let revindex: &BasicRevIndex = SourmashDiskRevIndex::as_rust(ptr);

    let coll = revindex.collection();

    let sigs: Vec<Signature> = coll
        .iter()
        .filter_map(|(_idx, record)| match coll.sig_from_record(record) {
            Ok(sig) => { Some(sig.into()) },
            Err(_) => None,
        })
        .collect();

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*mut SourmashSignature> = sigs
        .into_iter()
        .map(|x| Box::into_raw(Box::new(x)) as *mut SourmashSignature)
        .collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashSignature)
}
}


ffi_fn! {
unsafe fn disk_revindex_best_containment(
    db_ptr: *const SourmashDiskRevIndex,
    query_ptr: *const SourmashSignature,
    threshold_bp: u16
) -> Result<*mut SourmashSignature> {
    let revindex: &BasicRevIndex = SourmashDiskRevIndex::as_rust(db_ptr);
    let sig = SourmashSignature::as_rust(query_ptr);

    // extract KmerMinHash for query
    let query_mh: KmerMinHash = sig.clone()
        .try_into().expect("cannot get kmerminhash");
    let scaled = query_mh.scaled();
    let threshold = threshold_bp as u32 / scaled as u32;

    // do search & get first/best match
    let counter = revindex.counter_for_query(&query_mh);
    let (dataset_id, size) = counter.k_most_common_ordered(1)[0];

    if size as u32 >= threshold {
        // load into SigStore & convert to Signature.
        let match_sig = revindex.collection().sig_for_dataset(dataset_id)?;
        let match_sig: Signature = match_sig.into();

        Ok(SourmashSignature::from_rust(match_sig))
    } else {
        Ok(SourmashSignature::from_rust(Signature::default()))
    }
}
}
    
ffi_fn! {
unsafe fn disk_revindex_peek(
    db_ptr: *const SourmashDiskRevIndex,
    query_ptr: *const SourmashKmerMinHash,
) -> Result<*mut SourmashSignature> {
    let revindex: &BasicRevIndex = SourmashDiskRevIndex::as_rust(db_ptr);
    let query_mh = SourmashKmerMinHash::as_rust(query_ptr);
    
    // do search & get first/best match
    let counter = revindex.counter_for_query(&query_mh);
    let (dataset_id, _size) = counter.k_most_common_ordered(1)[0];

    // load into SigStore & convert to Signature.
    let match_sig = revindex.collection().sig_for_dataset(dataset_id)?;
    let match_sig: Signature = match_sig.into();

    Ok(SourmashSignature::from_rust(match_sig))
}
}
    
