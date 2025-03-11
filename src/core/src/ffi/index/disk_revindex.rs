use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::collection::{Collection, CollectionSet};
use crate::ffi::index::SourmashSearchResult;
use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::signature::SourmashSignature;
use crate::ffi::utils::ForeignObject;
use crate::index::revindex::disk_revindex::RevIndex as DDRevIndex;
use crate::index::revindex::RevIndex as BasicRevIndex;
use crate::index::revindex::{DatasetPicklist, RevIndexOps};
use std::collections::HashSet;
use std::ffi::CString;
use std::path::Path;
// use crate::index::Index;
// use crate::prelude::*;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
// use crate::sketch::Sketch;
// use crate::ScaledType;

pub struct SourmashDiskRevIndex;

impl ForeignObject for SourmashDiskRevIndex {
    type RustObject = BasicRevIndex;
}

pub struct SourmashDatasetPicklist;

impl ForeignObject for SourmashDatasetPicklist {
    type RustObject = DatasetPicklist;
}

unsafe fn retrieve_picklist(
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Option<DatasetPicklist> {
    if dataset_picklist_ptr.is_null() {
        None
    } else {
        let x = SourmashDatasetPicklist::as_rust(dataset_picklist_ptr);
        Some(x.clone())
    }
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

ffi_fn! {
unsafe fn disk_revindex_new_with_sigs( // @CTB rename to create
    sigs_ptr: *const *const SourmashSignature,
    insigs: usize,
    path_ptr: *const c_char,
) -> Result<()> {
    let sigs: Vec<Signature> = {
        assert!(!sigs_ptr.is_null());
        slice::from_raw_parts(sigs_ptr, insigs)
            .iter()
            .map(|sig| SourmashSignature::as_rust(*sig))
            .cloned()
            .collect()
    };

    let coll = Collection::from_sigs(sigs).expect("cannot create Collection");
    let cs: CollectionSet = coll.try_into().expect("cannot convert to CollectionSet");

    let rocksdb_path = {
        assert!(!path_ptr.is_null());
        CStr::from_ptr(path_ptr)
    }.to_str()?;

    let rocksdb_path = Path::new(rocksdb_path);

    let mut revindex = DDRevIndex::create(rocksdb_path, cs).expect("cannot create RocksDB");
    revindex.internalize_storage().expect("failed to internalize storage.");
    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn disk_revindex_free(ptr: *mut SourmashDiskRevIndex) {
    SourmashDiskRevIndex::drop(ptr);
}

ffi_fn! {
unsafe fn dataset_picklist_new_from_list(
    dataset_idxs_ptr: *const u32,
    insize: usize,
) -> Result<*const SourmashDatasetPicklist> {
    assert!(!dataset_idxs_ptr.is_null());
    let dids = HashSet::from_iter(
        slice::from_raw_parts(dataset_idxs_ptr as *mut u32, insize)
            .iter().copied()
    );

    let ds = DatasetPicklist {
        dataset_ids: dids
    };

    Ok(SourmashDatasetPicklist::from_rust(ds))
}
}

#[no_mangle]
pub unsafe extern "C" fn dataset_picklist_free(ptr: *mut SourmashDatasetPicklist) {
    SourmashDatasetPicklist::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn disk_revindex_len(ptr: *const SourmashDiskRevIndex) -> u64 {
    let revindex = SourmashDiskRevIndex::as_rust(ptr);
    revindex.collection().len() as u64
}

#[no_mangle]
pub unsafe extern "C" fn disk_revindex_ksize(ptr: *const SourmashDiskRevIndex) -> u32 {
    let revindex = SourmashDiskRevIndex::as_rust(ptr);

    revindex
        .collection()
        .manifest()
        .first()
        .expect("no records!?")
        .ksize()
}

#[no_mangle]
pub unsafe extern "C" fn disk_revindex_scaled(ptr: *const SourmashDiskRevIndex) -> u32 {
    let revindex = SourmashDiskRevIndex::as_rust(ptr);
    let (_, scaled) = revindex
        .collection()
        .min_max_scaled()
        .expect("no records!?");
    *scaled
}

#[no_mangle]
pub unsafe extern "C" fn disk_revindex_moltype(ptr: *const SourmashDiskRevIndex) -> *const c_char {
    let revindex = SourmashDiskRevIndex::as_rust(ptr);
    let moltype = revindex
        .collection()
        .manifest()
        .first()
        .expect("no records!?")
        .moltype();
    let moltype_str = moltype.to_string();
    let c_string = CString::new(moltype_str).expect("foo");
    c_string.as_ptr()
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
    threshold_bp: u16,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*mut SourmashSignature> {
    let revindex: &BasicRevIndex = SourmashDiskRevIndex::as_rust(db_ptr);
    let sig = SourmashSignature::as_rust(query_ptr);

    // extract KmerMinHash for query
    let query_mh: KmerMinHash = sig.clone()
        .try_into().expect("cannot get kmerminhash");
    let scaled = query_mh.scaled();
    let threshold = threshold_bp as u32 / scaled;

    // picklist?
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    // do search & get first/best match
    let counter = revindex.counter_for_query(&query_mh, dataset_picklist);
    let (dataset_id, size) = counter.k_most_common_ordered(1)[0];

    if size as u32 >= threshold {
        // load into SigStore & convert to Signature.
        let match_sig = revindex.collection().sig_for_dataset(dataset_id)?;
        let match_sig: Signature = match_sig.into();

        Ok(SourmashSignature::from_rust(match_sig))
    } else {
        Ok(SourmashSignature::from_rust(Signature::default())) // @CTB
    }
}
}

// implement prefetch/containment separately from search/jaccard

ffi_fn! {
unsafe fn disk_revindex_prefetch(
    db_ptr: *const SourmashDiskRevIndex,
    query_ptr: *const SourmashSignature,
    threshold_bp: u64,
    return_size: *mut usize,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*const *const SourmashSearchResult> {
    let revindex: &BasicRevIndex = SourmashDiskRevIndex::as_rust(db_ptr);
    let sig = SourmashSignature::as_rust(query_ptr);

    // extract KmerMinHash for query
    let query_mh: KmerMinHash = sig.clone()
        .try_into().expect("cannot get kmerminhash");
    let scaled = query_mh.scaled();
    let threshold_bp: u64 = threshold_bp as u64 / scaled as u64;

    // picklist?
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    // do search & get matches
    let counter = revindex.counter_for_query(&query_mh, dataset_picklist);

    // right now this iterates over all matches from 'counter.most_common()'.
    // we could probably truncate the search here in some way, yes?
    // but it would require changing this to a loop rather than using an
    // iterator I think.
    let results: Vec<(f64, Signature, String)> = counter
        .most_common()
        .into_iter()
        .filter_map(|(dataset_id, size)| {
            if size as u64 >= threshold_bp {
                let filename = "some rocksdb database"; // @CTB
                let sig: Signature = revindex
                    .collection()
                    .sig_for_dataset(dataset_id)
                    .expect("dataset not found")
                    .into();
                let f_cont = size as f64 / query_mh.size() as f64;

                Some((f_cont, sig, filename.to_owned()))
            } else {
                None
            }
        })
        .collect();

    // convert to ffi.
    let ptr_results: Vec<*const SourmashSearchResult> = results
        .into_iter()
        .map(|x| Box::into_raw(Box::new(x)) as *const SourmashSearchResult)
        .collect();

    let b = ptr_results.into_boxed_slice();
    *return_size = b.len();
    Ok(Box::into_raw(b) as *const *const SourmashSearchResult)
}
}

// implement search/jaccard separately from search/jaccard asdf @CTB

ffi_fn! {
unsafe fn disk_revindex_search_jaccard(
    db_ptr: *const SourmashDiskRevIndex,
    query_ptr: *const SourmashSignature,
    threshold: f64,
    return_size: *mut usize,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*const *const SourmashSearchResult> {
    let revindex: &BasicRevIndex = SourmashDiskRevIndex::as_rust(db_ptr);
    let sig = SourmashSignature::as_rust(query_ptr);

    // extract KmerMinHash for query
    let query_mh: KmerMinHash = sig.clone()
        .try_into().expect("cannot get kmerminhash");

    // picklist?
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    // do search
    let counter = revindex.counter_for_query(&query_mh, dataset_picklist);

    // retrieve/convert matches. I don't think there's a simple way to
    // truncate this without going through all the matches, so it's
    // potentially (much) more expensive than prefetch.
    let results: Vec<(f64, Signature, String)> = counter
        .most_common()
        .into_iter()
        .filter_map(|(dataset_id, _size)| {
            let filename = "some rocksdb database"; // @CTB
            let sig: Signature = revindex
                .collection()
                .sig_for_dataset(dataset_id)
                .expect("dataset not found")
                .into();

            let match_mh = sig.minhash().expect("cannot retrieve match");
            let f_match = query_mh.jaccard(match_mh).expect("cannot calculate Jaccard");

            if f_match >= threshold {
                Some((f_match, sig, filename.to_owned()))
            } else {
                None
            }
        })
        .collect();

    // convert to ffi.
    let ptr_results: Vec<*const SourmashSearchResult> = results
        .into_iter()
        .map(|x| Box::into_raw(Box::new(x)) as *const SourmashSearchResult)
        .collect();

    let b = ptr_results.into_boxed_slice();
    *return_size = b.len();
    Ok(Box::into_raw(b) as *const *const SourmashSearchResult)
}
}

// implement peek: used in 'gather' to retrieve best containment possible

ffi_fn! {
unsafe fn disk_revindex_peek(
    db_ptr: *const SourmashDiskRevIndex,
    query_ptr: *const SourmashKmerMinHash,
    threshold_bp: u64,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*mut SourmashSignature> {
    let revindex: &BasicRevIndex = SourmashDiskRevIndex::as_rust(db_ptr);
    let query_mh = SourmashKmerMinHash::as_rust(query_ptr);
    let scaled = query_mh.scaled();
    let threshold_bp: u64 = threshold_bp as u64 / scaled as u64;

    // picklist?
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    // do search & get first/best match
    let counter = revindex.counter_for_query(query_mh, dataset_picklist);
    let (dataset_id, size) = counter.k_most_common_ordered(1)[0];

    if size as u64 >= threshold_bp {
        // load into SigStore & convert to Signature.
        let match_sig = revindex.collection().sig_for_dataset(dataset_id)?;
        let match_sig: Signature = match_sig.into();

        Ok(SourmashSignature::from_rust(match_sig))
    } else {
        Ok(SourmashSignature::from_rust(Signature::default())) // @CTB
    }
}
}
