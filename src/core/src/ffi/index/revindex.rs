use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::collection::{Collection, CollectionSet};
use crate::encodings::*;
use crate::ffi::index::SourmashSearchResult;
use crate::ffi::index::SourmashStr;
use crate::ffi::manifest::SourmashManifest;
use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::signature::SourmashSignature;
use crate::ffi::utils::ForeignObject;
use crate::index::revindex::disk_revindex;
use crate::index::revindex::mem_revindex;
use crate::index::revindex::{self as module, CounterGather, DatasetPicklist, RevIndexOps};
use crate::manifest::Record;
use crate::prelude::*;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
use std::collections::HashSet;
use std::path::Path;

// FFI struct for base RevIndex struct & RevIndexOps trait

pub struct SourmashRevIndex;
impl ForeignObject for SourmashRevIndex {
    type RustObject = module::RevIndex;
}

// FFI struct for RevIndex-specific picklist of Idx

pub struct SourmashDatasetPicklist;
impl ForeignObject for SourmashDatasetPicklist {
    type RustObject = DatasetPicklist;
}

// FFI struct for CounterGather object to hold intermediate results for
// gather.

#[allow(non_camel_case_types)]
pub struct SourmashRevIndex_CounterGather;
impl ForeignObject for SourmashRevIndex_CounterGather {
    type RustObject = CounterGather;
}

pub unsafe fn retrieve_picklist(
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Option<DatasetPicklist> {
    if dataset_picklist_ptr.is_null() {
        None
    } else {
        let x = SourmashDatasetPicklist::as_rust(dataset_picklist_ptr);
        Some(x.clone())
    }
}

// Build new RevIndex struct from existing RocksDB/DiskRevIndex.

ffi_fn! {
unsafe fn revindex_new_from_rocksdb(
    path_ptr: *const c_char,
) -> Result<*mut SourmashRevIndex> {
    // FIXME use buffer + len instead of cstr
    let rocksdb_path = {
        assert!(!path_ptr.is_null());
        CStr::from_ptr(path_ptr)
    }.to_str()?;

    let rocksdb = disk_revindex::DiskRevIndex::open(
        rocksdb_path,
        true,
        None
    )?;

    Ok(SourmashRevIndex::from_rust(rocksdb))
}
}

// Create new DiskRevIndex from list of signatures.

ffi_fn! {
unsafe fn revindex_disk_create(
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

    let mut revindex = disk_revindex::DiskRevIndex::create(rocksdb_path, cs).expect("cannot create RocksDB");
    revindex.internalize_storage().expect("failed to internalize storage.");
    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn revindex_free(ptr: *mut SourmashRevIndex) {
    SourmashRevIndex::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn revindex_countergather_free(ptr: *mut SourmashRevIndex_CounterGather) {
    SourmashRevIndex_CounterGather::drop(ptr);
}

// create a DatasetPicklist from a collection of Idx (record references).

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
pub unsafe extern "C" fn revindex_len(
    ptr: *const SourmashRevIndex,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> u64 {
    let revindex = SourmashRevIndex::as_rust(ptr);
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    let coll = revindex.collection();

    // filter by picklist
    let records: Vec<(Idx, &Record)> = coll
        .iter()
        .filter_map(|(idx, record)| {
            if let Some(pl) = &dataset_picklist {
                if pl.dataset_ids.contains(&idx) {
                    Some((idx, record))
                } else {
                    None
                }
            } else {
                Some((idx, record))
            }
        })
        .collect();

    records.len() as u64
}

#[no_mangle]
pub unsafe extern "C" fn revindex_ksize(ptr: *const SourmashRevIndex) -> u32 {
    let revindex = SourmashRevIndex::as_rust(ptr);

    // note: here 'collection' is a CollectionSet, so all the same ksize.
    revindex
        .collection()
        .manifest()
        .first()
        .expect("no records!?")
        .ksize()
}

#[no_mangle]
pub unsafe extern "C" fn revindex_scaled(ptr: *const SourmashRevIndex) -> u32 {
    let revindex = SourmashRevIndex::as_rust(ptr);

    // note: here 'collection' is a CollectionSet, so all the same scaled.
    let (_, scaled) = revindex
        .collection()
        .min_max_scaled()
        .expect("no records!?");
    *scaled
}

#[no_mangle]
pub unsafe extern "C" fn revindex_moltype(ptr: *const SourmashRevIndex) -> SourmashStr {
    let revindex = SourmashRevIndex::as_rust(ptr);

    // note: here 'collection' is a CollectionSet, so all the same moltype.
    let moltype = revindex
        .collection()
        .manifest()
        .first()
        .expect("no records!?")
        .moltype();
    let moltype_str = moltype.to_string();
    moltype_str.into()
}

ffi_fn! {
unsafe fn revindex_manifest(ptr: *const SourmashRevIndex) -> Result<*mut SourmashManifest> {
    let revindex = SourmashRevIndex::as_rust(ptr);
    let mf = revindex.collection().manifest().clone();

    Ok(SourmashManifest::from_rust(mf))
}
}

ffi_fn! {
unsafe fn revindex_signatures(
    ptr: *const SourmashRevIndex,
    size: *mut usize,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*mut *mut SourmashSignature> {
    let revindex = &SourmashRevIndex::as_rust(ptr);
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    let coll = revindex.collection();

    // filter by picklist
    let records: Vec<(Idx, &Record)> = coll.iter()
        .filter_map(|(idx, record)| {
            if let Some(pl) = &dataset_picklist {
                if pl.dataset_ids.contains(&idx) {
                    Some((idx, record))
                } else {
                    None
                }
            } else {
                Some((idx, record))
            }
        })
        .collect();

    // load sigs
    let sigs: Vec<Signature> = records.iter()
        .filter_map(|(idx, record)| match coll.sig_from_record(record) {
            Ok(sig) => Some(sig.into()),
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

// prefetch/containment overlap -> all matches. Implement separately from
// Jaccard, as this can be done efficiently on RevIndexes.

ffi_fn! {
unsafe fn revindex_prefetch(
    db_ptr: *const SourmashRevIndex,
    query_ptr: *const SourmashSignature,
    threshold_bp: u64,
    return_size: *mut usize,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*const *const SourmashSearchResult> {
    let revindex = &SourmashRevIndex::as_rust(db_ptr);
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
    //
    // we could also adjust 'counter_for_query' to respect a specific
    // threshold...
    let filename = revindex.location();
    let results: Vec<(f64, Signature, String)> = counter
        .most_common()
        .into_iter()
        .filter_map(|(dataset_id, size)| {
            if size as u64 >= threshold_bp {
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

// implement jaccard search separately from containment analysis, since
// the latter can be done more efficiently on RevIndexes.

ffi_fn! {
unsafe fn revindex_search_jaccard(
    ptr: *const SourmashRevIndex,
    sig_ptr: *const SourmashSignature,
    threshold: f64,
    size: *mut usize,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*const *const SourmashSearchResult> {
    let revindex = SourmashRevIndex::as_rust(ptr);
    let sig = SourmashSignature::as_rust(sig_ptr);

    // picklist?
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

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
        .find_signatures(mh, threshold, dataset_picklist)?
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
// retrieve best match.

ffi_fn! {
unsafe fn revindex_best_containment(
    db_ptr: *const SourmashRevIndex,
    query_ptr: *const SourmashKmerMinHash,
    threshold_bp: u64,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*mut SourmashSignature> {
    let revindex = &SourmashRevIndex::as_rust(db_ptr);
    let query_mh = SourmashKmerMinHash::as_rust(query_ptr);
    let scaled = query_mh.scaled();
    let threshold_bp: u64 = threshold_bp as u64 / scaled as u64;

    // picklist?
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    // do search & get first/best match
    let counter = revindex.counter_for_query(query_mh, dataset_picklist);
    if counter.len() >= 1 {
        let (dataset_id, size) = counter.k_most_common_ordered(1)[0];

        if size as u64 >= threshold_bp {
            // load into SigStore & convert to Signature.
            let match_sig = revindex
                .collection()
                .sig_for_dataset(dataset_id)
                .expect("cannot load signature");
            let match_sig: Signature = match_sig.into();

            return Ok(SourmashSignature::from_rust(match_sig));
        }
    }

    Ok(SourmashSignature::from_rust(Signature::default()))
}
}

// return a CounterGather object with prefetch results

ffi_fn! {
unsafe fn revindex_prefetch_to_countergather(
    db_ptr: *const SourmashRevIndex,
    query_ptr: *const SourmashSignature,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*mut SourmashRevIndex_CounterGather> {
    let revindex = &SourmashRevIndex::as_rust(db_ptr);
    let sig = SourmashSignature::as_rust(query_ptr);

    // extract KmerMinHash for query
    let query_mh: KmerMinHash = sig.clone()
        .try_into().expect("cannot get kmerminhash");

    // picklist?
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    // do search & get matches
    let counter = revindex.prepare_gather_counters(&query_mh, dataset_picklist);

    Ok(SourmashRevIndex_CounterGather::from_rust(counter))
}
}

// decrement counters appropriately.

ffi_fn! {
unsafe fn revindex_countergather_consume(
    cg_ptr: *mut SourmashRevIndex_CounterGather,
    isect_ptr: *const SourmashKmerMinHash,
) -> Result<()> {
    let cg: &mut CounterGather = SourmashRevIndex_CounterGather::as_rust_mut(cg_ptr);
    let isect_mh = SourmashKmerMinHash::as_rust(isect_ptr);

    cg.consume(isect_mh);

    Ok(())
}
}

// retrieve top match.

ffi_fn! {
unsafe fn revindex_countergather_peek(
    cg_ptr: *const SourmashRevIndex_CounterGather,
    db_ptr: *const SourmashRevIndex,
    threshold_bp: u64,
) -> Result<*mut SourmashSignature> {
    let cg: &CounterGather = SourmashRevIndex_CounterGather::as_rust(cg_ptr);
    let revindex = &SourmashRevIndex::as_rust(db_ptr);

    let result = cg.peek(threshold_bp as usize);

    if let Some((dataset_id, _match_size)) = result {
        let match_sig = revindex
            .collection()
            .sig_for_dataset(dataset_id)
            .expect("cannot load signature");
        Ok(SourmashSignature::from_rust(match_sig.into()))
    } else {
        Ok(SourmashSignature::from_rust(Signature::default()))
    }
}
}

// retrieve all signatures for a CounterGather.

ffi_fn! {
unsafe fn revindex_countergather_signatures(
    cg_ptr: *const SourmashRevIndex_CounterGather,
    db_ptr: *const SourmashRevIndex,
    size: *mut usize,
) -> Result<*mut *mut SourmashSignature> {
    let cg: &CounterGather = SourmashRevIndex_CounterGather::as_rust(cg_ptr);
    let revindex = &SourmashRevIndex::as_rust(db_ptr);

    let coll = revindex.collection();
    let sigs: Vec<Signature> = cg
        .dataset_ids()
        .into_iter()
        .map(|idx| { coll
                     .sig_for_dataset(idx)
                     .expect("cannot retrieve sig!?")
                     .into()
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

// retrieve all hashes present in a CounterGather. Can be done efficiently.

ffi_fn! {
unsafe fn revindex_countergather_found_hashes(
    cg_ptr: *mut SourmashRevIndex_CounterGather,
    template_ptr: *const SourmashKmerMinHash,
) -> Result<*const SourmashKmerMinHash> {
    let cg: &CounterGather = SourmashRevIndex_CounterGather::as_rust_mut(cg_ptr);
    let template_mh = SourmashKmerMinHash::as_rust(template_ptr);

    let found_mh = cg.found_hashes(template_mh);
    Ok(SourmashKmerMinHash::from_rust(found_mh))
}
}

ffi_fn! {
unsafe fn revindex_countergather_len(
    cg_ptr: *mut SourmashRevIndex_CounterGather,
) -> Result<u64> {
    let cg: &CounterGather = SourmashRevIndex_CounterGather::as_rust_mut(cg_ptr);

    Ok(cg.len() as u64)
}
}

// convert a sketch template into a Selection, for use by the Rust layer.
// TODO: remove this when it is possible to pass Selection thru the FFI

pub fn from_template(template: &Sketch) -> Selection {
    let (num, scaled) = match template {
        Sketch::MinHash(mh) => (mh.num(), mh.scaled()),
        Sketch::LargeMinHash(mh) => (mh.num(), mh.scaled()),
        _ => unimplemented!(),
    };

    let (ksize, moltype) = match template {
        Sketch::MinHash(mh) => (mh.ksize() as u32, mh.hash_function()),
        Sketch::LargeMinHash(mh) => (mh.ksize() as u32, mh.hash_function()),
        _ => unimplemented!(),
    };

    let adj_ksize: u32 = match moltype {
        HashFunctions::Murmur64Dna => ksize,
        HashFunctions::Murmur64Protein => ksize / 3,
        HashFunctions::Murmur64Dayhoff => ksize / 3,
        HashFunctions::Murmur64Hp => ksize / 3,
        HashFunctions::Murmur64Skipm1n3 => ksize,
        HashFunctions::Murmur64Skipm2n3 => ksize,
        _ => ksize,
    };

    Selection::builder()
        .ksize(adj_ksize)
        .num(num)
        .scaled(scaled)
        .build()
}

// build a new MemRevIndex from a list of sigs.

ffi_fn! {
unsafe fn revindex_mem_new_with_sigs(
    search_sigs_ptr: *const *const SourmashSignature,
    insigs: usize,
    template_ptr: *const SourmashKmerMinHash,
) -> Result<*mut SourmashRevIndex> {
    let search_sigs: Vec<Signature> = {
        assert!(!search_sigs_ptr.is_null());
        slice::from_raw_parts(search_sigs_ptr, insigs)
            .iter()
            .map(|sig| SourmashSignature::as_rust(*sig))
            .cloned()
            .collect()
    };

    let template = {
        assert!(!template_ptr.is_null());
        //TODO: avoid clone here
        Sketch::MinHash(SourmashKmerMinHash::as_rust(template_ptr).clone())
    };

    let selection = from_template(&template);
    let revindex = mem_revindex::MemRevIndex::new_with_sigs(search_sigs, &selection, 0, None).expect("cannot create MemRevIndex");
    Ok(SourmashRevIndex::from_rust(revindex))
}
}
