use std::slice;

use crate::encodings::*;
use crate::ffi::index::disk_revindex::{
    retrieve_picklist, SourmashDatasetPicklist, SourmashRevIndex_CounterGather,
};
use crate::ffi::index::SourmashSearchResult;
use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::signature::SourmashSignature;
use crate::ffi::utils::ForeignObject;
use crate::index::revindex::mem_revindex;
// use crate::index::revindex::disk_revindex;
use crate::index::revindex::{CounterGather, RevIndexOps};
use crate::index::Index;
use crate::prelude::*;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
use crate::ScaledType;

pub struct SourmashMemRevIndex;

impl ForeignObject for SourmashMemRevIndex {
    type RustObject = mem_revindex::RevIndex;
}

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

/*
ffi_fn! {
unsafe fn revindex_new_with_paths(
    search_sigs_ptr: *const *const SourmashStr,
    insigs: usize,
    template_ptr: *const SourmashKmerMinHash,
    threshold: usize,
    queries_ptr: *const *const SourmashKmerMinHash,
    inqueries: usize,
    keep_sigs: bool,
) -> Result<*mut SourmashMemRevIndex> {
    let search_sigs: Vec<PathBuf> = {
        assert!(!search_sigs_ptr.is_null());
        slice::from_raw_parts(search_sigs_ptr, insigs)
            .iter()
            .map(|path| {
                let mut new_path = PathBuf::new();
                new_path.push(SourmashStr::as_rust(*path).as_str());
                new_path
            })
            .collect()
    };

    let template = {
        assert!(!template_ptr.is_null());
        //TODO: avoid clone here
        Sketch::MinHash(SourmashKmerMinHash::as_rust(template_ptr).clone())
    };

    let queries_vec: Vec<KmerMinHash>;
    let queries: Option<&[KmerMinHash]> = if queries_ptr.is_null() {
        None
    } else {
        queries_vec = slice::from_raw_parts(queries_ptr, inqueries)
            .iter()
            .map(|mh_ptr|
            // TODO: avoid this clone
          SourmashKmerMinHash::as_rust(*mh_ptr).clone())
            .collect();
        Some(queries_vec.as_ref())
    };

    let selection = from_template(&template);

    let revindex = mem_revindex::RevIndex::new(
        search_sigs.as_ref(),
        &selection,
        threshold,
        queries,
        keep_sigs,
    )?;
    Ok(SourmashMemRevIndex::from_rust(revindex))
}
}
*/

ffi_fn! {
unsafe fn revindex_new_with_sigs(
    search_sigs_ptr: *const *const SourmashSignature,
    insigs: usize,
    template_ptr: *const SourmashKmerMinHash,
) -> Result<*mut SourmashMemRevIndex> {
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
    let revindex = mem_revindex::RevIndex::new_with_sigs(search_sigs, &selection, 0, None)?;
    Ok(SourmashMemRevIndex::from_rust(revindex))
}
}

#[no_mangle]
pub unsafe extern "C" fn revindex_free(ptr: *mut SourmashMemRevIndex) {
    SourmashMemRevIndex::drop(ptr);
}

ffi_fn! {
unsafe fn revindex_search(
    ptr: *const SourmashMemRevIndex,
    sig_ptr: *const SourmashSignature,
    threshold: f64,
    do_containment: bool,
    _ignore_abundance: bool,
    size: *mut usize,
) -> Result<*const *const SourmashSearchResult> {
    let revindex = SourmashMemRevIndex::as_rust(ptr);
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
    ptr: *const SourmashMemRevIndex,
    sig_ptr: *const SourmashSignature,
    threshold: f64,
    _do_containment: bool,
    _ignore_abundance: bool,
    size: *mut usize,
) -> Result<*const *const SourmashSearchResult> {
    let revindex = SourmashMemRevIndex::as_rust(ptr);
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

    let mut cg = revindex.prepare_gather_counters(mh, None);

    let results: Vec<(f64, Signature, String)> = revindex
        .gather(&mut cg, threshold, mh, None)
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
pub unsafe extern "C" fn revindex_scaled(ptr: *const SourmashMemRevIndex) -> ScaledType {
    let revindex = SourmashMemRevIndex::as_rust(ptr);
    if let Sketch::MinHash(mh) = revindex.template() {
        mh.scaled()
    } else {
        unimplemented!()
    }
}

#[no_mangle]
pub unsafe extern "C" fn revindex_len(ptr: *const SourmashMemRevIndex) -> u64 {
    let revindex = SourmashMemRevIndex::as_rust(ptr);
    revindex.len() as u64
}

ffi_fn! {
unsafe fn revindex_signatures(
    ptr: *const SourmashMemRevIndex,
    size: *mut usize,
) -> Result<*mut *mut SourmashSignature> {
    let revindex = SourmashMemRevIndex::as_rust(ptr);

    let sigs = revindex.signatures();

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

// return a CounterGather object with prefetch results

ffi_fn! {
unsafe fn mem_revindex_prefetch_to_countergather(
    db_ptr: *const SourmashMemRevIndex,
    query_ptr: *const SourmashSignature,
    dataset_picklist_ptr: *const SourmashDatasetPicklist,
) -> Result<*mut SourmashRevIndex_CounterGather> {
    let revindex: &mem_revindex::RevIndex = SourmashMemRevIndex::as_rust(db_ptr);
    let sig = SourmashSignature::as_rust(query_ptr);

    // extract KmerMinHash for query
    let query_mh: KmerMinHash = sig.clone()
        .try_into().expect("cannot get kmerminhash");

    // picklist?
    let dataset_picklist = retrieve_picklist(dataset_picklist_ptr);

    // do search & get matches - @CTB picklist needed!
    let counter = revindex.prepare_gather_counters(&query_mh, dataset_picklist);

    Ok(SourmashRevIndex_CounterGather::from_rust(counter))
}
}

ffi_fn! {
unsafe fn mem_revindex_countergather_consume(
    cg_ptr: *mut SourmashRevIndex_CounterGather,
    isect_ptr: *const SourmashKmerMinHash,
) -> Result<()> {
    let cg: &mut CounterGather = SourmashRevIndex_CounterGather::as_rust_mut(cg_ptr);
    let isect_mh = SourmashKmerMinHash::as_rust(isect_ptr);

    cg.consume(isect_mh);

    Ok(())
}
}

ffi_fn! {
unsafe fn mem_revindex_countergather_peek(
    cg_ptr: *const SourmashRevIndex_CounterGather,
    db_ptr: *const SourmashMemRevIndex,
    threshold_bp: u64,
) -> Result<*mut SourmashSignature> {
    let cg: &CounterGather = SourmashRevIndex_CounterGather::as_rust(cg_ptr);
    let revindex: &mem_revindex::RevIndex = SourmashMemRevIndex::as_rust(db_ptr);

    let result = cg.peek(threshold_bp as usize);

    // if result.is_none() { // @CTB...
    // }

    let (dataset_id, _match_size) = result.unwrap();

    let match_sig = revindex.collection().sig_for_dataset(dataset_id)?;

    Ok(SourmashSignature::from_rust(match_sig.into()))
}
}

ffi_fn! {
unsafe fn mem_revindex_countergather_signatures(
    cg_ptr: *const SourmashRevIndex_CounterGather,
    db_ptr: *const SourmashMemRevIndex,
    size: *mut usize,
) -> Result<*mut *mut SourmashSignature> {
    let cg: &CounterGather = SourmashRevIndex_CounterGather::as_rust(cg_ptr);
    let revindex: &mem_revindex::RevIndex = SourmashMemRevIndex::as_rust(db_ptr);

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

ffi_fn! {
unsafe fn mem_revindex_countergather_found_hashes(
    cg_ptr: *mut SourmashRevIndex_CounterGather,
    template_ptr: *const SourmashKmerMinHash,
) -> Result<*const SourmashKmerMinHash> {
    let cg: &mut CounterGather = SourmashRevIndex_CounterGather::as_rust_mut(cg_ptr);
    let template_mh = SourmashKmerMinHash::as_rust(template_ptr);

    let found_mh = cg.found_hashes(template_mh);
    Ok(SourmashKmerMinHash::from_rust(found_mh))
}
}

ffi_fn! {
unsafe fn mem_revindex_countergather_len(
    cg_ptr: *mut SourmashRevIndex_CounterGather,
) -> Result<u64> {
    let cg: &mut CounterGather = SourmashRevIndex_CounterGather::as_rust_mut(cg_ptr);

    Ok(cg.len() as u64)
}
}
