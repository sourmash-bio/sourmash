use std::path::PathBuf;
use std::slice;
use std::sync::Arc;

use crate::index::revindex::{LinearRevIndex, RevIndex};
use crate::index::Index;
use crate::manifest::Manifest;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use crate::sketch::Sketch;
use crate::storage::Storage;

use crate::ffi::index::{
    SignatureIterator, SourmashSearchResult, SourmashSelection, SourmashSignatureIter,
};
use crate::ffi::manifest::SourmashManifest;
use crate::ffi::minhash::SourmashKmerMinHash;
use crate::ffi::signature::SourmashSignature;
use crate::ffi::storage::SourmashZipStorage;
use crate::ffi::utils::{ForeignObject, SourmashStr};

pub struct SourmashRevIndex;

impl ForeignObject for SourmashRevIndex {
    type RustObject = RevIndex;
}

ffi_fn! {
unsafe fn revindex_new_with_paths(
    search_sigs_ptr: *const *const SourmashStr,
    insigs: usize,
    template_ptr: *const SourmashKmerMinHash,
    threshold: usize,
    queries_ptr: *const *const SourmashKmerMinHash,
    inqueries: usize,
    keep_sigs: bool,
) -> Result<*mut SourmashRevIndex> {
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

    let revindex = RevIndex::new(
        search_sigs.as_ref(),
        &template,
        threshold,
        queries,
        keep_sigs,
    );
    Ok(SourmashRevIndex::from_rust(revindex))
}
}

ffi_fn! {
unsafe fn revindex_new_with_sigs(
    search_sigs_ptr: *const *const SourmashSignature,
    insigs: usize,
    template_ptr: *const SourmashKmerMinHash,
    threshold: usize,
    queries_ptr: *const *const SourmashKmerMinHash,
    inqueries: usize,
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
    let revindex = RevIndex::new_with_sigs(search_sigs, &template, threshold, queries);
    Ok(SourmashRevIndex::from_rust(revindex))
}
}

#[no_mangle]
pub unsafe extern "C" fn revindex_free(ptr: *mut SourmashRevIndex) {
    SourmashRevIndex::drop(ptr);
}

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
pub unsafe extern "C" fn revindex_scaled(ptr: *const SourmashRevIndex) -> u64 {
    let revindex = SourmashRevIndex::as_rust(ptr);
    if let Sketch::MinHash(mh) = revindex.template() {
        mh.scaled()
    } else {
        unimplemented!()
    }
}

#[no_mangle]
pub unsafe extern "C" fn revindex_len(ptr: *const SourmashRevIndex) -> u64 {
    let revindex = SourmashRevIndex::as_rust(ptr);
    revindex.len() as u64
}

ffi_fn! {
unsafe fn revindex_signatures(
    ptr: *const SourmashRevIndex,
    size: *mut usize,
) -> Result<*mut *mut SourmashSignature> {
    let revindex = SourmashRevIndex::as_rust(ptr);

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

//--------------------------------------------------

pub struct SourmashLinearIndex;

impl ForeignObject for SourmashLinearIndex {
    type RustObject = LinearRevIndex;
}

ffi_fn! {
unsafe fn linearindex_new(
    storage_ptr: *mut SourmashZipStorage,
    manifest_ptr: *mut SourmashManifest,
    selection_ptr: *mut SourmashSelection,
    use_manifest: bool,
) -> Result<*mut SourmashLinearIndex> {
    let storage = Arc::try_unwrap(*SourmashZipStorage::into_rust(storage_ptr)).ok().unwrap();

    let manifest = if manifest_ptr.is_null() {
        if use_manifest {
        // Load manifest from zipstorage
            Some(Manifest::from_reader(storage.load("SOURMASH-MANIFEST.csv")?.as_slice())?)
        } else {
            None
        }
    } else {
        Some(*SourmashManifest::into_rust(manifest_ptr))
    };

    let _selection = if !selection_ptr.is_null() {
        Some(SourmashSelection::into_rust(selection_ptr))
    } else {
        None
    };
    // TODO: how to extract a template? Probably from selection?
    let max_hash = max_hash_for_scaled(100);
    let template = Sketch::MinHash(
        KmerMinHash::builder()
            .num(0u32)
            .ksize(57)
            .hash_function(crate::encodings::HashFunctions::murmur64_protein)
            .max_hash(max_hash)
            .build(),
    );

    /*
    def __init__(self, storage, *, selection_dict=None,
                 traverse_yield_all=False, manifest=None, use_manifest=True):

        sig_files: Manifest,
        template: &Sketch,
        keep_sigs: bool,
        ref_sigs: Option<Vec<Signature>>,
        storage: Option<ZipStorage>,
    */

    let linear_index = LinearRevIndex::new(manifest, &template, false, None, Some(storage));

    Ok(SourmashLinearIndex::from_rust(linear_index))
}
}

#[no_mangle]
pub unsafe extern "C" fn linearindex_free(ptr: *mut SourmashLinearIndex) {
    SourmashLinearIndex::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn linearindex_manifest(
    ptr: *const SourmashLinearIndex,
) -> *const SourmashManifest {
    let index = SourmashLinearIndex::as_rust(ptr);
    SourmashManifest::from_rust(index.manifest())
}

ffi_fn! {
unsafe fn linearindex_set_manifest(
    ptr: *mut SourmashLinearIndex,
    manifest_ptr: *mut SourmashManifest,
) -> Result<()> {
    let index = SourmashLinearIndex::as_rust_mut(ptr);
    let manifest = SourmashManifest::into_rust(manifest_ptr);

    index.set_manifest(*manifest)?;
    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn linearindex_len(ptr: *const SourmashLinearIndex) -> u64 {
    let index = SourmashLinearIndex::as_rust(ptr);
    index.len() as u64
}

#[no_mangle]
pub unsafe extern "C" fn linearindex_location(ptr: *const SourmashLinearIndex) -> SourmashStr {
    let index = SourmashLinearIndex::as_rust(ptr);
    match index.location() {
        Some(x) => x,
        None => "".into(),
    }
    .into()
}

#[no_mangle]
pub unsafe extern "C" fn linearindex_storage(
    ptr: *const SourmashLinearIndex,
) -> *const SourmashZipStorage {
    let index = SourmashLinearIndex::as_rust(ptr);
    let storage = index.storage();

    match storage {
        Some(st) => SourmashZipStorage::from_rust(st),
        None => std::ptr::null::<SourmashZipStorage>(),
    }
}

#[no_mangle]
pub unsafe extern "C" fn linearindex_signatures(
    ptr: *const SourmashLinearIndex,
) -> *mut SourmashSignatureIter {
    let index = SourmashLinearIndex::as_rust(ptr);

    let iter = Box::new(index.signatures_iter());
    SourmashSignatureIter::from_rust(SignatureIterator { iter })
}

ffi_fn! {
unsafe fn linearindex_select(
    ptr: *mut SourmashLinearIndex,
    selection_ptr: *const SourmashSelection,
) -> Result<*mut SourmashLinearIndex> {
    let index = SourmashLinearIndex::into_rust(ptr);
    let selection = SourmashSelection::as_rust(selection_ptr);

    let new_index = index.select(selection)?;
    Ok(SourmashLinearIndex::from_rust(new_index))
}
}
