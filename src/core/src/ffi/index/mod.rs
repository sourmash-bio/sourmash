pub mod revindex;

use crate::index::{Selection, SigStore};
use crate::signature::Signature;

use crate::ffi::signature::SourmashSignature;
use crate::ffi::utils::{ForeignObject, SourmashStr};

pub struct SourmashSearchResult;

impl ForeignObject for SourmashSearchResult {
    type RustObject = (f64, Signature, String);
}

#[no_mangle]
pub unsafe extern "C" fn searchresult_free(ptr: *mut SourmashSearchResult) {
    SourmashSearchResult::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn searchresult_score(ptr: *const SourmashSearchResult) -> f64 {
    let result = SourmashSearchResult::as_rust(ptr);
    result.0
}

#[no_mangle]
pub unsafe extern "C" fn searchresult_filename(ptr: *const SourmashSearchResult) -> SourmashStr {
    let result = SourmashSearchResult::as_rust(ptr);
    (result.2).clone().into()
}

#[no_mangle]
pub unsafe extern "C" fn searchresult_signature(
    ptr: *const SourmashSearchResult,
) -> *mut SourmashSignature {
    let result = SourmashSearchResult::as_rust(ptr);
    SourmashSignature::from_rust((result.1).clone())
}

//================================================================

pub struct SourmashSelection;

impl ForeignObject for SourmashSelection {
    type RustObject = Selection;
}

#[no_mangle]
pub unsafe extern "C" fn selection_new() -> *mut SourmashSelection {
    SourmashSelection::from_rust(Selection::default())
}

#[no_mangle]
pub unsafe extern "C" fn selection_ksize(ptr: *const SourmashSelection) -> u32 {
    let sel = SourmashSelection::as_rust(ptr);
    if let Some(ksize) = sel.ksize() {
        ksize
    } else {
        todo!("empty ksize case not supported yet")
    }
}

#[no_mangle]
pub unsafe extern "C" fn selection_set_ksize(ptr: *mut SourmashSelection, new_ksize: u32) {
    let sel = SourmashSelection::as_rust_mut(ptr);
    sel.set_ksize(new_ksize);
}

//================================================================
//
pub struct SignatureIterator {
    iter: Box<dyn Iterator<Item = SigStore>>,
}

pub struct SourmashSignatureIter;

impl ForeignObject for SourmashSignatureIter {
    type RustObject = SignatureIterator;
}

#[no_mangle]
pub unsafe extern "C" fn signatures_iter_next(
    ptr: *mut SourmashSignatureIter,
) -> *const SourmashSignature {
    let iterator = SourmashSignatureIter::as_rust_mut(ptr);

    match iterator.iter.next() {
        Some(sig) => SourmashSignature::from_rust(sig.into()),
        None => std::ptr::null(),
    }
}
