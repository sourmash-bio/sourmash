pub mod revindex;

use crate::index::Selection;
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

pub struct SourmashSelection;

impl ForeignObject for SourmashSelection {
    type RustObject = Selection;
}

pub struct SignatureIterator {
    iter: Box<dyn Iterator<Item = Signature>>,
}

pub struct SourmashSignatureIter;

impl ForeignObject for SourmashSignatureIter {
    type RustObject = SignatureIterator;
}

#[no_mangle]
pub unsafe extern "C" fn signatures_iter_next(
    ptr: *mut SourmashSignatureIter,
) -> *const SourmashSignature {
    let mut iterator = SourmashSignatureIter::into_rust(ptr);

    match iterator.iter.next() {
        Some(sig) => SourmashSignature::from_rust(sig),
        None => std::ptr::null(),
    }
}
