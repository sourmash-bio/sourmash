use crate::index::{JaccardSearch, SearchType};
use crate::signature::Signature;

use crate::ffi::signature::SourmashSignature;
use crate::ffi::utils::{ForeignObject, SourmashStr};

pub struct SourmashSearchFn;

impl ForeignObject for SourmashSearchFn {
    type RustObject = JaccardSearch;
}

#[no_mangle]
pub unsafe extern "C" fn searchfn_free(ptr: *mut SourmashSearchFn) {
    SourmashSearchFn::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn searchfn_new(
    search_type: SearchType,
    threshold: f64,
) -> *mut SourmashSearchFn {
    SourmashSearchFn::from_rust(JaccardSearch::with_threshold(search_type, threshold))
}

/*
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
*/
