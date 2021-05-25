use crate::index::{JaccardSearch, SearchType};

use crate::ffi::utils::ForeignObject;

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
