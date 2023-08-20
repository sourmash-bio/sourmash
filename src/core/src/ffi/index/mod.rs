#[cfg(not(target_arch = "wasm32"))]
#[cfg(feature = "branchwater")]
pub mod revindex;

use crate::ffi::HashFunctions;
use crate::selection::Selection;
use crate::storage::SigStore;

use crate::signature::Signature;

use crate::ffi::picklist::SourmashPicklist;
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

#[no_mangle]
pub unsafe extern "C" fn selection_num(ptr: *const SourmashSelection) -> u32 {
    let sel = SourmashSelection::as_rust(ptr);
    if let Some(num) = sel.num() {
        num
    } else {
        todo!("empty num case not supported yet")
    }
}

#[no_mangle]
pub unsafe extern "C" fn selection_set_num(ptr: *mut SourmashSelection, new_num: u32) {
    let sel = SourmashSelection::as_rust_mut(ptr);
    sel.set_num(new_num);
}

#[no_mangle]
pub unsafe extern "C" fn selection_scaled(ptr: *const SourmashSelection) -> u32 {
    let sel = SourmashSelection::as_rust(ptr);
    if let Some(scaled) = sel.scaled() {
        scaled
    } else {
        todo!("empty scaled case not supported yet")
    }
}

#[no_mangle]
pub unsafe extern "C" fn selection_set_scaled(ptr: *mut SourmashSelection, new_scaled: u32) {
    let sel = SourmashSelection::as_rust_mut(ptr);
    sel.set_scaled(new_scaled);
}

#[no_mangle]
pub unsafe extern "C" fn selection_containment(ptr: *const SourmashSelection) -> bool {
    let sel = SourmashSelection::as_rust(ptr);
    if let Some(containment) = sel.containment() {
        containment
    } else {
        todo!("empty scaled case not supported yet")
    }
}

#[no_mangle]
pub unsafe extern "C" fn selection_set_containment(
    ptr: *mut SourmashSelection,
    new_containment: bool,
) {
    let sel = SourmashSelection::as_rust_mut(ptr);
    sel.set_containment(new_containment);
}

#[no_mangle]
pub unsafe extern "C" fn selection_abund(ptr: *const SourmashSelection) -> bool {
    let sel = SourmashSelection::as_rust(ptr);
    if let Some(abund) = sel.abund() {
        abund
    } else {
        todo!("empty abund case not supported yet")
    }
}

#[no_mangle]
pub unsafe extern "C" fn selection_set_abund(ptr: *mut SourmashSelection, new_abund: bool) {
    let sel = SourmashSelection::as_rust_mut(ptr);
    sel.set_abund(new_abund);
}

#[no_mangle]
pub unsafe extern "C" fn selection_moltype(ptr: *const SourmashSelection) -> HashFunctions {
    let sel = SourmashSelection::as_rust(ptr);
    if let Some(hash_function) = sel.moltype() {
        hash_function.into()
    } else {
        todo!("empty hash_function case not supported yet")
    }
}

#[no_mangle]
pub unsafe extern "C" fn selection_set_moltype(
    ptr: *mut SourmashSelection,
    new_moltype: HashFunctions,
) {
    let sel = SourmashSelection::as_rust_mut(ptr);
    sel.set_moltype(new_moltype.into());
}

#[no_mangle]
pub unsafe extern "C" fn selection_picklist(
    ptr: *const SourmashSelection,
) -> *const SourmashPicklist {
    let sel = SourmashSelection::as_rust(ptr);
    if let Some(picklist) = sel.picklist() {
        SourmashPicklist::from_rust(picklist)
    } else {
        todo!("empty picklist case not supported yet")
    }
}

#[no_mangle]
pub unsafe extern "C" fn selection_set_picklist(
    ptr: *mut SourmashSelection,
    new_picklist: *mut SourmashPicklist,
) {
    let sel = SourmashSelection::as_rust_mut(ptr);
    let pick = SourmashPicklist::into_rust(new_picklist);
    sel.set_picklist(*pick);
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
