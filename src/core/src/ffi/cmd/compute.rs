use std::slice;

use crate::cmd::ComputeParameters;

use crate::ffi::utils::ForeignObject;

pub struct SourmashComputeParameters;

impl ForeignObject for SourmashComputeParameters {
    type RustObject = ComputeParameters;
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_new() -> *mut SourmashComputeParameters {
    SourmashComputeParameters::from_rust(ComputeParameters::default())
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_free(ptr: *mut SourmashComputeParameters) {
    SourmashComputeParameters::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_seed(ptr: *const SourmashComputeParameters) -> u64 {
    let cp = SourmashComputeParameters::as_rust(ptr);
    cp.seed()
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_seed(
    ptr: *mut SourmashComputeParameters,
    new_seed: u64,
) {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);
    cp.set_seed(new_seed);
}

ffi_fn! {
unsafe fn computeparams_ksizes(ptr: *const SourmashComputeParameters, size: *mut usize) -> Result<*const u32> {
    let cp = SourmashComputeParameters::as_rust(ptr);
    let output = cp.ksizes().clone();
    *size = output.len();

    // FIXME use a SourmashSlice_u32?
    Ok(Box::into_raw(output.into_boxed_slice()) as *const u32)
}
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_ksizes_free(ptr: *mut u32, insize: usize) {
    // FIXME use a SourmashSlice_u32?
    if ptr.is_null() {
        return;
    }
    Vec::from_raw_parts(ptr as *mut u32, insize, insize);
}

ffi_fn! {
unsafe fn computeparams_set_ksizes(
    ptr: *mut SourmashComputeParameters,
    ksizes_ptr: *const u32,
    insize: usize,
  ) -> Result<()> {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);

    let ksizes = {
        assert!(!ksizes_ptr.is_null());
        slice::from_raw_parts(ksizes_ptr as *const u32, insize)
    };

    cp.set_ksizes(ksizes.into());

    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_protein(ptr: *const SourmashComputeParameters) -> bool {
    let cp = SourmashComputeParameters::as_rust(ptr);
    cp.protein()
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_protein(ptr: *mut SourmashComputeParameters, v: bool) {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);
    cp.set_protein(v);
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_dayhoff(ptr: *const SourmashComputeParameters) -> bool {
    let cp = SourmashComputeParameters::as_rust(ptr);
    cp.dayhoff()
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_dayhoff(ptr: *mut SourmashComputeParameters, v: bool) {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);
    cp.set_dayhoff(v);
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_hp(ptr: *const SourmashComputeParameters) -> bool {
    let cp = SourmashComputeParameters::as_rust(ptr);
    cp.hp()
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_hp(ptr: *mut SourmashComputeParameters, v: bool) {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);
    cp.set_hp(v);
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_dna(ptr: *const SourmashComputeParameters) -> bool {
    let cp = SourmashComputeParameters::as_rust(ptr);
    cp.dna()
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_dna(ptr: *mut SourmashComputeParameters, v: bool) {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);
    cp.set_dna(v);
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_track_abundance(
    ptr: *const SourmashComputeParameters,
) -> bool {
    let cp = SourmashComputeParameters::as_rust(ptr);
    cp.track_abundance()
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_track_abundance(
    ptr: *mut SourmashComputeParameters,
    v: bool,
) {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);
    cp.set_track_abundance(v);
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_num_hashes(ptr: *const SourmashComputeParameters) -> u32 {
    let cp = SourmashComputeParameters::as_rust(ptr);
    cp.num_hashes()
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_num_hashes(
    ptr: *mut SourmashComputeParameters,
    num: u32,
) {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);
    cp.set_num_hashes(num);
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_scaled(ptr: *const SourmashComputeParameters) -> u64 {
    let cp = SourmashComputeParameters::as_rust(ptr);
    cp.scaled()
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_scaled(
    ptr: *mut SourmashComputeParameters,
    scaled: u64,
) {
    let cp = SourmashComputeParameters::as_rust_mut(ptr);
    cp.set_scaled(scaled);
}
