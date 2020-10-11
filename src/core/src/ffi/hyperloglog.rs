use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::signature::SigsTrait;
use crate::sketch::hyperloglog::HyperLogLog;

use crate::ffi::utils::ForeignObject;

pub struct SourmashHyperLogLog;

impl ForeignObject for SourmashHyperLogLog {
    type RustObject = HyperLogLog;
}

#[no_mangle]
pub unsafe extern "C" fn hll_new() -> *mut SourmashHyperLogLog {
    SourmashHyperLogLog::from_rust(HyperLogLog::default())
}

#[no_mangle]
pub unsafe extern "C" fn hll_free(ptr: *mut SourmashHyperLogLog) {
    SourmashHyperLogLog::drop(ptr);
}

ffi_fn! {
unsafe fn hll_with_error_rate(
    error_rate: f64,
    ksize: usize,
) -> Result<*mut SourmashHyperLogLog> {
    let hll = HyperLogLog::with_error_rate(error_rate, ksize)?;
    Ok(SourmashHyperLogLog::from_rust(hll))
}
}

#[no_mangle]
pub unsafe extern "C" fn hll_ksize(ptr: *const SourmashHyperLogLog) -> usize {
    SourmashHyperLogLog::as_rust(ptr).ksize()
}

#[no_mangle]
pub unsafe extern "C" fn hll_cardinality(ptr: *const SourmashHyperLogLog) -> usize {
    SourmashHyperLogLog::as_rust(ptr).cardinality()
}

ffi_fn! {
unsafe fn hll_add_sequence(
  ptr: *mut SourmashHyperLogLog,
  sequence: *const c_char,
  insize: usize,
  force: bool
) -> Result<()> {

    let hll = SourmashHyperLogLog::as_rust_mut(ptr);

    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(sequence as *mut u8, insize)
    };

    hll.add_sequence(buf, force)
}
}

#[no_mangle]
pub unsafe extern "C" fn hll_add_hash(ptr: *mut SourmashHyperLogLog, hash: u64) {
    let hll = SourmashHyperLogLog::as_rust_mut(ptr);
    hll.add_hash(hash);
}
