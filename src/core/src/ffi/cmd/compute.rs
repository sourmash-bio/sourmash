use std::slice;

use crate::cmd::ComputeParameters;

#[no_mangle]
pub unsafe extern "C" fn computeparams_new() -> *mut ComputeParameters {
    Box::into_raw(Box::new(ComputeParameters::default())) as _
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_free(ptr: *mut ComputeParameters) {
    if ptr.is_null() {
        return;
    }
    Box::from_raw(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_seed(ptr: *mut ComputeParameters) -> u64 {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.seed
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_seed(ptr: *mut ComputeParameters, new_seed: u64) {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.seed = new_seed;
}

ffi_fn! {
unsafe fn computeparams_ksizes(ptr: *mut ComputeParameters, size: *mut usize) -> Result<*const u32> {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let output = cp.ksizes.clone();
    *size = output.len();

    Ok(Box::into_raw(output.into_boxed_slice()) as *const u32)
}
}

ffi_fn! {
unsafe fn computeparams_set_ksizes(
    ptr: *mut ComputeParameters,
    ksizes_ptr: *const u32,
    insize: usize,
  ) -> Result<()> {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let ksizes = {
        assert!(!ksizes_ptr.is_null());
        slice::from_raw_parts(ksizes_ptr as *const u32, insize)
    };

    cp.ksizes = ksizes.into();

    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_protein(ptr: *mut ComputeParameters) -> bool {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.protein
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_protein(ptr: *mut ComputeParameters, v: bool) {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.protein = v;
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_dayhoff(ptr: *mut ComputeParameters) -> bool {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.dayhoff
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_dayhoff(ptr: *mut ComputeParameters, v: bool) {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.dayhoff = v;
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_hp(ptr: *mut ComputeParameters) -> bool {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.hp
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_hp(ptr: *mut ComputeParameters, v: bool) {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.hp = v;
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_dna(ptr: *mut ComputeParameters) -> bool {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.dna
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_dna(ptr: *mut ComputeParameters, v: bool) {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.dna = v;
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_track_abundance(ptr: *mut ComputeParameters) -> bool {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.track_abundance
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_track_abundance(ptr: *mut ComputeParameters, v: bool) {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.track_abundance = v;
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_num_hashes(ptr: *mut ComputeParameters) -> u32 {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.num_hashes
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_num_hashes(ptr: *mut ComputeParameters, num: u32) {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.num_hashes = num;
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_scaled(ptr: *mut ComputeParameters) -> u64 {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.scaled
}

#[no_mangle]
pub unsafe extern "C" fn computeparams_set_scaled(ptr: *mut ComputeParameters, scaled: u64) {
    let cp = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    cp.scaled = scaled;
}
