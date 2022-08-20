use std::os::raw::c_char;
use std::slice;

use crate::ffi::utils::{ForeignObject, SourmashStr};
use crate::prelude::*;
use crate::storage::ZipStorage;

pub struct SourmashZipStorage;

impl ForeignObject for SourmashZipStorage {
    type RustObject = ZipStorage;
}

ffi_fn! {
unsafe fn zipstorage_new(ptr: *const c_char, insize: usize) -> Result<*mut SourmashZipStorage> {
    let path = {
        assert!(!ptr.is_null());
        let path = slice::from_raw_parts(ptr as *mut u8, insize);
        std::str::from_utf8(path)?
    };
    let zipstorage = ZipStorage::from_file(path)?;

    Ok(SourmashZipStorage::from_rust(zipstorage))
}
}

#[no_mangle]
pub unsafe extern "C" fn zipstorage_free(ptr: *mut SourmashZipStorage) {
    SourmashZipStorage::drop(ptr);
}

ffi_fn! {
unsafe fn zipstorage_load(ptr: *const SourmashZipStorage,
    path_ptr: *const c_char,
    insize: usize,
    size: *mut usize) -> Result<*const u8> {

    let storage = SourmashZipStorage::as_rust(ptr);

    let path = {
        assert!(!path_ptr.is_null());
        let path = slice::from_raw_parts(path_ptr as *mut u8, insize);
        std::str::from_utf8(path)?
    };

    let buffer = storage.load(path)?;

    let b = buffer.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *const u8)
}
}

ffi_fn! {
unsafe fn zipstorage_list_sbts(
    ptr: *const SourmashZipStorage,
    size: *mut usize,
) -> Result<*mut *mut SourmashStr> {
    let storage = SourmashZipStorage::as_rust(ptr);

    let sbts = storage.list_sbts()?;

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*mut SourmashStr> = sbts
        .into_iter()
        .map(|x| Box::into_raw(Box::new(SourmashStr::from_string(x))) as *mut SourmashStr)
        .collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashStr)
}
}

ffi_fn! {
unsafe fn zipstorage_filenames(
    ptr: *const SourmashZipStorage,
    size: *mut usize,
) -> Result<*mut *mut SourmashStr> {
    let storage = SourmashZipStorage::as_rust(ptr);

    let files = storage.filenames()?;

    // FIXME: use the ForeignObject trait, maybe define new method there...
    let ptr_sigs: Vec<*mut SourmashStr> = files
        .into_iter()
        .map(|x| Box::into_raw(Box::new(SourmashStr::from_string(x))) as *mut SourmashStr)
        .collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut SourmashStr)
}
}

ffi_fn! {
unsafe fn zipstorage_set_subdir(
    ptr: *mut SourmashZipStorage,
    path_ptr: *const c_char,
    insize: usize,
) -> Result<()> {
    let storage = SourmashZipStorage::as_rust_mut(ptr);

    let path = {
        assert!(!path_ptr.is_null());
        let path = slice::from_raw_parts(path_ptr as *mut u8, insize);
        std::str::from_utf8(path)?
    };

    storage.set_subdir(path.to_string());
    Ok(())
}
}

ffi_fn! {
unsafe fn zipstorage_path(ptr: *const SourmashZipStorage) -> Result<SourmashStr> {
    let storage = SourmashZipStorage::as_rust(ptr);

    if let Some(ref path) = storage.path() {
        Ok(path.clone().into())
    } else {
        Ok("".into())
    }
}
}

ffi_fn! {
unsafe fn zipstorage_subdir(ptr: *const SourmashZipStorage) -> Result<SourmashStr> {
    let storage = SourmashZipStorage::as_rust(ptr);

    if let Some(ref path) = storage.subdir() {
        Ok(path.clone().into())
    } else {
        Ok("".into())
    }
}
}
