use crate::manifest::{Manifest, Record};

use crate::ffi::utils::{ForeignObject, SourmashStr};

pub struct SourmashManifest;

impl ForeignObject for SourmashManifest {
    type RustObject = Manifest;
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn manifest_free(ptr: *mut SourmashManifest) {
    unsafe {
        SourmashManifest::drop(ptr);
    }
}

pub struct ManifestRowIterator {
    iter: Box<dyn Iterator<Item = &'static Record>>,
}

pub struct SourmashManifestRowIter;

impl ForeignObject for SourmashManifestRowIter {
    type RustObject = ManifestRowIterator;
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn manifest_rows_iter_next(
    ptr: *mut SourmashManifestRowIter,
) -> *const SourmashManifestRow {
    let iterator = unsafe { SourmashManifestRowIter::as_rust_mut(ptr) };

    match iterator.iter.next() {
        Some(row) => unsafe { SourmashManifestRow::from_rust(row.into()) },
        None => std::ptr::null(),
    }
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn manifest_rows(
    ptr: *const SourmashManifest,
) -> *mut SourmashManifestRowIter {
    let manifest = unsafe { SourmashManifest::as_rust(ptr) };

    let iter = Box::new(manifest.iter());
    unsafe { SourmashManifestRowIter::from_rust(ManifestRowIterator { iter }) }
}

#[repr(C)]
pub struct SourmashManifestRow {
    pub ksize: u32,
    pub with_abundance: bool,
    pub md5: SourmashStr,
    pub internal_location: SourmashStr,
    pub name: SourmashStr,
    pub moltype: SourmashStr,
    pub n_hashes: usize,
    pub num: u32,
    pub scaled: u32,
    pub filename: SourmashStr,
}

impl ForeignObject for SourmashManifestRow {
    type RustObject = SourmashManifestRow;
}

#[unsafe(no_mangle)]
pub unsafe extern "C" fn manifestrow_free(ptr: *mut SourmashManifestRow) {
    unsafe {
        SourmashManifestRow::drop(ptr);
    }
}

impl From<&Record> for SourmashManifestRow {
    fn from(record: &Record) -> SourmashManifestRow {
        Self {
            ksize: record.ksize(),
            with_abundance: record.with_abundance(),
            md5: record.md5().clone().into(),
            name: record.name().clone().into(),
            moltype: record.moltype().to_string().into(),
            internal_location: record.internal_location().to_string().into(),
            n_hashes: *record.n_hashes(),
            num: *record.num(),
            scaled: *record.scaled(),
            filename: record.filename().to_string().into(),
        }
    }
}
