use crate::manifest::{Manifest, Record};

use crate::ffi::utils::{ForeignObject, SourmashStr};

pub struct SourmashManifest;

impl ForeignObject for SourmashManifest {
    type RustObject = Manifest;
}

pub struct ManifestRowIterator {
    iter: Box<dyn Iterator<Item = &'static Record>>,
}

pub struct SourmashManifestRowIter;

impl ForeignObject for SourmashManifestRowIter {
    type RustObject = ManifestRowIterator;
}

#[no_mangle]
pub unsafe extern "C" fn manifest_rows_iter_next(
    ptr: *mut SourmashManifestRowIter,
) -> *const SourmashManifestRow {
    let iterator = SourmashManifestRowIter::as_rust_mut(ptr);

    match iterator.iter.next() {
        Some(row) => SourmashManifestRow::from_rust(row.into()),
        None => std::ptr::null(),
    }
}

#[no_mangle]
pub unsafe extern "C" fn manifest_rows(
    ptr: *const SourmashManifest,
) -> *mut SourmashManifestRowIter {
    let manifest = SourmashManifest::as_rust(ptr);

    let iter = Box::new(manifest.iter());
    SourmashManifestRowIter::from_rust(ManifestRowIterator { iter })
}

#[repr(C)]
pub struct SourmashManifestRow {
    pub ksize: u32,
    pub with_abundance: u8,
    pub md5: SourmashStr,
    pub internal_location: SourmashStr,
    pub name: SourmashStr,
    pub moltype: SourmashStr,
}

impl ForeignObject for SourmashManifestRow {
    type RustObject = SourmashManifestRow;
}

impl From<&Record> for SourmashManifestRow {
    fn from(record: &Record) -> SourmashManifestRow {
        Self {
            ksize: record.ksize(),
            with_abundance: record.with_abundance() as u8,
            md5: record.md5().into(),
            name: record.name().into(),
            moltype: record.moltype().into(),
            internal_location: record
                .internal_location()
                .to_str()
                .unwrap()
                .to_owned()
                .into(),
        }
    }
}
