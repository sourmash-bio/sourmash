use crate::manifest::Manifest;

use crate::ffi::utils::ForeignObject;

pub struct SourmashManifest;

impl ForeignObject for SourmashManifest {
    type RustObject = Manifest;
}
