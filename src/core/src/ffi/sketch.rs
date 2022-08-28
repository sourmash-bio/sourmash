use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::ffi::utils::ForeignObject;
use crate::signature::SigsTrait;
use crate::sketch::Sketch;

pub struct SourmashSketch;

impl ForeignObject for SourmashSketch {
    type RustObject = Sketch;
}
