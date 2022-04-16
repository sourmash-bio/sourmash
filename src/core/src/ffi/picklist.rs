use crate::picklist::Picklist;

use crate::ffi::utils::ForeignObject;

pub struct SourmashPicklist;

impl ForeignObject for SourmashPicklist{
    type RustObject = Picklist;
}
