use std::os::raw::c_char;
use std::slice;

use crate::picklist::Picklist;

use crate::ffi::utils::ForeignObject;

pub struct SourmashPicklist;

impl ForeignObject for SourmashPicklist {
    type RustObject = Picklist;
}

#[no_mangle]
pub unsafe extern "C" fn picklist_new() -> *mut SourmashPicklist {
    SourmashPicklist::from_rust(Picklist::default())
}

#[no_mangle]
pub unsafe extern "C" fn picklist_free(ptr: *mut SourmashPicklist) {
    SourmashPicklist::drop(ptr);
}

ffi_fn! {
unsafe fn picklist_set_coltype(
    ptr: *mut SourmashPicklist,
    coltype_ptr: *const c_char,
    insize: usize,
) -> Result<()> {
    let coltype = {
        assert!(!coltype_ptr.is_null());
        let coltype = slice::from_raw_parts(coltype_ptr as *mut u8, insize);
        std::str::from_utf8(coltype)?
    };
    let pl = SourmashPicklist::as_rust_mut(ptr);
    pl.set_coltype(coltype.to_string());

    Ok(())
}
}

ffi_fn! {
unsafe fn picklist_set_pickfile(
    ptr: *mut SourmashPicklist,
    prop_ptr: *const c_char,
    insize: usize,
) -> Result<()> {
    let prop = {
        assert!(!prop_ptr.is_null());
        let prop = slice::from_raw_parts(prop_ptr as *mut u8, insize);
        std::str::from_utf8(prop)?
    };
    let pl = SourmashPicklist::as_rust_mut(ptr);
    pl.set_pickfile(prop.to_string());

    Ok(())
}
}

ffi_fn! {
unsafe fn picklist_set_column_name(
    ptr: *mut SourmashPicklist,
    prop_ptr: *const c_char,
    insize: usize,
) -> Result<()> {
    let prop = {
        assert!(!prop_ptr.is_null());
        let prop = slice::from_raw_parts(prop_ptr as *mut u8, insize);
        std::str::from_utf8(prop)?
    };
    let pl = SourmashPicklist::as_rust_mut(ptr);
    pl.set_column_name(prop.to_string());

    Ok(())
}
}
