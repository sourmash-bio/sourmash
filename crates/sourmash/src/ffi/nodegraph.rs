use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use niffler::get_input;

use crate::sketch::minhash::KmerMinHash;
use crate::sketch::nodegraph::Nodegraph;

#[no_mangle]
pub unsafe extern "C" fn nodegraph_new() -> *mut Nodegraph {
    Box::into_raw(Box::new(Nodegraph::default())) as _
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_free(ptr: *mut Nodegraph) {
    if ptr.is_null() {
        return;
    }
    Box::from_raw(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_with_tables(
    ksize: usize,
    starting_size: usize,
    n_tables: usize,
) -> *mut Nodegraph {
    Box::into_raw(Box::new(Nodegraph::with_tables(
        starting_size,
        n_tables,
        ksize,
    ))) as _
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_count(ptr: *mut Nodegraph, h: u64) -> bool {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    ng.count(h)
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_get(ptr: *mut Nodegraph, h: u64) -> usize {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    ng.get(h)
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_expected_collisions(ptr: *mut Nodegraph) -> f64 {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    ng.expected_collisions()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_ksize(ptr: *mut Nodegraph) -> usize {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    ng.ksize()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_tablesize(ptr: *mut Nodegraph) -> usize {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    ng.tablesize()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_ntables(ptr: *mut Nodegraph) -> usize {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    ng.ntables()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_noccupied(ptr: *mut Nodegraph) -> usize {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    ng.noccupied()
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_matches(ptr: *mut Nodegraph, mh_ptr: *mut KmerMinHash) -> usize {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let mh = {
        assert!(!ptr.is_null());
        &mut *mh_ptr
    };

    ng.matches(mh)
}

#[no_mangle]
pub unsafe extern "C" fn nodegraph_update(ptr: *mut Nodegraph, optr: *mut Nodegraph) {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let ong = {
        assert!(!optr.is_null());
        &mut *optr
    };

    ng.update(ong);
}

ffi_fn! {
unsafe fn nodegraph_from_path(filename: *const c_char) -> Result<*mut Nodegraph> {
    let c_str = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };

    let (mut input, _) = get_input(c_str.to_str()?)?;
    let ng = Nodegraph::from_reader(&mut input)?;

    Ok(Box::into_raw(Box::new(ng)))
}
}

ffi_fn! {
unsafe fn nodegraph_from_buffer(ptr: *const c_char, insize: usize) -> Result<*mut Nodegraph> {
    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr as *mut u8, insize)
    };

    let ng = Nodegraph::from_reader(&mut &buf[..])?;

    Ok(Box::into_raw(Box::new(ng)))
}
}

ffi_fn! {
unsafe fn nodegraph_save(ptr: *mut Nodegraph, filename: *const c_char) -> Result<()> {
    let ng = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let c_str = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };

    ng.save(c_str.to_str()?)?;

    Ok(())
}
}
