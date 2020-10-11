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

/*
ffi_fn! {
unsafe fn nodegraph_from_path(filename: *const c_char) -> Result<*mut SourmashNodegraph> {
    // FIXME use buffer + len instead of c_str
    let c_str = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };

    let (mut input, _) = niffler::from_path(c_str.to_str()?)?;
    let ng = Nodegraph::from_reader(&mut input)?;

    Ok(SourmashNodegraph::from_rust(ng))
}
}

ffi_fn! {
unsafe fn nodegraph_from_buffer(ptr: *const c_char, insize: usize) -> Result<*mut SourmashNodegraph> {
    // FIXME use SourmashSlice_u8?
    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr as *mut u8, insize)
    };

    let ng = Nodegraph::from_reader(&mut &buf[..])?;

    Ok(SourmashNodegraph::from_rust(ng))
}
}

ffi_fn! {
unsafe fn nodegraph_save(ptr: *const SourmashNodegraph, filename: *const c_char) -> Result<()> {
    let ng = SourmashNodegraph::as_rust(ptr);

    // FIXME use buffer + len instead of c_str
    let c_str = {
        assert!(!filename.is_null());

        CStr::from_ptr(filename)
    };

    ng.save(c_str.to_str()?)?;

    Ok(())
}
}

ffi_fn! {
unsafe fn nodegraph_to_buffer(ptr: *const SourmashNodegraph, compression: u8, size: *mut usize) -> Result<*const u8> {
    let ng = SourmashNodegraph::as_rust(ptr);

    let mut buffer = vec![];
    {
      let mut writer = if compression > 0 {
          let level = match compression {
            1 => niffler::compression::Level::One,
            2 => niffler::compression::Level::Two,
            3 => niffler::compression::Level::Three,
            4 => niffler::compression::Level::Four,
            5 => niffler::compression::Level::Five,
            6 => niffler::compression::Level::Six,
            7 => niffler::compression::Level::Seven,
            8 => niffler::compression::Level::Eight,
            _ => niffler::compression::Level::Nine,
          };

          niffler::get_writer(Box::new(&mut buffer),
                              niffler::compression::Format::Gzip,
                              level)?
      } else {
          Box::new(&mut buffer)
      };
      ng.save_to_writer(&mut writer)?;
    }

    let b = buffer.into_boxed_slice();
    *size = b.len();

    // FIXME use SourmashSlice_u8?
    Ok(Box::into_raw(b) as *const u8)
}
}
*/
