use std::convert::TryInto;
use std::ffi::CStr;
use std::io;
use std::os::raw::c_char;
use std::slice;

use crate::cmd::ComputeParameters;
use crate::errors::SourmashError;
use crate::ffi::utils::SourmashStr;
use crate::signature::Signature;
use crate::sketch::minhash::{HashFunctions, KmerMinHash};
use crate::sketch::Sketch;

// Signature methods

#[no_mangle]
pub unsafe extern "C" fn signature_new() -> *mut Signature {
    Box::into_raw(Box::new(Signature::default())) as _
}

#[no_mangle]
pub unsafe extern "C" fn signature_from_params(ptr: *mut ComputeParameters) -> *mut Signature {
    let params = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    Box::into_raw(Box::new(Signature::from_params(params))) as _
}

#[no_mangle]
pub unsafe extern "C" fn signature_free(ptr: *mut Signature) {
    if ptr.is_null() {
        return;
    }
    Box::from_raw(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn signature_len(ptr: *mut Signature) -> usize {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    sig.signatures.len()
}

ffi_fn! {
unsafe fn signature_add_sequence(ptr: *mut Signature, sequence: *const c_char, force: bool) ->
    Result<()> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = {
        assert!(!sequence.is_null());

        CStr::from_ptr(sequence)
    };

    sig.add_sequence(c_str.to_bytes(), force)
}
}

ffi_fn! {
unsafe fn signature_add_protein(ptr: *mut Signature, sequence: *const c_char) ->
    Result<()> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = {
        assert!(!sequence.is_null());

        CStr::from_ptr(sequence)
    };

    sig.add_protein(c_str.to_bytes())
}
}

ffi_fn! {
unsafe fn signature_set_name(ptr: *mut Signature, name: *const c_char) ->
    Result<()> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = {
        assert!(!name.is_null());

        CStr::from_ptr(name)
    };

    if let Ok(name) = c_str.to_str() {
        sig.name = Some(name.to_string())
    }
    Ok(())
}
}

ffi_fn! {
unsafe fn signature_set_filename(ptr: *mut Signature, name: *const c_char) ->
    Result<()> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = {
        assert!(!name.is_null());

        CStr::from_ptr(name)
    };

    if let Ok(name) = c_str.to_str() {
        sig.filename = Some(name.to_string())
    }
    Ok(())
}
}

ffi_fn! {
unsafe fn signature_push_mh(ptr: *mut Signature, other: *const KmerMinHash) ->
    Result<()> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let mh = {
       assert!(!other.is_null());
       &*other
    };

    sig.signatures.push(Sketch::MinHash(mh.clone()));
    Ok(())
}
}

ffi_fn! {
unsafe fn signature_set_mh(ptr: *mut Signature, other: *const KmerMinHash) ->
    Result<()> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let mh = {
       assert!(!other.is_null());
       &*other
    };

    sig.signatures = vec![Sketch::MinHash(mh.clone())];
    Ok(())
}
}

ffi_fn! {
unsafe fn signature_get_name(ptr: *mut Signature) -> Result<SourmashStr> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    if let Some(ref name) = sig.name {
        Ok(SourmashStr::from_string(name.to_string()))
    } else {
        Ok(SourmashStr::from_string("".to_string()))
    }
}
}

ffi_fn! {
unsafe fn signature_get_filename(ptr: *mut Signature) -> Result<SourmashStr> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    if let Some(ref name) = sig.filename {
        Ok(SourmashStr::from_string(name.to_string()))
    } else {
        Ok(SourmashStr::from_string("".to_string()))
    }
}
}

ffi_fn! {
unsafe fn signature_get_license(ptr: *mut Signature) -> Result<SourmashStr> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    Ok(SourmashStr::from_string(sig.license.to_string()))
}
}

ffi_fn! {
unsafe fn signature_first_mh(ptr: *mut Signature) -> Result<*mut KmerMinHash> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    if let Some(item) = sig.signatures.get(0) {
        if let Sketch::MinHash(mh) = item {
          Ok(Box::into_raw(Box::new(mh.clone())) as _)
        } else {
          unimplemented!()
        }
    } else {
        // TODO: need to select the correct one
        unimplemented!()
    }
}
}

ffi_fn! {
unsafe fn signature_eq(ptr: *mut Signature, other: *mut Signature) -> Result<bool> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let other_sig = {
        assert!(!other.is_null());
        &mut *other
    };

    Ok(sig == other_sig)
}
}

ffi_fn! {
unsafe fn signature_save_json(ptr: *mut Signature) -> Result<SourmashStr> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let st = serde_json::to_string(sig)?;
    Ok(SourmashStr::from_string(st))
}
}

ffi_fn! {
unsafe fn signature_get_mhs(ptr: *mut Signature, size: *mut usize) -> Result<*mut *mut KmerMinHash> {
    let sig = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let output = sig.signatures.clone();

    let ptr_sigs: Vec<*mut Signature> = output.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut Signature
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut KmerMinHash)
}
}

ffi_fn! {
unsafe fn signatures_save_buffer(ptr: *mut *mut Signature, size: usize, compression: u8, osize: *mut usize) -> Result<*const u8> {
    let sigs = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr, size)
    };

    let rsigs: Vec<&Signature> = sigs.iter().map(|x| x.as_ref().unwrap()).collect();

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
      serde_json::to_writer(&mut writer, &rsigs)?;
    }

    let b = buffer.into_boxed_slice();
    *osize = b.len();

    Ok(Box::into_raw(b) as *const u8)
}
}

ffi_fn! {
unsafe fn signatures_load_path(ptr: *const c_char,
                               _ignore_md5sum: bool,
                               ksize: usize,
                               select_moltype: *const c_char,
                               size: *mut usize) -> Result<*mut *mut Signature> {
    let buf = {
        assert!(!ptr.is_null());
        CStr::from_ptr(ptr)
    };

    let moltype: Option<HashFunctions> = if select_moltype.is_null() {
          None
        } else {
          let mol = CStr::from_ptr(select_moltype).to_str()?;
          Some(mol.try_into()?)
    };

    // TODO: implement ignore_md5sum

    let k = match ksize {
      0 => None,
      x => Some(x)
    };

    let (mut input, _) = niffler::from_path(buf.to_str()?).map_err(|_| SourmashError::IOError)?;
    let filtered_sigs = Signature::load_signatures(&mut input, k, moltype, None)?;

    let ptr_sigs: Vec<*mut Signature> = filtered_sigs.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut Signature
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut Signature)
}
}

ffi_fn! {
unsafe fn signatures_load_buffer(ptr: *const c_char,
                                 insize: usize,
                                 _ignore_md5sum: bool,
                                 ksize: usize,
                                 select_moltype: *const c_char,
                                 size: *mut usize) -> Result<*mut *mut Signature> {
    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr as *mut u8, insize)
    };

    let moltype: Option<HashFunctions> = if select_moltype.is_null() {
          None
        } else {
          let mol = CStr::from_ptr(select_moltype).to_str()?;
          Some(mol.try_into()?)
    };

    let k = match ksize {
      0 => None,
      x => Some(x)
    };

    // TODO: implement ignore_md5sum

    let mut reader = io::BufReader::new(buf);
    let filtered_sigs = Signature::load_signatures(&mut reader, k, moltype, None)?;

    let ptr_sigs: Vec<*mut Signature> = filtered_sigs.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut Signature
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut Signature)
}
}
