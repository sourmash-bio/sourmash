use std::ffi::CStr;
use std::io;
use std::mem;
use std::os::raw::c_char;
use std::ptr;
use std::slice;

use ocf::get_input;
use serde_json;

use crate::utils::SourmashStr;
use crate::{KmerMinHash, Signature, _hash_murmur};

#[no_mangle]
pub extern "C" fn hash_murmur(kmer: *const c_char, seed: u64) -> u64 {
    let c_str = unsafe {
        assert!(!kmer.is_null());

        CStr::from_ptr(kmer)
    };

    _hash_murmur(c_str.to_bytes(), seed)
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_new(
    n: u32,
    k: u32,
    prot: bool,
    seed: u64,
    mx: u64,
    track_abundance: bool,
) -> *mut KmerMinHash {
    mem::transmute(Box::new(KmerMinHash::new(
        n,
        k,
        prot,
        seed,
        mx,
        track_abundance,
    )))
}

#[no_mangle]
pub extern "C" fn kmerminhash_free(ptr: *mut KmerMinHash) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        Box::from_raw(ptr);
    }
}

ffi_fn! {
unsafe fn kmerminhash_add_sequence(ptr: *mut KmerMinHash, sequence: *const c_char, force: bool) ->
    Result<()> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = {
        assert!(!sequence.is_null());

        CStr::from_ptr(sequence)
    };

    mh.add_sequence(c_str.to_bytes(), force)
}
}

#[no_mangle]
pub extern "C" fn kmerminhash_add_hash(ptr: *mut KmerMinHash, h: u64) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    mh.add_hash(h);
}

#[no_mangle]
pub extern "C" fn kmerminhash_add_word(ptr: *mut KmerMinHash, word: *const c_char) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = unsafe {
        assert!(!word.is_null());

        CStr::from_ptr(word)
    };

    mh.add_word(c_str.to_bytes());
}

ffi_fn! {
unsafe fn kmerminhash_get_mins(ptr: *mut KmerMinHash) -> Result<*const u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let output = mh.mins.clone();

    Ok(Box::into_raw(output.into_boxed_slice()) as *const u64)
}
}

ffi_fn! {
unsafe fn kmerminhash_get_abunds(ptr: *mut KmerMinHash) -> Result<*const u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    if let Some(ref abunds) = mh.abunds {
        let output = abunds.clone();
        Ok(Box::into_raw(output.into_boxed_slice()) as *const u64)
    } else {
        Ok(ptr::null())
    }
}
}

ffi_fn! {
unsafe fn kmerminhash_get_min_idx(ptr: *mut KmerMinHash, idx: u64) -> Result<u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    Ok(mh.mins[idx as usize])
}
}

#[no_mangle]
pub extern "C" fn kmerminhash_get_mins_size(ptr: *mut KmerMinHash) -> usize {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.mins.len()
}

#[no_mangle]
pub extern "C" fn kmerminhash_mins_push(ptr: *mut KmerMinHash, val: u64) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.mins.push(val)
}

ffi_fn! {
unsafe fn kmerminhash_get_abund_idx(ptr: *mut KmerMinHash, idx: u64) -> Result<u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    if let Some(ref mut abunds) = mh.abunds {
      Ok(abunds[idx as usize])
    } else {
      Ok(0)
    }
}
}

#[no_mangle]
pub extern "C" fn kmerminhash_get_abunds_size(ptr: *mut KmerMinHash) -> usize {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    if let Some(ref mut abunds) = mh.abunds {
        abunds.len()
    } else {
        0
    }
}

#[no_mangle]
pub extern "C" fn kmerminhash_abunds_push(ptr: *mut KmerMinHash, val: u64) {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    if let Some(ref mut abunds) = mh.abunds {
        abunds.push(val)
    }
}

#[no_mangle]
pub extern "C" fn kmerminhash_is_protein(ptr: *mut KmerMinHash) -> bool {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.is_protein
}

#[no_mangle]
pub extern "C" fn kmerminhash_seed(ptr: *mut KmerMinHash) -> u64 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.seed
}

#[no_mangle]
pub extern "C" fn kmerminhash_track_abundance(ptr: *mut KmerMinHash) -> bool {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.abunds.is_some()
}

#[no_mangle]
pub extern "C" fn kmerminhash_num(ptr: *mut KmerMinHash) -> u32 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.num
}

#[no_mangle]
pub extern "C" fn kmerminhash_ksize(ptr: *mut KmerMinHash) -> u32 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.ksize
}

#[no_mangle]
pub extern "C" fn kmerminhash_max_hash(ptr: *mut KmerMinHash) -> u64 {
    let mh = unsafe {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.max_hash
}

ffi_fn! {
unsafe fn kmerminhash_merge(ptr: *mut KmerMinHash, other: *const KmerMinHash) -> Result<()> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.merge(other_mh)?;
    Ok(())
}
}

ffi_fn! {
unsafe fn kmerminhash_add_from(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<()> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.add_from(other_mh)
}
}

ffi_fn! {
unsafe fn kmerminhash_count_common(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.count_common(other_mh)
}
}

ffi_fn! {
unsafe fn kmerminhash_intersection(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<u64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    if let Ok((_, size)) = mh.intersection(other_mh) {
        return Ok(size);
    }
    Ok(0)
}
}

ffi_fn! {
unsafe fn kmerminhash_compare(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<f64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.compare(other_mh)
}
}

// Signature methods

#[no_mangle]
pub unsafe extern "C" fn signature_new() -> *mut Signature {
    mem::transmute(Box::new(Signature::default()))
}

#[no_mangle]
pub extern "C" fn signature_free(ptr: *mut Signature) {
    if ptr.is_null() {
        return;
    }
    unsafe {
        Box::from_raw(ptr);
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

    sig.signatures.push(mh.clone());
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

    sig.signatures = vec![mh.clone()];
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

    if let Some(mh) = sig.signatures.get(0) {
        Ok(mem::transmute(Box::new(mh.clone())))
    } else {
        // TODO: this is totally wrong
        Ok(mem::transmute(Box::new(KmerMinHash::default())))
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
unsafe fn signatures_save_buffer(ptr: *mut *mut Signature, size: usize) -> Result<SourmashStr> {
    let sigs = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr, size)
    };

    let rsigs: Vec<&Signature> = sigs.iter().map(|x| x.as_ref().unwrap()).collect();
    let st = serde_json::to_string(&rsigs)?;
    Ok(SourmashStr::from_string(st))
}
}

ffi_fn! {
unsafe fn signatures_load_path(ptr: *const c_char,
                               ignore_md5sum: bool,
                               ksize: usize,
                               select_moltype: *const c_char,
                               size: *mut usize) -> Result<*mut *mut Signature> {
    let buf = {
        assert!(!ptr.is_null());
        CStr::from_ptr(ptr)
    };

    let moltype = {
        if select_moltype.is_null() {
          None
        } else {
          Some(CStr::from_ptr(select_moltype).to_str()?)
        }
    };

    // TODO: implement ignore_md5sum

    let (mut input, _) = get_input(buf.to_str()?)?;
    let filtered_sigs = Signature::load_signatures(&mut input, ksize, moltype, None)?;

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
                                 ignore_md5sum: bool,
                                 ksize: usize,
                                 select_moltype: *const c_char,
                                 size: *mut usize) -> Result<*mut *mut Signature> {
    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(ptr as *mut u8, insize)
    };

    let moltype = {
        if select_moltype.is_null() {
          None
        } else {
          Some(CStr::from_ptr(select_moltype).to_str()?)
        }
    };

    // TODO: implement ignore_md5sum

    let mut reader = io::BufReader::new(buf);
    let filtered_sigs = Signature::load_signatures(&mut reader, ksize, moltype, None)?;

    let ptr_sigs: Vec<*mut Signature> = filtered_sigs.into_iter().map(|x| {
      Box::into_raw(Box::new(x)) as *mut Signature
    }).collect();

    let b = ptr_sigs.into_boxed_slice();
    *size = b.len();

    Ok(Box::into_raw(b) as *mut *mut Signature)
}
}
