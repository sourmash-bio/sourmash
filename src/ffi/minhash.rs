use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::errors::SourmashError;
use crate::signature::SigsTrait;
use crate::sketch::minhash::{
    aa_to_dayhoff, aa_to_hp, translate_codon, HashFunctions, KmerMinHash,
};

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_new(
    n: u32,
    k: u32,
    prot: bool,
    dayhoff: bool,
    hp: bool,
    seed: u64,
    mx: u64,
    track_abundance: bool,
) -> *mut KmerMinHash {
    // TODO: at most one of (prot, dayhoff, hp) should be true

    let hash_function = if dayhoff {
        HashFunctions::murmur64_dayhoff
    } else if hp {
        HashFunctions::murmur64_hp
    } else if prot {
        HashFunctions::murmur64_protein
    } else {
        HashFunctions::murmur64_DNA
    };

    Box::into_raw(Box::new(KmerMinHash::new(
        n,
        k,
        hash_function,
        seed,
        mx,
        track_abundance,
    ))) as _
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_free(ptr: *mut KmerMinHash) {
    if ptr.is_null() {
        return;
    }
    Box::from_raw(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_slice_free(ptr: *mut u64, insize: usize) {
    if ptr.is_null() {
        return;
    }
    Vec::from_raw_parts(ptr as *mut u64, insize, insize);
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
pub unsafe extern "C" fn kmerminhash_add_hash(ptr: *mut KmerMinHash, h: u64) {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    mh.add_hash(h);
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_add_word(ptr: *mut KmerMinHash, word: *const c_char) {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let c_str = {
        assert!(!word.is_null());

        CStr::from_ptr(word)
    };

    mh.add_word(c_str.to_bytes());
}

ffi_fn! {
unsafe fn sourmash_translate_codon(codon: *const c_char) -> Result<c_char> {
    let c_str = {
        assert!(!codon.is_null());

        CStr::from_ptr(codon)
    };

    Ok(translate_codon(c_str.to_bytes())? as c_char)
}
}

#[no_mangle]
pub unsafe extern "C" fn sourmash_aa_to_dayhoff(aa: c_char) -> c_char {
    aa_to_dayhoff(aa as u8) as c_char
}

#[no_mangle]
pub unsafe extern "C" fn sourmash_aa_to_hp(aa: c_char) -> c_char {
    aa_to_hp(aa as u8) as c_char
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_remove_hash(ptr: *mut KmerMinHash, h: u64) {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    mh.remove_hash(h);
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_remove_many(
    ptr: *mut KmerMinHash,
    hashes_ptr: *const u64,
    insize: usize,
) {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    let hashes = {
        assert!(!hashes_ptr.is_null());
        slice::from_raw_parts(hashes_ptr as *mut u64, insize)
    };

    mh.remove_many(hashes).expect("Hash removal error");
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
        //throw error, can't get abund
        unimplemented!()
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
pub unsafe extern "C" fn kmerminhash_get_mins_size(ptr: *mut KmerMinHash) -> usize {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.mins.len()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_mins_push(ptr: *mut KmerMinHash, val: u64) {
    let mh = {
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
pub unsafe extern "C" fn kmerminhash_get_abunds_size(ptr: *mut KmerMinHash) -> usize {
    let mh = {
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
pub unsafe extern "C" fn kmerminhash_abunds_push(ptr: *mut KmerMinHash, val: u64) {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    if let Some(ref mut abunds) = mh.abunds {
        abunds.push(val)
    }
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_is_protein(ptr: *mut KmerMinHash) -> bool {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.is_protein()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_dayhoff(ptr: *mut KmerMinHash) -> bool {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.dayhoff()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_hp(ptr: *mut KmerMinHash) -> bool {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.hp()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_seed(ptr: *mut KmerMinHash) -> u64 {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.seed()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_track_abundance(ptr: *mut KmerMinHash) -> bool {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.abunds.is_some()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_disable_abundance(ptr: *mut KmerMinHash) {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.abunds = None;
}

ffi_fn! {
unsafe fn kmerminhash_enable_abundance(ptr: *mut KmerMinHash) -> Result<()> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    if !mh.mins.is_empty() {
      return Err(SourmashError::NonEmptyMinHash { message: "track_abundance=True".into()}.into());
    }

    mh.abunds = Some(vec![]);
    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_num(ptr: *mut KmerMinHash) -> u32 {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.num()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_ksize(ptr: *mut KmerMinHash) -> u32 {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.ksize() as u32
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_max_hash(ptr: *mut KmerMinHash) -> u64 {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.max_hash()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_hash_function(ptr: *mut KmerMinHash) -> HashFunctions {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    mh.hash_function()
}

ffi_fn! {
unsafe fn kmerminhash_hash_function_set(ptr: *mut KmerMinHash, hash_function: HashFunctions) -> Result<()> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };

    if !mh.mins.is_empty() {
      return Err(SourmashError::NonEmptyMinHash { message: "hash_function".into()}.into());
    }

    mh.hash_function = hash_function;
    Ok(())
}
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

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_is_compatible(
    ptr: *const KmerMinHash,
    other: *const KmerMinHash,
) -> bool {
    let mh = {
        assert!(!ptr.is_null());
        &*ptr
    };
    let other_mh = {
        assert!(!other.is_null());
        &*other
    };

    mh.check_compatible(other_mh).is_ok()
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

    if let Ok((_, size)) = mh.intersection_size(other_mh) {
        return Ok(size);
    }
    Ok(0)
}
}

ffi_fn! {
unsafe fn kmerminhash_containment_ignore_maxhash(ptr: *mut KmerMinHash, other: *const KmerMinHash)
    -> Result<f64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.containment_ignore_maxhash(&other_mh)
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

ffi_fn! {
unsafe fn kmerminhash_similarity(ptr: *mut KmerMinHash, other: *const KmerMinHash, ignore_abundance: bool)
    -> Result<f64> {
    let mh = {
        assert!(!ptr.is_null());
        &mut *ptr
    };
    let other_mh = {
       assert!(!other.is_null());
       &*other
    };

    mh.similarity(other_mh, ignore_abundance)
}
}
