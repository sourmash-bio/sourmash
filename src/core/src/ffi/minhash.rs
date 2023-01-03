use std::ffi::CStr;
use std::os::raw::c_char;
use std::slice;

use crate::encodings::{aa_to_dayhoff, aa_to_hp, translate_codon, HashFunctions};
use crate::ffi::utils::{ForeignObject, SourmashStr};
use crate::signature::SeqToHashes;
use crate::signature::SigsTrait;
use crate::sketch::minhash::KmerMinHash;

pub struct SourmashKmerMinHash;

impl ForeignObject for SourmashKmerMinHash {
    type RustObject = KmerMinHash;
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_new(
    scaled: u64,
    k: u32,
    hash_function: HashFunctions,
    seed: u64,
    track_abundance: bool,
    n: u32,
) -> *mut SourmashKmerMinHash {
    let mh = KmerMinHash::new(scaled, k, hash_function, seed, track_abundance, n);

    SourmashKmerMinHash::from_rust(mh)
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_free(ptr: *mut SourmashKmerMinHash) {
    SourmashKmerMinHash::drop(ptr);
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_slice_free(ptr: *mut u64, insize: usize) {
    // FIXME
    if ptr.is_null() {
        return;
    }
    Vec::from_raw_parts(ptr as *mut u64, insize, insize);
}

ffi_fn! {
unsafe fn kmerminhash_add_sequence(ptr: *mut SourmashKmerMinHash, sequence: *const c_char, force: bool) ->
    Result<()> {

    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    // FIXME: take buffer and len instead of c_char
    let c_str = {
        assert!(!sequence.is_null());

        CStr::from_ptr(sequence)
    };

    mh.add_sequence(c_str.to_bytes(), force)
}
}

ffi_fn! {
unsafe fn kmerminhash_seq_to_hashes(ptr: *mut SourmashKmerMinHash, sequence: *const c_char, insize: usize, force: bool, bad_kmers_as_zeroes: bool, is_protein: bool, size: *mut usize) ->
Result<*const u64> {

    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    let buf = {
        assert!(!ptr.is_null());
        slice::from_raw_parts(sequence as *const u8, insize)
    };

    let mut output: Vec<u64> = Vec::with_capacity(insize);

    if force && bad_kmers_as_zeroes{
        for hash_value in SeqToHashes::new(buf, mh.ksize(), force, is_protein, mh.hash_function(), mh.seed()){
            match hash_value{
                Ok(x) => output.push(x),
                Err(err) => return Err(err),
            }
        }
    }else{
        for hash_value in SeqToHashes::new(buf, mh.ksize(), force, is_protein, mh.hash_function(), mh.seed()){
            match hash_value{
                Ok(0) => continue,
                Ok(x) => output.push(x),
                Err(err) => return Err(err),
            }
        }
    }


    *size = output.len();

    // FIXME: make a SourmashSlice_u64 type?
    Ok(Box::into_raw(output.into_boxed_slice()) as *const u64)
}
}

ffi_fn! {
unsafe fn kmerminhash_add_protein(ptr: *mut SourmashKmerMinHash, sequence: *const c_char) ->
    Result<()> {

    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    // FIXME: take buffer and len instead of c_char
    let c_str = {
        assert!(!sequence.is_null());

        CStr::from_ptr(sequence)
    };

    mh.add_protein(c_str.to_bytes())
}
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_clear(ptr: *mut SourmashKmerMinHash) {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    mh.clear();
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_add_hash(ptr: *mut SourmashKmerMinHash, h: u64) {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    mh.add_hash(h);
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_add_hash_with_abundance(
    ptr: *mut SourmashKmerMinHash,
    h: u64,
    abundance: u64,
) {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    mh.add_hash_with_abundance(h, abundance);
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_add_word(ptr: *mut SourmashKmerMinHash, word: *const c_char) {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    // FIXME: take buffer and len instead of c_char
    let c_str = {
        assert!(!word.is_null());

        CStr::from_ptr(word)
    };

    mh.add_word(c_str.to_bytes());
}

ffi_fn! {
unsafe fn sourmash_translate_codon(codon: *const c_char) -> Result<c_char> {
    // FIXME: take buffer and len instead of c_char
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
pub unsafe extern "C" fn kmerminhash_remove_hash(ptr: *mut SourmashKmerMinHash, h: u64) {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    mh.remove_hash(h);
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_remove_many(
    ptr: *mut SourmashKmerMinHash,
    hashes_ptr: *const u64,
    insize: usize,
) {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    // FIXME: make a SourmashSlice_u64 type?
    let hashes = {
        assert!(!hashes_ptr.is_null());
        slice::from_raw_parts(hashes_ptr as *mut u64, insize)
    };

    // FIXME: proper exception here
    mh.remove_many(hashes).expect("Hash removal error");
}

ffi_fn! {
unsafe fn kmerminhash_get_mins(ptr: *const SourmashKmerMinHash, size: *mut usize) -> Result<*const u64> {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let output = mh.mins();
    *size = output.len();

    // FIXME: make a SourmashSlice_u64 type?
    Ok(Box::into_raw(output.into_boxed_slice()) as *const u64)
}
}

ffi_fn! {
unsafe fn kmerminhash_md5sum(ptr: *const SourmashKmerMinHash) -> Result<SourmashStr> {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let output = mh.md5sum();

    Ok(output.into())
}
}

ffi_fn! {
unsafe fn kmerminhash_add_many(
    ptr: *mut SourmashKmerMinHash,
    hashes_ptr: *const u64,
    insize: usize,
  ) -> Result<()> {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    // FIXME: make a SourmashSlice_u64 type?
    let hashes = {
        assert!(!hashes_ptr.is_null());
        slice::from_raw_parts(hashes_ptr as *const u64, insize)
    };

    for hash in hashes {
      mh.add_hash(*hash);
    }

    Ok(())
}
}

ffi_fn! {
unsafe fn kmerminhash_get_abunds(ptr: *mut SourmashKmerMinHash, size: *mut usize) -> Result<*const u64> {
    let mh = SourmashKmerMinHash::as_rust(ptr);

    if let Some(abunds) = mh.abunds() {
        *size = abunds.len();
        Ok(Box::into_raw(abunds.into_boxed_slice()) as *const u64)
    } else {
        //throw error, can't get abund
        unimplemented!()
    }
}
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_get_mins_size(ptr: *const SourmashKmerMinHash) -> usize {
    let mh = SourmashKmerMinHash::as_rust(ptr);

    mh.size()
}

ffi_fn! {
unsafe fn kmerminhash_set_abundances(
    ptr: *mut SourmashKmerMinHash,
    hashes_ptr: *const u64,
    abunds_ptr: *const u64,
    insize: usize,
    clear: bool,
) -> Result<()> {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);

    // FIXME: make a SourmashSlice_u64 type?
    let hashes = {
        assert!(!hashes_ptr.is_null());
        slice::from_raw_parts(hashes_ptr as *const u64, insize)
    };

    // FIXME: make a SourmashSlice_u64 type?
    let abunds = {
        assert!(!abunds_ptr.is_null());
        slice::from_raw_parts(abunds_ptr as *const u64, insize)
    };

    let mut pairs: Vec<_> = hashes.iter().cloned().zip(abunds.iter().cloned()).collect();
    pairs.sort_unstable();

    // Reset the minhash
    if clear {
        mh.clear();
    }

    mh.add_many_with_abund(&pairs)?;

    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_is_protein(ptr: *const SourmashKmerMinHash) -> bool {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.is_protein()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_dayhoff(ptr: *const SourmashKmerMinHash) -> bool {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.dayhoff()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_hp(ptr: *const SourmashKmerMinHash) -> bool {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.hp()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_seed(ptr: *const SourmashKmerMinHash) -> u64 {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.seed()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_track_abundance(ptr: *const SourmashKmerMinHash) -> bool {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.track_abundance()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_disable_abundance(ptr: *mut SourmashKmerMinHash) {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);
    mh.disable_abundance();
}

ffi_fn! {
unsafe fn kmerminhash_enable_abundance(ptr: *mut SourmashKmerMinHash) -> Result<()> {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);
    mh.enable_abundance()?;
    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_num(ptr: *const SourmashKmerMinHash) -> u32 {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.num()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_ksize(ptr: *const SourmashKmerMinHash) -> u32 {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.ksize() as u32
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_max_hash(ptr: *const SourmashKmerMinHash) -> u64 {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.max_hash()
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_hash_function(
    ptr: *const SourmashKmerMinHash,
) -> HashFunctions {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    mh.hash_function()
}

ffi_fn! {
unsafe fn kmerminhash_hash_function_set(ptr: *mut SourmashKmerMinHash, hash_function: HashFunctions) -> Result<()> {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);
    mh.set_hash_function(hash_function)
}
}

ffi_fn! {
unsafe fn kmerminhash_merge(ptr: *mut SourmashKmerMinHash, other: *const SourmashKmerMinHash) -> Result<()> {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);
    mh.merge(other_mh)?;
    Ok(())
}
}

#[no_mangle]
pub unsafe extern "C" fn kmerminhash_is_compatible(
    ptr: *const SourmashKmerMinHash,
    other: *const SourmashKmerMinHash,
) -> bool {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);
    mh.check_compatible(other_mh).is_ok()
}

ffi_fn! {
unsafe fn kmerminhash_add_from(ptr: *mut SourmashKmerMinHash, other: *const SourmashKmerMinHash)
    -> Result<()> {
    let mh = SourmashKmerMinHash::as_rust_mut(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);
    mh.add_from(other_mh)
}
}

ffi_fn! {
    unsafe fn kmerminhash_remove_from(ptr: *mut SourmashKmerMinHash, other: *const SourmashKmerMinHash)
    -> Result<()> {
        let mh = SourmashKmerMinHash::as_rust_mut(ptr);
        let other_mh = SourmashKmerMinHash::as_rust(other);
        mh.remove_from(other_mh)
    }
}

ffi_fn! {
unsafe fn kmerminhash_count_common(ptr: *const SourmashKmerMinHash, other: *const SourmashKmerMinHash, downsample: bool)
    -> Result<u64> {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);
    mh.count_common(other_mh, downsample)
}
}

ffi_fn! {
unsafe fn kmerminhash_intersection(ptr: *const SourmashKmerMinHash, other: *const SourmashKmerMinHash)
    -> Result<*mut SourmashKmerMinHash> {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);

    let isect = mh.intersection(other_mh)?;
    let mut new_mh = mh.clone();
    new_mh.clear();
    new_mh.add_many(&isect.0)?;

    Ok(SourmashKmerMinHash::from_rust(new_mh))
}
}

ffi_fn! {
unsafe fn kmerminhash_intersection_union_size(ptr: *const SourmashKmerMinHash, other: *const SourmashKmerMinHash, union_size: *mut u64)
    -> Result<u64> {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);

    if let Ok((common, union_s)) = mh.intersection_size(other_mh) {
        *union_size = union_s;
        return Ok(common);
    }

    *union_size = 0;
    Ok(0)
}
}

ffi_fn! {
unsafe fn kmerminhash_jaccard(ptr: *const SourmashKmerMinHash, other: *const SourmashKmerMinHash)
    -> Result<f64> {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);
    mh.jaccard(other_mh)
}
}

ffi_fn! {
unsafe fn kmerminhash_similarity(ptr: *const SourmashKmerMinHash, other: *const SourmashKmerMinHash, ignore_abundance: bool, downsample: bool)
    -> Result<f64> {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);
    mh.similarity(other_mh, ignore_abundance, downsample)
}
}
ffi_fn! {
unsafe fn kmerminhash_angular_similarity(ptr: *const SourmashKmerMinHash, other: *const SourmashKmerMinHash)
                                         -> Result<f64> {
    let mh = SourmashKmerMinHash::as_rust(ptr);
    let other_mh = SourmashKmerMinHash::as_rust(other);
    mh.angular_similarity(other_mh)
}
}
