//! # Foreign Function Interface for calling sourmash from a C API
//!
//! Primary client for now is the Python version, using CFFI and maturin.
#![allow(clippy::missing_safety_doc)]

#[macro_use]
pub mod utils;

pub mod cmd;
pub mod hyperloglog;
pub mod index;
pub mod manifest;
pub mod minhash;
pub mod nodegraph;
pub mod signature;
pub mod storage;

use std::ffi::CStr;
use std::os::raw::c_char;

use crate::_hash_murmur;

#[no_mangle]
pub unsafe extern "C" fn hash_murmur(kmer: *const c_char, seed: u64) -> u64 {
    let c_str = {
        assert!(!kmer.is_null());

        CStr::from_ptr(kmer)
    };

    _hash_murmur(c_str.to_bytes(), seed)
}

#[repr(u32)]
pub enum HashFunctions {
    Murmur64Dna = 1,
    Murmur64Protein = 2,
    Murmur64Dayhoff = 3,
    Murmur64Hp = 4,
    Murmur64Skipm1n3 = 5,
    Murmur64Skipm2n3 = 6,
}

impl From<HashFunctions> for crate::encodings::HashFunctions {
    fn from(v: HashFunctions) -> crate::encodings::HashFunctions {
        use crate::encodings::HashFunctions::{
            Murmur64Dayhoff, Murmur64Dna, Murmur64Hp, Murmur64Protein, Murmur64Skipm1n3,
            Murmur64Skipm2n3,
        };
        match v {
            HashFunctions::Murmur64Dna => Murmur64Dna,
            HashFunctions::Murmur64Protein => Murmur64Protein,
            HashFunctions::Murmur64Dayhoff => Murmur64Dayhoff,
            HashFunctions::Murmur64Hp => Murmur64Hp,
            HashFunctions::Murmur64Skipm1n3 => Murmur64Skipm1n3,
            HashFunctions::Murmur64Skipm2n3 => Murmur64Skipm2n3,
        }
    }
}

impl From<crate::encodings::HashFunctions> for HashFunctions {
    fn from(v: crate::encodings::HashFunctions) -> HashFunctions {
        use crate::encodings::HashFunctions::{
            Murmur64Dayhoff, Murmur64Dna, Murmur64Hp, Murmur64Protein, Murmur64Skipm1n3,
            Murmur64Skipm2n3,
        };
        match v {
            Murmur64Dna => HashFunctions::Murmur64Dna,
            Murmur64Protein => HashFunctions::Murmur64Protein,
            Murmur64Dayhoff => HashFunctions::Murmur64Dayhoff,
            Murmur64Hp => HashFunctions::Murmur64Hp,
            Murmur64Skipm1n3 => HashFunctions::Murmur64Skipm1n3,
            Murmur64Skipm2n3 => HashFunctions::Murmur64Skipm2n3,
            _ => todo!("Not supported, probably custom"),
        }
    }
}
