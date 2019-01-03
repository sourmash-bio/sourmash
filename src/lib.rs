pub mod errors;

#[macro_use]
pub mod utils;

pub mod index;

#[cfg(feature = "from-finch")]
pub mod from;

use serde::de::{Deserialize, Deserializer};
use serde::ser::{Serialize, SerializeStruct, Serializer};
use serde_derive::{Deserialize, Serialize};

use std::cmp::Ordering;
use std::collections::HashMap;
use std::fs::File;
use std::io;
use std::iter::{Iterator, Peekable};
use std::path::Path;
use std::str;

use cfg_if::cfg_if;
use failure::Error;
use lazy_static::lazy_static;
use murmurhash3::murmurhash3_x64_128;

use crate::errors::SourmashError;

cfg_if! {
    if #[cfg(target_arch = "wasm32")] {
        use wasm_bindgen::prelude::*;

        pub mod wasm;
    } else {
        pub mod ffi;
    }
}

pub fn _hash_murmur(kmer: &[u8], seed: u64) -> u64 {
    murmurhash3_x64_128(kmer, seed).0
}

#[cfg_attr(target_arch = "wasm32", wasm_bindgen)]
#[derive(Debug, Clone, PartialEq)]
pub struct KmerMinHash {
    pub num: u32,
    pub ksize: u32,
    pub is_protein: bool,
    pub seed: u64,
    pub max_hash: u64,
    mins: Vec<u64>,
    abunds: Option<Vec<u64>>,
}

impl Default for KmerMinHash {
    fn default() -> KmerMinHash {
        KmerMinHash {
            num: 1000,
            ksize: 21,
            is_protein: false,
            seed: 42,
            max_hash: 0,
            mins: Vec::with_capacity(1000),
            abunds: None,
        }
    }
}

impl Serialize for KmerMinHash {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = match &self.abunds {
            Some(_) => 8,
            _ => 7,
        };

        let mut md5_ctx = md5::Context::new();
        md5_ctx.consume(&self.ksize.to_string());
        self.mins
            .iter()
            .map(|x| md5_ctx.consume(x.to_string()))
            .count();

        let mut partial = serializer.serialize_struct("KmerMinHash", n_fields)?;
        partial.serialize_field("num", &self.num)?;
        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("seed", &self.seed)?;
        partial.serialize_field("max_hash", &self.max_hash)?;
        partial.serialize_field("mins", &self.mins)?;

        partial.serialize_field("md5sum", &format!("{:x}", md5_ctx.compute()))?;

        if let Some(abunds) = &self.abunds {
            partial.serialize_field("abundances", abunds)?;
        }

        partial.serialize_field(
            "molecule",
            match &self.is_protein {
                true => "protein",
                false => "DNA",
            },
        )?;

        partial.end()
    }
}

impl<'de> Deserialize<'de> for KmerMinHash {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct TempSig {
            num: u32,
            ksize: u32,
            seed: u64,
            max_hash: u64,
            md5sum: String,
            mins: Vec<u64>,
            abundances: Option<Vec<u64>>,
            molecule: String,
        }

        let tmpsig = TempSig::deserialize(deserializer)?;

        let num = if tmpsig.max_hash != 0 { 0 } else { tmpsig.num };

        Ok(KmerMinHash {
            num,
            ksize: tmpsig.ksize,
            seed: tmpsig.seed,
            max_hash: tmpsig.max_hash,
            mins: tmpsig.mins,
            abunds: tmpsig.abundances,
            is_protein: match tmpsig.molecule.as_ref() {
                "protein" => true,
                "DNA" => false,
                _ => false, // TODO: throw error
            },
        })
    }
}

impl KmerMinHash {
    pub fn new(
        num: u32,
        ksize: u32,
        is_protein: bool,
        seed: u64,
        max_hash: u64,
        track_abundance: bool,
    ) -> KmerMinHash {
        let mins: Vec<u64>;
        let abunds: Option<Vec<u64>>;

        if num > 0 {
            mins = Vec::with_capacity(num as usize);
        } else {
            mins = Vec::with_capacity(1000);
        }

        if track_abundance {
            abunds = Some(Vec::with_capacity(mins.capacity()));
        } else {
            abunds = None
        }

        KmerMinHash {
            num,
            ksize,
            is_protein,
            seed,
            max_hash,
            mins,
            abunds,
        }
    }

    pub fn check_compatible(&self, other: &KmerMinHash) -> Result<bool, Error> {
        if self.ksize != other.ksize {
            return Err(SourmashError::MismatchKSizes.into());
        }
        if self.is_protein != other.is_protein {
            return Err(SourmashError::MismatchDNAProt.into());
        }
        if self.max_hash != other.max_hash {
            return Err(SourmashError::MismatchMaxHash.into());
        }
        if self.seed != other.seed {
            return Err(SourmashError::MismatchSeed.into());
        }
        Ok(true)
    }

    pub fn add_hash(&mut self, hash: u64) {
        let current_max = match self.mins.last() {
            Some(&x) => x,
            None => u64::max_value(),
        };

        if hash <= self.max_hash || self.max_hash == 0 {
            // empty? add it, if within range / no range specified.
            if self.mins.is_empty() {
                self.mins.push(hash);
                if let Some(ref mut abunds) = self.abunds {
                    abunds.push(1);
                }
                return;
            } else if hash <= self.max_hash
                || current_max > hash
                || (self.mins.len() as u32) < self.num
            {
                // "good" hash - within range, smaller than current entry, or
                // still have space available
                let pos = match self.mins.binary_search(&hash) {
                    Ok(p) => p,
                    Err(p) => p,
                };

                if pos == self.mins.len() {
                    // at end - must still be growing, we know the list won't
                    // get too long
                    self.mins.push(hash);
                    if let Some(ref mut abunds) = self.abunds {
                        abunds.push(1);
                    }
                } else if self.mins[pos] != hash {
                    // didn't find hash in mins, so inserting somewhere
                    // in the middle; shrink list if needed.
                    self.mins.insert(pos, hash);
                    if let Some(ref mut abunds) = self.abunds {
                        abunds.insert(pos, 1);
                    }

                    // is it too big now?
                    if self.num != 0 && self.mins.len() > (self.num as usize) {
                        self.mins.pop();
                        if let Some(ref mut abunds) = self.abunds {
                            abunds.pop();
                        }
                    }
                } else if let Some(ref mut abunds) = self.abunds {
                    // pos == hash: hash value already in mins, inc count
                    abunds[pos] += 1;
                }
            }
        }
    }

    pub fn add_word(&mut self, word: &[u8]) {
        let hash = _hash_murmur(word, self.seed);
        self.add_hash(hash);
    }

    pub fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        let sequence: Vec<u8> = seq
            .iter()
            .map(|&x| (x as char).to_ascii_uppercase() as u8)
            .collect();
        if sequence.len() >= (self.ksize as usize) {
            if !self.is_protein {
                // dna
                for kmer in sequence.windows(self.ksize as usize) {
                    if _checkdna(kmer) {
                        let rc = revcomp(kmer);
                        if kmer < &rc {
                            self.add_word(kmer);
                        } else {
                            self.add_word(&rc);
                        }
                    } else if !force {
                        return Err(SourmashError::InvalidDNA {
                            message: String::from_utf8(kmer.to_vec()).unwrap(),
                        }
                        .into());
                    }
                }
            } else {
                // protein
                let rc = revcomp(&sequence);
                let aa_ksize = self.ksize / 3;

                for i in 0..3 {
                    let substr: Vec<u8> = sequence
                        .iter()
                        .cloned()
                        .skip(i)
                        .take(sequence.len() - i)
                        .collect();
                    let aa = to_aa(&substr);

                    aa.windows(aa_ksize as usize)
                        .map(|n| self.add_word(n))
                        .count();

                    let rc_substr: Vec<u8> =
                        rc.iter().cloned().skip(i).take(rc.len() - i).collect();
                    let aa_rc = to_aa(&rc_substr);

                    aa_rc
                        .windows(aa_ksize as usize)
                        .map(|n| self.add_word(n))
                        .count();
                }
            }
        }
        Ok(())
    }

    pub fn merge(&mut self, other: &KmerMinHash) -> Result<(), Error> {
        self.check_compatible(other)?;
        let max_size = self.mins.len() + other.mins.len();
        let mut merged: Vec<u64> = Vec::with_capacity(max_size);
        let mut merged_abunds: Vec<u64> = Vec::with_capacity(max_size);

        {
            let mut self_iter = self.mins.iter();
            let mut other_iter = other.mins.iter();

            let mut self_abunds_iter: Option<std::slice::Iter<'_, u64>>;
            if let Some(ref mut abunds) = self.abunds {
                self_abunds_iter = Some(abunds.iter());
            } else {
                self_abunds_iter = None;
            }

            let mut other_abunds_iter: Option<std::slice::Iter<'_, u64>>;
            if let Some(ref abunds) = other.abunds {
                other_abunds_iter = Some(abunds.iter());
            } else {
                other_abunds_iter = None;
            }

            let mut self_value = self_iter.next();
            let mut other_value = other_iter.next();
            while self_value.is_some() {
                let value = self_value.unwrap();
                match other_value {
                    None => {
                        merged.push(*value);
                        merged.extend(self_iter);
                        if let Some(sai) = self_abunds_iter {
                            merged_abunds.extend(sai);
                        }
                        break;
                    }
                    Some(x) if x < value => {
                        merged.push(*x);
                        other_value = other_iter.next();

                        if let Some(ref mut oai) = other_abunds_iter {
                            if let Some(v) = oai.next() {
                                merged_abunds.push(*v)
                            }
                        }
                    }
                    Some(x) if x == value => {
                        merged.push(*x);
                        other_value = other_iter.next();
                        self_value = self_iter.next();

                        if let Some(ref mut oai) = other_abunds_iter {
                            if let Some(v) = oai.next() {
                                if let Some(ref mut sai) = self_abunds_iter {
                                    if let Some(s) = sai.next() {
                                        merged_abunds.push(*v + *s)
                                    }
                                }
                            }
                        }
                    }
                    Some(x) if x > value => {
                        merged.push(*value);
                        self_value = self_iter.next();

                        if let Some(ref mut sai) = self_abunds_iter {
                            if let Some(v) = sai.next() {
                                merged_abunds.push(*v)
                            }
                        }
                    }
                    Some(_) => {}
                }
            }
            if let Some(value) = other_value {
                merged.push(*value);
            }
            merged.extend(other_iter);
            if let Some(oai) = other_abunds_iter {
                merged_abunds.extend(oai);
            }
        }

        if merged.len() < (self.num as usize) || (self.num as usize) == 0 {
            self.mins = merged;
            self.abunds = Some(merged_abunds);
        } else {
            self.mins = merged
                .iter()
                .map(|&x| x as u64)
                .take(self.num as usize)
                .collect();
            self.abunds = Some(merged_abunds) // TODO: reduce this one too
        }
        Ok(())
    }

    pub fn add_from(&mut self, other: &KmerMinHash) -> Result<(), Error> {
        for min in &other.mins {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many(&mut self, hashes: &[u64]) -> Result<(), Error> {
        for min in hashes {
            self.add_hash(*min);
        }
        Ok(())
    }

    pub fn add_many_with_abund(&mut self, hashes: &[(u64, u64)]) -> Result<(), Error> {
        for item in hashes {
            for _i in 0..item.1 {
                self.add_hash(item.0);
            }
        }
        Ok(())
    }

    pub fn count_common(&self, other: &KmerMinHash) -> Result<u64, Error> {
        self.check_compatible(other)?;
        let iter = Intersection {
            left: self.mins.iter().peekable(),
            right: other.mins.iter().peekable(),
        };

        Ok(iter.count() as u64)
    }

    pub fn intersection(&self, other: &KmerMinHash) -> Result<(Vec<u64>, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.is_protein,
            self.seed,
            self.max_hash,
            self.abunds.is_some(),
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection {
            left: self.mins.iter().peekable(),
            right: other.mins.iter().peekable(),
        };

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let it2 = Intersection {
            left: i1.iter().peekable(),
            right: combined_mh.mins.iter().peekable(),
        };

        let common: Vec<u64> = it2.cloned().collect();
        Ok((common, combined_mh.mins.len() as u64))
    }

    pub fn intersection_size(&self, other: &KmerMinHash) -> Result<(u64, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.is_protein,
            self.seed,
            self.max_hash,
            self.abunds.is_some(),
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection {
            left: self.mins.iter().peekable(),
            right: other.mins.iter().peekable(),
        };

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.into_iter().cloned().collect();
        let it2 = Intersection {
            left: i1.iter().peekable(),
            right: combined_mh.mins.iter().peekable(),
        };

        Ok((it2.count() as u64, combined_mh.mins.len() as u64))
    }

    pub fn compare(&self, other: &KmerMinHash) -> Result<f64, Error> {
        self.check_compatible(other)?;
        if let Ok((common, size)) = self.intersection_size(other) {
            return Ok(common as f64 / u64::max(1, size) as f64);
        } else {
            return Ok(0.0);
        }
    }

    pub fn size(&self) -> usize {
        self.mins.len()
    }

    pub fn to_vec(&self) -> Vec<u64> {
        self.mins.clone()
    }
}

struct Intersection<T, I: Iterator<Item = T>> {
    left: Peekable<I>,
    right: Peekable<I>,
}

impl<T: Ord, I: Iterator<Item = T>> Iterator for Intersection<T, I> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            let res = match (self.left.peek(), self.right.peek()) {
                (Some(ref left_key), Some(ref right_key)) => left_key.cmp(right_key),
                _ => return None,
            };

            match res {
                Ordering::Less => {
                    self.left.next();
                }
                Ordering::Greater => {
                    self.right.next();
                }
                Ordering::Equal => {
                    self.right.next();
                    return self.left.next();
                }
            }
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Signature {
    #[serde(default = "default_class")]
    pub class: String,

    #[serde(default)]
    pub email: String,
    pub hash_function: String,

    pub filename: Option<String>,
    pub name: Option<String>,

    #[serde(default = "default_license")]
    pub license: String,

    pub signatures: Vec<KmerMinHash>,

    #[serde(default = "default_version")]
    pub version: f64,
}

fn default_license() -> String {
    "CC0".to_string()
}

fn default_class() -> String {
    "sourmash_signature".to_string()
}

fn default_version() -> f64 {
    0.4
}

impl Signature {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Vec<Signature>, Error> {
        let mut reader = io::BufReader::new(File::open(path)?);
        Ok(Signature::from_reader(&mut reader)?)
    }

    pub fn from_reader<R>(rdr: &mut R) -> Result<Vec<Signature>, Error>
    where
        R: io::Read,
    {
        let sigs: Vec<Signature> = serde_json::from_reader(rdr)?;
        Ok(sigs)
    }

    pub fn load_signatures<R>(
        buf: &mut R,
        ksize: usize,
        moltype: Option<&str>,
        scaled: Option<u64>,
    ) -> Result<Vec<Signature>, Error>
    where
        R: io::Read,
    {
        let orig_sigs = Signature::from_reader(buf)?;

        let flat_sigs = orig_sigs.into_iter().flat_map(|s| {
            s.signatures
                .iter()
                .map(|mh| {
                    let mut new_s = s.clone();
                    new_s.signatures = vec![mh.clone()];
                    new_s
                })
                .collect::<Vec<Signature>>()
        });

        let filtered_sigs = flat_sigs.filter_map(|mut sig| {
            let good_mhs: Vec<KmerMinHash> = sig
                .signatures
                .into_iter()
                .filter(|mh| {
                    if ksize == 0 || ksize == mh.ksize as usize {
                        match moltype {
                            Some(x) => {
                                if (x.to_lowercase() == "dna" && !mh.is_protein)
                                    || (x.to_lowercase() == "protein" && mh.is_protein)
                                {
                                    return true;
                                }
                            }
                            None => return true,
                        };
                    };
                    false
                })
                .collect();

            if good_mhs.is_empty() {
                return None;
            };

            sig.signatures = good_mhs;
            Some(sig)
        });

        Ok(filtered_sigs.collect())
    }
}

impl Default for Signature {
    fn default() -> Signature {
        Signature {
            class: default_class(),
            email: "".to_string(),
            hash_function: "0.murmur64".to_string(),
            license: default_license(),
            filename: None,
            name: None,
            signatures: Vec::<KmerMinHash>::new(),
            version: default_version(),
        }
    }
}

impl PartialEq for Signature {
    fn eq(&self, other: &Signature) -> bool {
        let metadata = self.class == other.class
            && self.email == other.email
            && self.hash_function == other.hash_function
            && self.filename == other.filename
            && self.name == other.name;

        let mh = &self.signatures[0];
        let other_mh = &other.signatures[0];
        metadata && (mh == other_mh)
    }
}

#[inline]
fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|n| match *n as char {
            'A' | 'a' => 'T',
            'T' | 't' => 'A',
            'C' | 'c' => 'G',
            'G' | 'g' => 'C',
            x => x,
        } as u8) // TODO: error?
        .collect()
}

lazy_static! {
    static ref CODONTABLE: HashMap<&'static str, u8> = {
        let mut m = HashMap::new();

        m.insert("TTT", b'F');
        m.insert("TTC", b'F');
        m.insert("TTA", b'L');
        m.insert("TTG", b'L');

        m.insert("TCT", b'S');
        m.insert("TCC", b'S');
        m.insert("TCA", b'S');
        m.insert("TCG", b'S');

        m.insert("TAT", b'Y');
        m.insert("TAC", b'Y');
        m.insert("TAA", b'*');
        m.insert("TAG", b'*');

        m.insert("TGT", b'C');
        m.insert("TGC", b'C');
        m.insert("TGA", b'*');
        m.insert("TGG", b'W');

        m.insert("CTT", b'L');
        m.insert("CTC", b'L');
        m.insert("CTA", b'L');
        m.insert("CTG", b'L');

        m.insert("CCT", b'P');
        m.insert("CCC", b'P');
        m.insert("CCA", b'P');
        m.insert("CCG", b'P');

        m.insert("CAT", b'H');
        m.insert("CAC", b'H');
        m.insert("CAA", b'Q');
        m.insert("CAG", b'Q');

        m.insert("CGT", b'R');
        m.insert("CGC", b'R');
        m.insert("CGA", b'R');
        m.insert("CGG", b'R');

        m.insert("ATT", b'I');
        m.insert("ATC", b'I');
        m.insert("ATA", b'I');
        m.insert("ATG", b'M');

        m.insert("ACT", b'T');
        m.insert("ACC", b'T');
        m.insert("ACA", b'T');
        m.insert("ACG", b'T');

        m.insert("AAT", b'N');
        m.insert("AAC", b'N');
        m.insert("AAA", b'K');
        m.insert("AAG", b'K');

        m.insert("AGT", b'S');
        m.insert("AGC", b'S');
        m.insert("AGA", b'R');
        m.insert("AGG", b'R');

        m.insert("GTT", b'V');
        m.insert("GTC", b'V');
        m.insert("GTA", b'V');
        m.insert("GTG", b'V');

        m.insert("GCT", b'A');
        m.insert("GCC", b'A');
        m.insert("GCA", b'A');
        m.insert("GCG", b'A');

        m.insert("GAT", b'D');
        m.insert("GAC", b'D');
        m.insert("GAA", b'E');
        m.insert("GAG", b'E');

        m.insert("GGT", b'G');
        m.insert("GGC", b'G');
        m.insert("GGA", b'G');
        m.insert("GGG", b'G');

        m
    };
}

#[inline]
fn to_aa(seq: &[u8]) -> Vec<u8> {
    let mut converted: Vec<u8> = Vec::with_capacity(seq.len() / 3);

    for chunk in seq.chunks(3) {
        if chunk.len() != 3 {
            break;
        }
        if let Some(codon) = CODONTABLE.get(str::from_utf8(chunk).unwrap()) {
            converted.push(*codon);
        }
    }

    converted
}

#[inline]
fn _checkdna(seq: &[u8]) -> bool {
    for n in seq.iter() {
        match *n as char {
            'A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'T' | 't' => (),
            _ => return false,
        }
    }
    true
}
