use std::cmp::Ordering;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::f64::consts::PI;
use std::iter::{Iterator, Peekable};
use std::str;

use failure::Error;
use lazy_static::lazy_static;
use serde::de::{Deserialize, Deserializer};
use serde::ser::{Serialize, SerializeStruct, Serializer};
use serde_derive::Deserialize;

use crate::_hash_murmur;
use crate::errors::SourmashError;
use crate::signature::SigsTrait;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq)]
#[repr(u32)]
pub enum HashFunctions {
    murmur64_DNA = 1,
    murmur64_protein = 2,
    murmur64_dayhoff = 3,
    murmur64_hp = 4,
}

impl TryFrom<&str> for HashFunctions {
    type Error = Error;

    fn try_from(moltype: &str) -> Result<Self, Self::Error> {
        match moltype.to_lowercase().as_ref() {
            "dna" => Ok(HashFunctions::murmur64_DNA),
            "dayhoff" => Ok(HashFunctions::murmur64_dayhoff),
            "hp" => Ok(HashFunctions::murmur64_hp),
            "protein" => Ok(HashFunctions::murmur64_protein),
            _ => unimplemented!(),
        }
    }
}

#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(Debug, Clone, PartialEq)]
pub struct KmerMinHash {
    num: u32,
    ksize: u32,
    pub(crate) hash_function: HashFunctions,
    seed: u64,
    max_hash: u64,
    pub(crate) mins: Vec<u64>,
    pub(crate) abunds: Option<Vec<u64>>,
}

impl Default for KmerMinHash {
    fn default() -> KmerMinHash {
        KmerMinHash {
            num: 1000,
            ksize: 21,
            hash_function: HashFunctions::murmur64_DNA,
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

        let mut partial = serializer.serialize_struct("KmerMinHash", n_fields)?;
        partial.serialize_field("num", &self.num)?;
        partial.serialize_field("ksize", &self.ksize)?;
        partial.serialize_field("seed", &self.seed)?;
        partial.serialize_field("max_hash", &self.max_hash)?;
        partial.serialize_field("mins", &self.mins)?;
        partial.serialize_field("md5sum", &self.md5sum())?;

        if let Some(abunds) = &self.abunds {
            partial.serialize_field("abundances", abunds)?;
        }

        partial.serialize_field(
            "molecule",
            match &self.is_protein() {
                true => {
                    if self.dayhoff() {
                        "dayhoff"
                    } else if self.hp() {
                        "hp"
                    } else {
                        "protein"
                    }
                }
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
            //md5sum: String,
            mins: Vec<u64>,
            abundances: Option<Vec<u64>>,
            molecule: String,
        }

        let tmpsig = TempSig::deserialize(deserializer)?;

        let num = if tmpsig.max_hash != 0 { 0 } else { tmpsig.num };
        let hash_function = match tmpsig.molecule.to_lowercase().as_ref() {
            "protein" => HashFunctions::murmur64_protein,
            "dayhoff" => HashFunctions::murmur64_dayhoff,
            "dna" => HashFunctions::murmur64_DNA,
            _ => unimplemented!(), // TODO: throw error here
        };

        Ok(KmerMinHash {
            num,
            ksize: tmpsig.ksize,
            seed: tmpsig.seed,
            max_hash: tmpsig.max_hash,
            mins: tmpsig.mins,
            abunds: tmpsig.abundances,
            hash_function,
        })
    }
}

impl KmerMinHash {
    pub fn new(
        num: u32,
        ksize: u32,
        hash_function: HashFunctions,
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
            hash_function,
            seed,
            max_hash,
            mins,
            abunds,
        }
    }

    pub fn num(&self) -> u32 {
        self.num
    }

    pub fn is_protein(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_protein
    }

    fn is_dna(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_DNA
    }

    pub fn seed(&self) -> u64 {
        self.seed
    }

    pub fn max_hash(&self) -> u64 {
        self.max_hash
    }

    pub fn md5sum(&self) -> String {
        let mut md5_ctx = md5::Context::new();
        md5_ctx.consume(self.ksize().to_string());
        self.mins
            .iter()
            .for_each(|x| md5_ctx.consume(x.to_string()));
        format!("{:x}", md5_ctx.compute())
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

    pub fn remove_hash(&mut self, hash: u64) {
        if let Ok(pos) = self.mins.binary_search(&hash) {
            if self.mins[pos] == hash {
                self.mins.remove(pos);
                if let Some(ref mut abunds) = self.abunds {
                    abunds.remove(pos);
                }
            }
        };
    }

    pub fn remove_many(&mut self, hashes: &[u64]) -> Result<(), Error> {
        for min in hashes {
            self.remove_hash(*min);
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
            self.abunds = if merged_abunds.is_empty() {
                None
            } else {
                Some(merged_abunds)
            };
        } else {
            self.mins = merged.into_iter().take(self.num as usize).collect();
            self.abunds = if merged_abunds.is_empty() {
                None
            } else {
                Some(merged_abunds.into_iter().take(self.num as usize).collect())
            }
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
        let iter = Intersection::new(self.mins.iter(), other.mins.iter());

        Ok(iter.count() as u64)
    }

    pub fn intersection(&self, other: &KmerMinHash) -> Result<(Vec<u64>, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.hash_function,
            self.seed,
            self.max_hash,
            self.abunds.is_some(),
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection::new(self.mins.iter(), other.mins.iter());

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let it2 = Intersection::new(i1.iter(), combined_mh.mins.iter());

        let common: Vec<u64> = it2.cloned().collect();
        Ok((common, combined_mh.mins.len() as u64))
    }

    pub fn intersection_size(&self, other: &KmerMinHash) -> Result<(u64, u64), Error> {
        self.check_compatible(other)?;

        let mut combined_mh = KmerMinHash::new(
            self.num,
            self.ksize,
            self.hash_function,
            self.seed,
            self.max_hash,
            self.abunds.is_some(),
        );

        combined_mh.merge(&self)?;
        combined_mh.merge(&other)?;

        let it1 = Intersection::new(self.mins.iter(), other.mins.iter());

        // TODO: there is probably a way to avoid this Vec here,
        // and pass the it1 as left in it2.
        let i1: Vec<u64> = it1.cloned().collect();
        let it2 = Intersection::new(i1.iter(), combined_mh.mins.iter());

        Ok((it2.count() as u64, combined_mh.mins.len() as u64))
    }

    pub fn compare(&self, other: &KmerMinHash) -> Result<f64, Error> {
        self.check_compatible(other)?;
        if let Ok((common, size)) = self.intersection_size(other) {
            Ok(common as f64 / u64::max(1, size) as f64)
        } else {
            Ok(0.0)
        }
    }

    pub fn similarity(&self, other: &KmerMinHash, ignore_abundance: bool) -> Result<f64, Error> {
        self.check_compatible(other)?;

        if ignore_abundance {
            if let Ok((common, size)) = self.intersection_size(other) {
                Ok(common as f64 / u64::max(1, size) as f64)
            } else {
                Ok(0.0)
            }
        } else {
            if self.abunds.is_none() || other.abunds.is_none() {
                // TODO: throw error, we need abundance for this
                unimplemented!()
            }

            let abunds = self.abunds.as_ref().unwrap();
            let other_abunds = other.abunds.as_ref().unwrap();

            let mut prod = 0;
            let mut other_iter = other.mins.iter().enumerate();
            let mut next_hash = other_iter.next();
            let a_sq: u64 = abunds.iter().map(|a| (a * a)).sum();
            let b_sq: u64 = other_abunds.iter().map(|a| (a * a)).sum();

            for (i, hash) in self.mins.iter().enumerate() {
                while let Some((j, k)) = next_hash {
                    match k.cmp(hash) {
                        Ordering::Less => next_hash = other_iter.next(),
                        Ordering::Equal => {
                            // Calling `get_unchecked` here is safe since
                            // both `i` and `j` are valid indices
                            // (`i` and `j` came from valid iterator calls)
                            unsafe {
                                prod += abunds.get_unchecked(i) * other_abunds.get_unchecked(j);
                            }
                            break;
                        }
                        Ordering::Greater => break,
                    }
                }
            }

            let norm_a = (a_sq as f64).sqrt();
            let norm_b = (b_sq as f64).sqrt();

            if norm_a == 0. || norm_b == 0. {
                return Ok(0.0);
            }

            let prod = f64::min(prod as f64 / (norm_a * norm_b), 1.);
            let distance = 2. * prod.acos() / PI;
            Ok(1. - distance)
        }
    }

    pub fn containment_ignore_maxhash(&self, other: &KmerMinHash) -> Result<f64, Error> {
        let it = Intersection::new(self.mins.iter(), other.mins.iter());

        Ok(it.count() as f64 / self.size() as f64)
    }

    pub fn dayhoff(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_dayhoff
    }

    pub fn hp(&self) -> bool {
        self.hash_function == HashFunctions::murmur64_hp
    }

    pub fn hash_function(&self) -> HashFunctions {
        self.hash_function
    }

    pub fn mins(&self) -> Vec<u64> {
        self.mins.clone()
    }
}

impl SigsTrait for KmerMinHash {
    fn size(&self) -> usize {
        self.mins.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.mins.clone()
    }

    fn ksize(&self) -> usize {
        self.ksize as usize
    }

    fn check_compatible(&self, other: &KmerMinHash) -> Result<(), Error> {
        /*
        if self.num != other.num {
            return Err(SourmashError::MismatchNum {
                n1: self.num,
                n2: other.num,
            }
            .into());
        }
        */
        if self.ksize != other.ksize {
            return Err(SourmashError::MismatchKSizes.into());
        }
        if self.hash_function != other.hash_function {
            // TODO: fix this error
            return Err(SourmashError::MismatchDNAProt.into());
        }
        if self.max_hash != other.max_hash {
            return Err(SourmashError::MismatchMaxHash.into());
        }
        if self.seed != other.seed {
            return Err(SourmashError::MismatchSeed.into());
        }
        Ok(())
    }

    fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        let sequence: Vec<u8> = seq
            .iter()
            .map(|&x| (x as char).to_ascii_uppercase() as u8)
            .collect();
        if sequence.len() >= (self.ksize as usize) {
            if self.is_dna() {
                for kmer in sequence.windows(self.ksize as usize) {
                    if _checkdna(kmer) {
                        let rc = revcomp(kmer);
                        if kmer < rc.as_slice() {
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
                    let aa = to_aa(&substr, self.dayhoff(), self.hp())?;

                    aa.windows(aa_ksize as usize).for_each(|n| self.add_word(n));

                    let rc_substr: Vec<u8> =
                        rc.iter().cloned().skip(i).take(rc.len() - i).collect();
                    let aa_rc = to_aa(&rc_substr, self.dayhoff(), self.hp())?;

                    aa_rc
                        .windows(aa_ksize as usize)
                        .for_each(|n| self.add_word(n));
                }
            }
        }
        Ok(())
    }
}

struct Intersection<T, I: Iterator<Item = T>> {
    iter: Peekable<I>,
    other: Peekable<I>,
}

impl<T, I: Iterator<Item = T>> Intersection<T, I> {
    pub fn new(left: I, right: I) -> Self {
        Intersection {
            iter: left.peekable(),
            other: right.peekable(),
        }
    }
}

impl<T: Ord, I: Iterator<Item = T>> Iterator for Intersection<T, I> {
    type Item = T;

    fn next(&mut self) -> Option<T> {
        loop {
            let res = match (self.iter.peek(), self.other.peek()) {
                (Some(ref left_key), Some(ref right_key)) => left_key.cmp(right_key),
                _ => return None,
            };

            match res {
                Ordering::Less => {
                    self.iter.next();
                }
                Ordering::Greater => {
                    self.other.next();
                }
                Ordering::Equal => {
                    self.other.next();
                    return self.iter.next();
                }
            }
        }
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
      [
        // F
        ("TTT", b'F'), ("TTC", b'F'),
        // L
        ("TTA", b'L'), ("TTG", b'L'),

        // S
        ("TCT", b'S'), ("TCC", b'S'), ("TCA", b'S'), ("TCG", b'S'), ("TCN", b'S'),

        // Y
        ("TAT", b'Y'), ("TAC", b'Y'),
        // *
        ("TAA", b'*'), ("TAG", b'*'),

        // *
        ("TGA", b'*'),
        // C
        ("TGT", b'C'), ("TGC", b'C'),
        // W
        ("TGG", b'W'),

        // L
        ("CTT", b'L'), ("CTC", b'L'), ("CTA", b'L'), ("CTG", b'L'), ("CTN", b'L'),

        // P
        ("CCT", b'P'), ("CCC", b'P'), ("CCA", b'P'), ("CCG", b'P'), ("CCN", b'P'),

        // H
        ("CAT", b'H'), ("CAC", b'H'),
        // Q
        ("CAA", b'Q'), ("CAG", b'Q'),

        // R
        ("CGT", b'R'), ("CGC", b'R'), ("CGA", b'R'), ("CGG", b'R'), ("CGN", b'R'),

        // I
        ("ATT", b'I'), ("ATC", b'I'), ("ATA", b'I'),
        // M
        ("ATG", b'M'),

        // T
        ("ACT", b'T'), ("ACC", b'T'), ("ACA", b'T'), ("ACG", b'T'), ("ACN", b'T'),

        // N
        ("AAT", b'N'), ("AAC", b'N'),
        // K
        ("AAA", b'K'), ("AAG", b'K'),

        // S
        ("AGT", b'S'), ("AGC", b'S'),
        // R
        ("AGA", b'R'), ("AGG", b'R'),

        // V
        ("GTT", b'V'), ("GTC", b'V'), ("GTA", b'V'), ("GTG", b'V'), ("GTN", b'V'),

        // A
        ("GCT", b'A'), ("GCC", b'A'), ("GCA", b'A'), ("GCG", b'A'), ("GCN", b'A'),

        // D
        ("GAT", b'D'), ("GAC", b'D'),
        // E
        ("GAA", b'E'), ("GAG", b'E'),

        // G
        ("GGT", b'G'), ("GGC", b'G'), ("GGA", b'G'), ("GGG", b'G'), ("GGN", b'G'),
        ].iter().cloned().collect()
    };
}

// Dayhoff table from
// Peris, P., López, D., & Campos, M. (2008).
// IgTM: An algorithm to predict transmembrane domains and topology in
// proteins. BMC Bioinformatics, 9(1), 1029–11.
// http://doi.org/10.1186/1471-2105-9-367
//
// Original source:
// Dayhoff M. O., Schwartz R. M., Orcutt B. C. (1978).
// A model of evolutionary change in proteins,
// in Atlas of Protein Sequence and Structure,
// ed Dayhoff M. O., editor.
// (Washington, DC: National Biomedical Research Foundation; ), 345–352.
//
// | Amino acid    | Property              | Dayhoff |
// |---------------|-----------------------|---------|
// | C             | Sulfur polymerization | a       |
// | A, G, P, S, T | Small                 | b       |
// | D, E, N, Q    | Acid and amide        | c       |
// | H, K, R       | Basic                 | d       |
// | I, L, M, V    | Hydrophobic           | e       |
// | F, W, Y       | Aromatic              | f       |
lazy_static! {
    static ref DAYHOFFTABLE: HashMap<u8, u8> = {
      [
        // a
        (b'C', b'a'),

        // b
        (b'A', b'b'), (b'G', b'b'), (b'P', b'b'), (b'S', b'b'), (b'T', b'b'),

        // c
        (b'D', b'c'), (b'E', b'c'), (b'N', b'c'), (b'Q', b'c'),

        // d
        (b'H', b'd'), (b'K', b'd'), (b'R', b'd'),

        // e
        (b'I', b'e'), (b'L', b'e'), (b'M', b'e'), (b'V', b'e'),

        // e
        (b'F', b'f'), (b'W', b'f'), (b'Y', b'f'),
        ].iter().cloned().collect()
    };
}

// HP Hydrophobic/hydrophilic mapping
// From: Phillips, R., Kondev, J., Theriot, J. (2008).
// Physical Biology of the Cell. New York: Garland Science, Taylor & Francis Group. ISBN: 978-0815341635

//
// | Amino acid                            | HP
// |---------------------------------------|---------|
// | A, F, G, I, L, M, P, V, W, Y          | h       |
// | N, C, S, T, D, E, R, H, K, Q          | p       |
lazy_static! {
    static ref HPTABLE: HashMap<u8, u8> = {
        [
            // h
            (b'A', b'h'), (b'F', b'h'), (b'G', b'h'), (b'I', b'h'), (b'L', b'h'),
            (b'M', b'h'), (b'P', b'h'), (b'V', b'h'), (b'W', b'h'), (b'Y', b'h'),

            // p
            (b'N', b'p'), (b'C', b'p'), (b'S', b'p'), (b'T', b'p'), (b'D', b'p'),
            (b'E', b'p'), (b'R', b'p'), (b'H', b'p'), (b'K', b'p'), (b'Q', b'p'),
        ]
        .iter()
        .cloned()
        .collect()
    };
}

#[inline]
pub(crate) fn translate_codon(codon: &[u8]) -> Result<u8, Error> {
    if codon.len() == 1 {
        return Ok(b'X');
    }

    if codon.len() == 2 {
        let mut v = codon.to_vec();
        v.push(b'N');
        match CODONTABLE.get(str::from_utf8(v.as_slice()).unwrap()) {
            Some(aa) => return Ok(*aa),
            None => return Ok(b'X'),
        }
    }

    if codon.len() == 3 {
        match CODONTABLE.get(str::from_utf8(codon).unwrap()) {
            Some(aa) => return Ok(*aa),
            None => return Ok(b'X'),
        }
    }

    Err(SourmashError::InvalidCodonLength {
        message: format!("{}", codon.len()),
    }
    .into())
}

#[inline]
pub(crate) fn aa_to_dayhoff(aa: u8) -> char {
    match DAYHOFFTABLE.get(&aa) {
        Some(letter) => *letter as char,
        None => 'X',
    }
}

pub(crate) fn aa_to_hp(aa: u8) -> char {
    match HPTABLE.get(&aa) {
        Some(letter) => *letter as char,
        None => 'X',
    }
}

#[inline]
fn to_aa(seq: &[u8], dayhoff: bool, hp: bool) -> Result<Vec<u8>, Error> {
    let mut converted: Vec<u8> = Vec::with_capacity(seq.len() / 3);

    for chunk in seq.chunks(3) {
        if chunk.len() < 3 {
            break;
        }

        let residue = translate_codon(chunk)?;
        if dayhoff {
            converted.push(aa_to_dayhoff(residue) as u8);
        } else if hp {
            converted.push(aa_to_hp(residue) as u8);
        } else {
            converted.push(residue);
        }
    }

    Ok(converted)
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
