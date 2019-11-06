use serde::de::{Deserialize, Deserializer};
use serde::ser::{Serialize, SerializeStruct, Serializer};
use serde_derive::Deserialize;

use std::cmp::Ordering;
use std::collections::HashMap;
use std::iter::{Iterator, Peekable};
use std::str;

use failure::Error;
use lazy_static::lazy_static;

use crate::_hash_murmur;
use crate::errors::SourmashError;
use crate::signature::SigsTrait;

#[cfg(all(target_arch = "wasm32", target_vendor = "unknown"))]
use wasm_bindgen::prelude::*;

#[cfg_attr(all(target_arch = "wasm32", target_vendor = "unknown"), wasm_bindgen)]
#[derive(Debug, Clone, PartialEq)]
pub struct KmerMinHash {
    num: u32,
    ksize: u32,
    is_protein: bool,
    dayhoff: bool,
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
            is_protein: false,
            dayhoff: false,
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
            match &self.is_protein {
                true => {
                    if self.dayhoff {
                        "dayhoff"
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
        let molecule = tmpsig.molecule.to_lowercase();

        Ok(KmerMinHash {
            num,
            ksize: tmpsig.ksize,
            seed: tmpsig.seed,
            max_hash: tmpsig.max_hash,
            mins: tmpsig.mins,
            abunds: tmpsig.abundances,
            is_protein: match molecule.as_ref() {
                "protein" => true,
                "dayhoff" => true,
                "dna" => false,
                _ => unimplemented!(),
            },
            dayhoff: molecule == "dayhoff",
        })
    }
}

impl KmerMinHash {
    pub fn new(
        num: u32,
        ksize: u32,
        is_protein: bool,
        dayhoff: bool,
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
            dayhoff,
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
        self.is_protein
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
            .map(|x| md5_ctx.consume(x.to_string()))
            .count();
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
            self.dayhoff,
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
            self.dayhoff,
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

    pub fn dayhoff(&self) -> bool {
        self.dayhoff
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
        if self.ksize != other.ksize {
            return Err(SourmashError::MismatchKSizes.into());
        }
        if self.is_protein != other.is_protein {
            return Err(SourmashError::MismatchDNAProt.into());
        }
        if self.dayhoff != other.dayhoff {
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
            if !self.is_protein {
                // dna
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
                    let aa = to_aa(&substr, self.dayhoff)?;

                    aa.windows(aa_ksize as usize)
                        .map(|n| self.add_word(n))
                        .count();

                    let rc_substr: Vec<u8> =
                        rc.iter().cloned().skip(i).take(rc.len() - i).collect();
                    let aa_rc = to_aa(&rc_substr, self.dayhoff)?;

                    aa_rc
                        .windows(aa_ksize as usize)
                        .map(|n| self.add_word(n))
                        .count();
                }
            }
        }
        Ok(())
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
        ].into_iter().cloned().collect()
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
        ].into_iter().cloned().collect()
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

#[inline]
fn to_aa(seq: &[u8], dayhoff: bool) -> Result<Vec<u8>, Error> {
    let mut converted: Vec<u8> = Vec::with_capacity(seq.len() / 3);

    for chunk in seq.chunks(3) {
        if chunk.len() < 3 {
            break;
        }

        let residue = translate_codon(chunk)?;
        if dayhoff {
            converted.push(aa_to_dayhoff(residue) as u8);
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
