//! # Compressed representations of genomic data
//!
//! A signature is a collection of sketches for a genomic dataset.

use core::iter::FusedIterator;

use std::fs::File;
use std::io;
use std::path::Path;
use std::str;

use cfg_if::cfg_if;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::encodings::{aa_to_dayhoff, aa_to_hp, revcomp, to_aa, HashFunctions, VALID};
use crate::prelude::*;
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
use crate::Error;
use crate::HashIntoType;

// TODO: this is the behavior expected from Sketch, but that name is already
// used. Sketchable?
pub trait SigsTrait {
    fn size(&self) -> usize;
    fn to_vec(&self) -> Vec<u64>;
    fn ksize(&self) -> usize;
    fn check_compatible(&self, other: &Self) -> Result<(), Error>;
    fn seed(&self) -> u64;

    fn hash_function(&self) -> HashFunctions;

    fn add_hash(&mut self, hash: HashIntoType);

    fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        let ready_hashes = SeqToHashes::new(
            seq,
            self.ksize(),
            force,
            false,
            self.hash_function(),
            self.seed(),
        );

        for hash_value in ready_hashes {
            match hash_value {
                Ok(0) => continue,
                Ok(x) => self.add_hash(x),
                Err(err) => return Err(err),
            }
        }

        // Should be always ok
        Ok(())
    }

    fn add_protein(&mut self, seq: &[u8]) -> Result<(), Error> {
        let ready_hashes = SeqToHashes::new(
            seq,
            self.ksize(),
            false,
            true,
            self.hash_function(),
            self.seed(),
        );

        for hash_value in ready_hashes {
            match hash_value {
                Ok(0) => continue,
                Ok(x) => self.add_hash(x),
                Err(err) => return Err(err),
            }
        }

        // Should be always ok
        Ok(())
    }
}

impl SigsTrait for Sketch {
    fn size(&self) -> usize {
        match *self {
            Sketch::MinHash(ref mh) => mh.size(),
            Sketch::LargeMinHash(ref mh) => mh.size(),
            Sketch::HyperLogLog(ref hll) => hll.size(),
        }
    }

    fn to_vec(&self) -> Vec<u64> {
        match *self {
            Sketch::MinHash(ref mh) => mh.to_vec(),
            Sketch::LargeMinHash(ref mh) => mh.to_vec(),
            Sketch::HyperLogLog(ref hll) => hll.to_vec(),
        }
    }

    fn ksize(&self) -> usize {
        match *self {
            Sketch::MinHash(ref mh) => mh.ksize(),
            Sketch::LargeMinHash(ref mh) => mh.ksize(),
            Sketch::HyperLogLog(ref hll) => hll.ksize(),
        }
    }

    fn seed(&self) -> u64 {
        match *self {
            Sketch::MinHash(ref mh) => mh.seed(),
            Sketch::LargeMinHash(ref mh) => mh.seed(),
            Sketch::HyperLogLog(ref hll) => hll.seed(),
        }
    }

    fn hash_function(&self) -> HashFunctions {
        match *self {
            Sketch::MinHash(ref mh) => mh.hash_function(),
            Sketch::LargeMinHash(ref mh) => mh.hash_function(),
            Sketch::HyperLogLog(ref hll) => hll.hash_function(),
        }
    }

    fn add_hash(&mut self, hash: HashIntoType) {
        match *self {
            Sketch::MinHash(ref mut mh) => mh.add_hash(hash),
            Sketch::LargeMinHash(ref mut mh) => mh.add_hash(hash),
            Sketch::HyperLogLog(ref mut hll) => hll.add_hash(hash),
        }
    }

    fn check_compatible(&self, other: &Self) -> Result<(), Error> {
        match *self {
            Sketch::MinHash(ref mh) => match other {
                Sketch::MinHash(ref ot) => mh.check_compatible(ot),
                _ => Err(Error::MismatchSignatureType),
            },
            Sketch::LargeMinHash(ref mh) => match other {
                Sketch::LargeMinHash(ref ot) => mh.check_compatible(ot),
                _ => Err(Error::MismatchSignatureType),
            },
            Sketch::HyperLogLog(ref hll) => match other {
                Sketch::HyperLogLog(ref ot) => hll.check_compatible(ot),
                _ => Err(Error::MismatchSignatureType),
            },
        }
    }

    fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        match *self {
            Sketch::MinHash(ref mut mh) => mh.add_sequence(seq, force),
            Sketch::LargeMinHash(ref mut mh) => mh.add_sequence(seq, force),
            Sketch::HyperLogLog(_) => unimplemented!(),
        }
    }

    fn add_protein(&mut self, seq: &[u8]) -> Result<(), Error> {
        match *self {
            Sketch::MinHash(ref mut mh) => mh.add_protein(seq),
            Sketch::LargeMinHash(ref mut mh) => mh.add_protein(seq),
            Sketch::HyperLogLog(_) => unimplemented!(),
        }
    }
}

#[derive(Debug, Clone)]
pub enum ReadingFrame {
    DNA {
        fw: Vec<u8>,
        rc: Vec<u8>,
        len: usize, // len gives max_index for kmer iterator
    },
    Protein {
        fw: Vec<u8>, // Only forward frame
        len: usize,
    },
}

impl std::fmt::Display for ReadingFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadingFrame::DNA { fw, rc, len } => {
                let fw_str = String::from_utf8(fw.clone()).expect("Invalid UTF-8 sequence in fw");
                let rc_str = String::from_utf8(rc.clone()).expect("Invalid UTF-8 sequence in rc");
                write!(
                    f,
                    "Type: DNA ({}bp), Forward: {}, Reverse Complement: {}",
                    len, fw_str, rc_str
                )
            }
            ReadingFrame::Protein { fw, len } => {
                let fw_str = String::from_utf8(fw.clone()).expect("Invalid UTF-8 sequence in fw");
                write!(f, "Type: Protein ({}aa), Forward: {}", len, fw_str)
            }
        }
    }
}

impl ReadingFrame {
    pub fn new_dna(sequence: &[u8]) -> Self {
        let fw = sequence.to_vec();
        let rc = revcomp(sequence);
        let len = sequence.len();
        ReadingFrame::DNA { fw, rc, len }
    }

    pub fn new_protein(sequence: &[u8], dayhoff: bool, hp: bool) -> Self {
        let fw: Vec<u8> = if dayhoff {
            sequence.iter().map(|&aa| aa_to_dayhoff(aa)).collect()
        } else if hp {
            sequence.iter().map(|&aa| aa_to_hp(aa)).collect()
        } else {
            sequence.to_vec() // protein, as-is.
        };

        let len = fw.len();
        ReadingFrame::Protein { fw, len }
    }

    pub fn new_skipmer(seq: &[u8], start: usize, m: usize, n: usize) -> Self {
        if start >= n {
            panic!("Skipmer frame number must be < n ({})", n);
        }
        // Generate forward skipmer frame
        let fw: Vec<u8> = seq
            .iter()
            .skip(start)
            .enumerate()
            .filter_map(|(i, &base)| if i % n < m { Some(base) } else { None })
            .collect();

        let len = fw.len();
        let rc = revcomp(&fw);
        ReadingFrame::DNA { fw, rc, len }
    }

    pub fn new_translated(sequence: &[u8], frame_number: usize, dayhoff: bool, hp: bool) -> Self {
        if frame_number > 2 {
            panic!("Frame number must be 0, 1, or 2");
        }

        // translate sequence
        let fw: Vec<u8> = sequence
            .iter()
            .cloned()
            .skip(frame_number) // skip the initial bases for the frame
            .take(sequence.len() - frame_number) // adjust length based on skipped bases
            .collect::<Vec<u8>>() // collect the DNA subsequence
            .chunks(3) // group into codons (triplets)
            .filter_map(|codon| to_aa(codon, dayhoff, hp).ok()) // translate each codon
            .flatten() // flatten the nested results into a single sequence
            .collect();

        let len = fw.len();

        // return protein reading frame
        ReadingFrame::Protein { fw, len }
    }

    pub fn get_fw(&self) -> Option<&[u8]> {
        match self {
            ReadingFrame::DNA { fw, .. } => Some(fw),
            ReadingFrame::Protein { fw, .. } => Some(fw),
        }
    }

    pub fn get_rc(&self) -> Option<&[u8]> {
        match self {
            ReadingFrame::DNA { rc, .. } => Some(rc),
            ReadingFrame::Protein { .. } => None,
        }
    }

    pub fn kmer_count(&self, k_size: usize) -> usize {
        match self {
            ReadingFrame::DNA { len, .. } | ReadingFrame::Protein { len, .. } => {
                if *len >= k_size {
                    len - k_size + 1
                } else {
                    0 // No k-mers possible if len is smaller than k_size
                }
            }
        }
    }

    pub fn kmer_iter(&self, ksize: usize, seed: u64, force: bool) -> KmerIterator {
        KmerIterator::new(self, ksize, seed, force)
    }

}

#[derive(Debug, Clone)]
pub struct KmerIterator<'a> {
    frame: &'a ReadingFrame, // Reference to the ReadingFrame
    ksize: usize,
    index: usize,
    seed: u64,
    force: bool,
}

impl<'a> KmerIterator<'a> {
    pub fn new(frame: &'a ReadingFrame, ksize: usize, seed: u64, force: bool) -> Self {
        Self {
            frame,
            ksize,
            index: 0,
            seed,
            force,
        }
    }

    fn out_of_bounds(&self, length: usize) -> bool {
        self.index + self.ksize > length
    }

    fn validate_dna_kmer(&self, kmer: &[u8]) -> Result<(), Error> {
        for &nt in kmer {
            if !VALID[nt as usize] {
                return Err(Error::InvalidDNA {
                    message: String::from_utf8_lossy(kmer).to_string(),
                });
            }
        }
        Ok(())
    }
}

impl<'a> Iterator for KmerIterator<'a> {
    type Item = Result<u64, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.frame {
            ReadingFrame::DNA { fw, rc, len, .. } => {
                if self.out_of_bounds(*len) {
                    return None;
                }

                let kmer = &fw[self.index..self.index + self.ksize];
                if !self.force {
                    if let Err(e) = self.validate_dna_kmer(kmer) {
                        self.index += 1;
                        return Some(Err(e));
                    }
                }

                let krc = &rc[rc.len() - self.ksize - self.index..rc.len() - self.index];
                let hash = crate::_hash_murmur(std::cmp::min(kmer, krc), self.seed);
                // NTP TESTING
                eprintln!(
                    "Forward DNA k-mer: {}, Reverse Complement k-mer: {}, hash: {}",
                    String::from_utf8_lossy(kmer),
                    String::from_utf8_lossy(krc),
                    hash,
                    );
                self.index += 1;
                Some(Ok(hash))
            }
            ReadingFrame::Protein { fw, len, .. } => {
                if self.out_of_bounds(*len) {
                    return None;
                }
                let kmer = &fw[self.index..self.index + self.ksize];
                let hash = crate::_hash_murmur(kmer, self.seed);
                // NTP TESTING
                eprintln!(
                    "Protein k-mer: {}, hash: {}",
                    String::from_utf8_lossy(kmer),
                    hash
                );
                self.index += 1;
                Some(Ok(hash))
            }
        }
    }
}

pub struct SeqToHashes<'a> {
    k_size: usize,
    force: bool,
    seed: u64,
    frames: Vec<ReadingFrame>,
    frame_index: usize, // Index of the current frame
    current_kmer_iter: Option<KmerIterator<'a>>,
}

impl<'a> SeqToHashes<'a> {
    pub fn new(
        seq: &[u8],
        k_size: usize,
        force: bool,
        is_protein: bool,
        hash_function: HashFunctions,
        seed: u64,
    ) -> Self {
        let mut ksize: usize = k_size;

        // Adjust kmer size for protein-based hash functions
        if is_protein || hash_function.protein() || hash_function.dayhoff() || hash_function.hp() {
            ksize = k_size / 3;
        }

        // uppercase the sequence
        let sequence = seq.to_ascii_uppercase();

        // Generate frames based on sequence type and hash function
        let frames = if is_protein {
            Self::protein_frames(&sequence, &hash_function)
        } else if hash_function.protein() || hash_function.dayhoff() || hash_function.hp() {
            Self::translated_frames(&sequence, &hash_function)
        } else if hash_function.skipm1n3() || hash_function.skipm2n3() {
            Self::skipmer_frames(&sequence, &hash_function)
        } else {
            Self::dna_frames(&sequence)
        };

        SeqToHashes {
            k_size: ksize,
            force,
            seed,
            frames,
            frame_index: 0,
            current_kmer_iter: None,
        }
    }

    /// generate frames from DNA: 1 DNA frame (fw+rc)
    fn dna_frames(seq: &[u8]) -> Vec<ReadingFrame> {
        vec![ReadingFrame::new_dna(&seq)]
    }

    /// generate frames from protein: 1 protein frame
    fn protein_frames(seq: &[u8], hash_function: &HashFunctions) -> Vec<ReadingFrame> {
        vec![ReadingFrame::new_protein(
            &seq,
            hash_function.dayhoff(),
            hash_function.hp(),
        )]
    }

    /// generate translated frames: 6 protein frames
    fn translated_frames(seq: &[u8], hash_function: &HashFunctions) -> Vec<ReadingFrame> {
        let revcomp_sequence = revcomp(&seq);
        (0..3)
            .flat_map(|frame_number| {
                vec![
                    ReadingFrame::new_translated(
                        &seq,
                        frame_number,
                        hash_function.dayhoff(),
                        hash_function.hp(),
                    ),
                    ReadingFrame::new_translated(
                        &revcomp_sequence,
                        frame_number,
                        hash_function.dayhoff(),
                        hash_function.hp(),
                    ),
                ]
            })
            .collect()
    }

    /// generate skipmer frames: 3 DNA frames (each with fw+rc)
    fn skipmer_frames(seq: &[u8], hash_function: &HashFunctions) -> Vec<ReadingFrame> {
        let (m, n) = if hash_function.skipm1n3() {
            (1, 3)
        } else {
            (2, 3)
        };
        (0..3)
            .flat_map(|frame_number| vec![ReadingFrame::new_skipmer(&seq, frame_number, m, n)])
            .collect()
    }
}

impl<'a> Iterator for SeqToHashes<'a> {
    type Item = Result<u64, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.frame_index < self.frames.len() {
            // Initialize the kmer_iter for the current frame if it is None
            if self.current_kmer_iter.is_none() {
                let frame = &self.frames[self.frame_index];
                self.current_kmer_iter = Some(frame.kmer_iter(self.k_size, self.seed, self.force));
            }

            // Attempt to get the next k-mer from the current iterator
            if let Some(ref mut kmer_iter) = self.current_kmer_iter {
                if let Some(hash_result) = kmer_iter.next() {
                    return Some(hash_result);
                }
            }

            // If the current iterator is exhausted, move to the next frame
            self.current_kmer_iter = None;
            self.frame_index += 1;
        }

        // All frames and iterators are exhausted
        None
    }
}

// impl<'a> Iterator for SeqToHashes<'a> {
//     type Item = Result<u64, Error>;

//     fn next(&mut self) -> Option<Self::Item> {
//         while self.frame_index < self.frames.len() {
//             // Initialize kmer_iter for the current frame if it is None
//             if self.current_kmer_iter.is_none() {
//                 let frame = &self.frames[self.frame_index];
//                 self.current_kmer_iter = Some(frame.kmer_iter(self.k_size, self.seed, self.force));
//             }

//             // Attempt to get the next hash from the current iterator
//             if let Some(ref mut kmer_iter) = self.current_kmer_iter {
//                 if let Some(hash_result) = kmer_iter.next() {
//                     return Some(hash_result);
//                 }
//             }

//             // If the current iterator is exhausted, move to the next frame
//             self.current_kmer_iter = None;
//             self.frame_index += 1;
//         }

//         // All frames and iterators are exhausted
//         None
//     }
// }
// impl Iterator for SeqToHashes {
//     type Item = Result<u64, Error>;

//     fn next(&mut self) -> Option<Self::Item> {
//         // Iterate over the frames using frame_index
//         while self.frame_index < self.frames.len() {
//             let frame = &self.frames[self.frame_index];
//             // Create a KmerIterator for the current frame
//             let mut kmer_iter = frame.kmer_iter(self.k_size, self.seed, self.force);

//             // Process k-mers in the current frame
//             if let Some(hash_result) = kmer_iter.next() {
//                 return Some(hash_result); // Return the next hash
//             }

//             // Move to the next frame if the current one is exhausted
//             self.frame_index += 1;
//         }

//         // All frames exhausted
//         None
//     }
// }

#[derive(Serialize, Deserialize, Debug, Clone, TypedBuilder)]
#[cfg_attr(
    feature = "rkyv",
    derive(rkyv::Serialize, rkyv::Deserialize, rkyv::Archive)
)]
pub struct Signature {
    #[serde(default = "default_class")]
    #[builder(default = default_class())]
    class: String,

    #[serde(default)]
    #[builder(default)]
    email: String,

    #[builder(setter(into))]
    hash_function: String,

    #[builder(default)]
    filename: Option<String>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub(crate) name: Option<String>,

    #[serde(default = "default_license")]
    #[builder(default = default_license())]
    license: String,

    pub(crate) signatures: Vec<Sketch>,

    #[serde(default = "default_version")]
    #[builder(default = default_version())]
    version: f64,
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
    pub fn name(&self) -> String {
        if let Some(name) = &self.name {
            name.clone()
        } else if let Some(filename) = &self.filename {
            filename.clone()
        } else {
            self.md5sum()
        }
    }

    pub fn set_name(&mut self, name: &str) {
        self.name = Some(name.into())
    }

    pub fn filename(&self) -> String {
        if let Some(filename) = &self.filename {
            filename.clone()
        } else {
            "".into()
        }
    }

    pub fn set_filename(&mut self, name: &str) {
        self.filename = Some(name.into())
    }

    pub fn size(&self) -> usize {
        self.signatures.len()
    }

    pub fn sketches(&self) -> Vec<Sketch> {
        self.signatures.clone()
    }

    pub fn reset_sketches(&mut self) {
        self.signatures = vec![];
    }

    pub fn push(&mut self, sketch: Sketch) {
        self.signatures.push(sketch);
    }

    pub fn license(&self) -> String {
        self.license.clone()
    }

    pub fn class(&self) -> String {
        self.class.clone()
    }

    pub fn hash_function(&self) -> String {
        self.hash_function.clone()
    }

    pub fn email(&self) -> String {
        self.email.clone()
    }

    pub fn md5sum(&self) -> String {
        if self.signatures.len() == 1 {
            match &self.signatures[0] {
                Sketch::MinHash(mh) => mh.md5sum(),
                Sketch::LargeMinHash(mh) => mh.md5sum(),
                Sketch::HyperLogLog(_) => unimplemented!(),
            }
        } else {
            // TODO: select the correct signature
            unimplemented!()
        }
    }

    pub fn select_sketch(&self, sketch: &Sketch) -> Option<&Sketch> {
        if let Sketch::MinHash(template) = sketch {
            for sk in &self.signatures {
                if let Sketch::MinHash(mh) = sk {
                    if mh.check_compatible(template).is_ok() {
                        return Some(sk);
                    }
                } else {
                    unimplemented!()
                }
            }
        } else {
            unimplemented!()
        }
        None
    }

    // return single corresponding sketch
    pub fn get_sketch(&self) -> Option<&Sketch> {
        if self.signatures.len() != 1 {
            if self.signatures.len() > 1 {
                todo!("Multiple sketches found! Please run select first.");
            }
            return None;
        }
        self.signatures.iter().find(|sk| {
            matches!(
                sk,
                Sketch::MinHash(_) | Sketch::LargeMinHash(_) | Sketch::HyperLogLog(_)
            )
        })
    }

    // return minhash directly
    pub fn minhash(&self) -> Option<&KmerMinHash> {
        if self.signatures.len() != 1 {
            if self.signatures.len() > 1 {
                todo!("Multiple sketches found! Please run select first.");
            }
            return None;
        }
        self.signatures.iter().find_map(|sk| {
            if let Sketch::MinHash(mh) = sk {
                Some(mh)
            } else {
                None
            }
        })
    }

    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Vec<Signature>, Error> {
        let mut reader = io::BufReader::new(File::open(path)?);
        Signature::from_reader(&mut reader)
    }

    pub fn from_reader<R>(rdr: R) -> Result<Vec<Signature>, Error>
    where
        R: io::Read,
    {
        let (rdr, _format) = niffler::get_reader(Box::new(rdr))?;

        let sigs: Vec<Signature> = serde_json::from_reader(rdr)?;
        Ok(sigs)
    }

    pub fn load_signatures<R>(
        buf: R,
        ksize: Option<usize>,
        moltype: Option<HashFunctions>,
        _scaled: Option<u64>,
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
            let good_mhs: Vec<Sketch> = sig
                .signatures
                .into_iter()
                .filter(|sig| {
                    match sig {
                        Sketch::MinHash(mh) => {
                            if let Some(k) = ksize {
                                if k != mh.ksize() {
                                    return false;
                                }
                            };

                            match &moltype {
                                Some(x) => {
                                    if mh.hash_function() == *x {
                                        return true;
                                    }
                                }
                                None => return true, // TODO: match previous behavior
                            };
                        }
                        Sketch::LargeMinHash(mh) => {
                            if let Some(k) = ksize {
                                if k != mh.ksize() {
                                    return false;
                                }
                            };

                            match &moltype {
                                Some(x) => {
                                    if mh.hash_function() == *x {
                                        return true;
                                    }
                                }
                                None => return true, // TODO: match previous behavior
                            };
                        }
                        Sketch::HyperLogLog(_) => unimplemented!(),
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

    pub fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        cfg_if! {
        if #[cfg(feature = "parallel")] {
            self.signatures
                .par_iter_mut()
                .try_for_each(|sketch| {
                    sketch.add_sequence(seq, force) }
                )?;
        } else {
            for sketch in self.signatures.iter_mut(){
                sketch.add_sequence(seq, force)?;
            }
        }
        }

        Ok(())
    }

    pub fn add_protein(&mut self, seq: &[u8]) -> Result<(), Error> {
        cfg_if! {
        if #[cfg(feature = "parallel")] {
            self.signatures
                .par_iter_mut()
                .try_for_each(|sketch| {
                    sketch.add_protein(seq) }
                )?;
        } else {
            self.signatures
                .iter_mut()
                .try_for_each(|sketch| {
                    sketch.add_protein(seq) }
                )?;
        }
        }

        Ok(())
    }

    pub fn iter_mut(&mut self) -> IterMut<'_> {
        let length = self.signatures.len();
        IterMut {
            iter: self.signatures.iter_mut(),
            length,
        }
    }

    pub fn iter(&self) -> Iter<'_> {
        let length = self.signatures.len();
        Iter {
            iter: self.signatures.iter(),
            length,
        }
    }
}

pub struct IterMut<'a> {
    iter: std::slice::IterMut<'a, Sketch>,
    length: usize,
}

impl<'a> IntoIterator for &'a mut Signature {
    type Item = &'a mut Sketch;
    type IntoIter = IterMut<'a>;

    fn into_iter(self) -> IterMut<'a> {
        self.iter_mut()
    }
}

impl<'a> Iterator for IterMut<'a> {
    type Item = &'a mut Sketch;

    fn next(&mut self) -> Option<&'a mut Sketch> {
        if self.length == 0 {
            None
        } else {
            self.length -= 1;
            self.iter.next()
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.length, Some(self.length))
    }
}

pub struct Iter<'a> {
    iter: std::slice::Iter<'a, Sketch>,
    length: usize,
}

impl<'a> Iterator for Iter<'a> {
    type Item = &'a Sketch;

    fn next(&mut self) -> Option<&'a Sketch> {
        if self.length == 0 {
            None
        } else {
            self.length -= 1;
            self.iter.next()
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        (self.length, Some(self.length))
    }
}

impl FusedIterator for Iter<'_> {}

impl ExactSizeIterator for Iter<'_> {
    fn len(&self) -> usize {
        self.length
    }
}

impl Clone for Iter<'_> {
    fn clone(&self) -> Self {
        Iter {
            iter: self.iter.clone(),
            length: self.length,
        }
    }
}

impl ToWriter for Signature {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: io::Write,
    {
        serde_json::to_writer(writer, &vec![&self])?;
        Ok(())
    }
}

impl Select for Signature {
    fn select(mut self, selection: &Selection) -> Result<Self, Error> {
        self.signatures.retain(|s| {
            let mut valid = true;
            valid = if let Some(ksize) = selection.ksize() {
                let k = s.ksize() as u32;
                let adjusted_ksize = match s.hash_function() {
                    HashFunctions::Murmur64Protein
                    | HashFunctions::Murmur64Dayhoff
                    | HashFunctions::Murmur64Hp => ksize * 3,
                    _ => ksize,
                };
                k == adjusted_ksize
            } else {
                valid
            };
            // keep compatible scaled if applicable
            valid = if let Some(sel_scaled) = selection.scaled() {
                match s {
                    Sketch::MinHash(mh) => valid && mh.scaled() <= sel_scaled,
                    // TODO: test LargeMinHash
                    // Sketch::LargeMinHash(lmh) => valid && lmh.scaled() <= sel_scaled as u64,
                    _ => valid, // other sketch types or invalid cases
                }
            } else {
                valid // if selection.scaled() is None, keep prior valid
            };
            /*
            valid = if let Some(abund) = selection.abund() {
                valid && *s.with_abundance() == abund
            } else {
                valid
            };
            valid = if let Some(moltype) = selection.moltype() {
                valid && s.moltype() == moltype
            } else {
                valid
            };
            */

            valid
        });

        // downsample the retained sketches if needed.
        if let Some(sel_scaled) = selection.scaled() {
            for sketch in self.signatures.iter_mut() {
                // TODO: also account for LargeMinHash
                if let Sketch::MinHash(mh) = sketch {
                    if mh.scaled() < sel_scaled {
                        *sketch = Sketch::MinHash(mh.clone().downsample_scaled(sel_scaled)?);
                    }
                }
            }
        }
        Ok(self)
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
            signatures: Vec::<Sketch>::new(),
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

        // TODO: find the right signature
        // as long as we have a matching
        if let Sketch::MinHash(mh) = &self.signatures[0] {
            if let Sketch::MinHash(other_mh) = &other.signatures[0] {
                return metadata && (mh == other_mh);
            }
        } else {
            unimplemented!()
        }
        metadata
    }
}

impl TryInto<KmerMinHash> for Signature {
    type Error = Error;

    fn try_into(self) -> Result<KmerMinHash, Error> {
        match self.signatures.len() {
            1 => self
                .signatures
                .into_iter()
                .find_map(|sk| {
                    if let Sketch::MinHash(mh) = sk {
                        Some(mh)
                    } else {
                        None
                    }
                })
                .ok_or(Error::NoMinHashFound),
            0 => Err(Error::EmptySignature),
            _ => Err(Error::MultipleSketchesFound),
        }
    }
}

#[cfg(test)]
mod test {

    use super::*;
    use std::fs::File;
    use std::io::{BufReader, Read};
    use std::path::PathBuf;

    use needletail::parse_fastx_reader;

    use crate::cmd::ComputeParameters;
    use crate::encodings::HashFunctions;
    use crate::signature::SeqToHashes;
    use crate::signature::SigsTrait;

    use super::Signature;

    use crate::prelude::Select;
    use crate::selection::Selection;
    use crate::sketch::Sketch;

    #[test]
    fn load_sig() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/.sbt.v3/60f7e23c24a8d94791cc7a8680c493f9");

        let mut reader = BufReader::new(File::open(filename).unwrap());
        let sigs = Signature::load_signatures(
            &mut reader,
            Some(31),
            Some("DNA".try_into().unwrap()),
            None,
        )
        .unwrap();
        let _sig_data = sigs[0].clone();
        // TODO: check sig_data
    }

    #[test]
    fn load_signature() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/genome-s10+s11.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        assert_eq!(sigs.len(), 4);

        let sig = sigs.get(0).unwrap();
        assert_eq!(sig.class, "sourmash_signature");
        assert_eq!(sig.email, "");
        if let Some(ref filename) = sig.filename {
            assert_eq!(filename, "-");
        }
        assert_eq!(sig.hash_function, "0.murmur64");
        if let Some(ref name) = sig.name {
            assert_eq!(name, "genome-s10+s11");
        }
        assert_eq!(sig.signatures.len(), 1);
    }

    #[test]
    fn signature_from_computeparams() {
        let params = ComputeParameters::builder()
            .ksizes(vec![2, 3, 4])
            .num_hashes(3u32)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGC", false).unwrap();

        assert_eq!(sig.signatures.len(), 3);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 2);
        assert_eq!(sig.signatures[2].size(), 1);
    }

    #[test]
    fn signature_slow_path() {
        let params = ComputeParameters::builder()
            .ksizes(vec![2, 3, 4, 5])
            .num_hashes(3u32)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCTN", true).unwrap();

        assert_eq!(sig.signatures.len(), 4);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 3);
        assert_eq!(sig.signatures[2].size(), 2);
        assert_eq!(sig.signatures[3].size(), 1);
    }

    #[test]
    fn signature_add_sequence_protein() {
        let params = ComputeParameters::builder()
            .ksizes(vec![3, 6])
            .num_hashes(3u32)
            .protein(true)
            .dna(false)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCAT", false).unwrap();

        assert_eq!(sig.signatures.len(), 2);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 1);
    }

    #[test]
    fn signature_add_protein() {
        let params = ComputeParameters::builder()
            .ksizes(vec![3, 6])
            .num_hashes(3u32)
            .protein(true)
            .dna(false)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_protein(b"AGY").unwrap();

        assert_eq!(sig.signatures.len(), 2);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 2);
    }

    #[test]
    fn signature_skipm2n3_add_sequence() {
        let params = ComputeParameters::builder()
            .ksizes(vec![3, 4, 5, 6])
            .num_hashes(3u32)
            .dna(false)
            .skipm2n3(true)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCATGA", false).unwrap();

        assert_eq!(sig.signatures.len(), 4);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 3);
        assert_eq!(sig.signatures[2].size(), 2);
        assert_eq!(sig.signatures[3].size(), 1);
    }

    #[test]
    fn signature_skipm1n3_add_sequence() {
        let params = ComputeParameters::builder()
            .ksizes(vec![3, 4, 5, 6])
            .num_hashes(3u32)
            .dna(false)
            .skipm1n3(true)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCATGA", false).unwrap();

        assert_eq!(sig.signatures.len(), 4);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 3);
        assert_eq!(sig.signatures[1].size(), 3);
        assert_eq!(sig.signatures[2].size(), 2);
        assert_eq!(sig.signatures[3].size(), 1);
    }

    #[test]
    #[should_panic(expected = "not implemented")]
    fn signature_skipm2n3_add_sequence_too_small() {
        let params = ComputeParameters::builder()
            .ksizes(vec![2])
            .num_hashes(3u32)
            .dna(false)
            .skipm2n3(true)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCATGA", false).unwrap();
    }

    #[test]
    #[should_panic(expected = "not implemented")]
    fn signature_skipm1n3_add_sequence_too_small() {
        let params = ComputeParameters::builder()
            .ksizes(vec![2])
            .num_hashes(3u32)
            .dna(false)
            .skipm1n3(true)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCATGA", false).unwrap();
    }

    #[test]
    fn signature_add_sequence_cp() {
        let mut cp = ComputeParameters::default();
        cp.set_dayhoff(true);
        cp.set_protein(true);
        cp.set_hp(true);
        cp.set_dna(true);

        let mut sig = Signature::from_params(&cp);

        let mut data: Vec<u8> = vec![];
        let mut f = File::open("../../tests/test-data/ecoli.genes.fna").unwrap();
        let _ = f.read_to_end(&mut data);

        let mut parser = parse_fastx_reader(&data[..]).unwrap();
        while let Some(record) = parser.next() {
            let record = record.unwrap();
            sig.add_sequence(&record.seq(), false).unwrap();
        }

        assert_eq!(sig.size(), 12);
        for sk in sig.iter() {
            assert_eq!(sk.size(), 500);
        }
    }

    #[test]
    fn load_minhash_from_signature() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47.fa.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        assert_eq!(sigs.len(), 1);

        let sig = sigs.get(0).unwrap();
        let mh = sig.minhash().unwrap();
        assert_eq!(mh.scaled(), 1000);
    }

    #[test]
    fn load_single_sketch_from_signature() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47.fa.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        assert_eq!(sigs.len(), 1);

        let sig = sigs.get(0).unwrap();
        let mhdirect = sig.minhash().unwrap();
        let sketch = sig.get_sketch().unwrap();
        if let Sketch::MinHash(mh) = sketch {
            assert_eq!(mh.scaled(), 1000);
            assert_eq!(mhdirect, mh); // should be the same
        } else {
            // error
            assert!(false);
        }
    }

    #[test]
    #[should_panic]
    fn get_sketch_multisketch_panic() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47.fa.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        assert_eq!(sigs.len(), 1);

        let sig = sigs.get(0).unwrap();
        let mut mhdirect = sig.minhash().unwrap().clone();
        // change slightly and push into new_sig
        mhdirect.add_sequence(b"ATGGA", false).unwrap();
        let new_sketch = Sketch::MinHash(mhdirect.clone());
        let mut new_sig = sig.clone();
        new_sig.push(new_sketch);
        // check there are now two sketches in new_sig
        assert_eq!(new_sig.signatures.len(), 2);

        let _ = new_sig.get_sketch();
    }

    #[test]
    #[should_panic]
    fn load_minhash_multisketch_panic() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47.fa.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        assert_eq!(sigs.len(), 1);

        let sig = sigs.get(0).unwrap();
        let mut mhdirect = sig.minhash().unwrap().clone();
        // change slightly and push into new_sig
        mhdirect.add_sequence(b"ATGGA", false).unwrap();
        let new_sketch = Sketch::MinHash(mhdirect.clone());
        let mut new_sig = sig.clone();
        new_sig.push(new_sketch);
        // check there are now two sketches in new_sig
        assert_eq!(new_sig.signatures.len(), 2);

        let _ = new_sig.minhash();
    }

    #[test]
    fn selection_with_downsample() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47+63-multisig.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        // create Selection object
        let mut selection = Selection::default();
        selection.set_scaled(2000);
        // iterate and check scaled
        for sig in &sigs {
            let modified_sig = sig.clone().select(&selection).unwrap();
            for sketch in modified_sig.iter() {
                if let Sketch::MinHash(mh) = sketch {
                    dbg!("scaled: {:?}", mh.scaled());
                    assert_eq!(mh.scaled(), 2000);
                }
            }
        }
    }

    #[test]
    fn selection_protein() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push(
            "../../tests/test-data/prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig",
        );

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        // create Selection object
        let mut selection = Selection::default();
        let prot_ksize = 19;
        selection.set_ksize(prot_ksize);
        let selected_sig = sigs[0].clone().select(&selection).unwrap();
        let mh = selected_sig.minhash().unwrap();
        assert_eq!(mh.ksize(), prot_ksize as usize * 3);
    }

    #[test]
    fn selection_dayhoff() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push(
            "../../tests/test-data/prot/dayhoff/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig",
        );

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        // create Selection object
        let mut selection = Selection::default();
        let prot_ksize = 19;
        selection.set_ksize(prot_ksize);
        selection.set_moltype(crate::encodings::HashFunctions::Murmur64Dayhoff);
        let selected_sig = sigs[0].clone().select(&selection).unwrap();
        let mh = selected_sig.minhash().unwrap();
        assert_eq!(mh.ksize(), prot_ksize as usize * 3);
    }

    #[test]
    fn selection_hp() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename
            .push("../../tests/test-data/prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        // create Selection object
        let mut selection = Selection::default();
        let prot_ksize = 19;
        selection.set_ksize(prot_ksize);
        selection.set_moltype(crate::encodings::HashFunctions::Murmur64Hp);
        let selected_sig = sigs[0].clone().select(&selection).unwrap();
        let mh = selected_sig.minhash().unwrap();
        assert_eq!(mh.ksize(), prot_ksize as usize * 3);
    }

    #[test]
    fn selection_protein2() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push(
            "../../tests/test-data/prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig",
        );

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        // create Selection object
        let mut selection = Selection::default();
        let prot_ksize = 19;
        selection.set_ksize(prot_ksize * 3);
        let selected_sig = sigs[0].clone().select(&selection).unwrap();
        let mh = selected_sig.minhash();
        assert!(mh.is_none());
    }

    #[test]
    fn selection_scaled_too_low() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47+63-multisig.sig");

        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");

        // create Selection object
        let mut selection = Selection::default();
        selection.set_scaled(100);
        // iterate and check no sigs are returned (original scaled is 1000)
        for sig in &sigs {
            let modified_sig = sig.clone().select(&selection).unwrap();
            assert_eq!(modified_sig.size(), 0);
        }
    }


    #[test]
    fn test_seqtohashes_skipm2n3() {
        let sequence = b"AGTCGTCA";
        // let rc_seq = b"TGACGACT";
        let k_size = 5;
        let seed = 42;
        let force = true; // Force skip over invalid bases if needed
        let is_protein = false;

        // Initialize SeqToHashes iterator using the new constructor
        let mut seq_to_hashes = SeqToHashes::new(
            sequence,
            k_size,
            force,
            is_protein,
            HashFunctions::Murmur64Skipm2n3,
            seed,
        );

        // Define expected hashes for the skipmer configuration.
        let expected_kmers = ["AGCGC", "GTGTA"];
        // rc of the k-mer, not of the sequence, then skipmerized. Correct?
        let expected_krc = ["GCGCT", "TACAC"];

        // Compute expected hashes by hashing each k-mer with its reverse complement
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .zip(expected_krc.iter())
            .map(|(kmer, krc)| {
                // Convert both kmer and krc to byte slices and pass to _hash_murmur
                crate::_hash_murmur(std::cmp::min(kmer.as_bytes(), krc.as_bytes()), seed)
            })
            .collect();

        // Compare each produced hash from the iterator with the expected hash
        for expected_hash in expected_hashes {
            let hash = seq_to_hashes.next().unwrap().ok().unwrap();
            assert_eq!(hash, expected_hash, "Mismatch in skipmer hash");
        }
    }

    #[test]
    fn test_reading_frame_new_dna() {
        let sequence = b"AGTCGT";
        let hash_function = HashFunctions::Murmur64Dna;

        let frames = ReadingFrame::new_dna(sequence);

        assert_eq!(frames.get_fw(), Some(sequence.as_slice()));
        assert_eq!(frames.get_rc(), Some(b"ACGACT".as_slice()));
    }

    #[test]
    fn test_reading_frames_new_is_protein() {
        let sequence = b"MVLSPADKTNVKAAW";

        let frames = ReadingFrame::new_protein(sequence, false, false);

        assert_eq!(frames.get_fw(), Some(sequence.as_slice()));
        assert_eq!(frames.get_rc(), None); // No reverse complement for protein
    }

    #[test]
    fn test_reading_frame_translate() {
        let sequence = b"AGTCGTCGAGCT";
        let revcomp = revcomp(sequence);
        let mut frames = Vec::new();
        for i in 0..3 {
            let frame_fw = ReadingFrame::new_translated(sequence, i, false, false);
            eprintln!("frame: {}", frame_fw);
            frames.push(frame_fw);
            let frame_rc = ReadingFrame::new_translated(&revcomp, i, false, false);
            eprintln!("frame: {}", frame_rc);
            frames.push(frame_rc);
        }

        assert_eq!(frames[0].get_fw(), Some(b"SRRA".as_slice()));
        assert_eq!(frames[1].get_fw(), Some(b"SSTT".as_slice()));
        assert_eq!(frames[2].get_fw(), Some(b"VVE".as_slice()));
        assert_eq!(frames[3].get_fw(), Some(b"ARR".as_slice()));
        assert_eq!(frames[4].get_fw(), Some(b"SSS".as_slice()));
        assert_eq!(frames[5].get_fw(), Some(b"LDD".as_slice()));
    }

    #[test]
    fn test_reading_frame_skipmer_m1n3() {
        let sequence = b"AGTCGTCGAGCT";
        let m = 1;
        let n = 3;

        let mut frames = Vec::new();
        for start in 0..3 {
            let frame = ReadingFrame::new_skipmer(sequence, start, m, n);
            eprintln!("frame: {}", frame);
            frames.push(frame);
        }

        assert_eq!(frames.len(), 3); // Three skipmer frames

        // Expected skipmer sequences for m=1, n=3 (keep-1, skip-2)
        assert_eq!(frames[0].get_fw(), Some(b"ACCG".as_slice()));
        assert_eq!(frames[0].get_rc(), Some(b"CGGT".as_slice()));

        assert_eq!(frames[1].get_fw(), Some(b"GGGC".as_slice()));
        assert_eq!(frames[1].get_rc(), Some(b"GCCC".as_slice()));

        assert_eq!(frames[2].get_fw(), Some(b"TTAT".as_slice()));
        assert_eq!(frames[2].get_rc(), Some(b"ATAA".as_slice()));
    }

    #[test]
    fn test_reading_frame_skipmer_m2n3() {
        let sequence = b"AGTCGTCGAGCT";
        let m = 2;
        let n = 3;

        let mut frames = Vec::new();
        for start in 0..3 {
            let frame = ReadingFrame::new_skipmer(sequence, start, m, n);
            eprintln!("frame: {}", frame);
            frames.push(frame);
        }

        assert_eq!(frames.len(), 3); // Three skipmer frames

        // Expected skipmer sequences for m=1, n=3 (keep-1, skip-2)
        assert_eq!(frames[0].get_fw(), Some(b"AGCGCGGC".as_slice()));
        assert_eq!(frames[0].get_rc(), Some(b"GCCGCGCT".as_slice()));

        assert_eq!(frames[1].get_fw(), Some(b"GTGTGACT".as_slice()));
        assert_eq!(frames[1].get_rc(), Some(b"AGTCACAC".as_slice()));

        assert_eq!(frames[2].get_fw(), Some(b"TCTCAGT".as_slice()));
        assert_eq!(frames[2].get_rc(), Some(b"ACTGAGA".as_slice()));
    }

    #[test]
    fn test_reading_frame_kmer_iter() {
        let sequence = b"AGTCGT";
        let frame = ReadingFrame::new_dna(sequence);
        assert_eq!(frame.get_fw(), Some(sequence.as_slice()));
        assert_eq!(frame.get_rc(), Some(b"ACGACT".as_slice()));

        // Create a KmerIterator
        let mut kmer_iter = KmerIterator::new(&frame, 3, 42, false);

        // Expected k-mers from the forward and reverse complement sequence
        let expected_kmers = vec![
            (b"AGT".to_vec(), b"ACT".to_vec()),
            (b"GTC".to_vec(), b"GAC".to_vec()),
            (b"TCG".to_vec(), b"CGA".to_vec()),
            (b"CGT".to_vec(), b"ACG".to_vec()),
        ];

        // Compute expected hashes
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .map(|(fw_kmer, rc_kmer)| crate::_hash_murmur(std::cmp::min(fw_kmer, rc_kmer), 42))
            .collect();

        // Collect hashes produced by the kmer_iter
        let mut produced_hashes = Vec::new();

        while let Some(result) = kmer_iter.next() {
            match result {
                Ok(hash) => produced_hashes.push(hash),
                Err(e) => panic!("Error encountered during k-mer iteration: {:?}", e),
            }
        }

        // Assert that produced hashes match expected hashes
        assert_eq!(
            produced_hashes, expected_hashes,
            "Hashes do not match in order"
        );

        // Debugging output for verification
        eprintln!(
            "Expected hashes: {:?}\nProduced hashes: {:?}",
            expected_hashes, produced_hashes
        );
    }

    #[test]
    fn test_kmer_iter_is_protein() {
        let sequence = b"MVLSPADKTNVKAAW";
        let hash_function = HashFunctions::Murmur64Protein;

        let frame =
            ReadingFrame::new_protein(sequence, hash_function.dayhoff(), hash_function.hp());

        let kmer_iter = frame.kmer_iter(3, 42, false);

        // Expected k-mers for protein sequence
        let expected_kmers = vec![
            b"MVL".to_vec(),
            b"VLS".to_vec(),
            b"LSP".to_vec(),
            b"SPA".to_vec(),
            b"PAD".to_vec(),
            b"ADK".to_vec(),
            b"DKT".to_vec(),
            b"KTN".to_vec(),
            b"TNV".to_vec(),
            b"NVK".to_vec(),
            b"VKA".to_vec(),
            b"KAA".to_vec(),
            b"AAW".to_vec(),
        ];

        // Compute hashes for expected k-mers
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .map(|fw_kmer| crate::_hash_murmur(fw_kmer, 42))
            .collect();

        // Collect hashes produced by the kmer_iter
        let produced_hashes: Vec<u64> = kmer_iter.map(|result| result.unwrap()).collect();

        // Check that produced hashes match expected hashes in order
        assert_eq!(
            produced_hashes, expected_hashes,
            "Hashes do not match in order"
        );
    }

    #[test]
    fn test_kmer_iter_translate_frames() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Protein;
        let k_size =3;
        let seed = 42;
        let force = false;
        let is_protein = false;

        let sth = SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed);
        let frames = sth.frames;

        assert_eq!(frames[0].get_fw(), Some(b"SRRA".as_slice()));
        assert_eq!(frames[1].get_fw(), Some(b"SSTT".as_slice()));
        assert_eq!(frames[2].get_fw(), Some(b"VVE".as_slice()));
        assert_eq!(frames[3].get_fw(), Some(b"ARR".as_slice()));
        assert_eq!(frames[4].get_fw(), Some(b"SSS".as_slice()));
        assert_eq!(frames[5].get_fw(), Some(b"LDD".as_slice()));

        // Six translated frames
        assert_eq!(frames.len(), 6);

        // Expected k-mers for translated frames
        let f1_kmers =  vec![b"SRR".as_slice(), b"RRA".as_slice()];
        let f2_kmers = vec![b"SST".as_slice(), b"STT".as_slice()];
        let f3_kmers = vec![b"VVE".as_slice()];
        let f4_kmers = vec![b"ARR".as_slice()];
        let f5_kmers = vec![b"SSS".as_slice()];
        let f6_kmers = vec![b"LDD".as_slice()];
        let expected_kmers = vec![f1_kmers, f2_kmers, f3_kmers, f4_kmers, f5_kmers, f6_kmers];

        for (frame, expected_frame_kmers) in frames.iter().zip(expected_kmers.iter()) {
            let kmer_iter = frame.kmer_iter(k_size, seed, false);
            // Compute hashes for expected k-mers
            let expected_hashes: Vec<u64> = expected_frame_kmers
                .iter()
                .map(|fw_kmer| crate::_hash_murmur(fw_kmer, 42))
                .collect();

            // Collect hashes produced by the kmer_iter
            let produced_hashes: Vec<u64> = kmer_iter.map(|result| result.unwrap()).collect();

            // Check that produced hashes match expected hashes in order
            assert_eq!(
                produced_hashes, expected_hashes,
                "Hashes do not match in order for frame"
            );
        }
    }

    #[test]
    fn test_seqtohashes_kmer_iter_skipmer_m1n3() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Skipm1n3;
        let k_size = 3;
        let is_protein = false;
        let seed = 42;
        let force = false;

        let sth = SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed);
        // get hashes from sth
        let frames = sth.frames.clone();
        let sth_hashes: Vec<u64> = sth.map(|result| result.unwrap()).collect();
        eprintln!("sth_hashes: {:?}", sth_hashes);

        // Three skipmer frames
        assert_eq!(frames.len(), 3);
        assert_eq!(frames[0].get_fw(), Some(b"ACCG".as_slice()));
        assert_eq!(frames[0].get_rc(), Some(b"CGGT".as_slice()));
        let f1_kmers =  vec![(b"ACC".as_slice(), b"GGT".as_slice()), (b"CCG".as_slice(), b"CGG".as_slice())];

        assert_eq!(frames[1].get_fw(), Some(b"GGGC".as_slice()));
        assert_eq!(frames[1].get_rc(), Some(b"GCCC".as_slice()));
        let f2_kmers =  vec![(b"GGG".as_slice(), b"CCC".as_slice()), (b"GGC".as_slice(), b"GCC".as_slice())];

        assert_eq!(frames[2].get_fw(), Some(b"TTAT".as_slice()));
        assert_eq!(frames[2].get_rc(), Some(b"ATAA".as_slice()));
        let f3_kmers =  vec![(b"TTA".as_slice(), b"TAA".as_slice()), (b"TAT".as_slice(), b"ATA".as_slice())];

        // Expected k-mers for skipmer (m=1, n=3)
        let expected_kmers = vec![f1_kmers, f2_kmers, f3_kmers];

        let mut all_expected = Vec::new();
        for (frame, expected_frame_kmers) in frames.iter().zip(expected_kmers.iter()) {
            let kmer_iter = frame.kmer_iter(3, 42, false);

            // Compute hashes for expected k-mers
            let expected_hashes: Vec<u64> = expected_frame_kmers
                .iter()
                .map(|(fw_kmer, rc_kmer)| crate::_hash_murmur(std::cmp::min(fw_kmer, rc_kmer), 42))
                .collect();

            // Collect hashes produced by the kmer_iter
            let produced_hashes: Vec<u64> = kmer_iter.map(|result| result.unwrap()).collect();

            // Check that produced hashes match expected hashes in order
            assert_eq!(
                produced_hashes, expected_hashes,
                "Hashes do not match in order for frame"
            );
            // keep track of all expected hashes
            all_expected.extend(expected_hashes);
        }

        // Check that sth hashes match expected hashes in order
        // assert_eq!(
        //     sth_hashes, all_expected,
        //     "Hashes do not match in order for SeqToHashes"
        // );
    }


    #[test]
    fn test_kmer_iter_skipmer_m2n3() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Skipm2n3;
        let k_size = 7;
        let is_protein = false;
        let seed = 42;
        let force = false;

        let sth = SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed);
        let frames = sth.frames;

        // Three skipmer frames
        assert_eq!(frames.len(), 3);
        assert_eq!(frames[0].get_fw(), Some(b"AGCGCGGC".as_slice()));
        assert_eq!(frames[0].get_rc(), Some(b"GCCGCGCT".as_slice()));
        let f1_kmers =  vec![(b"AGCGCGG".as_slice(), b"CCGCGCT".as_slice()),
                                                  (b"GCGCGGC".as_slice(), b"GCCGCGC".as_slice())];

        assert_eq!(frames[1].get_fw(), Some(b"GTGTGACT".as_slice()));
        assert_eq!(frames[1].get_rc(), Some(b"AGTCACAC".as_slice()));
        let f2_kmers =  vec![(b"GTGTGAC".as_slice(), b"GTCACAC".as_slice()),
                                                  (b"TGTGACT".as_slice(), b"AGTCACA".as_slice())];

        assert_eq!(frames[2].get_fw(), Some(b"TCTCAGT".as_slice()));
        assert_eq!(frames[2].get_rc(), Some(b"ACTGAGA".as_slice()));
        let f3_kmers =  vec![(b"TCTCAGT".as_slice(), b"ACTGAGA".as_slice())];


        // Expected k-mers for skipmer (m=2, n=3)
        let expected_kmers = vec![f1_kmers, f2_kmers, f3_kmers];
       
        for (frame, expected_frame_kmers) in frames.iter().zip(expected_kmers.iter()) {
            // Compute hashes for expected k-mers
            let expected_hashes: Vec<u64> = expected_frame_kmers
                .iter()
                .map(|(fw_kmer, rc_kmer)| crate::_hash_murmur(std::cmp::min(fw_kmer, rc_kmer), 42))
                .collect();

            // Collect hashes produced by the kmer_iter
            let kmer_iter = frame.kmer_iter(k_size, seed, force);
            let produced_hashes: Vec<u64> = kmer_iter.map(|result| result.unwrap()).collect();

            // Check that produced hashes match expected hashes in order
            eprintln!("expected: {:?}, produced: {:?}", expected_hashes, produced_hashes);
            assert_eq!(
                produced_hashes, expected_hashes,
                "Hashes do not match in order for frame"
            );
        }
    }

#[test]
    fn test_seqtohashes_dna() {
        let sequence = b"AGTCGTCA";
        let k_size = 7;
        let seed = 42;
        let force = true; // Force skip over invalid bases if needed
        let is_protein = false;
        // Initialize SeqToHashes iterator using the new constructor
        let mut seq_to_hashes = SeqToHashes::new(
            sequence,
            k_size,
            force,
            is_protein,
            HashFunctions::Murmur64Dna,
            seed,
        );

        // Define expected hashes for the kmer configuration.
        let expected_kmers = ["AGTCGTC", "GTCGTCA"];
        let expected_krc = ["GACGACT", "TGACGAC"];

        // Compute expected hashes by hashing each k-mer with its reverse complement
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .zip(expected_krc.iter())
            .map(|(kmer, krc)| {
                // Convert both kmer and krc to byte slices and pass to _hash_murmur
                crate::_hash_murmur(std::cmp::min(kmer.as_bytes(), krc.as_bytes()), seed)
            })
            .collect();

        // Compare each produced hash from the iterator with the expected hash
        for expected_hash in expected_hashes {
            let hash = seq_to_hashes.next().unwrap().ok().unwrap();
            assert_eq!(hash, expected_hash, "Mismatch in DNA hash");
        }
    }
}