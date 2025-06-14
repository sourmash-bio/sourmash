//! # Compressed representations of genomic data
//!
//! A signature is a collection of sketches for a genomic dataset.

use core::iter::FusedIterator;

use std::fs::File;
use std::io;
use std::path::Path;
use std::str;

use cfg_if::cfg_if;
use itertools::Itertools;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::encodings::{aa_to_dayhoff, aa_to_hp, revcomp, to_aa, HashFunctions, VALID};
use crate::errors::SourmashError;
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
        )?;

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
        )?;

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
        fw: Vec<u8>,
        len: usize,
    },
}

impl std::fmt::Display for ReadingFrame {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ReadingFrame::DNA { fw, rc, len } => {
                let fw_str = String::from_utf8_lossy(fw).to_string();
                let rc_str = String::from_utf8_lossy(rc).to_string();
                write!(
                    f,
                    "Type: DNA ({len}bp), Forward: {fw_str}, Reverse Complement: {rc_str}"
                )
            }
            ReadingFrame::Protein { fw, len } => {
                let fw_str = String::from_utf8_lossy(fw).to_string();
                write!(f, "Type: Protein ({len}aa), Forward: {fw_str}")
            }
        }
    }
}

impl ReadingFrame {
    pub fn new_dna(sequence: &[u8]) -> Self {
        let fw = sequence.to_ascii_uppercase();
        let rc = revcomp(&fw);
        let len = sequence.len();
        ReadingFrame::DNA { fw, rc, len }
    }

    pub fn new_protein(sequence: &[u8], dayhoff: bool, hp: bool) -> Self {
        let seq = sequence.to_ascii_uppercase();
        let fw: Vec<u8> = if dayhoff {
            seq.iter().map(|&aa| aa_to_dayhoff(aa)).collect()
        } else if hp {
            seq.iter().map(|&aa| aa_to_hp(aa)).collect()
        } else {
            seq
        };

        let len = fw.len();
        ReadingFrame::Protein { fw, len }
    }

    pub fn new_skipmer(
        sequence: &[u8],
        start: usize,
        m: usize,
        n: usize,
    ) -> Result<Self, SourmashError> {
        let seq = sequence.to_ascii_uppercase();
        if start >= n {
            return Err(SourmashError::InvalidSkipmerFrame { start, n });
        }
        // do we need to round up? (+1)
        let mut fw = Vec::with_capacity(((seq.len() * m) + 1) / n);
        seq.iter().skip(start).enumerate().for_each(|(i, &base)| {
            if i % n < m {
                fw.push(base.to_ascii_uppercase());
            }
        });

        let len = fw.len();
        let rc = revcomp(&fw);
        Ok(ReadingFrame::DNA { fw, rc, len })
    }

    // this is the only one that doesn't uppercase in here b/c more efficient to uppercase externally :/
    pub fn new_translated(
        sequence: &[u8],
        frame_number: usize,
        dayhoff: bool,
        hp: bool,
    ) -> Result<Self, SourmashError> {
        if frame_number > 2 {
            return Err(SourmashError::InvalidTranslateFrame { frame_number });
        }

        // Translate sequence into amino acids
        let mut fw = Vec::with_capacity(sequence.len() / 3);
        // NOTE: b/c of chunks(3), we only process full codons and ignore leftover bases (e.g. 1 or 2 at end of frame)
        sequence
            .iter()
            .skip(frame_number) // Skip the initial bases for the frame
            .take(sequence.len() - frame_number) // Adjust length based on skipped bases
            .chunks(3) // Group into codons (triplets) using itertools
            .into_iter()
            .filter_map(|chunk| {
                let codon: Vec<u8> = chunk.cloned().collect(); // Collect the chunk into a Vec<u8>
                to_aa(&codon, dayhoff, hp).ok() // Translate the codon
            })
            .for_each(|aa| fw.extend(aa)); // Extend `fw` with amino acids

        let len = fw.len();

        // return protein reading frame
        Ok(ReadingFrame::Protein { fw, len })
    }

    /// Get the forward sequence.
    #[inline]
    pub fn fw(&self) -> &[u8] {
        match self {
            ReadingFrame::DNA { fw, .. } => fw,
            ReadingFrame::Protein { fw, .. } => fw,
        }
    }

    /// Get the reverse complement sequence (if DNA).
    #[inline]
    pub fn rc(&self) -> &[u8] {
        match self {
            ReadingFrame::DNA { rc, .. } => rc,
            _ => panic!("Reverse complement is only available for DNA frames"),
        }
    }

    #[inline]
    pub fn length(&self) -> usize {
        match self {
            ReadingFrame::DNA { len, .. } => *len,
            ReadingFrame::Protein { len, .. } => *len,
        }
    }

    /// Get the type of the frame as a string.
    pub fn frame_type(&self) -> &'static str {
        match self {
            ReadingFrame::DNA { .. } => "DNA",
            ReadingFrame::Protein { .. } => "Protein",
        }
    }
}

pub struct SeqToHashes {
    k_size: usize,
    force: bool,
    seed: u64,
    frames: Vec<ReadingFrame>,
    frame_index: usize,         // Index of the current frame
    kmer_index: usize,          // Current k-mer index within the frame
    last_position_check: usize, // Index of last base we validated
}

impl SeqToHashes {
    pub fn new(
        seq: &[u8],
        k_size: usize,
        force: bool,
        is_protein: bool,
        hash_function: HashFunctions,
        seed: u64,
    ) -> Result<Self, SourmashError> {
        let mut ksize: usize = k_size;

        // Adjust kmer size for protein-based hash functions
        if is_protein || hash_function.protein() || hash_function.dayhoff() || hash_function.hp() {
            ksize = k_size / 3;
        }

        // Generate frames based on sequence type and hash function
        let frames = if hash_function.dna() {
            Self::dna_frames(seq)
        } else if is_protein {
            Self::protein_frames(seq, &hash_function)
        } else if hash_function.protein() || hash_function.dayhoff() || hash_function.hp() {
            Self::translated_frames(seq, &hash_function)?
        } else if hash_function.skipm1n3() || hash_function.skipm2n3() {
            Self::skipmer_frames(seq, &hash_function, ksize)?
        } else {
            return Err(SourmashError::InvalidHashFunction {
                function: format!("{hash_function:?}"),
            });
        };

        Ok(SeqToHashes {
            k_size: ksize,
            force,
            seed,
            frames,
            frame_index: 0,
            kmer_index: 0,
            last_position_check: 0,
        })
    }

    /// generate frames from DNA: 1 DNA frame (fw+rc)
    fn dna_frames(seq: &[u8]) -> Vec<ReadingFrame> {
        vec![ReadingFrame::new_dna(seq)]
    }

    /// generate frames from protein: 1 protein frame
    fn protein_frames(seq: &[u8], hash_function: &HashFunctions) -> Vec<ReadingFrame> {
        vec![ReadingFrame::new_protein(
            seq,
            hash_function.dayhoff(),
            hash_function.hp(),
        )]
    }

    /// generate translated frames: 6 protein frames
    fn translated_frames(
        seq: &[u8],
        hash_function: &HashFunctions,
    ) -> Result<Vec<ReadingFrame>, SourmashError> {
        // since we need to revcomp BEFORE making ReadingFrames, uppercase the sequence here
        let sequence = seq.to_ascii_uppercase();
        let revcomp_sequence = revcomp(&sequence);
        let frames = (0..3)
            .flat_map(|frame_number| {
                vec![
                    ReadingFrame::new_translated(
                        &sequence,
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
            .collect::<Result<Vec<_>, _>>()?;

        Ok(frames)
    }

    /// generate skipmer frames: 3 DNA frames (each with fw+rc)
    fn skipmer_frames(
        seq: &[u8],
        hash_function: &HashFunctions,
        ksize: usize,
    ) -> Result<Vec<ReadingFrame>, SourmashError> {
        let (m, n) = if hash_function.skipm1n3() {
            (1, 3)
        } else {
            (2, 3)
        };
        if ksize < n {
            return Err(SourmashError::InvalidSkipmerSize { ksize, n });
        }
        let frames = (0..3)
            .flat_map(|frame_number| vec![ReadingFrame::new_skipmer(seq, frame_number, m, n)])
            .collect::<Result<Vec<_>, _>>()?;

        Ok(frames)
    }

    fn out_of_bounds(&self, frame: &ReadingFrame) -> bool {
        self.kmer_index + self.k_size > frame.length()
    }
}

impl Iterator for SeqToHashes {
    type Item = Result<u64, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        while self.frame_index < self.frames.len() {
            let frame = &self.frames[self.frame_index];

            // Do we need to move to the next frame?
            if self.out_of_bounds(frame) {
                self.frame_index += 1;
                self.kmer_index = 0; // Reset for the next frame
                self.last_position_check = 0;
                continue;
            }

            let result = match frame {
                ReadingFrame::DNA { .. } => {
                    let kmer = &frame.fw()[self.kmer_index..self.kmer_index + self.k_size];
                    let rc = frame.rc();

                    // Validate k-mer bases
                    for j in std::cmp::max(self.kmer_index, self.last_position_check)
                        ..self.kmer_index + self.k_size
                    {
                        if !VALID[frame.fw()[j] as usize] {
                            if !self.force {
                                // Return an error if force is false
                                return Some(Err(Error::InvalidDNA {
                                    message: String::from_utf8(kmer.to_vec()).unwrap(),
                                }));
                            } else {
                                // Skip the invalid k-mer
                                self.kmer_index += 1;
                                return Some(Ok(0));
                            }
                        }
                        self.last_position_check += 1;
                    }

                    // Compute canonical hash
                    //   For a ksize = 3, and a sequence AGTCGT (len = 6):
                    //                   +-+---------+---------------+-------+
                    //   seq      RC     |i|i + ksize|len - ksize - i|len - i|
                    //  AGTCGT   ACGACT  +-+---------+---------------+-------+
                    //  +->         +->  |0|    2    |       3       |   6   |
                    //   +->       +->   |1|    3    |       2       |   5   |
                    //    +->     +->    |2|    4    |       1       |   4   |
                    //     +->   +->     |3|    5    |       0       |   3   |
                    //                   +-+---------+---------------+-------+
                    // (leaving this table here because I had to draw to
                    //  get the indices correctly)
                    let krc = &rc[frame.length() - self.k_size - self.kmer_index
                        ..frame.length() - self.kmer_index];
                    let hash = crate::_hash_murmur(std::cmp::min(kmer, krc), self.seed);
                    Ok(hash)
                }
                ReadingFrame::Protein { .. } => {
                    let kmer = &frame.fw()[self.kmer_index..self.kmer_index + self.k_size];
                    Ok(crate::_hash_murmur(kmer, self.seed))
                }
            };

            self.kmer_index += 1; // Advance k-mer index for valid k-mers
            return Some(result);
        }
        None // No more frames or k-mers
    }
}

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
    pub fn name(&self) -> Option<String> {
        self.name.clone()
    }

    /// return name, if not None; or "" if None.
    pub fn name_str(&self) -> String {
        self.name().unwrap_or("".into())
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

impl ToWriter for Vec<&Signature> {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: io::Write,
    {
        serde_json::to_writer(writer, &self)?;
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

    use std::fs::File;
    use std::io::{BufReader, Read};
    use std::path::PathBuf;

    use needletail::parse_fastx_reader;

    use crate::cmd::ComputeParameters;
    use crate::encodings::HashFunctions;
    use crate::signature::{ReadingFrame, SeqToHashes, SigsTrait};

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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

        assert_eq!(sig.name_str(), "");
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
        eprintln!("{:?}", sig.signatures[2]);
        assert_eq!(sig.signatures[2].size(), 3);
        assert_eq!(sig.signatures[3].size(), 1);
    }

    #[test]
    fn signature_skipm1n3_add_sequence() {
        let params = ComputeParameters::builder()
            .ksizes(vec![3, 4, 5, 6])
            .num_hashes(10u32)
            .dna(false)
            .skipm1n3(true)
            .build();

        let mut sig = Signature::from_params(&params);
        sig.add_sequence(b"ATGCATGAATGAC", false).unwrap();

        assert_eq!(sig.signatures.len(), 4);
        dbg!(&sig.signatures);
        assert_eq!(sig.signatures[0].size(), 5);
        assert_eq!(sig.signatures[1].size(), 4);
        assert_eq!(sig.signatures[2].size(), 1);
        assert_eq!(sig.signatures[3].size(), 0);
    }

    #[test]
    fn signature_skipm2n3_add_sequence_too_small() {
        let ksize = 2;
        let params = ComputeParameters::builder()
            .ksizes(vec![ksize])
            .num_hashes(10u32)
            .dna(false)
            .skipm2n3(true)
            .build();

        let mut sig = Signature::from_params(&params);
        let result = sig.add_sequence(b"ATGCATGA", false);

        match result {
            Err(error) => {
                // Convert the error to a string and check the message
                let error_message = format!("{}", error);
                assert_eq!(
                    error_message,
                    "Skipmer ksize must be >= n (3), but got ksize: 2"
                );
            }
            _ => panic!("Expected SourmashError::InvalidSkipmerSize"),
        }
    }

    #[test]
    fn signature_skipm1n3_add_sequence_too_small() {
        let params = ComputeParameters::builder()
            .ksizes(vec![2])
            .num_hashes(10u32)
            .dna(false)
            .skipm1n3(true)
            .build();

        let mut sig = Signature::from_params(&params);
        let result = sig.add_sequence(b"ATGCATGA", false);

        match result {
            Err(error) => {
                // Convert the error to a string and check the message
                let error_message = format!("{}", error);
                assert_eq!(
                    error_message,
                    "Skipmer ksize must be >= n (3), but got ksize: 2"
                );
            }
            _ => panic!("Expected SourmashError::InvalidSkipmerSize"),
        }
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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
        let sigs = Signature::from_reader(reader).expect("Loading error");

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
    fn test_readingframe_dna() {
        let sequence = b"AGTCGT";
        let frame = ReadingFrame::new_dna(sequence);

        assert_eq!(frame.fw(), sequence.as_slice());
        assert_eq!(frame.rc(), b"ACGACT".as_slice());
    }

    #[test]
    fn test_fw_dna() {
        let dna_frame = ReadingFrame::DNA {
            fw: b"ATCG".to_vec(),
            rc: b"CGAT".to_vec(),
            len: 4,
        };
        assert_eq!(dna_frame.fw(), b"ATCG");
    }

    #[test]
    fn test_rc_dna() {
        let dna_frame = ReadingFrame::DNA {
            fw: b"ATCG".to_vec(),
            rc: b"CGAT".to_vec(),
            len: 4,
        };
        assert_eq!(dna_frame.rc(), b"CGAT");
    }

    #[test]
    fn test_length_dna() {
        let dna_frame = ReadingFrame::DNA {
            fw: b"ATCG".to_vec(),
            rc: b"CGAT".to_vec(),
            len: 4,
        };
        assert_eq!(dna_frame.length(), 4);
    }

    #[test]
    fn test_frame_type_dna() {
        let dna_frame = ReadingFrame::DNA {
            fw: b"ATCG".to_vec(),
            rc: b"CGAT".to_vec(),
            len: 4,
        };
        assert_eq!(dna_frame.frame_type(), "DNA");
    }

    #[test]
    fn test_fw_protein() {
        let protein_frame = ReadingFrame::Protein {
            fw: b"MVHL".to_vec(),
            len: 4,
        };
        assert_eq!(protein_frame.fw(), b"MVHL");
    }

    #[test]
    #[should_panic(expected = "Reverse complement is only available for DNA frames")]
    fn test_rc_protein_panics() {
        let protein_frame = ReadingFrame::Protein {
            fw: b"MVHL".to_vec(),
            len: 4,
        };
        protein_frame.rc();
    }

    #[test]
    fn test_length_protein() {
        let protein_frame = ReadingFrame::Protein {
            fw: b"MVHL".to_vec(),
            len: 4,
        };
        assert_eq!(protein_frame.length(), 4);
    }

    #[test]
    fn test_frame_type_protein() {
        let protein_frame = ReadingFrame::Protein {
            fw: b"MVHL".to_vec(),
            len: 4,
        };
        assert_eq!(protein_frame.frame_type(), "Protein");
    }

    #[test]
    fn test_readingframe_display_protein() {
        // Create a Protein ReadingFrame
        let protein_frame = ReadingFrame::Protein {
            fw: b"MVHLK".to_vec(),
            len: 5,
        };

        let output = format!("{}", protein_frame);
        // Assert the output matches the expected format
        assert_eq!(output, "Type: Protein (5aa), Forward: MVHLK");
    }

    #[test]
    fn test_seqtohashes_frames_dna() {
        let sequence = b"AGTCGT";
        let hash_function = HashFunctions::Murmur64Dna;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames.clone();

        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].fw(), sequence.as_slice());
        assert_eq!(frames[0].rc(), b"ACGACT".as_slice());
    }

    #[test]
    fn test_seqtohashes_frames_is_protein() {
        let sequence = b"MVLSPADKTNVKAAW";
        let hash_function = HashFunctions::Murmur64Protein;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = true;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames.clone();

        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].fw(), sequence.as_slice());
    }

    #[test]
    fn test_readingframe_protein() {
        let sequence = b"MVLSPADKTNVKAAW";
        let hash_function = HashFunctions::Murmur64Protein;
        let frame =
            ReadingFrame::new_protein(sequence, hash_function.dayhoff(), hash_function.hp());

        assert_eq!(frame.fw(), sequence.as_slice());
    }

    #[test]
    #[should_panic]
    fn test_seqtohashes_frames_is_protein_try_access_rc() {
        // test panic if trying to access rc
        let sequence = b"MVLSPADKTNVKAAW";
        let hash_function = HashFunctions::Murmur64Protein;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = true;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames.clone();

        // protein frame doesn't have rc; this should panic
        eprintln!("{:?}", frames[0].rc());
    }

    #[test]
    fn test_seqtohashes_frames_is_protein_dayhoff() {
        let sequence = b"MVLSPADKTNVKAAW";
        let dayhoff_seq = b"eeebbbcdbcedbbf";
        let hash_function = HashFunctions::Murmur64Dayhoff;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = true;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames.clone();

        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].fw(), dayhoff_seq.as_slice());
    }

    #[test]
    fn test_seqtohashes_frames_is_protein_hp() {
        let sequence = b"MVLSPADKTNVKAAW";
        let hp_seq = b"hhhphhpppphphhh";
        let hash_function = HashFunctions::Murmur64Hp;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = true;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames.clone();

        assert_eq!(frames.len(), 1);
        assert_eq!(frames[0].fw(), hp_seq.as_slice());
    }

    #[test]
    fn test_seqtohashes_frames_translate_protein() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Protein;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames.clone();

        assert_eq!(frames[0].fw(), b"SRRA".as_slice());
        assert_eq!(frames[1].fw(), b"SSTT".as_slice());
        assert_eq!(frames[2].fw(), b"VVE".as_slice());
        assert_eq!(frames[3].fw(), b"ARR".as_slice());
        assert_eq!(frames[4].fw(), b"SSS".as_slice());
        assert_eq!(frames[5].fw(), b"LDD".as_slice());
    }

    #[test]
    fn test_readingframe_translate() {
        let sequence = b"AGTCGT";
        let frame_start = 3; // four frames but translate can only

        let result = ReadingFrame::new_translated(sequence, frame_start, false, false);

        match result {
            Err(error) => {
                // Convert the error to a string and check the message
                let error_message = format!("{}", error);
                assert_eq!(error_message, "Frame number must be 0, 1, or 2, but got 3");
            }
            _ => panic!("Expected SourmashError::InvalidTranslateFrame"),
        }
    }

    #[test]
    fn test_readingframe_skipmer() {
        let sequence = b"AGTCGT";
        let m = 2;
        let n = 3;
        let num_frames = 4; // four frames but n is only 3

        let result = ReadingFrame::new_skipmer(sequence, num_frames, m, n);

        match result {
            Err(error) => {
                // Convert the error to a string and check the message
                let error_message = format!("{}", error);
                assert_eq!(
                    error_message,
                    "Skipmer frame number must be < n (3), but got start: 4"
                );
            }
            _ => panic!("Expected SourmashError::InvalidSkipmerFrame"),
        }
    }

    #[test]
    fn test_seqtohashes_frames_skipmer_m1n3() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Skipm1n3; // Represents m=1, n=3
        let k_size = 3; // K-mer size is not directly relevant for skipmer frame validation
        let seed = 42; // Seed is also irrelevant for frame structure
        let force = false;
        let is_protein = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames.clone();

        eprintln!("Frames: {:?}", frames);

        assert_eq!(frames.len(), 3); // Three skipmer frames

        // Expected skipmer sequences for m=1, n=3 (keep-1, skip-2)
        assert_eq!(frames[0].fw(), b"ACCG".as_slice());
        assert_eq!(frames[0].rc(), b"CGGT".as_slice());

        assert_eq!(frames[1].fw(), b"GGGC".as_slice());
        assert_eq!(frames[1].rc(), b"GCCC".as_slice());

        assert_eq!(frames[2].fw(), b"TTAT".as_slice());
        assert_eq!(frames[2].rc(), b"ATAA".as_slice());
    }

    #[test]
    fn test_seqtohashes_frames_skipmer_m2n3() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Skipm2n3;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames;
        eprintln!("Frames: {:?}", frames);

        assert_eq!(frames.len(), 3); // Three skipmer frames

        // Expected skipmer sequences for m=1, n=3 (keep-1, skip-2)
        assert_eq!(frames[0].fw(), b"AGCGCGGC".as_slice());
        assert_eq!(frames[0].rc(), b"GCCGCGCT".as_slice());

        assert_eq!(frames[1].fw(), b"GTGTGACT".as_slice());
        assert_eq!(frames[1].rc(), b"AGTCACAC".as_slice());

        assert_eq!(frames[2].fw(), b"TCTCAGT".as_slice());
        assert_eq!(frames[2].rc(), b"ACTGAGA".as_slice());
    }

    #[test]
    fn test_seqtohashes_dna() {
        let sequence = b"AGTCGT";
        let hash_function = HashFunctions::Murmur64Dna;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();

        // Expected k-mers from the forward and reverse complement sequence
        let expected_kmers = vec![
            (b"AGT".to_vec(), b"ACT".to_vec()),
            (b"GTC".to_vec(), b"GAC".to_vec()),
            (b"TCG".to_vec(), b"CGA".to_vec()),
            (b"CGT".to_vec(), b"ACG".to_vec()),
        ];

        // Compute expected hashes from expected kmers
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .map(|(fw_kmer, rc_kmer)| crate::_hash_murmur(std::cmp::min(fw_kmer, rc_kmer), seed))
            .collect();

        // Collect hashes from SeqToHashes
        let sth_hashes: Vec<u64> = sth.map(|result| result.unwrap()).collect();
        eprintln!("SeqToHashes hashes: {:?}", sth_hashes);

        // Check that SeqToHashes matches expected hashes in order
        assert_eq!(
            sth_hashes, expected_hashes,
            "Hashes do not match in order for SeqToHashes"
        );
    }

    #[test]
    fn test_seqtohashes_dna_2() {
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
        )
        .unwrap();

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

    #[test]
    fn test_seqtohashes_is_protein() {
        let sequence = b"MVLSPADKTNVKAAW";
        let hash_function = HashFunctions::Murmur64Protein;
        let k_size = 3;
        let seed = 42;
        let force = false;
        let is_protein = true;

        let sth =
            SeqToHashes::new(sequence, k_size * 3, force, is_protein, hash_function, seed).unwrap();

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

        // Collect hashes from SeqToHashes
        let sth_hashes: Vec<u64> = sth.map(|result| result.unwrap()).collect();
        eprintln!("SeqToHashes hashes: {:?}", sth_hashes);

        // Check that SeqToHashes matches expected hashes in order
        assert_eq!(sth_hashes, expected_hashes, "Hashes do not match in order");
    }

    #[test]
    fn test_seqtohashes_translate() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Protein;
        let k_size = 9; // needs to be *3 for protein
        let seed = 42;
        let force = false;
        let is_protein = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();

        let expected_kmers = vec![
            b"SRR".as_slice(),
            b"RRA".as_slice(),
            b"SST".as_slice(),
            b"STT".as_slice(),
            b"VVE".as_slice(),
            b"ARR".as_slice(),
            b"SSS".as_slice(),
            b"LDD".as_slice(),
        ];

        // Compute expected hashes
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .map(|fw_kmer| crate::_hash_murmur(fw_kmer, seed))
            .collect();

        // Collect hashes from SeqToHashes
        let sth_hashes: Vec<u64> = sth.map(|result| result.unwrap()).collect();
        eprintln!("SeqToHashes hashes: {:?}", sth_hashes);

        // Check that SeqToHashes matches expected hashes in order
        assert_eq!(
            sth_hashes, expected_hashes,
            "Hashes do not match in order for SeqToHashes"
        );
    }

    #[test]
    fn test_seqtohashes_skipm1n3() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Skipm1n3;
        let k_size = 3;
        let is_protein = false;
        let seed = 42;
        let force = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        // Expected k-mers for skipmer (m=1, n=3) across all frames
        let expected_kmers = vec![
            (b"ACC".as_slice(), b"GGT".as_slice()),
            (b"CCG".as_slice(), b"CGG".as_slice()),
            (b"GGG".as_slice(), b"CCC".as_slice()),
            (b"GGC".as_slice(), b"GCC".as_slice()),
            (b"TTA".as_slice(), b"TAA".as_slice()),
            (b"TAT".as_slice(), b"ATA".as_slice()),
        ];

        // Compute expected hashes
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .map(|(fw_kmer, rc_kmer)| crate::_hash_murmur(std::cmp::min(fw_kmer, rc_kmer), seed))
            .collect();

        // Collect hashes from SeqToHashes
        let sth_hashes: Vec<u64> = sth.map(|result| result.unwrap()).collect();
        eprintln!("SeqToHashes hashes: {:?}", sth_hashes);

        // Check that SeqToHashes matches expected hashes in order
        assert_eq!(
            sth_hashes, expected_hashes,
            "Hashes do not match in order for SeqToHashes"
        );
    }

    #[test]
    fn test_seq2hashes_skipm2n3() {
        let sequence = b"AGTCGTCGAGCT";
        let hash_function = HashFunctions::Murmur64Skipm2n3;
        let k_size = 7;
        let is_protein = false;
        let seed = 42;
        let force = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();

        // Expected k-mers for skipmer (m=2, n=3)
        let expected_kmers = vec![
            (b"AGCGCGG".as_slice(), b"CCGCGCT".as_slice()),
            (b"GCGCGGC".as_slice(), b"GCCGCGC".as_slice()),
            (b"GTGTGAC".as_slice(), b"GTCACAC".as_slice()),
            (b"TGTGACT".as_slice(), b"AGTCACA".as_slice()),
            (b"TCTCAGT".as_slice(), b"ACTGAGA".as_slice()),
        ];

        // Compute expected hashes
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .map(|(fw_kmer, rc_kmer)| crate::_hash_murmur(std::cmp::min(fw_kmer, rc_kmer), seed))
            .collect();

        // Collect hashes from SeqToHashes
        let sth_hashes: Vec<u64> = sth.map(|result| result.unwrap()).collect();
        eprintln!("SeqToHashes hashes: {:?}", sth_hashes);

        // Check that SeqToHashes matches expected hashes in order
        assert_eq!(
            sth_hashes, expected_hashes,
            "Hashes do not match in order for SeqToHashes"
        );
    }

    #[test]
    fn test_seqtohashes_skipm2n3_2() {
        let sequence = b"AGTCGTCA";
        let hash_function = HashFunctions::Murmur64Skipm2n3;
        let k_size = 5;
        let seed = 42;
        let force = true;
        let is_protein = false;

        let sth =
            SeqToHashes::new(sequence, k_size, force, is_protein, hash_function, seed).unwrap();
        let frames = sth.frames.clone();
        for fr in frames {
            eprintln!("{}", fr);
        }

        let expected_kmers = vec![
            (b"AGCGC".as_slice(), b"GCGCT".as_slice()),
            (b"GCGCA".as_slice(), b"TGCGC".as_slice()),
            (b"GTGTA".as_slice(), b"TACAC".as_slice()),
        ];

        // Compute expected hashes
        let expected_hashes: Vec<u64> = expected_kmers
            .iter()
            .map(|(fw_kmer, rc_kmer)| crate::_hash_murmur(std::cmp::min(fw_kmer, rc_kmer), seed))
            .collect();

        // Collect hashes from SeqToHashes
        let sth_hashes: Vec<u64> = sth.map(|result| result.unwrap()).collect();
        eprintln!("SeqToHashes hashes: {:?}", sth_hashes);

        // Check that SeqToHashes matches expected hashes in order
        assert_eq!(
            sth_hashes, expected_hashes,
            "Hashes do not match in order for SeqToHashes"
        );
    }
}
