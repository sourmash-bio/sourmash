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

// Iterator for converting sequence to hashes
pub struct SeqToHashes {
    sequence: Vec<u8>,
    kmer_index: usize,
    k_size: usize,
    max_index: usize,
    force: bool,
    is_protein: bool,
    hash_function: HashFunctions,
    seed: u64,
    hashes_buffer: Vec<u64>,

    dna_configured: bool,
    dna_rc: Vec<u8>,
    dna_ksize: usize,
    dna_len: usize,
    dna_last_position_check: usize,

    prot_configured: bool,
    aa_seq: Vec<u8>,
    translate_iter_step: usize,
}

impl SeqToHashes {
    pub fn new(
        seq: &[u8],
        k_size: usize,
        force: bool,
        is_protein: bool,
        hash_function: HashFunctions,
        seed: u64,
    ) -> SeqToHashes {
        let mut ksize: usize = k_size;

        // Divide the kmer size by 3 if protein
        if is_protein || !hash_function.dna() {
            ksize = k_size / 3;
        }

        // By setting _max_index to 0, the iterator will return None and exit
        let _max_index = if seq.len() >= ksize {
            seq.len() - ksize + 1
        } else {
            0
        };

        SeqToHashes {
            // Here we convert the sequence to upper case
            sequence: seq.to_ascii_uppercase(),
            k_size: ksize,
            kmer_index: 0,
            max_index: _max_index,
            force,
            is_protein,
            hash_function,
            seed,
            hashes_buffer: Vec::with_capacity(1000),
            dna_configured: false,
            dna_rc: Vec::with_capacity(1000),
            dna_ksize: 0,
            dna_len: 0,
            dna_last_position_check: 0,
            prot_configured: false,
            aa_seq: Vec::new(),
            translate_iter_step: 0,
        }
    }
}

/*
Iterator that return a kmer hash for all modes except translate.
In translate mode:
    - all the frames are processed at once and converted to hashes.
    - all the hashes are stored in `hashes_buffer`
    - after processing all the kmers, `translate_iter_step` is incremented
      per iteration to iterate over all the indeces of the `hashes_buffer`.
    - the iterator will die once `translate_iter_step` == length(hashes_buffer)
More info https://github.com/sourmash-bio/sourmash/pull/1946
*/

impl Iterator for SeqToHashes {
    type Item = Result<u64, Error>;

    fn next(&mut self) -> Option<Self::Item> {
        if (self.kmer_index < self.max_index) || !self.hashes_buffer.is_empty() {
            // Processing DNA or Translated DNA
            if !self.is_protein {
                // Setting the parameters only in the first iteration
                if !self.dna_configured {
                    self.dna_ksize = self.k_size;
                    self.dna_len = self.sequence.len();
                    if self.dna_len < self.dna_ksize
                        || (!self.hash_function.dna() && self.dna_len < self.k_size * 3)
                    {
                        return None;
                    }
                    // pre-calculate the reverse complement for the full sequence...
                    self.dna_rc = revcomp(&self.sequence);
                    self.dna_configured = true;
                }

                // Processing DNA
                if self.hash_function.dna() {
                    let kmer = &self.sequence[self.kmer_index..self.kmer_index + self.dna_ksize];

                    for j in std::cmp::max(self.kmer_index, self.dna_last_position_check)
                        ..self.kmer_index + self.dna_ksize
                    {
                        if !VALID[self.sequence[j] as usize] {
                            if !self.force {
                                return Some(Err(Error::InvalidDNA {
                                    message: String::from_utf8(kmer.to_vec()).unwrap(),
                                }));
                            } else {
                                self.kmer_index += 1;
                                // Move the iterator to the next step
                                return Some(Ok(0));
                            }
                        }
                        self.dna_last_position_check += 1;
                    }

                    // ... and then while moving the k-mer window forward for the sequence
                    // we move another window backwards for the RC.
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

                    let krc = &self.dna_rc[self.dna_len - self.dna_ksize - self.kmer_index
                        ..self.dna_len - self.kmer_index];
                    let hash = crate::_hash_murmur(std::cmp::min(kmer, krc), self.seed);
                    self.kmer_index += 1;
                    Some(Ok(hash))
                } else if self.hashes_buffer.is_empty() && self.translate_iter_step == 0 {
                    // Processing protein by translating DNA
                    // TODO: Implement iterator over frames instead of hashes_buffer.

                    for frame_number in 0..3 {
                        let substr: Vec<u8> = self
                            .sequence
                            .iter()
                            .cloned()
                            .skip(frame_number)
                            .take(self.sequence.len() - frame_number)
                            .collect();

                        let aa = to_aa(
                            &substr,
                            self.hash_function.dayhoff(),
                            self.hash_function.hp(),
                        )
                        .unwrap();

                        aa.windows(self.k_size).for_each(|n| {
                            let hash = crate::_hash_murmur(n, self.seed);
                            self.hashes_buffer.push(hash);
                        });

                        let rc_substr: Vec<u8> = self
                            .dna_rc
                            .iter()
                            .cloned()
                            .skip(frame_number)
                            .take(self.dna_rc.len() - frame_number)
                            .collect();
                        let aa_rc = to_aa(
                            &rc_substr,
                            self.hash_function.dayhoff(),
                            self.hash_function.hp(),
                        )
                        .unwrap();

                        aa_rc.windows(self.k_size).for_each(|n| {
                            let hash = crate::_hash_murmur(n, self.seed);
                            self.hashes_buffer.push(hash);
                        });
                    }
                    Some(Ok(0))
                } else {
                    if self.translate_iter_step == self.hashes_buffer.len() {
                        self.hashes_buffer.clear();
                        self.kmer_index = self.max_index;
                        return Some(Ok(0));
                    }
                    let curr_idx = self.translate_iter_step;
                    self.translate_iter_step += 1;
                    Some(Ok(self.hashes_buffer[curr_idx]))
                }
            } else {
                // Processing protein
                // The kmer size is already divided by 3

                if self.hash_function.protein() {
                    let aa_kmer = &self.sequence[self.kmer_index..self.kmer_index + self.k_size];
                    let hash = crate::_hash_murmur(aa_kmer, self.seed);
                    self.kmer_index += 1;
                    Some(Ok(hash))
                } else {
                    if !self.prot_configured {
                        self.aa_seq = match &self.hash_function {
                            HashFunctions::Murmur64Dayhoff => {
                                self.sequence.iter().cloned().map(aa_to_dayhoff).collect()
                            }
                            HashFunctions::Murmur64Hp => {
                                self.sequence.iter().cloned().map(aa_to_hp).collect()
                            }
                            invalid => {
                                return Some(Err(Error::InvalidHashFunction {
                                    function: format!("{}", invalid),
                                }));
                            }
                        };
                    }

                    let aa_kmer = &self.aa_seq[self.kmer_index..self.kmer_index + self.k_size];
                    let hash = crate::_hash_murmur(aa_kmer, self.seed);
                    self.kmer_index += 1;
                    Some(Ok(hash))
                }
            }
        } else {
            // End the iterator
            None
        }
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
                    Sketch::MinHash(mh) => valid && mh.scaled() <= sel_scaled as u64,
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
                    if (mh.scaled() as u32) < sel_scaled {
                        *sketch = Sketch::MinHash(mh.downsample_scaled(sel_scaled as u64)?);
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

#[cfg(test)]
mod test {
    use std::fs::File;
    use std::io::{BufReader, Read};
    use std::path::PathBuf;

    use needletail::parse_fastx_reader;

    use crate::cmd::ComputeParameters;
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
}
