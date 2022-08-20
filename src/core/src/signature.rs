//! # Compressed representations of genomic data
//!
//! A signature is a collection of sketches for a genomic dataset.

use std::fs::File;
use std::io;
use std::iter::Iterator;
use std::path::Path;
use std::str;

use cfg_if::cfg_if;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::encodings::{aa_to_dayhoff, aa_to_hp, revcomp, to_aa, HashFunctions, VALID};
use crate::prelude::*;
use crate::sketch::Sketch;
use crate::Error;
use crate::HashIntoType;

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
            max_index: _max_index as usize,
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
                    self.dna_ksize = self.k_size as usize;
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

                        aa.windows(self.k_size as usize).for_each(|n| {
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

                        aa_rc.windows(self.k_size as usize).for_each(|n| {
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
                        self.aa_seq = match self.hash_function {
                            HashFunctions::murmur64_dayhoff => {
                                self.sequence.iter().cloned().map(aa_to_dayhoff).collect()
                            }
                            HashFunctions::murmur64_hp => {
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
                                if k != mh.ksize() as usize {
                                    return false;
                                }
                            };

                            match moltype {
                                Some(x) => {
                                    if mh.hash_function() == x {
                                        return true;
                                    }
                                }
                                None => return true, // TODO: match previous behavior
                            };
                        }
                        Sketch::LargeMinHash(mh) => {
                            if let Some(k) = ksize {
                                if k != mh.ksize() as usize {
                                    return false;
                                }
                            };

                            match moltype {
                                Some(x) => {
                                    if mh.hash_function() == x {
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
    use std::convert::TryInto;
    use std::fs::File;
    use std::io::{BufReader, Read};
    use std::path::PathBuf;

    use needletail::parse_fastx_reader;

    use crate::cmd::ComputeParameters;
    use crate::signature::SigsTrait;

    use super::Signature;

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

        let sketches = sig.sketches();
        assert_eq!(sketches.len(), 12);
        for sk in sketches {
            assert_eq!(sk.size(), 500);
        }
    }
}
