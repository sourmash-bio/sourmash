//! # Compressed representations of genomic data
//!
//! A signature is a collection of sketches for a genomic dataset.

use std::fs::File;
use std::io;
use std::iter::Iterator;
use std::path::Path;
use std::str;

use failure::Error;
use serde_derive::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

use crate::errors::SourmashError;
use crate::index::storage::ToWriter;
use crate::sketch::minhash::HashFunctions;
use crate::sketch::Sketch;

pub trait SigsTrait {
    fn size(&self) -> usize;
    fn to_vec(&self) -> Vec<u64>;
    fn check_compatible(&self, other: &Self) -> Result<(), Error>;
    fn add_sequence(&mut self, seq: &[u8], _force: bool) -> Result<(), Error>;
    fn ksize(&self) -> usize;
}

impl SigsTrait for Sketch {
    fn size(&self) -> usize {
        match *self {
            Sketch::UKHS(ref ukhs) => ukhs.size(),
            Sketch::MinHash(ref mh) => mh.size(),
        }
    }

    fn to_vec(&self) -> Vec<u64> {
        match *self {
            Sketch::UKHS(ref ukhs) => ukhs.to_vec(),
            Sketch::MinHash(ref mh) => mh.to_vec(),
        }
    }

    fn ksize(&self) -> usize {
        match *self {
            Sketch::UKHS(ref ukhs) => ukhs.ksize(),
            Sketch::MinHash(ref mh) => mh.ksize(),
        }
    }

    fn check_compatible(&self, other: &Self) -> Result<(), Error> {
        match *self {
            Sketch::UKHS(ref ukhs) => match other {
                Sketch::UKHS(ref ot) => ukhs.check_compatible(ot),
                _ => Err(SourmashError::MismatchSignatureType.into()),
            },
            Sketch::MinHash(ref mh) => match other {
                Sketch::MinHash(ref ot) => mh.check_compatible(ot),
                _ => Err(SourmashError::MismatchSignatureType.into()),
            },
        }
    }

    fn add_sequence(&mut self, seq: &[u8], force: bool) -> Result<(), Error> {
        match *self {
            Sketch::UKHS(ref mut ukhs) => ukhs.add_sequence(seq, force),
            Sketch::MinHash(ref mut mh) => mh.add_sequence(seq, force),
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, TypedBuilder)]
pub struct Signature {
    #[serde(default = "default_class")]
    #[builder(default_code = "default_class()")]
    pub class: String,

    #[serde(default)]
    #[builder(default)]
    pub email: String,

    pub hash_function: String,

    #[builder(default)]
    pub filename: Option<String>,

    #[serde(skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,

    #[serde(default = "default_license")]
    #[builder(default_code = "default_license()")]
    pub license: String,

    pub signatures: Vec<Sketch>,

    #[serde(default = "default_version")]
    #[builder(default_code = "default_version()")]
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
    pub fn name(&self) -> String {
        if let Some(name) = &self.name {
            name.clone()
        } else if let Some(filename) = &self.filename {
            filename.clone()
        } else {
            // TODO md5sum case
            unimplemented!()
        }
    }

    pub fn filename(&self) -> String {
        if let Some(filename) = &self.filename {
            filename.clone()
        } else {
            unimplemented!()
        }
    }

    pub fn md5sum(&self) -> String {
        if self.signatures.len() == 1 {
            match &self.signatures[0] {
                Sketch::MinHash(mh) => mh.md5sum(),
                Sketch::UKHS(hs) => hs.md5sum(),
            }
        } else {
            // TODO: select the correct signature
            unimplemented!()
        }
    }

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
                        Sketch::UKHS(hs) => {
                            if let Some(k) = ksize {
                                if k != hs.ksize() as usize {
                                    return false;
                                }
                            };

                            match moltype {
                                Some(x) => {
                                    if x == HashFunctions::murmur64_DNA {
                                        return true;
                                    } else {
                                        // TODO: draff only supports dna for now
                                        unimplemented!()
                                    }
                                }
                                None => unimplemented!(),
                            };
                        }
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

impl ToWriter for Signature {
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: io::Write,
    {
        match serde_json::to_writer(writer, &vec![&self]) {
            Ok(_) => Ok(()),
            Err(_) => Err(SourmashError::SerdeError.into()),
        }
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
    use std::io::BufReader;
    use std::path::PathBuf;

    use super::Signature;

    #[test]
    fn load_sig() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/.sbt.v3/60f7e23c24a8d94791cc7a8680c493f9");

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
}
