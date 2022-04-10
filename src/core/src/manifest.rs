use std::io::Read;
use std::ops::Deref;
use std::path::PathBuf;

use serde::{Deserialize, Serialize};

use crate::Error;

#[derive(Debug, Serialize, Deserialize)]
pub struct Record {
    internal_location: String,
    /*
    md5: String,
    md5short: String,
    ksize: String,
    moltype: String,
    num: String,
    scaled: String, n_hashes: String,
    with_abundance: String,
    name: String,
    filename: String,
    */
}

#[derive(Debug, Default, Serialize, Deserialize)]
pub struct Manifest {
    records: Vec<Record>,
}

impl Record {
    pub fn internal_location(&self) -> PathBuf {
        self.internal_location.clone().into()
    }
}

impl Manifest {
    pub fn from_reader<R: Read>(rdr: R) -> Result<Self, Error> {
        let mut records = vec![];

        let mut rdr = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .from_reader(rdr);
        for result in rdr.deserialize() {
            let record: Record = result?;
            records.push(record);
        }
        Ok(Manifest { records })
    }

    pub fn internal_locations(&self) -> impl Iterator<Item = &str> {
        self.records.iter().map(|r| r.internal_location.as_str())
    }
}

impl From<&[PathBuf]> for Manifest {
    fn from(v: &[PathBuf]) -> Self {
        Manifest {
            records: v
                .into_iter()
                .map(|p| Record {
                    internal_location: p.to_str().unwrap().into(),
                })
                .collect(),
        }
    }
}

impl Deref for Manifest {
    type Target = Vec<Record>;

    fn deref(&self) -> &Self::Target {
        &self.records
    }
}
