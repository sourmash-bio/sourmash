use std::io::Read;

use serde::Deserialize;

use crate::Error;

#[derive(Debug, Deserialize)]
struct Record {
    internal_location: String,
    /*
    md5: String,
    md5short: String,
    ksize: String,
    moltype: String,
    num: String,
    scaled: String,
    n_hashes: String,
    with_abundance: String,
    name: String,
    filename: String,
    */
}

#[derive(Debug)]
pub struct Manifest {
    records: Vec<Record>,
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
