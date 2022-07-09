use std::convert::TryInto;
use std::io::{Read, Write};
use std::ops::Deref;

use camino::Utf8PathBuf as PathBuf;
use getset::{CopyGetters, Getters, Setters};
use serde::de;
use serde::{Deserialize, Serialize};

use crate::encodings::HashFunctions;
use crate::index::Selection;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::Sketch;
use crate::Result;

#[derive(Debug, Serialize, Deserialize, Clone, CopyGetters, Getters, Setters)]
pub struct Record {
    #[getset(get = "pub", set = "pub")]
    internal_location: PathBuf,

    #[getset(get = "pub", set = "pub")]
    ksize: u32,

    #[getset(get = "pub", set = "pub")]
    #[serde(deserialize_with = "to_bool")]
    with_abundance: bool,

    #[getset(get = "pub", set = "pub")]
    md5: String,

    #[getset(get = "pub", set = "pub")]
    name: String,

    moltype: String,
    /*
    md5short: String,
    num: String,
    scaled: String,
    n_hashes: String,
    filename: String,
    */
}

fn to_bool<'de, D>(deserializer: D) -> std::result::Result<bool, D::Error>
where
    D: de::Deserializer<'de>,
{
    match String::deserialize(deserializer)?
        .to_ascii_lowercase()
        .as_ref()
    {
        "0" | "false" => Ok(false),
        "1" | "true" => Ok(true),
        other => Err(de::Error::invalid_value(
            de::Unexpected::Str(other),
            &"0/1 or true/false are the only supported values",
        )),
    }
}

#[derive(Debug, Default, Serialize, Deserialize, Clone)]
pub struct Manifest {
    records: Vec<Record>,
}

impl Record {
    pub fn from_sig(sig: &Signature, path: &str) -> Vec<Self> {
        sig.iter()
            .map(|sketch| {
                let (ksize, md5, with_abundance, moltype) = match sketch {
                    Sketch::MinHash(mh) => (
                        mh.ksize() as u32,
                        mh.md5sum(),
                        mh.track_abundance(),
                        mh.hash_function(),
                    ),
                    Sketch::LargeMinHash(mh) => (
                        mh.ksize() as u32,
                        mh.md5sum(),
                        mh.track_abundance(),
                        mh.hash_function(),
                    ),
                    _ => unimplemented!(),
                };

                Self {
                    internal_location: path.into(),
                    moltype: moltype.to_string(),
                    name: sig.name(),
                    ksize,
                    md5,
                    with_abundance,
                }
            })
            .collect()
    }

    pub fn moltype(&self) -> HashFunctions {
        self.moltype.as_str().try_into().unwrap()
    }

    pub fn check_compatible(&self, other: &Record) -> Result<()> {
        /*
        if self.num != other.num {
            return Err(Error::MismatchNum {
                n1: self.num,
                n2: other.num,
            }
            .into());
        }
        */
        use crate::Error;

        if self.ksize() != other.ksize() {
            return Err(Error::MismatchKSizes);
        }
        if self.moltype() != other.moltype() {
            // TODO: fix this error
            return Err(Error::MismatchDNAProt);
        }
        /*
        if self.scaled() < other.scaled() {
            return Err(Error::MismatchScaled);
        }
        if self.seed() != other.seed() {
            return Err(Error::MismatchSeed);
        }
        */
        Ok(())
    }
}

impl Manifest {
    pub fn from_reader<R: Read>(rdr: R) -> Result<Self> {
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

    pub fn to_writer<W: Write>(&self, mut wtr: W) -> Result<()> {
        wtr.write_all(b"# SOURMASH-MANIFEST-VERSION: 1.0\n")?;

        let mut wtr = csv::Writer::from_writer(wtr);

        for record in &self.records {
            wtr.serialize(record)?;
        }

        Ok(())
    }

    pub fn internal_locations(&self) -> impl Iterator<Item = &str> {
        self.records.iter().map(|r| r.internal_location.as_str())
    }

    pub fn iter(&self) -> impl Iterator<Item = &Record> {
        self.records.iter()
    }

    pub fn select_to_manifest(&self, selection: &Selection) -> Result<Self> {
        let rows = self.records.iter().filter(|row| {
            let mut valid = true;
            valid = if let Some(ksize) = selection.ksize() {
                row.ksize == ksize
            } else {
                valid
            };
            valid = if let Some(abund) = selection.abund() {
                valid && *row.with_abundance() == abund
            } else {
                valid
            };
            valid = if let Some(moltype) = selection.moltype() {
                valid && row.moltype() == moltype
            } else {
                valid
            };
            valid
        });

        Ok(Manifest {
            records: rows.cloned().collect(),
        })

        /*
        matching_rows = self.rows
        if ksize:
            matching_rows = ( row for row in matching_rows
                              if row['ksize'] == ksize )
        if moltype:
            matching_rows = ( row for row in matching_rows
                              if row['moltype'] == moltype )
        if scaled or containment:
            if containment and not scaled:
                raise ValueError("'containment' requires 'scaled' in Index.select'")

            matching_rows = ( row for row in matching_rows
                              if row['scaled'] and not row['num'] )
        if num:
            matching_rows = ( row for row in matching_rows
                              if row['num'] and not row['scaled'] )

        if abund:
            # only need to concern ourselves if abundance is _required_
            matching_rows = ( row for row in matching_rows
                              if row['with_abundance'] )

        if picklist:
            matching_rows = ( row for row in matching_rows
                              if picklist.matches_manifest_row(row) )

        # return only the internal filenames!
        for row in matching_rows:
            yield row
        */
    }
}

impl From<Vec<Record>> for Manifest {
    fn from(records: Vec<Record>) -> Self {
        Manifest { records }
    }
}

impl From<&[PathBuf]> for Manifest {
    fn from(v: &[PathBuf]) -> Self {
        Manifest {
            records: v
                .iter()
                .map(|p| Record {
                    internal_location: p.clone(),
                    ksize: 0,              // FIXME
                    with_abundance: false, // FIXME
                    md5: "".into(),        // FIXME
                    name: "".into(),       // FIXME
                    moltype: "".into(),    // FIXME
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
