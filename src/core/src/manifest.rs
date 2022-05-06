use std::convert::TryInto;
use std::io::Read;
use std::ops::Deref;
use std::path::PathBuf;

use serde::de;
use serde::{Deserialize, Serialize};

use crate::encodings::HashFunctions;
use crate::index::Selection;
use crate::Error;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Record {
    internal_location: String,
    ksize: u32,

    #[serde(deserialize_with = "to_bool")]
    with_abundance: bool,

    md5: String,
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

fn to_bool<'de, D>(deserializer: D) -> Result<bool, D::Error>
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
    pub fn internal_location(&self) -> PathBuf {
        self.internal_location.clone().into()
    }

    pub fn ksize(&self) -> u32 {
        self.ksize
    }

    pub fn with_abundance(&self) -> bool {
        self.with_abundance
    }

    pub fn md5(&self) -> &str {
        self.md5.as_ref()
    }

    pub fn name(&self) -> &str {
        self.name.as_ref()
    }

    pub fn moltype(&self) -> HashFunctions {
        self.moltype.as_str().try_into().unwrap()
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

    pub fn iter(&self) -> impl Iterator<Item = &Record> {
        self.records.iter()
    }

    pub fn select_to_manifest(&self, selection: &Selection) -> Result<Self, Error> {
        let rows = self.records.iter().filter(|row| {
            let mut valid = true;
            valid = if let Some(ksize) = selection.ksize() {
                row.ksize == ksize
            } else {
                valid
            };
            valid = if let Some(abund) = selection.abund() {
                valid && row.with_abundance() == abund
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

impl From<&[PathBuf]> for Manifest {
    fn from(v: &[PathBuf]) -> Self {
        Manifest {
            records: v
                .iter()
                .map(|p| Record {
                    internal_location: p.to_str().unwrap().into(),
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
