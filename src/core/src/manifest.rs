use std::convert::TryInto;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::ops::Deref;

use camino::Utf8PathBuf as PathBuf;
use getset::{CopyGetters, Getters, Setters};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::de;
use serde::{Deserialize, Serialize};

use crate::encodings::HashFunctions;
use crate::prelude::*;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::Sketch;
use crate::Result;

#[derive(Debug, Serialize, Deserialize, Clone, CopyGetters, Getters, Setters, PartialEq)]
pub struct Record {
    #[getset(get = "pub", set = "pub")]
    internal_location: PathBuf,

    #[getset(get = "pub", set = "pub")]
    md5: String,

    md5short: String,

    #[getset(get = "pub", set = "pub")]
    ksize: u32,

    moltype: String,

    num: u32,
    scaled: u64,
    n_hashes: usize,

    #[getset(get = "pub", set = "pub")]
    #[serde(serialize_with = "intbool", deserialize_with = "to_bool")]
    with_abundance: bool,

    #[getset(get = "pub", set = "pub")]
    name: String,

    filename: String,
}

fn intbool<S>(x: &bool, s: S) -> std::result::Result<S::Ok, S::Error>
where
    S: serde::Serializer,
{
    if *x {
        s.serialize_i32(1)
    } else {
        s.serialize_i32(0)
    }
}

fn to_bool<'de, D>(deserializer: D) -> std::result::Result<bool, D::Error>
where
    D: de::Deserializer<'de>,
{
    match String::deserialize(deserializer)?
        .to_ascii_lowercase()
        .as_ref()
    {
        "0" | "false" | "False" => Ok(false),
        "1" | "true" | "True" => Ok(true),
        other => Err(de::Error::invalid_value(
            de::Unexpected::Str(other),
            &"0/1, true/false, True/False are the only supported values",
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
                let (ksize, md5, with_abundance, moltype, n_hashes, num, scaled) = match sketch {
                    Sketch::MinHash(mh) => (
                        mh.ksize() as u32,
                        mh.md5sum(),
                        mh.track_abundance(),
                        mh.hash_function(),
                        mh.size(),
                        mh.num(),
                        mh.scaled(),
                    ),
                    Sketch::LargeMinHash(mh) => (
                        mh.ksize() as u32,
                        mh.md5sum(),
                        mh.track_abundance(),
                        mh.hash_function(),
                        mh.size(),
                        mh.num(),
                        mh.scaled(),
                    ),
                    _ => unimplemented!(),
                };

                let md5short = md5[0..8].into();

                Self {
                    internal_location: path.into(),
                    moltype: moltype.to_string(),
                    name: sig.name(),
                    ksize,
                    md5,
                    md5short,
                    with_abundance,
                    filename: sig.filename(),
                    n_hashes,
                    num,
                    scaled,
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
}

impl Select for Manifest {
    fn select(self, selection: &Selection) -> Result<Self> {
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
            valid = if let Some(scaled) = selection.scaled() {
                // num sigs have row.scaled = 0, don't include them
                valid && row.scaled != 0 && row.scaled <= scaled as u64
            } else {
                valid
            };
            valid = if let Some(num) = selection.num() {
                valid && row.num == num
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
    fn from(paths: &[PathBuf]) -> Self {
        #[cfg(feature = "parallel")]
        let iter = paths.par_iter();

        #[cfg(not(feature = "parallel"))]
        let iter = paths.iter();

        let records: Vec<Record> = iter
            .flat_map(|p| {
                let recs: Vec<Record> = Signature::from_path(p)
                    .unwrap_or_else(|_| panic!("Error processing {:?}", p))
                    .into_iter()
                    .flat_map(|v| Record::from_sig(&v, p.as_str()))
                    .collect();
                recs
            })
            .collect();

        Manifest { records }
    }
}

impl From<&PathBuf> for Manifest {
    fn from(pathlist: &PathBuf) -> Self {
        let file = File::open(pathlist).unwrap_or_else(|_| panic!("Failed to open {:?}", pathlist));
        let reader = BufReader::new(file);

        let paths: Vec<PathBuf> = reader
            .lines()
            .map(|line| line.unwrap_or_else(|_| panic!("Failed to read line from {:?}", pathlist)))
            .map(PathBuf::from)
            .collect();

        paths.as_slice().into()
    }
}

impl Deref for Manifest {
    type Target = Vec<Record>;

    fn deref(&self) -> &Self::Target {
        &self.records
    }
}

#[cfg(test)]
mod test {
    use camino::Utf8PathBuf as PathBuf;
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

    use super::Manifest;

    #[test]
    fn manifest_from_pathlist() {
        let temp_dir = TempDir::new().unwrap();
        let utf8_output = PathBuf::from_path_buf(temp_dir.path().to_path_buf())
            .expect("Path should be valid UTF-8");
        let mut filename = utf8_output.join("sig-pathlist.txt");
        //convert to camino utf8pathbuf
        filename = PathBuf::from(filename);
        // build sig filenames
        let base_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let test_sigs = vec![
            "../../tests/test-data/47.fa.sig",
            "../../tests/test-data/63.fa.sig",
        ];

        let full_paths: Vec<_> = test_sigs
            .into_iter()
            .map(|sig| base_path.join(sig))
            .collect();

        // write a file in test directory with a filename on each line
        let mut pathfile = File::create(&filename).unwrap();
        for sigfile in &full_paths {
            writeln!(pathfile, "{}", sigfile).unwrap();
        }

        // load into manifest
        let manifest = Manifest::from(&filename);
        assert_eq!(manifest.len(), 2);
    }

    #[test]
    #[should_panic(expected = "Failed to open \"no-exist\"")]
    fn manifest_from_pathlist_nonexistent_file() {
        let filename = PathBuf::from("no-exist");
        let _manifest = Manifest::from(&filename);
    }

    #[test]
    #[should_panic]
    fn manifest_from_pathlist_badfile() {
        let temp_dir = TempDir::new().unwrap();
        let utf8_output = PathBuf::from_path_buf(temp_dir.path().to_path_buf())
            .expect("Path should be valid UTF-8");
        let mut filename = utf8_output.join("sig-pathlist.txt");
        //convert to camino utf8pathbuf
        filename = PathBuf::from(filename);

        let mut pathfile = File::create(&filename).unwrap();
        write!(pathfile, "Valid line\n").unwrap();
        pathfile.write_all(&[0xED, 0xA0, 0x80]).unwrap(); // invalid UTF-8

        // load into manifest
        let _manifest = Manifest::from(&filename);
    }

    #[test]
    #[should_panic]
    fn manifest_from_paths_badpath() {
        let base_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let test_sigs = vec![
            PathBuf::from("no-exist"),
            PathBuf::from("../../tests/test-data/63.fa.sig"),
        ];

        let full_paths: Vec<PathBuf> = test_sigs
            .into_iter()
            .map(|sig| base_path.join(sig))
            .collect();

        // load into manifest
        let _manifest = Manifest::from(&full_paths[..]); // pass full_paths as a slice
    }

    #[test]
    fn test_manifest_to_writer_bools() {
        let base_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

        let test_sigs = vec![
            PathBuf::from("../../tests/test-data/47.fa.sig"),
            PathBuf::from("../../tests/test-data/track_abund/63.fa.sig"),
        ];

        let full_paths: Vec<PathBuf> = test_sigs
            .into_iter()
            .map(|sig| base_path.join(sig))
            .collect();

        let manifest = Manifest::from(&full_paths[..]); // pass full_paths as a slice

        let temp_dir = TempDir::new().unwrap();
        let utf8_output = PathBuf::from_path_buf(temp_dir.path().to_path_buf())
            .expect("Path should be valid UTF-8");

        let filename = utf8_output.join("sigs.manifest.csv");
        let mut wtr = File::create(&filename).expect("Failed to create file");

        manifest.to_writer(&mut wtr).unwrap();

        // check that we can reopen the file as a manifest + properly check abund
        let infile = File::open(&filename).expect("Failed to open file");
        let m2 = Manifest::from_reader(&infile).unwrap();
        for record in m2.iter() {
            eprintln!("{:?}", record.name());
            if record.name().contains("OS185") {
                assert_eq!(record.with_abundance(), &false)
            } else {
                assert_eq!(record.with_abundance(), &true)
            }
        }
    }
}
