use std::collections::HashMap;
use std::convert::TryInto;
use std::fmt::Write as FmtWrite;
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::ops::Deref;

use camino::Utf8PathBuf as PathBuf;
use getset::{CopyGetters, Getters, Setters};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use serde::de;
use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize, Serializer};

use crate::encodings::HashFunctions;
use crate::prelude::*;
use crate::signature::SigsTrait;
use crate::sketch::Sketch;
use crate::Result;

#[derive(Debug, Serialize, Deserialize, Clone, CopyGetters, Getters, Setters, Eq)]
pub struct Record {
    #[getset(get = "pub", set = "pub")]
    internal_location: PathBuf,

    #[getset(get = "pub", set = "pub")]
    md5: String,

    md5short: String,

    #[getset(get_copy = "pub", set = "pub")]
    ksize: u32,

    moltype: String,

    num: u32,
    scaled: u64,
    n_hashes: usize,

    #[getset(get_copy = "pub", set = "pub")]
    #[serde(serialize_with = "intbool", deserialize_with = "to_bool")]
    with_abundance: bool,

    #[getset(get = "pub", set = "pub")]
    name: String,

    #[getset(get = "pub", set = "pub")]
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

impl PartialEq for Record {
    fn eq(&self, other: &Self) -> bool {
        self.ksize == other.ksize
            && self.moltype == other.moltype
            && self.num == other.num
            && self.scaled == other.scaled
            && self.with_abundance == other.with_abundance
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
                let (mut ksize, md5, with_abundance, moltype, n_hashes, num, scaled, hash_function) = match sketch
                {
                    Sketch::MinHash(mh) => (
                        mh.ksize() as u32,
                        mh.md5sum(),
                        mh.track_abundance(),
                        mh.hash_function(),
                        mh.size(),
                        mh.num(),
                        mh.scaled(),
                        mh.hash_function(),
                    ),
                    Sketch::LargeMinHash(mh) => (
                        mh.ksize() as u32,
                        mh.md5sum(),
                        mh.track_abundance(),
                        mh.hash_function(),
                        mh.size(),
                        mh.num(),
                        mh.scaled(),
                        mh.hash_function(),
                    ),
                    _ => unimplemented!(),
                };

                let md5short = md5[0..8].into();

                ksize = match hash_function {
                    HashFunctions::Murmur64Protein | HashFunctions::Murmur64Dayhoff | HashFunctions::Murmur64Hp => ksize / 3,
                    _ => ksize,
                };

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

    pub fn summarize(&self) -> Vec<RecordSummary> {
        let mut summary_map: HashMap<(u32, String, u32, u64, bool), RecordSummary> = HashMap::new();

        for record in self.iter() {
            let key = (
                record.ksize(),
                record.moltype.clone(),
                record.num,
                record.scaled,
                record.with_abundance(),
            );
            let entry = summary_map.entry(key).or_insert_with(|| RecordSummary {
                ksize: record.ksize(),
                moltype: record.moltype.clone(),
                num: record.num,
                scaled: record.scaled,
                with_abundance: record.with_abundance(),
                count: 0,
                total_n_hashes: 0,
            });

            entry.count += 1;
            entry.total_n_hashes += record.n_hashes;
        }

        summary_map.into_values().collect()
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
                valid && row.with_abundance() == abund
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

#[derive(Debug, Hash, PartialEq, Eq, Clone, Deserialize)]
pub struct RecordSummary {
    ksize: u32,
    moltype: String,
    num: u32,
    scaled: u64,
    with_abundance: bool,
    // Fields for counting and aggregation
    count: usize,
    total_n_hashes: usize,
}

impl Serialize for RecordSummary {
    fn serialize<S>(&self, serializer: S) -> std::result::Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("RecordSummary", 7)?;
        state.serialize_field("ksize", &self.ksize)?;
        state.serialize_field("moltype", &self.moltype)?;
        state.serialize_field("num", &self.num)?;
        state.serialize_field("scaled", &self.scaled)?;
        state.serialize_field("with_abundance", &self.with_abundance)?;
        state.serialize_field("count", &self.count)?;
        state.serialize_field("total_n_hashes", &self.total_n_hashes)?;
        state.end()
    }
}

impl std::fmt::Display for RecordSummary {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{} sketches with {}, k={}, scaled={}  {} total hashes",
            self.count, self.moltype, self.ksize, self.scaled, self.total_n_hashes
        )
    }
}

pub fn write_summary(summaries: &[RecordSummary]) -> String {
    let mut output = String::new();
    FmtWrite::write_str(&mut output, "num signatures: ").unwrap();
    writeln!(
        output,
        "{}",
        summaries.iter().map(|s| s.count).sum::<usize>()
    )
    .unwrap();
    FmtWrite::write_str(&mut output, "** examining manifest...\n").unwrap();
    writeln!(
        output,
        "total hashes: {}",
        summaries.iter().map(|s| s.total_n_hashes).sum::<usize>()
    )
    .unwrap();
    FmtWrite::write_str(&mut output, "summary of sketches:\n").unwrap();
    for summary in summaries {
        writeln!(output, "   {}", summary).unwrap();
    }

    output
}

#[cfg(test)]
mod test {
    use super::{write_summary, Manifest};
    use crate::collection::Collection;
    use crate::encodings::HashFunctions;
    use crate::selection::{Select, Selection};
    use camino::Utf8PathBuf as PathBuf;
    use std::fs::File;
    use std::io::Write;
    use tempfile::TempDir;

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
        writeln!(pathfile, "Valid line").unwrap();
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
    fn manifest_to_writer_bools() {
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
                assert_eq!(record.with_abundance(), false)
            } else {
                assert_eq!(record.with_abundance(), true)
            }
        }
    }

    #[test]
    fn manifest_to_writer_moltype_dna() {
        let base_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

        let test_sigs = vec![PathBuf::from("../../tests/test-data/47.fa.sig")];

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
            eprintln!("{:?} {}", record.name(), record.moltype());
            assert_eq!(record.moltype().to_string(), "DNA");
        }
    }

    #[test]
    fn manifest_selection() {
        let base_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

        let test_sigs = vec![PathBuf::from("../../tests/test-data/prot/all.zip")];

        let full_paths: Vec<PathBuf> = test_sigs
            .into_iter()
            .map(|sig| base_path.join(sig))
            .collect();

        let collection = Collection::from_zipfile(&full_paths[0]).unwrap();
        let manifest = collection.manifest().clone();

        // check selection on manifest works
        let mut selection = Selection::default();
        selection.set_ksize(19);
        let prot_collect = manifest.select(&selection).unwrap();
        // eprintln!("{}", &prot_collect);
        assert_eq!(prot_collect.len(), 6);
        selection.set_moltype(HashFunctions::Murmur64Protein);

        let manifest = collection.manifest().clone();
        let protein_only = manifest.select(&selection).unwrap();
        assert_eq!(protein_only.len(), 2);

        let manifest = collection.manifest().clone();
        selection = Selection::default();
        selection.set_scaled(100);
        let scaled100 = manifest.select(&selection).unwrap();
        assert_eq!(scaled100.len(), 6);
    }

    #[test]
    fn manifest_summarize() {
        let base_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

        let test_sigs = vec![
            PathBuf::from(
                "../../tests/test-data/prot/protein/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig",
            ),
            PathBuf::from(
                "../../tests/test-data/prot/protein/GCA_001593935.1_ASM159393v1_protein.faa.gz.sig",
            ),
        ];

        let full_paths: Vec<PathBuf> = test_sigs
            .into_iter()
            .map(|sig| base_path.join(sig))
            .collect();

        let manifest = Manifest::from(&full_paths[..]); // pass full_paths as a slice

        let summaries = &manifest.summarize();
        // let summaries = &manifest.summarize();
        let serialized_summaries = serde_json::to_string(&summaries).unwrap();
        let output = write_summary(&summaries);
        let expected_output = "num signatures: 2\n** examining manifest...\ntotal hashes: 8214\nsummary of sketches:\n   2 sketches with protein, k=19, scaled=100  8214 total hashes\n";
        let expected_serialized = "[{\"ksize\":19,\"moltype\":\"protein\",\"num\":0,\"scaled\":100,\"with_abundance\":false,\"count\":2,\"total_n_hashes\":8214}]";
        assert_eq!(output, expected_output);
        assert_eq!(serialized_summaries.trim(), expected_serialized);
    }
}
