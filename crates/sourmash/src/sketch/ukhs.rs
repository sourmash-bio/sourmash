use failure::Error;
use serde_derive::{Deserialize, Serialize};

use crate::signature::SigsTrait;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlatUKHS {}

impl FlatUKHS {
    pub fn md5sum(&self) -> String {
        unimplemented!()
    }
}

impl SigsTrait for FlatUKHS {
    fn size(&self) -> usize {
        unimplemented!()
    }

    fn to_vec(&self) -> Vec<u64> {
        unimplemented!()
    }

    fn ksize(&self) -> usize {
        unimplemented!()
    }

    fn check_compatible(&self, _other: &Self) -> Result<(), Error> {
        unimplemented!()
    }

    fn add_sequence(&mut self, _seq: &[u8], _force: bool) -> Result<(), Error> {
        unimplemented!()
    }
}

/* FIXME bring back after succint-rs changes

use std::f64::consts::PI;
use std::fs::File;
use std::hash::BuildHasherDefault;
use std::io::{BufReader, BufWriter, Read, Write};
use std::mem;
use std::path::Path;

use itertools::Itertools;
use pdatastructs::hyperloglog::HyperLogLog;
use ukhs;

use crate::errors::SourmashError;
use crate::index::sbt::NoHashHasher;
use crate::index::storage::ToWriter;
use crate::sketch::nodegraph::Nodegraph;

#[derive(Clone)]
pub struct UKHS<T> {
    ukhs: ukhs::UKHS,
    buckets: Vec<T>,
}

impl<'a, T> UKHS<T> {
    pub fn buckets(&'a self) -> impl Iterator<Item = &'a T> {
        self.buckets.iter()
    }
}

impl UKHS<u64> {
    pub fn md5sum(&self) -> String {
        // TODO: review this!
        let mut md5_ctx = md5::Context::new();
        md5_ctx.consume(self.ukhs.k().to_string());
        self.buckets
            .iter()
            .for_each(|x| md5_ctx.consume(x.to_string()));
        format!("{:x}", md5_ctx.compute())
    }
}

impl<T> std::fmt::Debug for UKHS<T>
where
    T: std::marker::Sync + std::fmt::Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "UKHS [W={}, K={}, buckets: {:?}]",
            self.ukhs.w(),
            self.ukhs.k(),
            self.buckets
        )
    }
}

impl<T> Default for UKHS<T>
where
    T: Sync + Default,
    UKHS<T>: UKHSTrait<Storage = T>,
{
    fn default() -> Self {
        UKHS::new(7, 21).unwrap()
    }
}

pub type HLL = HyperLogLog<u64, BuildHasherDefault<NoHashHasher>>;
pub type MemberUKHS = UKHS<Nodegraph>;
pub type FlatUKHS = UKHS<u64>;
pub type UniqueUKHS = UKHS<HLL>;

pub trait UKHSTrait: SigsTrait {
    type Storage;

    fn new(ksize: usize, wsize: usize) -> Result<UKHS<Self::Storage>, Error>;

    fn reset(&mut self);

    fn merge(&mut self, other: &Self);

    fn distance(&self, other: &Self) -> f64;

    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write;

    fn save<P: AsRef<Path>>(&self, path: P, _name: &str) -> Result<(), Error> {
        let file = File::open(&path)?;
        let mut writer = BufWriter::new(file);

        self.to_writer(&mut writer)
    }

    fn load<P: AsRef<Path>>(path: P) -> Result<FlatUKHS, Error> {
        let file = File::open(&path)?;
        let reader = BufReader::new(file);

        let ukhs = FlatUKHS::from_reader(reader)?;
        Ok(ukhs)
    }

    fn from_reader<R>(rdr: R) -> Result<FlatUKHS, Error>
    where
        R: Read,
    {
        let ukhs = serde_json::from_reader(rdr)?;
        Ok(ukhs)
    }
}

impl<T> ToWriter for UKHS<T>
where
    UKHS<T>: UKHSTrait<Storage = T>,
{
    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write,
    {
        <UKHS<T> as UKHSTrait>::to_writer(self, writer)
    }
}

impl UKHSTrait for UKHS<u64> {
    type Storage = u64;

    fn new(ksize: usize, wsize: usize) -> Result<Self, Error> {
        let wk_ukhs = ukhs::UKHS::new(ksize, wsize)?;
        let len = wk_ukhs.len();

        Ok(UKHS {
            ukhs: wk_ukhs,
            buckets: vec![0; len],
        })
    }

    fn reset(&mut self) {
        self.buckets = vec![0; self.ukhs.len()];
    }

    fn merge(&mut self, other: &Self) {
        let max_buckets = self
            .buckets
            .iter()
            .zip(other.buckets.iter())
            //.map(|(c1, c2)| (*c1 / 2) + (*c2 / 2))
            .map(|(c1, c2)| u64::max(*c1, *c2))
            //.map(|(c1, c2)| u64::min(*c1, *c2))
            .collect();
        mem::replace(&mut self.buckets, max_buckets);
    }

    fn distance(&self, other: &Self) -> f64 {
        // TODO: don't iterate twice...
        let prod: f64 = self
            .buckets
            .iter()
            .zip(other.buckets.iter())
            .map(|(a, b)| (a * b) as f64)
            .sum();
        let a_sq: f64 = self.buckets.iter().map(|a| (a * a) as f64).sum();
        let b_sq: f64 = other.buckets.iter().map(|a| (a * a) as f64).sum();

        if a_sq == 0. || b_sq == 0. {
            return 0.;
        }

        let d = f64::min(prod / (a_sq.sqrt() * b_sq.sqrt()), 1.);

        // TODO: which distance?
        //
        // this is the angular distance, which is a proper metric.
        2. * d.acos() / PI
        //
        // this is the cosine distance as defined by scipy
        //1. - d

        // This is the weighted Jaccard distance
        // TODO: don't iterate twice...
        //let mins: u64 = self
        //    .buckets
        //    .iter()
        //    .zip(other.buckets.iter())
        //    .map(|(a, b)| u64::min(*a, *b))
        //    .sum();
        //let maxs: u64 = self
        //    .buckets
        //    .iter()
        //    .zip(other.buckets.iter())
        //    .map(|(a, b)| u64::max(*a, *b))
        //    .sum();
        //
        //1. - (mins as f64 / maxs as f64)
    }

    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write,
    {
        match serde_json::to_writer(writer, &self) {
            Ok(_) => Ok(()),
            Err(_) => Err(SourmashError::SerdeError.into()),
        }
    }
}

impl SigsTrait for UKHS<u64> {
    fn size(&self) -> usize {
        self.buckets.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.buckets.clone()
    }

    fn ksize(&self) -> usize {
        // TODO: return k or w here?
        self.ukhs.w()
    }

    fn check_compatible(&self, _other: &Self) -> Result<(), Error> {
        unimplemented!()
    }

    fn add_sequence(&mut self, seq: &[u8], _force: bool) -> Result<(), Error> {
        // TODO: is seq.len() > W?
        let it: Vec<(u64, u64)> = self.ukhs.hash_iter_sequence(seq)?.collect();

        // This one update every unikmer bucket with w_hash
        //it.into_iter()
        //    .map(|(_, k_hash)| {
        //        self.buckets[self.ukhs.query_bucket(k_hash).unwrap()] += 1;
        //    })
        //    .count();
        //

        // Only update the bucket for the minimum unikmer found
        for (_, group) in &it.into_iter().group_by(|(w, _)| *w) {
            let (_, unikmer) = group.min().unwrap();
            self.buckets[self.ukhs.query_bucket(unikmer).unwrap()] += 1;
        }

        Ok(())
    }
}

impl UKHSTrait for UKHS<Nodegraph> {
    type Storage = Nodegraph;

    fn new(ksize: usize, wsize: usize) -> Result<Self, Error> {
        let wk_ukhs = ukhs::UKHS::new(ksize, wsize)?;
        let len = wk_ukhs.len();

        Ok(UKHS {
            ukhs: wk_ukhs,
            buckets: vec![Nodegraph::with_tables(100_000, 4, wsize); len],
        })
    }

    fn reset(&mut self) {
        self.buckets = vec![Nodegraph::with_tables(100_000, 4, self.ukhs.w()); self.ukhs.len()];
    }

    fn merge(&mut self, _other: &Self) {
        unimplemented!()
    }

    fn distance(&self, _other: &Self) -> f64 {
        unimplemented!()
    }

    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write,
    {
        // TODO: avoid cloning?
        let flat: FlatUKHS = self.into();

        match serde_json::to_writer(writer, &flat) {
            Ok(_) => Ok(()),
            Err(_) => Err(SourmashError::SerdeError.into()),
        }
    }
}

impl SigsTrait for UKHS<Nodegraph> {
    fn size(&self) -> usize {
        self.buckets.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.buckets
            .iter()
            .map(|b| b.unique_kmers() as u64)
            .collect()
    }

    fn ksize(&self) -> usize {
        // TODO: return k or w here?
        self.ukhs.w()
    }

    fn check_compatible(&self, _other: &Self) -> Result<(), Error> {
        unimplemented!()
    }

    fn add_sequence(&mut self, seq: &[u8], _force: bool) -> Result<(), Error> {
        let it: Vec<(u64, u64)> = self.ukhs.hash_iter_sequence(seq)?.collect();

        // This one update every unikmer bucket with w_hash
        //it.into_iter()
        //    .map(|(w_hash, k_hash)| {
        //        self.buckets[self.ukhs.query_bucket(k_hash).unwrap()].count(w_hash);
        //    })
        //    .count();

        // Only update the bucket for the minimum unikmer found
        for (w_hash, group) in &it.into_iter().group_by(|(w, _)| *w) {
            let (_, unikmer) = group.min().unwrap();
            self.buckets[self.ukhs.query_bucket(unikmer).unwrap()].count(w_hash);
        }

        Ok(())
    }
}

impl From<MemberUKHS> for FlatUKHS {
    fn from(other: MemberUKHS) -> Self {
        let buckets = other.to_vec(); // TODO: implement into_vec?
        let ukhs = other.ukhs;

        FlatUKHS { ukhs, buckets }
    }
}

impl From<&MemberUKHS> for FlatUKHS {
    fn from(other: &MemberUKHS) -> Self {
        FlatUKHS {
            ukhs: other.ukhs.clone(),
            buckets: other.to_vec(), // TODO: also implement into_vec?
        }
    }
}

impl UKHSTrait for UKHS<HLL> {
    type Storage = HLL;

    fn new(ksize: usize, wsize: usize) -> Result<Self, Error> {
        let wk_ukhs = ukhs::UKHS::new(ksize, wsize)?;
        let len = wk_ukhs.len();

        let bh = BuildHasherDefault::<NoHashHasher>::default();

        Ok(UKHS {
            ukhs: wk_ukhs,
            buckets: vec![HLL::with_hash(14, bh); len], // TODO: space usage is 2^b
        })
    }

    fn reset(&mut self) {
        let bh = BuildHasherDefault::<NoHashHasher>::default();
        self.buckets = vec![HLL::with_hash(14, bh); self.ukhs.len()];
    }

    fn merge(&mut self, _other: &Self) {
        unimplemented!()
    }

    fn distance(&self, _other: &Self) -> f64 {
        unimplemented!()
    }

    fn to_writer<W>(&self, writer: &mut W) -> Result<(), Error>
    where
        W: Write,
    {
        // TODO: avoid cloning?
        let flat: FlatUKHS = self.into();

        match serde_json::to_writer(writer, &flat) {
            Ok(_) => Ok(()),
            Err(_) => Err(SourmashError::SerdeError.into()),
        }
    }
}

impl SigsTrait for UKHS<HLL> {
    fn size(&self) -> usize {
        self.buckets.len()
    }

    fn to_vec(&self) -> Vec<u64> {
        self.buckets.iter().map(|b| b.count() as u64).collect()
    }

    fn ksize(&self) -> usize {
        // TODO: return k or w here?
        self.ukhs.w()
    }

    fn check_compatible(&self, _other: &Self) -> Result<(), Error> {
        unimplemented!()
    }

    fn add_sequence(&mut self, seq: &[u8], _force: bool) -> Result<(), Error> {
        let it: Vec<(u64, u64)> = self.ukhs.hash_iter_sequence(seq)?.collect();

        // This one update every unikmer bucket with w_hash
        //it.into_iter()
        //    .map(|(w_hash, k_hash)| {
        //        self.buckets[self.ukhs.query_bucket(k_hash).unwrap()].add(&w_hash);
        //    })
        //    .count();

        // Only update the bucket for the minimum unikmer found
        for (w_hash, group) in &it.into_iter().group_by(|(w, _)| *w) {
            let (_, unikmer) = group.min().unwrap();
            self.buckets[self.ukhs.query_bucket(unikmer).unwrap()].add(&w_hash);
        }

        Ok(())
    }
}

impl From<UniqueUKHS> for FlatUKHS {
    fn from(other: UniqueUKHS) -> Self {
        let buckets = other.to_vec(); // TODO: implement into_vec?
        let ukhs = other.ukhs;

        FlatUKHS { ukhs, buckets }
    }
}

impl From<&UniqueUKHS> for FlatUKHS {
    fn from(other: &UniqueUKHS) -> Self {
        FlatUKHS {
            ukhs: other.ukhs.clone(),
            buckets: other.to_vec(), // TODO: also implement into_vec?
        }
    }
}

impl Serialize for UKHS<u64> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let n_fields = 5;

        let buckets = self.buckets.to_vec();

        let mut partial = serializer.serialize_struct("UKHS", n_fields)?;
        partial.serialize_field("signature", &buckets)?;
        partial.serialize_field("W", &self.ukhs.w())?;
        partial.serialize_field("K", &self.ukhs.k())?;
        partial.serialize_field("size", &self.buckets.len())?;
        partial.serialize_field("name", "")?;

        // TODO: properly set name
        //partial.serialize_field("name", &self.name)?;

        partial.end()
    }
}

impl<'de> Deserialize<'de> for UKHS<u64> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        #[derive(Deserialize)]
        struct TempUKHS {
            signature: Vec<u64>,
            #[serde(rename = "K")]
            k: usize,
            #[serde(rename = "W")]
            w: usize,
            //size: usize,
            //name: String,
        }

        let tmpukhs = TempUKHS::deserialize(deserializer)?;

        //TODO: remove this unwrap, need to map Failure error to serde error?
        let mut u = UKHS::<u64>::new(tmpukhs.k, tmpukhs.w).unwrap();

        u.buckets = tmpukhs.signature;

        //TODO: what to do with name?

        Ok(u)
    }
}

impl<T> PartialEq for UKHS<T>
where
    T: PartialEq,
{
    fn eq(&self, other: &UKHS<T>) -> bool {
        self.buckets
            .iter()
            .zip(other.buckets.iter())
            .all(|(b1, b2)| b1 == b2)
            && self.ukhs.w() == other.ukhs.w()
            && self.ukhs.k() == self.ukhs.k()
    }
}

// Removed this for now, because calling .into() in these doesn't
// transfer all the important information...
//impl From<FlatUKHS> for Dataset<Signature> {
//    fn from(other: FlatUKHS) -> Dataset<Signature> {
//        let data = Lazy::new();
//        data.get_or_create(|| other.into());
//
//        Dataset::builder()
//            .data(Rc::new(data))
//            .filename("")
//            .name("")
//            .metadata("")
//            .storage(None)
//            .build()
//    }
//}
//
//impl From<FlatUKHS> for Signature {
//    fn from(other: FlatUKHS) -> Signature {
//        Signature::builder()
//            .hash_function("nthash") // TODO: spec!
//            .class("draff_signature") // TODO: spec!
//            .name(Some("draff_file".into())) // TODO: spec!
//            .signatures(vec![Sketch::UKHS(other)])
//            .build()
//    }
//}

#[cfg(test)]
mod test {
    use std::path::PathBuf;

    use needletail::parse_sequence_path;

    use super::{FlatUKHS, MemberUKHS, UKHSTrait};
    use crate::signature::SigsTrait;

    #[test]
    fn ukhs_add_sequence() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/ecoli.genes.fna");

        let mut ukhs = MemberUKHS::new(9, 21).unwrap();

        parse_sequence_path(
            filename,
            |_| {},
            |record| {
                ukhs.add_sequence(&record.seq, false).unwrap();
            },
        )
        .expect("error parsing");

        // TODO: find test case...
        //assert_eq!(ukhs.to_vec(), [1, 2, 3]);
    }

    #[test]
    fn ukhs_writer_reader() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("tests/test-data/ecoli.genes.fna");

        let mut ukhs = FlatUKHS::new(9, 21).unwrap();

        parse_sequence_path(
            filename,
            |_| {},
            |record| {
                ukhs.add_sequence(&record.seq, false).unwrap();
            },
        )
        .expect("error parsing");

        let mut buffer = Vec::new();
        ukhs.to_writer(&mut buffer).unwrap();

        match FlatUKHS::from_reader(&buffer[..]) {
            Ok(new_ukhs) => {
                assert_eq!(ukhs.buckets, new_ukhs.buckets);
            }
            Err(e) => {
                dbg!(e);
                assert_eq!(1, 0);
            }
        }
    }
}
*/
