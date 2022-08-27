pub mod mem_revindex;
pub mod revindex;

use std::collections::{BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::iter::FromIterator;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use byteorder::{LittleEndian, WriteBytesExt};
use rkyv::{Archive, Deserialize, Serialize};

use crate::index::revindex::mem_revindex::GatherResult;
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use crate::sketch::Sketch;

use crate::encodings::Color;

//type DB = rocksdb::DBWithThreadMode<rocksdb::SingleThreaded>;
type DB = rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>;

type DatasetID = u64;
type SigCounter = counter::Counter<DatasetID>;
type QueryColors = HashMap<Color, Datasets>;
type HashToColor = HashMap<DatasetID, Color>;

const HASHES: &str = "hashes";
const SIGS: &str = "signatures";
const COLORS: &str = "colors";

pub enum RevIndex {
    //Color(color_revindex::ColorRevIndex),
    Plain(revindex::RevIndex),
}

impl RevIndex {
    /* TODO: need the repair_cf variant, not available in rocksdb-rust yet
        pub fn repair(index: &Path, colors: bool) {
            if colors {
                color_revindex::repair(index);
            } else {
                revindex::repair(index);
            }
        }
    */

    pub fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter {
        match self {
            //Self::Color(db) => db.counter_for_query(query),
            Self::Plain(db) => db.counter_for_query(query),
        }
    }

    pub fn matches_from_counter(
        &self,
        counter: SigCounter,
        threshold: usize,
    ) -> Vec<(String, usize)> {
        match self {
            //Self::Color(db) => todo!(), //db.matches_from_counter(counter, threshold),
            Self::Plain(db) => db.matches_from_counter(counter, threshold),
        }
    }

    pub fn prepare_gather_counters(
        &self,
        query: &KmerMinHash,
    ) -> (SigCounter, QueryColors, HashToColor) {
        match self {
            //Self::Color(_db) => todo!(), //db.prepare_gather_counters(query),
            Self::Plain(db) => db.prepare_gather_counters(query),
        }
    }

    pub fn index(
        &self,
        index_sigs: Vec<PathBuf>,
        template: &Sketch,
        threshold: f64,
        save_paths: bool,
    ) {
        match self {
            //Self::Color(db) => db.index(index_sigs, template, threshold, save_paths),
            Self::Plain(db) => db.index(index_sigs, template, threshold, save_paths),
        }
    }

    pub fn compact(&self) {
        match self {
            //Self::Color(db) => db.compact(),
            Self::Plain(db) => db.compact(),
        };
    }

    pub fn flush(&self) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            //Self::Color(db) => db.flush(),
            Self::Plain(db) => db.flush(),
        }
    }

    pub fn convert(&self, output_db: RevIndex) -> Result<(), Box<dyn std::error::Error>> {
        match self {
            //Self::Color(_db) => todo!(),
            Self::Plain(db) => db.convert(output_db),
        }
    }

    pub fn check(&self, quick: bool) {
        match self {
            //Self::Color(db) => db.check(quick),
            Self::Plain(db) => db.check(quick),
        }
    }

    pub fn create(index: &Path, colors: bool) -> Self {
        if colors {
            todo!() //color_revindex::ColorRevIndex::create(index)
        } else {
            revindex::RevIndex::create(index)
        }
    }

    pub fn open(index: &Path, read_only: bool) -> Self {
        let opts = Self::db_options();
        let cfs = DB::list_cf(&opts, index).unwrap();

        if cfs.into_iter().any(|c| c == COLORS) {
            // TODO: ColorRevIndex can't be read-only for now,
            //       due to pending unmerged colors
            todo!() //color_revindex::ColorRevIndex::open(index, false)
        } else {
            revindex::RevIndex::open(index, read_only)
        }
    }

    fn db_options() -> rocksdb::Options {
        let mut opts = rocksdb::Options::default();
        opts.set_max_open_files(500);

        // Updated defaults from
        // https://github.com/facebook/rocksdb/wiki/Setup-Options-and-Basic-Tuning#other-general-options
        opts.set_bytes_per_sync(1048576);
        let mut block_opts = rocksdb::BlockBasedOptions::default();
        block_opts.set_block_size(16 * 1024);
        block_opts.set_cache_index_and_filter_blocks(true);
        block_opts.set_pin_l0_filter_and_index_blocks_in_cache(true);
        block_opts.set_format_version(5);
        opts.set_block_based_table_factory(&block_opts);
        // End of updated defaults

        opts
    }

    pub fn gather(
        &self,
        counter: SigCounter,
        query_colors: QueryColors,
        hash_to_color: HashToColor,
        threshold: usize,
        query: &KmerMinHash,
        template: &Sketch,
    ) -> Result<Vec<GatherResult>, Box<dyn std::error::Error>> {
        match self {
            //Self::Color(_db) => todo!(),
            Self::Plain(db) => db.gather(
                counter,
                query_colors,
                hash_to_color,
                threshold,
                query,
                template,
            ),
        }
    }
}

#[derive(Debug, PartialEq, Clone, Archive, Serialize, Deserialize)]
enum SignatureData {
    Empty,
    Internal(Signature),
    External(String),
}

impl Default for SignatureData {
    fn default() -> Self {
        SignatureData::Empty
    }
}

impl SignatureData {
    fn from_slice(slice: &[u8]) -> Option<Self> {
        // TODO: avoid the aligned vec allocation here
        let mut vec = rkyv::AlignedVec::new();
        vec.extend_from_slice(slice);
        let archived_value = unsafe { rkyv::archived_root::<Self>(vec.as_ref()) };
        let inner = archived_value.deserialize(&mut rkyv::Infallible).unwrap();
        Some(inner)
    }

    fn as_bytes(&self) -> Option<Vec<u8>> {
        let bytes = rkyv::to_bytes::<_, 256>(self).unwrap();
        Some(bytes.into_vec())

        /*
        let mut serializer = DefaultSerializer::default();
        let v = serializer.serialize_value(self).unwrap();
        debug_assert_eq!(v, 0);
        let buf = serializer.into_serializer().into_inner();
        debug_assert!(Datasets::from_slice(&buf.to_vec()).is_some());
        Some(buf.to_vec())
        */
    }
}

fn check_compatible_downsample(me: &KmerMinHash, other: &KmerMinHash) -> Result<(), crate::Error> {
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

    if me.ksize() != other.ksize() {
        return Err(Error::MismatchKSizes);
    }
    if me.hash_function() != other.hash_function() {
        // TODO: fix this error
        return Err(Error::MismatchDNAProt);
    }
    if me.max_hash() < other.max_hash() {
        return Err(Error::MismatchScaled);
    }
    if me.seed() != other.seed() {
        return Err(Error::MismatchSeed);
    }
    Ok(())
}

pub fn prepare_query(search_sig: &Signature, template: &Sketch) -> Option<KmerMinHash> {
    let mut search_mh = None;
    if let Some(Sketch::MinHash(mh)) = search_sig.select_sketch(template) {
        search_mh = Some(mh.clone());
    } else {
        // try to find one that can be downsampled
        if let Sketch::MinHash(template_mh) = template {
            for sketch in search_sig.sketches() {
                if let Sketch::MinHash(ref_mh) = sketch {
                    if check_compatible_downsample(&ref_mh, template_mh).is_ok() {
                        let max_hash = max_hash_for_scaled(template_mh.scaled());
                        let mh = ref_mh.downsample_max_hash(max_hash).unwrap();
                        search_mh = Some(mh);
                    }
                }
            }
        }
    }
    search_mh
}

#[derive(Debug, PartialEq, Clone, Archive, Serialize, Deserialize, Hash)]
pub enum Datasets {
    Empty,
    Unique(DatasetID),
    Many(BTreeSet<DatasetID>),
}

impl IntoIterator for Datasets {
    type Item = DatasetID;
    type IntoIter = Box<dyn Iterator<Item = Self::Item>>;

    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::Empty => Box::new(std::iter::empty()),
            Self::Unique(v) => Box::new(std::iter::once(v)),
            Self::Many(v) => Box::new(v.into_iter()),
        }
    }
}

impl Default for Datasets {
    fn default() -> Self {
        Datasets::Empty
    }
}

impl Extend<DatasetID> for Datasets {
    fn extend<T>(&mut self, iter: T)
    where
        T: IntoIterator<Item = DatasetID>,
    {
        if let Self::Many(v) = self {
            v.extend(iter);
            return;
        }

        let mut it = iter.into_iter();
        while let Some(value) = it.next() {
            match self {
                Self::Empty => *self = Datasets::Unique(value),
                Self::Unique(v) => {
                    if *v != value {
                        *self = Self::Many([*v, value].iter().copied().collect());
                    }
                }
                Self::Many(v) => {
                    v.extend(it);
                    return;
                }
            }
        }
    }
}

impl Datasets {
    fn new(vals: &[DatasetID]) -> Self {
        if vals.is_empty() {
            Self::Empty
        } else if vals.len() == 1 {
            Self::Unique(vals[0])
        } else {
            Self::Many(BTreeSet::from_iter(vals.iter().cloned()))
        }
    }

    fn from_slice(slice: &[u8]) -> Option<Self> {
        // TODO: avoid the aligned vec allocation here
        let mut vec = rkyv::AlignedVec::new();
        vec.extend_from_slice(slice);
        let archived_value = unsafe { rkyv::archived_root::<Datasets>(vec.as_ref()) };
        let inner = archived_value.deserialize(&mut rkyv::Infallible).unwrap();
        Some(inner)
    }

    fn as_bytes(&self) -> Option<Vec<u8>> {
        let bytes = rkyv::to_bytes::<_, 256>(self).unwrap();
        Some(bytes.into_vec())

        /*
        let mut serializer = DefaultSerializer::default();
        let v = serializer.serialize_value(self).unwrap();
        debug_assert_eq!(v, 0);
        let buf = serializer.into_serializer().into_inner();
        debug_assert!(Datasets::from_slice(&buf.to_vec()).is_some());
        Some(buf.to_vec())
        */
    }

    fn union(&mut self, other: Datasets) {
        match self {
            Datasets::Empty => match other {
                Datasets::Empty => (),
                Datasets::Unique(_) | Datasets::Many(_) => *self = other,
            },
            Datasets::Unique(v) => match other {
                Datasets::Empty => (),
                Datasets::Unique(o) => {
                    if *v != o {
                        *self = Datasets::Many([*v, o].iter().copied().collect())
                    }
                }
                Datasets::Many(o) => {
                    let mut new_hashset: BTreeSet<DatasetID> = [*v].iter().copied().collect();
                    new_hashset.extend(o.into_iter());
                    *self = Datasets::Many(new_hashset);
                }
            },
            Datasets::Many(ref mut v) => v.extend(other.into_iter()),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Unique(_) => 1,
            Self::Many(ref v) => v.len(),
        }
    }

    fn contains(&self, value: &DatasetID) -> bool {
        match self {
            Self::Empty => false,
            Self::Unique(v) => v == value,
            Self::Many(ref v) => v.contains(value),
        }
    }
}

fn sig_save_to_db(
    db: Arc<DB>,
    mut search_sig: Signature,
    search_mh: KmerMinHash,
    size: u64,
    threshold: f64,
    save_paths: bool,
    filename: &Path,
    dataset_id: u64,
) {
    // Save signature to DB
    let sig = if search_mh.is_empty() || size < threshold as u64 {
        SignatureData::Empty
    } else if save_paths {
        SignatureData::External(filename.to_str().unwrap().to_string())
    } else {
        search_sig.reset_sketches();
        search_sig.push(Sketch::MinHash(search_mh));
        SignatureData::Internal(search_sig)
    };

    let sig_bytes = sig.as_bytes().unwrap();
    let cf_sigs = db.cf_handle(SIGS).unwrap();
    let mut hash_bytes = [0u8; 8];
    (&mut hash_bytes[..])
        .write_u64::<LittleEndian>(dataset_id)
        .expect("error writing bytes");
    db.put_cf(&cf_sigs, &hash_bytes[..], sig_bytes.as_slice())
        .expect("error saving sig");
}

fn stats_for_cf(db: Arc<DB>, cf_name: &str, deep_check: bool, quick: bool) {
    todo!()
    /*
    use byteorder::ReadBytesExt;
    use histogram::Histogram;
    use numsep::{separate, Locale};
    use log::info;

    let cf = db.cf_handle(cf_name).unwrap();

    let iter = db.iterator_cf(&cf, rocksdb::IteratorMode::Start);
    let mut kcount = 0;
    let mut vcount = 0;
    let mut vcounts = Histogram::new();
    let mut datasets: Datasets = Default::default();

    for (key, value) in iter {
        let _k = (&key[..]).read_u64::<LittleEndian>().unwrap();
        kcount += key.len();

        //println!("Saw {} {:?}", k, Datasets::from_slice(&value));
        vcount += value.len();

        if !quick && deep_check {
            let v = Datasets::from_slice(&value).expect("Error with value");
            vcounts.increment(v.len() as u64).unwrap();
            datasets.union(v);
        }
        //println!("Saw {} {:?}", k, value);
    }

    info!("*** {} ***", cf_name);
    use size::Size;
    let ksize = Size::from_bytes(kcount);
    let vsize = Size::from_bytes(vcount);
    if !quick && cf_name == COLORS {
        info!(
            "total datasets: {}",
            separate(datasets.len(), Locale::English)
        );
    }
    info!("total keys: {}", separate(kcount / 8, Locale::English));

    info!("k: {}", ksize.to_string());
    info!("v: {}", vsize.to_string());

    if !quick && kcount > 0 && deep_check {
        info!("max v: {}", vcounts.maximum().unwrap());
        info!("mean v: {}", vcounts.mean().unwrap());
        info!("stddev: {}", vcounts.stddev().unwrap());
        info!("median v: {}", vcounts.percentile(50.0).unwrap());
        info!("p25 v: {}", vcounts.percentile(25.0).unwrap());
        info!("p75 v: {}", vcounts.percentile(75.0).unwrap());
    }
    */
}
