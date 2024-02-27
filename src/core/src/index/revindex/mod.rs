pub mod disk_revindex;
pub mod mem_revindex;

use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::path::Path;
use std::sync::Arc;

use byteorder::{LittleEndian, WriteBytesExt};
use enum_dispatch::enum_dispatch;
use getset::{Getters, Setters};
use nohash_hasher::BuildNoHashHasher;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};

use crate::collection::CollectionSet;
use crate::encodings::{Color, Colors, Idx};
use crate::index::{GatherResult, SigCounter};
use crate::prelude::*;
use crate::signature::Signature;
use crate::sketch::minhash::KmerMinHash;
use crate::sketch::Sketch;
use crate::HashIntoType;
use crate::Result;

type DB = rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>;

type QueryColors = HashMap<Color, Datasets>;
type HashToColorT = HashMap<HashIntoType, Color, BuildNoHashHasher<HashIntoType>>;
#[derive(Serialize, Deserialize)]
pub struct HashToColor(HashToColorT);

// Column families
const HASHES: &str = "hashes";
const COLORS: &str = "colors";
const METADATA: &str = "metadata";

// DB metadata saved in the METADATA column family
const MANIFEST: &str = "manifest";
const STORAGE_SPEC: &str = "storage_spec";
const VERSION: &str = "version";

#[enum_dispatch(RevIndexOps)]
pub enum RevIndex {
    //Color(color_revindex::ColorRevIndex),
    Plain(disk_revindex::RevIndex),
    //Mem(mem_revindex::RevIndex),
}

#[enum_dispatch]
pub trait RevIndexOps {
    /* TODO: need the repair_cf variant, not available in rocksdb-rust yet
        pub fn repair(index: &Path, colors: bool);
    */

    fn counter_for_query(&self, query: &KmerMinHash) -> SigCounter;

    fn matches_from_counter(&self, counter: SigCounter, threshold: usize) -> Vec<(String, usize)>;

    fn prepare_gather_counters(
        &self,
        query: &KmerMinHash,
    ) -> (SigCounter, QueryColors, HashToColor);

    fn update(self, collection: CollectionSet) -> Result<RevIndex>
    where
        Self: Sized;

    fn compact(&self);

    fn flush(&self) -> Result<()>;

    fn convert(&self, output_db: RevIndex) -> Result<()>;

    fn check(&self, quick: bool) -> DbStats;

    fn gather(
        &self,
        counter: SigCounter,
        query_colors: QueryColors,
        hash_to_color: HashToColor,
        threshold: usize,
        query: &KmerMinHash,
        selection: Option<Selection>,
    ) -> Result<Vec<GatherResult>>;
}

impl HashToColor {
    fn new() -> Self {
        HashToColor(HashMap::<
            HashIntoType,
            Color,
            BuildNoHashHasher<HashIntoType>,
        >::with_hasher(BuildNoHashHasher::default()))
    }

    fn get(&self, hash: &HashIntoType) -> Option<&Color> {
        self.0.get(hash)
    }

    fn len(&self) -> usize {
        self.0.len()
    }

    fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    fn add_to(&mut self, colors: &mut Colors, dataset_id: Idx, matched_hashes: Vec<u64>) {
        let mut color = None;

        matched_hashes.into_iter().for_each(|hash| {
            color = Some(colors.update(color, &[dataset_id]).unwrap());
            self.0.insert(hash, color.unwrap());
        });
    }

    fn reduce_hashes_colors(
        a: (HashToColor, Colors),
        b: (HashToColor, Colors),
    ) -> (HashToColor, Colors) {
        let ((small_hashes, small_colors), (mut large_hashes, mut large_colors)) =
            if a.0.len() > b.0.len() {
                (b, a)
            } else {
                (a, b)
            };

        small_hashes.0.into_iter().for_each(|(hash, color)| {
            large_hashes
                .0
                .entry(hash)
                .and_modify(|entry| {
                    // Hash is already present.
                    // Update the current color by adding the indices from
                    // small_colors.
                    let ids = small_colors.indices(&color);
                    let new_color = large_colors.update(Some(*entry), ids).unwrap();
                    *entry = new_color;
                })
                .or_insert_with(|| {
                    // In this case, the hash was not present yet.
                    // we need to create the same color from small_colors
                    // into large_colors.
                    let ids = small_colors.indices(&color);
                    let new_color = large_colors.update(None, ids).unwrap();
                    assert_eq!(new_color, color);
                    new_color
                });
        });

        (large_hashes, large_colors)
    }
}

impl FromIterator<(HashIntoType, Color)> for HashToColor {
    fn from_iter<T>(iter: T) -> Self
    where
        T: IntoIterator<Item = (HashIntoType, Color)>,
    {
        HashToColor(HashToColorT::from_iter(iter))
    }
}

impl RevIndex {
    /* TODO: need the repair_cf variant, not available in rocksdb-rust yet
        pub fn repair(index: &Path, colors: bool) {
            if colors {
                color_revindex::repair(index);
            } else {
                disk_revindex::repair(index);
            }
        }
    */

    pub fn create<P: AsRef<Path>>(
        index: P,
        collection: CollectionSet,
        colors: bool,
    ) -> Result<Self> {
        if colors {
            todo!() //color_revindex::ColorRevIndex::create(index)
        } else {
            disk_revindex::RevIndex::create(index.as_ref(), collection)
        }
    }

    pub fn open<P: AsRef<Path>>(index: P, read_only: bool, spec: Option<&str>) -> Result<Self> {
        let opts = Self::db_options();
        let cfs = DB::list_cf(&opts, index.as_ref()).unwrap();

        if cfs.into_iter().any(|c| c == COLORS) {
            // TODO: ColorRevIndex can't be read-only for now,
            //       due to pending unmerged colors
            todo!() //color_revindex::ColorRevIndex::open(index, false)
        } else {
            disk_revindex::RevIndex::open(index, read_only, spec)
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
        block_opts.set_format_version(6);
        opts.set_block_based_table_factory(&block_opts);
        // End of updated defaults

        opts.increase_parallelism(rayon::current_num_threads() as i32);
        //opts.max_background_jobs = 6;
        // opts.optimize_level_style_compaction();
        // opts.optimize_universal_style_compaction();

        opts
    }
}

pub fn prepare_query(search_sig: Signature, selection: &Selection) -> Option<KmerMinHash> {
    let sig = search_sig.select(selection).ok();

    sig.and_then(|sig| {
        if let Sketch::MinHash(mh) = sig.sketches().swap_remove(0) {
            Some(mh)
        } else {
            None
        }
    })
}

#[derive(Debug, Default, Clone)]
pub enum Datasets {
    #[default]
    Empty,
    Unique(Idx),
    Many(RoaringBitmap),
}

impl Hash for Datasets {
    fn hash<H>(&self, state: &mut H)
    where
        H: Hasher,
    {
        match self {
            Self::Empty => todo!(),
            Self::Unique(v) => v.hash(state),
            Self::Many(v) => {
                for value in v.iter() {
                    value.hash(state);
                }
            }
        }
    }
}

impl IntoIterator for Datasets {
    type Item = Idx;
    type IntoIter = Box<dyn Iterator<Item = Self::Item>>;

    fn into_iter(self) -> Self::IntoIter {
        match self {
            Self::Empty => Box::new(std::iter::empty()),
            Self::Unique(v) => Box::new(std::iter::once(v)),
            Self::Many(v) => Box::new(v.into_iter()),
        }
    }
}

impl Extend<Idx> for Datasets {
    fn extend<T>(&mut self, iter: T)
    where
        T: IntoIterator<Item = Idx>,
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
    fn new(vals: &[Idx]) -> Self {
        if vals.is_empty() {
            Self::Empty
        } else if vals.len() == 1 {
            Self::Unique(vals[0])
        } else {
            Self::Many(RoaringBitmap::from_sorted_iter(vals.iter().copied()).unwrap())
        }
    }

    fn from_slice(slice: &[u8]) -> Option<Self> {
        use byteorder::ReadBytesExt;

        if slice.len() == 8 {
            // Unique
            Some(Self::Unique(
                (&slice[..]).read_u32::<LittleEndian>().unwrap(),
            ))
        } else if slice.len() == 1 {
            // Empty
            Some(Self::Empty)
        } else {
            // Many
            Some(Self::Many(RoaringBitmap::deserialize_from(slice).unwrap()))
        }
    }

    fn as_bytes(&self) -> Option<Vec<u8>> {
        match self {
            Self::Empty => Some(vec![42_u8]),
            Self::Unique(v) => {
                let mut buf = vec![0u8; 8];
                (&mut buf[..])
                    .write_u32::<LittleEndian>(*v)
                    .expect("error writing bytes");
                Some(buf)
            }
            Self::Many(v) => {
                let mut buf = vec![];
                v.serialize_into(&mut buf).unwrap();
                Some(buf)
            }
        }
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
                Datasets::Many(mut o) => {
                    o.extend([*v]);
                    *self = Datasets::Many(o);
                }
            },
            Datasets::Many(ref mut v) => v.extend(other),
        }
    }

    fn len(&self) -> usize {
        match self {
            Self::Empty => 0,
            Self::Unique(_) => 1,
            Self::Many(ref v) => v.len() as usize,
        }
    }

    /*
    fn contains(&self, value: &Idx) -> bool {
        match self {
            Self::Empty => false,
            Self::Unique(v) => v == value,
            Self::Many(ref v) => v.contains(*value),
        }
    }
    */
}

#[derive(Getters, Setters, Debug)]
pub struct DbStats {
    #[getset(get = "pub")]
    total_datasets: usize,

    #[getset(get = "pub")]
    total_keys: usize,

    #[getset(get = "pub")]
    kcount: usize,

    #[getset(get = "pub")]
    vcount: usize,

    #[getset(get = "pub")]
    vcounts: histogram::Histogram,
}

fn stats_for_cf(db: Arc<DB>, cf_name: &str, deep_check: bool, quick: bool) -> DbStats {
    use byteorder::ReadBytesExt;
    use histogram::Histogram;

    let cf = db.cf_handle(cf_name).unwrap();

    let iter = db.iterator_cf(&cf, rocksdb::IteratorMode::Start);
    let mut kcount = 0;
    let mut vcount = 0;
    // Using power values from https://docs.rs/histogram/0.8.3/histogram/struct.Config.html#resulting-size
    let mut vcounts = Histogram::new(12, 64).expect("Error initializing histogram");
    let mut datasets: Datasets = Default::default();

    for result in iter {
        let (key, value) = result.unwrap();
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

    DbStats {
        total_datasets: datasets.len(),
        total_keys: kcount / 8,
        kcount,
        vcount,
        vcounts,
    }
}

#[cfg(test)]
mod test {

    use camino::Utf8PathBuf as PathBuf;
    use tempfile::TempDir;

    use crate::collection::Collection;
    use crate::prelude::*;
    use crate::selection::Selection;
    use crate::Result;

    use super::{prepare_query, RevIndex, RevIndexOps};

    #[test]
    fn revindex_index() -> Result<()> {
        let mut basedir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        basedir.push("../../tests/test-data/scaled/");

        let siglist: Vec<_> = (10..=12)
            .map(|i| {
                let mut filename = basedir.clone();
                filename.push(format!("genome-s{}.fa.gz.sig", i));
                filename
            })
            .collect();

        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let output = TempDir::new()?;

        let mut query = None;
        let query_sig = Signature::from_path(&siglist[0])?
            .swap_remove(0)
            .select(&selection)?;
        if let Some(q) = prepare_query(query_sig, &selection) {
            query = Some(q);
        }
        let query = query.unwrap();

        let collection = Collection::from_paths(&siglist)?.select(&selection)?;
        let index = RevIndex::create(output.path(), collection.try_into()?, false)?;

        let counter = index.counter_for_query(&query);
        let matches = index.matches_from_counter(counter, 0);

        assert_eq!(matches, [("../genome-s10.fa.gz".into(), 48)]);

        Ok(())
    }

    #[test]
    fn revindex_update() -> Result<()> {
        let mut basedir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        basedir.push("../../tests/test-data/scaled/");

        let siglist: Vec<_> = (10..=11)
            .map(|i| {
                let mut filename = basedir.clone();
                filename.push(format!("genome-s{}.fa.gz.sig", i));
                filename
            })
            .collect();

        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let output = TempDir::new()?;

        let mut new_siglist = siglist.clone();
        {
            let collection = Collection::from_paths(&siglist)?.select(&selection)?;
            RevIndex::create(output.path(), collection.try_into()?, false)?;
        }

        let mut filename = basedir.clone();
        filename.push("genome-s12.fa.gz.sig");
        new_siglist.push(filename);

        let mut query = None;
        let query_sig = Signature::from_path(&new_siglist[2])?
            .swap_remove(0)
            .select(&selection)?;
        if let Some(q) = prepare_query(query_sig, &selection) {
            query = Some(q);
        }
        let query = query.unwrap();

        let new_collection = Collection::from_paths(&new_siglist)?.select(&selection)?;
        let index =
            RevIndex::open(output.path(), false, None)?.update(new_collection.try_into()?)?;

        let counter = index.counter_for_query(&query);
        let matches = index.matches_from_counter(counter, 0);

        assert!(matches[0].0.ends_with("/genome-s12.fa.gz"));
        assert_eq!(matches[0].1, 45);

        Ok(())
    }

    #[test]
    fn revindex_load_and_gather() -> Result<()> {
        let mut basedir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        basedir.push("../../tests/test-data/scaled/");

        let siglist: Vec<_> = (10..=12)
            .map(|i| {
                let mut filename = basedir.clone();
                filename.push(format!("genome-s{}.fa.gz.sig", i));
                filename
            })
            .collect();

        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let output = TempDir::new()?;

        let mut query = None;
        let query_sig = Signature::from_path(&siglist[0])?
            .swap_remove(0)
            .select(&selection)?;
        if let Some(q) = prepare_query(query_sig, &selection) {
            query = Some(q);
        }
        let query = query.unwrap();

        {
            let collection = Collection::from_paths(&siglist)?.select(&selection)?;
            let _index = RevIndex::create(output.path(), collection.try_into()?, false);
        }

        let index = RevIndex::open(output.path(), true, None)?;

        let (counter, query_colors, hash_to_color) = index.prepare_gather_counters(&query);

        let matches = index.gather(
            counter,
            query_colors,
            hash_to_color,
            0,
            &query,
            Some(selection),
        )?;

        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].name(), "../genome-s10.fa.gz");
        assert_eq!(matches[0].f_match(), 1.0);

        Ok(())
    }

    #[test]
    fn revindex_move() -> Result<()> {
        let basedir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

        let mut zip_collection = basedir.clone();
        zip_collection.push("../../tests/test-data/track_abund/track_abund.zip");

        let outdir = TempDir::new()?;

        let zip_copy = PathBuf::from(
            outdir
                .path()
                .join("sigs.zip")
                .into_os_string()
                .into_string()
                .unwrap(),
        );
        std::fs::copy(zip_collection, zip_copy.as_path())?;

        let selection = Selection::builder().ksize(31).scaled(10000).build();
        let collection = Collection::from_zipfile(zip_copy.as_path())?.select(&selection)?;
        let output = outdir.path().join("index");

        let query = prepare_query(collection.sig_for_dataset(0)?.into(), &selection).unwrap();

        {
            RevIndex::create(output.as_path(), collection.try_into()?, false)?;
        }

        {
            let index = RevIndex::open(output.as_path(), false, None)?;

            let counter = index.counter_for_query(&query);
            let matches = index.matches_from_counter(counter, 0);

            assert!(matches[0].0.starts_with("NC_009665.1"));
            assert_eq!(matches[0].1, 514);
        }

        let new_zip = outdir
            .path()
            .join("new_sigs.zip")
            .into_os_string()
            .into_string()
            .unwrap();
        std::fs::rename(zip_copy, &new_zip)?;

        // RevIndex can't know where the new sigs are
        assert!(RevIndex::open(output.as_path(), false, None).is_err());

        let index = RevIndex::open(output.as_path(), false, Some(&format!("zip://{}", new_zip)))?;

        let counter = index.counter_for_query(&query);
        let matches = index.matches_from_counter(counter, 0);

        assert!(matches[0].0.starts_with("NC_009665.1"));
        assert_eq!(matches[0].1, 514);

        Ok(())
    }
}
