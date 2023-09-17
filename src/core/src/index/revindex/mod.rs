pub mod disk_revindex;
pub mod mem_revindex;

use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::path::Path;
use std::sync::Arc;

use byteorder::{LittleEndian, WriteBytesExt};
use enum_dispatch::enum_dispatch;
use nohash_hasher::BuildNoHashHasher;
use roaring::RoaringBitmap;
use serde::{Deserialize, Serialize};

use crate::collection::CollectionSet;
use crate::encodings::{Color, Colors, Idx};
use crate::index::{GatherResult, SigCounter};
use crate::signature::{Signature, SigsTrait};
use crate::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use crate::sketch::Sketch;
use crate::HashIntoType;
use crate::Result;

type DB = rocksdb::DBWithThreadMode<rocksdb::MultiThreaded>;

type QueryColors = HashMap<Color, Datasets>;
type HashToColorT = HashMap<HashIntoType, Color, BuildNoHashHasher<HashIntoType>>;
#[derive(Serialize, Deserialize)]
pub struct HashToColor(HashToColorT);

const HASHES: &str = "hashes";
const COLORS: &str = "colors";
const METADATA: &str = "metadata";
const MANIFEST: &str = "manifest";
const STORAGE_SPEC: &str = "storage_spec";

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

    fn check(&self, quick: bool);

    fn gather(
        &self,
        counter: SigCounter,
        query_colors: QueryColors,
        hash_to_color: HashToColor,
        threshold: usize,
        query: &KmerMinHash,
        template: &Sketch,
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

    pub fn open<P: AsRef<Path>>(index: P, read_only: bool) -> Result<Self> {
        let opts = Self::db_options();
        let cfs = DB::list_cf(&opts, index.as_ref()).unwrap();

        if cfs.into_iter().any(|c| c == COLORS) {
            // TODO: ColorRevIndex can't be read-only for now,
            //       due to pending unmerged colors
            todo!() //color_revindex::ColorRevIndex::open(index, false)
        } else {
            disk_revindex::RevIndex::open(index, read_only)
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

        opts.increase_parallelism(8);
        // opts.optimize_level_style_compaction();
        // opts.optimize_universal_style_compaction();

        opts
    }
}

fn check_compatible_downsample(me: &KmerMinHash, other: &KmerMinHash) -> Result<()> {
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

#[derive(Debug, Default, PartialEq, Clone)]
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

fn stats_for_cf(db: Arc<DB>, cf_name: &str, deep_check: bool, quick: bool) {
    use byteorder::ReadBytesExt;
    use histogram::Histogram;
    use log::info;
    use numsep::{separate, Locale};

    let cf = db.cf_handle(cf_name).unwrap();

    let iter = db.iterator_cf(&cf, rocksdb::IteratorMode::Start);
    let mut kcount = 0;
    let mut vcount = 0;
    let mut vcounts = Histogram::new();
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
}

#[cfg(test)]
mod test {

    use camino::Utf8PathBuf as PathBuf;
    use tempfile::TempDir;

    use crate::collection::Collection;
    use crate::prelude::*;
    use crate::selection::Selection;
    use crate::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
    use crate::sketch::Sketch;
    use crate::Result;

    use super::{prepare_query, RevIndex, RevIndexOps};

    fn build_template(ksize: u8, scaled: usize) -> Sketch {
        let max_hash = max_hash_for_scaled(scaled as u64);
        let template_mh = KmerMinHash::builder()
            .num(0u32)
            .ksize(ksize as u32)
            .max_hash(max_hash)
            .build();
        Sketch::MinHash(template_mh)
    }

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

        let template = build_template(31, 10000);
        let output = TempDir::new()?;

        let query_sig = Signature::from_path(&siglist[0])?;
        let mut query = None;
        for sig in &query_sig {
            if let Some(q) = prepare_query(sig, &template) {
                query = Some(q);
            }
        }
        let query = query.unwrap();

        let collection =
            Collection::from_paths(&siglist)?.select(&Selection::from_template(&template))?;
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

        let template = build_template(31, 10000);
        let output = TempDir::new()?;

        let mut new_siglist = siglist.clone();
        {
            let collection =
                Collection::from_paths(&siglist)?.select(&Selection::from_template(&template))?;
            RevIndex::create(output.path(), collection.try_into()?, false)?;
        }

        let mut filename = basedir.clone();
        filename.push("genome-s12.fa.gz.sig");
        new_siglist.push(filename);

        let query_sig = Signature::from_path(&new_siglist[2])?;
        let mut query = None;
        for sig in &query_sig {
            if let Some(q) = prepare_query(sig, &template) {
                query = Some(q);
            }
        }
        let query = query.unwrap();

        let new_collection =
            Collection::from_paths(&new_siglist)?.select(&Selection::from_template(&template))?;
        let index = RevIndex::open(output.path(), false)?.update(new_collection.try_into()?)?;

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

        let template = build_template(31, 10000);
        let output = TempDir::new()?;

        let query_sig = Signature::from_path(&siglist[0])?;
        let mut query = None;
        for sig in &query_sig {
            if let Some(q) = prepare_query(sig, &template) {
                query = Some(q);
            }
        }
        let query = query.unwrap();

        {
            let collection =
                Collection::from_paths(&siglist)?.select(&Selection::from_template(&template))?;
            let _index = RevIndex::create(output.path(), collection.try_into()?, false);
        }

        let index = RevIndex::open(output.path(), true)?;

        let counter = index.counter_for_query(&query);
        let matches = index.matches_from_counter(counter, 0);

        assert_eq!(matches, [("../genome-s10.fa.gz".into(), 48)]);

        Ok(())
    }
}
