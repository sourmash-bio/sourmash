use std::ops::{Deref, DerefMut};

use camino::Utf8Path as Path;
use camino::Utf8PathBuf as PathBuf;

use crate::encodings::Idx;
use crate::manifest::{Manifest, Record};
use crate::prelude::*;
use crate::storage::{FSStorage, InnerStorage, MemStorage, SigStore, ZipStorage};
use crate::{Error, Result, ScaledType};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

/// a Manifest and Storage, combined. Can contain any collection of signatures.

#[derive(Clone)]
pub struct Collection {
    manifest: Manifest,
    storage: InnerStorage,
}

/// A consistent collection of signatures. Can be created using `select`.

#[derive(Clone)]
pub struct CollectionSet {
    collection: Collection,
}

impl Deref for CollectionSet {
    type Target = Collection;

    fn deref(&self) -> &Self::Target {
        &self.collection
    }
}

impl DerefMut for CollectionSet {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.collection
    }
}

impl TryFrom<Collection> for CollectionSet {
    type Error = crate::Error;

    fn try_from(collection: Collection) -> Result<Self> {
        let first = if let Some(first) = collection.manifest.first() {
            first
        } else {
            // empty collection is consistent ¯\_(ツ)_/¯
            return Ok(Self { collection });
        };

        let (min_scaled, max_scaled) = collection.min_max_scaled().expect("empty collection!?");
        if min_scaled != max_scaled {
            return Err(Error::MismatchScaled);
        }

        collection
            .manifest
            .iter()
            .skip(1)
            .try_for_each(|c| first.check_compatible(c))?;

        Ok(Self { collection })
    }
}

impl CollectionSet {
    pub fn into_inner(self) -> Collection {
        self.collection
    }

    pub fn selection(&self) -> Selection {
        todo!("Extract selection from first sig")
    }

    /// Replace the storage with a new one.
    ///
    /// # Safety
    ///
    /// This method doesn't check if the manifest matches what is in the
    /// storage (which can be expensive). It is up to the caller to
    /// guarantee the manifest and storage are in sync.
    pub unsafe fn set_storage_unchecked(&mut self, storage: InnerStorage) {
        self.storage = storage;
    }
}

impl Collection {
    pub fn new(manifest: Manifest, storage: InnerStorage) -> Self {
        Self { manifest, storage }
    }

    pub fn iter(&self) -> impl Iterator<Item = (Idx, &Record)> {
        self.manifest.iter().enumerate().map(|(i, r)| (i as Idx, r))
    }

    #[cfg(feature = "parallel")]
    pub fn par_iter(&self) -> impl IndexedParallelIterator<Item = (Idx, &Record)> {
        self.manifest
            .par_iter()
            .enumerate()
            .map(|(i, r)| (i as Idx, r))
    }

    pub fn len(&self) -> usize {
        self.manifest.len()
    }

    pub fn is_empty(&self) -> bool {
        self.manifest.len() == 0
    }

    pub fn manifest(&self) -> &Manifest {
        &self.manifest
    }

    pub fn storage(&self) -> &InnerStorage {
        &self.storage
    }

    pub fn check_superset(&self, other: &Collection) -> Result<usize> {
        self.iter()
            .zip(other.iter())
            .all(|((id1, rec1), (id2, rec2))| id1 == id2 && rec1 == rec2)
            .then(|| self.len())
            // TODO: right error here
            .ok_or(Error::MismatchKSizes)
    }

    pub fn from_zipfile<P: AsRef<Path>>(zipfile: P) -> Result<Self> {
        let storage = ZipStorage::from_file(zipfile)?;
        // Load manifest from standard location in zipstorage
        let manifest = Manifest::from_reader(storage.load("SOURMASH-MANIFEST.csv")?.as_slice())?;
        Ok(Self {
            manifest,
            storage: InnerStorage::new(storage),
        })
    }

    #[cfg(all(feature = "branchwater", not(target_arch = "wasm32")))]
    pub fn from_rocksdb<P: AsRef<Path>>(dirname: P) -> Result<Self> {
        use crate::index::revindex::{RevIndex, RevIndexOps};

        let path = dirname.as_ref().as_str().to_string();
        let index = RevIndex::open(path, true, None)?;
        let collection: Collection = index.collection().clone().into_inner();

        Ok(collection)
    }

    pub fn from_sigs(sigs: Vec<Signature>) -> Result<Self> {
        let storage = MemStorage::new();

        #[cfg(feature = "parallel")]
        let iter = sigs.into_par_iter();

        #[cfg(not(feature = "parallel"))]
        let iter = sigs.into_iter();

        let records: Vec<_> = iter
            .enumerate()
            .flat_map(|(i, sig)| {
                let path = format!("{i}");
                let mut record = Record::from_sig(&sig, &path);
                let path = storage.save_sig(&path, sig).expect("Error saving sig");
                record.iter_mut().for_each(|rec| {
                    rec.set_internal_location(path.clone().into());
                });
                record
            })
            .collect();

        Ok(Self {
            manifest: records.into(),
            storage: InnerStorage::new(storage),
        })
    }

    pub fn from_paths(paths: &[PathBuf]) -> Result<Self> {
        // TODO:
        // - figure out if there is a common path between sigs for FSStorage?

        Ok(Self {
            manifest: paths.into(),
            storage: InnerStorage::new(
                FSStorage::builder()
                    .fullpath("".into())
                    .subdir("".into())
                    .build(),
            ),
        })
    }

    pub fn record_for_dataset(&self, dataset_id: Idx) -> Result<&Record> {
        Ok(&self.manifest[dataset_id as usize])
    }

    pub fn sig_for_dataset(&self, dataset_id: Idx) -> Result<SigStore> {
        let match_path = if self.manifest.is_empty() {
            ""
        } else {
            self.manifest[dataset_id as usize]
                .internal_location()
                .as_str()
        };

        let selection = Selection::from_record(&self.manifest[dataset_id as usize])?;
        let sig = self.storage.load_sig(match_path)?.select(&selection)?;
        assert_eq!(sig.signatures.len(), 1);
        Ok(sig)
    }

    pub fn sig_from_record(&self, record: &Record) -> Result<SigStore> {
        let match_path = record.internal_location().as_str();
        let selection = Selection::from_record(record)?;
        let sig = self.storage.load_sig(match_path)?.select(&selection)?;
        assert_eq!(sig.signatures.len(), 1);
        Ok(sig)
    }

    pub fn intersect_manifest(&mut self, mf: &Manifest) {
        self.manifest = self.manifest.intersect_manifest(mf);
    }

    // CTB: question, should we do something about num here?
    pub fn min_max_scaled(&self) -> Option<(&ScaledType, &ScaledType)> {
        self.manifest.first().map(|first| {
            self.manifest
                .iter()
                .fold((first.scaled(), first.scaled()), |f, r| {
                    (f.0.min(r.scaled()), f.1.max(r.scaled()))
                })
        })
    }
}

impl Select for Collection {
    fn select(mut self, selection: &Selection) -> Result<Self> {
        self.manifest = self.manifest.select(selection)?;
        Ok(self)
    }
}

#[cfg(test)]
mod test {
    use camino::Utf8PathBuf as PathBuf;
    use std::fs::File;
    use std::io::BufReader;

    use super::Collection;

    use crate::encodings::HashFunctions;
    use crate::manifest::Manifest;
    use crate::prelude::Select;
    use crate::selection::Selection;
    use crate::signature::Signature;
    #[cfg(all(feature = "branchwater", not(target_arch = "wasm32")))]
    use crate::Result;

    #[test]
    fn sigstore_selection_with_downsample() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47+63-multisig.sig");
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs = Signature::from_reader(reader).expect("Loading error");
        // create Selection object
        let mut selection = Selection::default();
        selection.set_scaled(2000);
        // load sigs into collection + select compatible signatures
        let cl = Collection::from_sigs(sigs)
            .unwrap()
            .select(&selection)
            .unwrap();
        // count collection length
        assert_eq!(cl.len(), 6);
        for (idx, _rec) in cl.iter() {
            // need to pass select again here so we actually downsample
            let this_sig = cl.sig_for_dataset(idx).unwrap().select(&selection).unwrap();
            let this_mh = this_sig.minhash().unwrap();
            assert_eq!(this_mh.scaled(), 2000);
        }
    }

    #[test]
    fn sigstore_selection_with_downsample_too_low() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47+63-multisig.sig");
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs = Signature::from_reader(reader).expect("Loading error");
        // create Selection object
        let mut selection = Selection::default();
        selection.set_scaled(500);
        // load sigs into collection + select compatible signatures
        let cl = Collection::from_sigs(sigs)
            .unwrap()
            .select(&selection)
            .unwrap();
        // no sigs should remain
        assert_eq!(cl.len(), 0);
    }

    #[test]
    fn sigstore_selection_scaled_handle_num_sig() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        // four num=500 sigs
        filename.push("../../tests/test-data/genome-s11.fa.gz.sig");
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs = Signature::from_reader(reader).expect("Loading error");
        assert_eq!(sigs.len(), 4);
        // create Selection object
        let mut selection = Selection::default();
        selection.set_scaled(1000);
        // load sigs into collection + select compatible signatures
        let cl = Collection::from_sigs(sigs)
            .unwrap()
            .select(&selection)
            .unwrap();
        // no sigs should remain
        assert_eq!(cl.len(), 0);
    }

    #[test]
    fn sigstore_selection_num() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        // four num=500 sigs
        filename.push("../../tests/test-data/genome-s11.fa.gz.sig");
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs = Signature::from_reader(reader).expect("Loading error");
        let sigs_copy = sigs.clone();
        assert_eq!(sigs.len(), 4);
        // create Selection object
        let mut selection = Selection::default();
        selection.set_num(500);
        // load sigs into collection + select compatible signatures
        let cl = Collection::from_sigs(sigs)
            .unwrap()
            .select(&selection)
            .unwrap();
        // all sigs should remain
        assert_eq!(cl.len(), 4);
        //now select diff num and none should remain
        selection.set_num(100);
        let cl2 = Collection::from_sigs(sigs_copy)
            .unwrap()
            .select(&selection)
            .unwrap();
        assert_eq!(cl2.len(), 0);
    }

    #[test]
    fn sigstore_selection_num_handle_scaled_sig() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        // four num=500 sigs
        filename.push("../../tests/test-data/47+63-multisig.sig");
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs = Signature::from_reader(reader).expect("Loading error");
        assert_eq!(sigs.len(), 6);
        // create Selection object
        let mut selection = Selection::default();
        selection.set_num(500);
        // load sigs into collection + select compatible signatures
        let cl = Collection::from_sigs(sigs)
            .unwrap()
            .select(&selection)
            .unwrap();
        // no sigs should remain
        assert_eq!(cl.len(), 0);
    }

    #[test]
    fn collection_intersect_manifest() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        // four num=500 sigs
        filename.push("../../tests/test-data/genome-s11.fa.gz.sig");
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs = Signature::from_reader(reader).expect("Loading error");
        assert_eq!(sigs.len(), 4);
        // load sigs into collection + select compatible signatures
        let mut cl = Collection::from_sigs(sigs).unwrap();
        // all sigs should remain
        assert_eq!(cl.len(), 4);

        // grab first record
        let manifest = cl.manifest();
        let record = manifest.iter().next().unwrap().clone();
        let vr = vec![record];

        // now intersect:
        let manifest2 = Manifest::from(vr);
        cl.intersect_manifest(&manifest2);
        assert_eq!(cl.len(), 1);
    }

    #[test]
    fn sigstore_sig_from_record() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/47+63-multisig.sig");
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs = Signature::from_reader(reader).expect("Loading error");
        // create Selection object
        let mut selection = Selection::default();
        selection.set_scaled(2000);
        // load sigs into collection + select compatible signatures
        let cl = Collection::from_sigs(sigs)
            .unwrap()
            .select(&selection)
            .unwrap();
        // no sigs should remain
        assert_eq!(cl.len(), 6);
        for (_idx, rec) in cl.iter() {
            dbg!("record scaled is: {}", rec.scaled());
            let this_sig = cl.sig_from_record(rec).unwrap();
            let this_mh = this_sig.minhash().unwrap();
            assert_eq!(this_mh.scaled(), 2000);
        }
    }

    #[test]
    #[should_panic] // for now...
    fn sigstore_sig_from_record_2() {
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/short.sig.gz");
        let v = [filename];
        let collection = Collection::from_paths(&v).expect("no sigs!?");

        // pull off first record
        let v: Vec<_> = collection.iter().collect();
        let (_idx, rec) = v.first().expect("no records in collection?!");

        // this will panic with "unimplemented" because there are two
        // sketches and that is not supported.
        let _first_sig = collection.sig_from_record(rec).expect("no sig!?");
    }

    #[test]
    fn sigstore_selection_moltype_zip() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename.push("../../tests/test-data/prot/hp.zip");
        // create Selection object
        let mut selection = Selection::default();
        selection.set_scaled(200);
        selection.set_moltype(HashFunctions::Murmur64Hp);
        // load sigs into collection + select compatible signatures
        let cl = Collection::from_zipfile(&filename)
            .unwrap()
            .select(&selection)
            .unwrap();
        // count collection length
        assert_eq!(cl.len(), 2);
        for (idx, _rec) in cl.iter() {
            let this_sig = cl.sig_for_dataset(idx).unwrap();
            let this_mh = this_sig.minhash().unwrap();
            assert_eq!(this_mh.scaled(), 200);
        }
    }

    #[test]
    fn sigstore_selection_moltype_sig() {
        // load test sigs
        let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        filename
            .push("../../tests/test-data/prot/hp/GCA_001593925.1_ASM159392v1_protein.faa.gz.sig");
        let file = File::open(filename).unwrap();
        let reader = BufReader::new(file);
        let sigs = Signature::from_reader(reader).expect("Loading error");
        // create Selection object
        let mut selection = Selection::default();
        selection.set_moltype(HashFunctions::Murmur64Hp);
        // load sigs into collection + select compatible signatures
        let cl = Collection::from_sigs(sigs)
            .unwrap()
            .select(&selection)
            .unwrap();
        // count collection length
        assert_eq!(cl.len(), 1);
        for (idx, _rec) in cl.iter() {
            // need to pass select again here so we actually downsample
            let this_sig = cl.sig_for_dataset(idx).unwrap().select(&selection).unwrap();
            let this_mh = this_sig.minhash().unwrap();
            assert_eq!(this_mh.scaled(), 100);
        }
    }

    #[test]
    fn collection_from_collectionset() -> () {
        use crate::collection::CollectionSet;

        let base_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

        let test_sigs = vec![PathBuf::from("../../tests/test-data/prot/all.zip")];

        let full_paths: Vec<PathBuf> = test_sigs
            .into_iter()
            .map(|sig| base_path.join(sig))
            .collect();

        let collection = Collection::from_zipfile(&full_paths[0]).unwrap();

        let mut selection = Selection::default();
        selection.set_moltype(HashFunctions::Murmur64Protein);
        selection.set_scaled(200);

        let collection = collection.select(&selection).expect("should pass");
        let (min_scaled, max_scaled) = collection.min_max_scaled().expect("not empty");
        assert_eq!(*min_scaled, *max_scaled);
        assert_eq!(*min_scaled, 200);
        let _cs: CollectionSet = collection.try_into().expect("should pass");
    }

    #[test]
    #[should_panic]
    fn collection_from_collectionset_fail() -> () {
        use crate::collection::CollectionSet;

        let base_path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

        let test_sigs = vec![PathBuf::from("../../tests/test-data/prot/all.zip")];

        let full_paths: Vec<PathBuf> = test_sigs
            .into_iter()
            .map(|sig| base_path.join(sig))
            .collect();

        let collection = Collection::from_zipfile(&full_paths[0]).unwrap();
        let _cs: CollectionSet = collection.try_into().expect("should fail");
    }

    #[test]
    #[cfg(all(feature = "branchwater", not(target_arch = "wasm32")))]
    fn collection_from_rocksdb_storage() -> Result<()> {
        use crate::index::revindex::{RevIndex, RevIndexOps};
        use camino::Utf8PathBuf as PathBuf;
        use tempfile::TempDir;

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
        let output: PathBuf = outdir.path().join("index").try_into().unwrap();

        // Step 1: create an index
        let index = RevIndex::create(output.as_path(), collection.clone().try_into()?)?;

        // Step 2: internalize the storage for the index
        {
            let mut index = index;
            index
                .internalize_storage()
                .expect("Error internalizing storage");
        }

        // Step 3: Create a new collection from rocksdb
        let new_collection = Collection::from_rocksdb(output.as_path())?;

        // Step 4: assert all content is the same
        for (a, b) in collection.iter().zip(new_collection.iter()) {
            assert_eq!(a, b);
        }

        Ok(())
    }
}
