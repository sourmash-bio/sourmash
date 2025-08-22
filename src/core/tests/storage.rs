use std::path::PathBuf;

use tempfile::TempDir;

use sourmash::prelude::Select;
use sourmash::selection::Selection;
use sourmash::signature::Signature;
use sourmash::storage::{FSStorage, InnerStorage, Storage, StorageArgs, ZipStorage};

#[test]
fn zipstorage_load_file() -> Result<(), Box<dyn std::error::Error>> {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/v6.sbt.zip");

    let zs = ZipStorage::from_file(filename.to_str().unwrap())?;

    let data = zs.load("v6.sbt.json")?;

    let description: serde_json::Value = serde_json::from_slice(&data[..])?;
    assert_eq!(description["version"], 6);

    Ok(())
}

#[test]
fn zipstorage_load_manifest() -> Result<(), Box<dyn std::error::Error>> {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/prot/protein.sbt.zip");

    let zs = ZipStorage::from_file(filename.to_str().unwrap())?;

    let _data = zs.load("protein.manifest.csv").expect("error loading file");

    Ok(())
}

#[test]
fn zipstorage_list_sbts() -> Result<(), Box<dyn std::error::Error>> {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/v6.sbt.zip");

    let zs = ZipStorage::from_file(filename.to_str().unwrap())?;

    let sbts = zs.list_sbts()?;

    assert_eq!(sbts.len(), 1);

    Ok(())
}

#[cfg(feature = "parallel")]
#[test]
fn zipstorage_parallel_access() -> Result<(), Box<dyn std::error::Error>> {
    use rayon::prelude::*;
    use sourmash::signature::SigsTrait;

    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/v6.sbt.zip");

    let zs = ZipStorage::from_file(filename.to_str().unwrap())?;

    let total_hashes: usize = [
        ".sbt.v3/f71e78178af9e45e6f1d87a0c53c465c",
        ".sbt.v3/f0c834bc306651d2b9321fb21d3e8d8f",
        ".sbt.v3/4e94e60265e04f0763142e20b52c0da1",
        ".sbt.v3/6d6e87e1154e95b279e5e7db414bc37b",
        ".sbt.v3/0107d767a345eff67ecdaed2ee5cd7ba",
        ".sbt.v3/b59473c94ff2889eca5d7165936e64b3",
        ".sbt.v3/60f7e23c24a8d94791cc7a8680c493f9",
    ]
    .par_iter()
    .map(|path| {
        let data = zs.load(path).unwrap();
        let sigs = Signature::from_reader(&data[..]).expect("Loading error");
        sigs.iter()
            .map(|v| v.sketches().iter().map(|mh| mh.size()).sum::<usize>())
            .sum::<usize>()
    })
    .sum();

    assert_eq!(total_hashes, 3500);

    Ok(())
}

#[test]
fn innerstorage_save_sig() -> Result<(), Box<dyn std::error::Error>> {
    let output = TempDir::new()?;

    let fst = FSStorage::new("".into(), output.path().as_os_str().to_str().unwrap());

    let instorage = InnerStorage::new(fst);

    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/genome-s10.fa.gz.sig");

    let sig = Signature::from_path(filename)?.swap_remove(0);
    let new_path = instorage.save_sig("test", sig.clone())?;
    dbg!(new_path);

    let loaded_sig = instorage.load_sig("test")?;

    assert_eq!(sig.name_str(), loaded_sig.name());
    assert_eq!(sig.md5sum(), loaded_sig.md5sum());

    Ok(())
}

#[test]
fn innerstorage_load() -> Result<(), Box<dyn std::error::Error>> {
    let output = TempDir::new()?;

    let fst = FSStorage::new("".into(), output.path().as_os_str().to_str().unwrap());

    let instorage = InnerStorage::new(fst);

    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/genome-s10.fa.gz.sig");

    let sig = Signature::from_path(filename)?.swap_remove(0);
    let new_path = instorage.save_sig("test", sig.clone())?;
    dbg!(new_path);

    let raw_data = instorage.load("test")?;
    let loaded_sig = Signature::from_reader(raw_data.as_slice())?.swap_remove(0);

    assert_eq!(sig.name(), loaded_sig.name());
    assert_eq!(sig.md5sum(), loaded_sig.md5sum());

    Ok(())
}

#[test]
fn innerstorage_args() -> Result<(), Box<dyn std::error::Error>> {
    let output = TempDir::new()?;
    let path = output.path().as_os_str().to_str().unwrap();

    let fst = FSStorage::new("".into(), path);

    let instorage = InnerStorage::new(fst);

    let args = instorage.args();

    assert!(matches!(args, StorageArgs::FSStorage { .. }));
    let StorageArgs::FSStorage { path: p } = args;
    assert_eq!(p, path);

    Ok(())
}

#[test]
fn innerstorage_from_args() -> Result<(), Box<dyn std::error::Error>> {
    let output = TempDir::new()?;
    let path = output.path().as_os_str().to_str().unwrap();

    let fst = FSStorage::new("".into(), path);
    let args = fst.args();

    let instorage = InnerStorage::new(FSStorage::from(&args));
    let inargs = instorage.args();

    assert!(matches!(inargs, StorageArgs::FSStorage { .. }));
    let StorageArgs::FSStorage { path: p1 } = inargs;
    assert_eq!(p1, path);

    assert!(matches!(args, StorageArgs::FSStorage { .. }));
    let StorageArgs::FSStorage { path: p2 } = args;
    assert_eq!(p2, path);

    Ok(())
}

#[test]
fn wortstorage_genomes() -> Result<(), Box<dyn std::error::Error>> {
    let storage = InnerStorage::from_spec("wort://".to_string())?;

    // Using the following signature from wort:
    // https://wort.sourmash.bio/view/genomes/GCA_000250945.2/
    // which can be downloaded from the API with the URL
    // https://wort.sourmash.bio/v1/view/genomes/GCA_000250945.2

    let mut selection = Selection::default();
    selection.set_ksize(31);

    let raw_data = storage.load("genomes/GCA_000250945.2")?;
    let loaded_sig = Signature::from_reader(raw_data.as_slice())?
        .swap_remove(0)
        .select(&selection)?;

    assert_eq!(
        loaded_sig.name(),
        Some("GCA_000250945.2 Enterococcus faecium Aus0004 strain=Aus0004, ASM25094v2".to_string())
    );
    assert_eq!(loaded_sig.md5sum(), "f6b8b19547211001f87ef397bf4ac1e1");

    Ok(())
}

#[test]
fn wortstorage_collection() -> Result<(), Box<dyn std::error::Error>> {
    use sourmash::collection::{Collection, CollectionSet};
    use sourmash::manifest::Manifest;

    let storage = InnerStorage::from_spec("wort://".to_string())?;

    let manifest = r#"# SOURMASH-MANIFEST-VERSION: 1.0
internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename
genomes/GCA_027604085.1,4bbf422430fe90c3b4d63032d604af19,4bbf4224,21,DNA,0,1000,27639,True,"GCA_027604085.1 Chytriomyces hyalinus strain=JEL0345, ASM2760408v1",/dev/fd/63
genomes/GCA_027604165.1,ff1f063c431f819b08be2955ee54ade0,ff1f063c,21,DNA,0,1000,28704,True,"GCA_027604165.1 Chytriomyces hyalinus strain=ARG085, ASM2760416v1",/dev/fd/63
genomes/GCA_027604745.1,0b626d3a644f585ca09d4b34875ecd6a,0b626d3a,21,DNA,0,1000,28490,True,"GCA_027604745.1 Chytriomyces hyalinus strain=JEL0176, ASM2760474v1",/dev/fd/63
genomes/GCA_900079185.1,9862ba9f9c58681ca28b4e7dd3f6a53c,9862ba9f,21,DNA,0,1000,37958,True,"GCA_900079185.1 Absidia glauca strain=CBS 101.48 substr. RVII-324 met-, AG_v1",/dev/fd/63
genomes/GCA_027604105.1,b8bfefac8cff5a1e774e553ea6d53ea4,b8bfefac,21,DNA,0,1000,28568,True,"GCA_027604105.1 Chytriomyces hyalinus strain=ARG121, ASM2760410v1",/dev/fd/63"#;

    let manifest = Manifest::from_reader(manifest.as_bytes())?;
    let collection: CollectionSet = Collection::new(manifest, storage).try_into()?;

    assert_eq!(collection.len(), 5);

    Ok(())
}

#[test]
#[cfg(all(feature = "branchwater", not(target_arch = "wasm32")))]
fn wort_collection_to_rocksdb() -> sourmash::Result<()> {
    use camino::Utf8PathBuf as PathBuf;
    use tempfile::TempDir;

    use sourmash::collection::{Collection, CollectionSet};
    use sourmash::index::revindex::{prepare_query, RevIndex, RevIndexOps};
    use sourmash::manifest::Manifest;

    let storage = InnerStorage::from_spec("wort://".to_string())?;

    let manifest = r#"# SOURMASH-MANIFEST-VERSION: 1.0
internal_location,md5,md5short,ksize,moltype,num,scaled,n_hashes,with_abundance,name,filename
genomes/GCA_027604085.1,4bbf422430fe90c3b4d63032d604af19,4bbf4224,21,DNA,0,1000,27639,True,"GCA_027604085.1 Chytriomyces hyalinus strain=JEL0345, ASM2760408v1",/dev/fd/63
genomes/GCA_027604165.1,ff1f063c431f819b08be2955ee54ade0,ff1f063c,21,DNA,0,1000,28704,True,"GCA_027604165.1 Chytriomyces hyalinus strain=ARG085, ASM2760416v1",/dev/fd/63
genomes/GCA_027604745.1,0b626d3a644f585ca09d4b34875ecd6a,0b626d3a,21,DNA,0,1000,28490,True,"GCA_027604745.1 Chytriomyces hyalinus strain=JEL0176, ASM2760474v1",/dev/fd/63
genomes/GCA_900079185.1,9862ba9f9c58681ca28b4e7dd3f6a53c,9862ba9f,21,DNA,0,1000,37958,True,"GCA_900079185.1 Absidia glauca strain=CBS 101.48 substr. RVII-324 met-, AG_v1",/dev/fd/63
genomes/GCA_027604105.1,b8bfefac8cff5a1e774e553ea6d53ea4,b8bfefac,21,DNA,0,1000,28568,True,"GCA_027604105.1 Chytriomyces hyalinus strain=ARG121, ASM2760410v1",/dev/fd/63"#;

    let manifest = Manifest::from_reader(manifest.as_bytes())?;
    let collection: CollectionSet = Collection::new(manifest, storage).try_into()?;

    let basedir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let outdir = TempDir::new()?;
    let output: PathBuf = outdir.path().join("index").try_into().unwrap();

    // Step 1: create an index
    let index = RevIndex::create(output.as_path(), collection.clone())?;

    // Step 2: internalize the storage for the index
    {
        let mut index = index;
        index
            .internalize_storage()
            .expect("Error internalizing storage");
    }

    let index = RevIndex::open(output.as_path(), true, None)?;

    // Step 3: Create a new collection from rocksdb
    let new_collection: CollectionSet = Collection::from_rocksdb(output.as_path())?.try_into()?;

    // Step 4: assert all content is the same
    for (a, b) in collection.iter().zip(new_collection.iter()) {
        assert_eq!(a, b);
    }

    // Step 5: can we search and get results from the new index?
    let query_sig = new_collection.sig_for_dataset(0)?;
    let selection = new_collection.selection();

    let query_mh =
        prepare_query(query_sig.into(), &selection).expect("can't get compatible MinHash");

    let results = index.find_signatures(&query_mh, 1.0, None)?;
    assert_eq!(results.len(), 1);

    Ok(())
}
