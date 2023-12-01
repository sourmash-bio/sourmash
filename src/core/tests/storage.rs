use std::path::PathBuf;

use tempfile::TempDir;

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
        let sigs: Vec<Signature> = serde_json::from_reader(&data[..]).expect("Loading error");
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

    assert_eq!(sig.name(), loaded_sig.name());
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
