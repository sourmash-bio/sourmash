use std::path::PathBuf;

use sourmash::storage::{Storage, ZipStorage};

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
