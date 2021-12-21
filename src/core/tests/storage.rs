use std::fs::File;
use std::path::PathBuf;

use sourmash::storage::{Storage, ZipStorage};

#[test]
fn zipstorage_load_file() -> Result<(), Box<dyn std::error::Error>> {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/v6.sbt.zip");

    let zs = ZipStorage::new(filename.to_str().unwrap())?;

    let data = zs.load("v6.sbt.json")?;

    let description: serde_json::Value = serde_json::from_slice(&data[..])?;
    assert_eq!(description["version"], 6);

    Ok(())
}

#[test]
fn zipstorage_load_slice() -> Result<(), Box<dyn std::error::Error>> {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/v6.sbt.zip");

    let zip_file = File::open(filename)?;
    let mapping = unsafe { memmap2::Mmap::map(&zip_file)? };

    let zs = ZipStorage::from_slice(&mapping)?;

    let data = zs.load("v6.sbt.json")?;

    let description: serde_json::Value = serde_json::from_slice(&data[..])?;
    assert_eq!(description["version"], 6);

    Ok(())
}
