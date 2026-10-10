//! Scan a disk RevIndex for dataset IDs outside the manifest and for values that do not parse.
//! Usage: cargo run --release --features branchwater --example scan_revindex_ids -- <index path>

use byteorder::{LittleEndian, ReadBytesExt};

use sourmash::index::revindex::{Datasets, RevIndex, RevIndexOps};

// Name of the hashes column family in `sourmash::storage::rocksdb`.
const HASHES: &str = "hashes";

fn hex(bytes: &[u8]) -> String {
    bytes
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let path = std::env::args()
        .nth(1)
        .ok_or("usage: scan_revindex_ids <index path>")?;

    let index = RevIndex::open(&path, true, None)?;
    let n_datasets = index.collection().len() as u32;
    eprintln!("manifest has {n_datasets} datasets");

    let RevIndex::Disk(disk) = index else {
        return Err("not a disk RevIndex".into());
    };
    // SAFETY: the index is open read-only and the scan only reads.
    let db = unsafe { disk.db() };
    let cf = db.cf_handle(HASHES).ok_or("no hashes column family")?;

    let mut n_keys = 0u64;
    let mut n_bad = 0u64;
    for item in db.iterator_cf(&cf, rocksdb::IteratorMode::Start) {
        let (key, value) = item?;
        n_keys += 1;
        let hash = (&key[..]).read_u64::<LittleEndian>()?;

        match Datasets::from_slice(&value) {
            Err(e) => {
                n_bad += 1;
                println!(
                    "hash={hash} len={} parse_error={e} bytes=[{}]",
                    value.len(),
                    hex(&value)
                );
            }
            Ok(datasets) => {
                let out_of_range: Vec<u32> =
                    datasets.into_iter().filter(|&i| i >= n_datasets).collect();
                if !out_of_range.is_empty() {
                    n_bad += 1;
                    let bytes = if value.len() <= 64 {
                        hex(&value)
                    } else {
                        format!("{} ...", hex(&value[..64]))
                    };
                    println!(
                        "hash={hash} len={} out_of_range={out_of_range:?} bytes=[{bytes}]",
                        value.len()
                    );
                }
            }
        }
    }

    eprintln!("scanned {n_keys} keys, {n_bad} bad values");
    Ok(())
}
