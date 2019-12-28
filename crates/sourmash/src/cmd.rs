use failure::Error;

use crate::index::MHBT;

/* FIXME bring back after succint-rs changes
pub fn count_unique(index_path: &str) -> Result<(), Error> {
    let index = MHBT::from_path(index_path)?;

    info!("Loaded index: {}", index_path);

    let mut hll = pdatastructs::hyperloglog::HyperLogLog::new(16);

    let mut total_hashes = 0u64;

    for (n, sig) in index.signatures().iter().enumerate() {
        if n % 1000 == 0 {
            info!("Processed {} signatures", n);
            info!("Unique hashes in {}: {}", index_path, hll.count());
            info!("Total hashes in {}: {}", index_path, total_hashes);
        };
        if let Sketch::MinHash(mh) = &sig.signatures[0] {
            for hash in mh.mins() {
                hll.add(&hash);
                total_hashes += 1;
            }
        }
    }

    info!("Unique k-mers in {}: {}", index_path, hll.count());
    info!("Total hashes in {}: {}", index_path, total_hashes);

    Ok(())
}
*/

pub fn prepare(index_path: &str) -> Result<(), Error> {
    let mut index = MHBT::from_path(index_path)?;

    // TODO equivalent to fill_internal in python
    //unimplemented!();

    index.save_file(index_path, None)?;

    Ok(())
}
