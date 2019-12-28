use std::path::{Path, PathBuf};
use std::rc::Rc;

use failure::Error;
use log::info;
use needletail::parse_sequence_path;

use sourmash::index::linear::LinearIndex;
use sourmash::index::storage::{FSStorage, Storage};
use sourmash::index::{Comparable, Index, MHBT};
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::Sketch;

use crate::index::UKHSTree;
use crate::sketch::ukhs::{FlatUKHS, UKHSTrait, UniqueUKHS};

pub fn draff_index(sig_files: Vec<&str>, outfile: &str) -> Result<(), Error> {
    let storage: Rc<dyn Storage> = Rc::new(
        FSStorage::new(".", ".draff"), // TODO: use outfile
    );

    //let mut index = UKHSTree::builder().storage(Rc::clone(&storage)).build();
    let mut index = LinearIndex::<Signature>::builder()
        .storage(Rc::clone(&storage))
        .build();

    for filename in sig_files {
        // TODO: check for stdin? can also use get_input()?

        let ukhs_sig = FlatUKHS::load(&filename)?;

        let path = Path::new(filename);
        let basename = path.file_name().unwrap().to_str().unwrap().to_owned();

        let sig = Signature::builder()
            .hash_function("nthash") // TODO: spec!
            .class("draff_signature") // TODO: spec!
            .name(Some(basename.to_owned()))
            .filename(Some(basename.to_owned()))
            .signatures(vec![Sketch::UKHS(ukhs_sig)])
            .build();

        index.insert(sig)?;
    }

    // TODO: implement to_writer and use this?
    //let mut output = get_output(outfile, CompressionFormat::No)?;
    //index.to_writer(&mut output)?

    index.save_file(outfile, None)
}

pub fn draff_compare(sigs: Vec<&str>) -> Result<(), Error> {
    let mut dists = vec![vec![0.; sigs.len()]; sigs.len()];
    let loaded_sigs: Vec<FlatUKHS> = sigs.iter().map(|s| FlatUKHS::load(s).unwrap()).collect();

    for (i, sig1) in loaded_sigs.iter().enumerate() {
        for (j, sig2) in loaded_sigs.iter().enumerate() {
            if i >= j {
                dists[i][j] = 1. - sig1.distance(sig2);
                dists[j][i] = dists[i][j];
            }
        }
    }

    for row in dists {
        println!("{:.2?}", row);
    }

    Ok(())
}

pub fn draff_search(index: &str, query: &str) -> Result<(), Error> {
    let index = UKHSTree::from_path(index)?;

    let ukhs_sig = FlatUKHS::load(&query)?;

    let sig = Signature::builder()
        .hash_function("nthash") // TODO: spec!
        .class("draff_signature") // TODO: spec!
        .name(Some(query.to_owned()))
        .filename(Some(query.to_owned()))
        .signatures(vec![Sketch::UKHS(ukhs_sig)])
        .build();

    for found in index.search(&sig, 0.9, false)? {
        println!("{:.2}: {:?}", sig.similarity(found), found);
    }

    Ok(())
}

pub fn draff_signature(files: Vec<&str>, k: usize, w: usize) -> Result<(), Error> {
    for filename in files {
        // TODO: check for stdin?

        let mut ukhs = UniqueUKHS::new(k, w)?;

        info!("Build signature for {} with W={}, K={}...", filename, w, k);

        parse_sequence_path(
            filename,
            |_| {},
            |record| {
                // if there is anything other than ACGT in sequence,
                // it is replaced with A.
                // This matches khmer and screed behavior
                //
                // NOTE: sourmash is different! It uses the force flag to drop
                // k-mers that are not ACGT
                let seq: Vec<u8> = record
                    .seq
                    .iter()
                    .map(|&x| match x as char {
                        'A' | 'C' | 'G' | 'T' => x,
                        'a' | 'c' | 'g' | 't' => x.to_ascii_uppercase(),
                        _ => b'A',
                    })
                    .collect();

                ukhs.add_sequence(&seq, false)
                    .expect("Error adding sequence");
            },
        )?;

        let mut outfile = PathBuf::from(filename);
        outfile.set_extension("sig");

        /*
        let mut output = get_output(outfile.to_str().unwrap(), CompressionFormat::No)?;

        let flat: FlatUKHS = ukhs.into();
        flat.to_writer(&mut output)?
        */
    }
    info!("Done.");

    Ok(())
}
