use std::fs::File;
use std::io;
use std::path::Path;
use std::rc::Rc;

use clap::{load_yaml, App};
use exitfailure::ExitFailure;
use failure::Error;
//use human_panic::setup_panic;
use lazy_init::Lazy;
use log::{info, LevelFilter};
use ocf::{get_output, CompressionFormat};
use serde::ser::SerializeStruct;
use serde::{Serialize, Serializer};

use sourmash::cmd::{
    count_unique, draff_compare, draff_index, draff_search, draff_signature, prepare,
};
use sourmash::index::linear::LinearIndex;
use sourmash::index::sbt::scaffold;
use sourmash::index::search::{
    search_minhashes, search_minhashes_containment, search_minhashes_find_best,
};
use sourmash::index::{Comparable, Dataset, Index, MHBT};
use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::Sketch;

struct Query<T> {
    data: T,
}

impl Query<Signature> {
    fn ksize(&self) -> u64 {
        // TODO: select the correct signature
        self.data.signatures[0].ksize() as u64
    }

    fn moltype(&self) -> String {
        // TODO: this might panic
        match &self.data.signatures[0] {
            Sketch::MinHash(mh) => {
                if mh.is_protein() {
                    "protein".into()
                } else {
                    "DNA".into()
                }
            }
            Sketch::UKHS(_) => {
                // TODO: draff only supports dna for now
                "DNA".into()
            }
        }
    }

    fn name(&self) -> String {
        self.data.name().clone()
    }
}

impl From<Query<Signature>> for Dataset<Signature> {
    fn from(other: Query<Signature>) -> Dataset<Signature> {
        let data = Lazy::new();
        data.get_or_create(|| other.data);

        Dataset::builder()
            .data(Rc::new(data))
            .filename("")
            .name("")
            .metadata("")
            .storage(None)
            .build()
    }
}

fn load_query_signature(
    query: &str,
    ksize: usize,
    moltype: Option<&str>,
    scaled: Option<u64>,
) -> Result<Query<Signature>, Error> {
    let mut reader = io::BufReader::new(File::open(query)?);
    let sigs = Signature::load_signatures(&mut reader, ksize, moltype, scaled)?;

    //dbg!(&sigs);
    // TODO: what if we have more than one left?
    let data = sigs[0].clone();

    Ok(Query { data })
}

struct Database {
    data: Indices,
    path: String,
}

enum Indices {
    MHBT(MHBT),
    LinearIndex(LinearIndex<Dataset<Signature>>),
}

impl Index for Database {
    type Item = Dataset<Signature>;

    fn find<F>(
        &self,
        search_fn: F,
        sig: &Self::Item,
        threshold: f64,
    ) -> Result<Vec<&Self::Item>, Error>
    where
        F: Fn(&dyn Comparable<Self::Item>, &Self::Item, f64) -> bool,
    {
        match &self.data {
            Indices::MHBT(data) => data.find(search_fn, sig, threshold),
            Indices::LinearIndex(data) => data.find(search_fn, sig, threshold),
        }
    }

    fn insert(&mut self, node: &Self::Item) -> Result<(), Error> {
        match &mut self.data {
            Indices::MHBT(data) => data.insert(node),
            Indices::LinearIndex(data) => data.insert(node),
        }
    }

    fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), Error> {
        match &self.data {
            Indices::MHBT(data) => data.save(path),
            Indices::LinearIndex(data) => data.save(path),
        }
    }

    fn load<P: AsRef<Path>>(_path: P) -> Result<(), Error> {
        unimplemented!();
    }

    fn datasets(&self) -> Vec<Self::Item> {
        match &self.data {
            Indices::MHBT(data) => data.datasets(),
            Indices::LinearIndex(data) => data.datasets(),
        }
    }
}

fn load_sbts_and_sigs(
    filenames: &[&str],
    query: &Query<Signature>,
    _containment: bool,
    traverse: bool,
) -> Result<Vec<Database>, Error> {
    let mut dbs = Vec::default();

    let _ksize = query.ksize();
    let _moltype = query.moltype();

    let n_signatures = 0;
    let mut n_databases = 0;

    for path in filenames {
        if traverse && Path::new(path).is_dir() {
            continue;
        }

        if let Ok(data) = MHBT::from_path(path) {
            // TODO: check compatible
            dbs.push(Database {
                data: Indices::MHBT(data),
                path: String::from(*path),
            });
            info!("loaded SBT {}", path);
            n_databases += 1;
            continue;
        } else if let Ok(data) = LinearIndex::<Dataset<Signature>>::from_path(path) {
            // TODO: check compatible
            dbs.push(Database {
                data: Indices::LinearIndex(data),
                path: String::from(*path),
            });
            info!("loaded LinearIndex {}", path);
            n_databases += 1;
            continue;
        }

        // TODO: load sig, need to change Database
        // IDEA: put sig into a LinearIndex, and replace Database with a Box<dyn Index>?
    }

    if n_signatures > 0 && n_databases > 0 {
        info!(
            "loaded {} signatures and {} databases total.",
            n_signatures, n_databases
        );
    } else if n_signatures > 0 {
        info!("loaded {} signatures.", n_signatures);
    } else if n_databases > 0 {
        info!("loaded {} databases.", n_databases);
    } else {
        return Err(failure::err_msg("Couldn't load any databases"));
    }

    Ok(dbs)
}

struct Results {
    similarity: f64,
    match_sig: Signature,
    db: String,
}

impl Serialize for Results {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut partial = serializer.serialize_struct("Results", 4)?;
        partial.serialize_field("similarity", &self.similarity)?;
        partial.serialize_field("name", &self.match_sig.name())?;
        partial.serialize_field("filename", &self.db)?;
        partial.serialize_field("md5", &self.match_sig.md5sum())?;
        partial.end()
    }
}

fn search_databases(
    query: Query<Signature>,
    databases: &[Database],
    threshold: f64,
    containment: bool,
    best_only: bool,
    _ignore_abundance: bool,
) -> Result<Vec<Results>, Error> {
    let mut results = Vec::default();

    let search_fn = if best_only {
        search_minhashes_find_best()
    } else if containment {
        search_minhashes_containment
    } else {
        search_minhashes
    };
    let query_leaf = query.into();

    // TODO: set up scaled for DB and query

    for db in databases {
        let matches = db.find(search_fn, &query_leaf, threshold).unwrap();
        for dataset in matches.into_iter() {
            let similarity = query_leaf.similarity(dataset);

            // should always be true, but... better safe than sorry.
            if similarity >= threshold {
                results.push(Results {
                    similarity,
                    match_sig: dataset.clone().into(),
                    db: db.path.clone(),
                })
            }
        }
    }

    results.sort_by(|a, b| b.similarity.partial_cmp(&a.similarity).unwrap());
    Ok(results)
}

fn main() -> Result<(), ExitFailure> {
    //setup_panic!();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let yml = load_yaml!("sourmash.yml");
    let m = App::from_yaml(yml).get_matches();

    match m.subcommand_name() {
        Some("draff") => {
            let cmd = m.subcommand_matches("draff").unwrap();
            let inputs = cmd
                .values_of("inputs")
                .map(|vals| vals.collect::<Vec<_>>())
                .unwrap();

            let ksize: usize = cmd.value_of("ksize").unwrap().parse().unwrap();
            let wsize: usize = cmd.value_of("wsize").unwrap().parse().unwrap();

            draff_signature(inputs, ksize, wsize)?;
        }
        Some("draff_search") => {
            let cmd = m.subcommand_matches("draff_search").unwrap();

            let index: &str = cmd.value_of("index").unwrap();
            let query: &str = cmd.value_of("query").unwrap();

            draff_search(index, query)?;
        }
        Some("draff_compare") => {
            let cmd = m.subcommand_matches("draff_compare").unwrap();
            let inputs = cmd
                .values_of("inputs")
                .map(|vals| vals.collect::<Vec<_>>())
                .unwrap();

            draff_compare(inputs)?;
        }
        Some("prepare") => {
            let cmd = m.subcommand_matches("prepare").unwrap();
            let index: &str = cmd.value_of("index").unwrap();

            prepare(index)?;
        }
        Some("index") => {
            let cmd = m.subcommand_matches("index").unwrap();
            let inputs = cmd
                .values_of("inputs")
                .map(|vals| vals.collect::<Vec<_>>())
                .unwrap();

            let output: &str = cmd.value_of("output").unwrap();

            draff_index(inputs, output)?;
        }
        Some("scaffold") => {
            let cmd = m.subcommand_matches("scaffold").unwrap();
            let sbt_file = cmd.value_of("current_sbt").unwrap();

            let sbt = MHBT::from_path(sbt_file)?;
            let mut new_sbt: MHBT = scaffold(sbt.datasets(), sbt.storage());

            new_sbt.save_file("test", None)?;

            assert_eq!(new_sbt.datasets().len(), sbt.datasets().len());
        }
        Some("count_unique") => {
            let cmd = m.subcommand_matches("count_unique").unwrap();

            let index: &str = cmd.value_of("index").unwrap();

            count_unique(index)?;
        }
        Some("search") => {
            let cmd = m.subcommand_matches("search").unwrap();

            if cmd.is_present("quiet") {
                log::set_max_level(LevelFilter::Warn);
            }

            let query = load_query_signature(
                cmd.value_of("query").unwrap(),
                if cmd.is_present("ksize") {
                    cmd.value_of("ksize").unwrap().parse().unwrap()
                } else {
                    // TODO default k
                    unimplemented!()
                },
                Some("dna"), // TODO: select moltype,
                if cmd.is_present("scaled") {
                    Some(cmd.value_of("scaled").unwrap().parse().unwrap())
                } else {
                    None
                },
            )?;

            info!(
                "loaded query: {}... (k={}, {})",
                query.name(),
                query.ksize(),
                query.moltype()
            );

            // TODO: check args.scaled and downsample

            let containment = cmd.is_present("containment");
            let traverse_directory = cmd.is_present("traverse-directory");
            let databases = load_sbts_and_sigs(
                &cmd.values_of("databases")
                    .map(|vals| vals.collect::<Vec<_>>())
                    .unwrap(),
                &query,
                containment,
                traverse_directory,
            )?;

            if databases.is_empty() {
                return Err(failure::err_msg("Nothing found to search!").into());
            }

            let best_only = cmd.is_present("best-only");
            let threshold = cmd.value_of("threshold").unwrap().parse().unwrap();
            let ignore_abundance = cmd.is_present("ignore-abundance");
            let results = search_databases(
                query,
                &databases,
                threshold,
                containment,
                best_only,
                ignore_abundance,
            )?;

            let num_results = if best_only {
                1
            } else {
                cmd.value_of("num-results").unwrap().parse().unwrap()
            };

            let n_matches = if num_results == 0 || results.len() <= num_results {
                println!("{} matches:", results.len());
                results.len()
            } else {
                println!("{} matches; showing first {}:", results.len(), num_results);
                num_results
            };

            println!("similarity   match");
            println!("----------   -----");
            for sr in &results[..n_matches] {
                println!(
                    "{:>5.1}%       {:60}",
                    sr.similarity * 100.,
                    sr.match_sig.name()
                );
            }

            if best_only {
                info!("** reporting only one match because --best-only was set")
            }

            if let Some(output) = cmd.value_of("output") {
                let mut wrt = csv::Writer::from_path(output)?;

                for sr in &results[..n_matches] {
                    wrt.serialize(sr)?;
                }
                wrt.flush()?;
            };

            if let Some(outname) = cmd.value_of("save-matches") {
                let writer = get_output(outname, CompressionFormat::No)?;

                info!("saving all matched signatures to \"{}\"", outname);

                let sigs: Vec<Signature> = results.into_iter().map(|sr| sr.match_sig).collect();
                serde_json::to_writer(writer, &sigs)?;
            }
        }
        _ => {
            println!("{:?}", m);
        }
    }
    Ok(())
}
