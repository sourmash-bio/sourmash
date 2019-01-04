use std::fs::File;
use std::io;
use std::path::Path;
use std::rc::Rc;

use clap::{load_yaml, App};
use exitfailure::ExitFailure;
use failure::{Error, ResultExt};
use human_panic::setup_panic;
use log::{debug, error, info, LevelFilter};

use sourmash::index::nodegraph::Nodegraph;
use sourmash::index::sbt::{scaffold, Node, MHBT, SBT};
use sourmash::index::search::search_minhashes;
use sourmash::index::{Index, Leaf, LeafBuilder};
use sourmash::Signature;

struct Query<T> {
    data: T,
}

impl Query<Signature> {
    fn ksize(&self) -> u64 {
        // TODO: this might panic
        self.data.signatures[0].ksize as u64
    }

    fn moltype(&self) -> String {
        // TODO: this might panic
        if self.data.signatures[0].is_protein {
            "protein".into()
        } else {
            "DNA".into()
        }
    }

    fn name(&self) -> String {
        self.data.name.clone().unwrap()
    }
}

impl From<Query<Signature>> for Leaf<Signature> {
    fn from(other: Query<Signature>) -> Leaf<Signature> {
        let leaf = LeafBuilder::default().build().unwrap();
        //leaf.data.get_or_create(|| data.query);
        leaf
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

    debug!("{:?}", sigs);
    // TODO: what if we have more than one left?
    let data = sigs[0].clone();

    Ok(Query { data })
}

struct Database {
    data: MHBT,
    path: String,
    is_index: bool,
}

fn load_sbts_and_sigs(
    filenames: &[&str],
    query: &Query<Signature>,
    containment: bool,
    traverse: bool,
) -> Result<Vec<Database>, Error> {
    let mut dbs = Vec::default();

    let ksize = query.ksize();
    let moltype = query.moltype();

    let mut n_signatures = 0;
    let mut n_databases = 0;

    for path in filenames {
        if traverse && Path::new(path).is_dir() {
            continue;
        }

        if let Ok(data) = MHBT::from_path(path) {
            // TODO: check compatible
            dbs.push(Database {
                data,
                path: String::from(*path),
                is_index: true,
            });
            info!("loaded SBT {}", path);
            n_databases += 1;
            continue;
        }

        // TODO: load sig, need to change Database
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
}

fn search_databases(
    query: Query<Signature>,
    databases: &[Database],
    threshold: f64,
    containment: bool,
    best_only: bool,
    ignore_abundance: bool,
) -> Result<Vec<Results>, Error> {
    let mut results = Vec::default();

    let search_fn = search_minhashes;
    let query_leaf = query.into();

    for db in databases {
        for dataset in db.data.find(search_fn, &query_leaf, threshold) {}
    }

    Ok(results)
}

fn main() -> Result<(), ExitFailure> {
    //setup_panic!();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let yml = load_yaml!("sourmash.yml");
    let m = App::from_yaml(yml).get_matches();

    match m.subcommand_name() {
        Some("scaffold") => {
            let cmd = m.subcommand_matches("scaffold").unwrap();
            let sbt_file = cmd.value_of("current_sbt").unwrap();

            let sbt = MHBT::from_path(sbt_file)?;
            let new_sbt: MHBT = scaffold(sbt.leaves());

            assert_eq!(new_sbt.leaves().len(), 100);
            Ok(())
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
                    0
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

            if databases.len() == 0 {
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

            let n_matches = if num_results == 0 || results.len() < num_results {
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
                    "{:>6.1}%       {:60}",
                    sr.similarity * 100.,
                    sr.match_sig.name.clone().unwrap()
                );
            }

            if best_only {
                info!("** reporting only one match because --best-only was set")
            }

            /*
            if args.output:
                fieldnames = ['similarity', 'name', 'filename', 'md5']
                w = csv.DictWriter(args.output, fieldnames=fieldnames)
                w.writeheader()

                for sr in &results:
                    d = dict(sr._asdict())
                    del d['match_sig']
                    w.writerow(d)

            if args.save_matches:
                outname = args.save_matches.name
                info!("saving all matched signatures to \"{}\"", outname)
                Signature::save_signatures([sr.match_sig for sr in results], args.save_matches)
            */

            Ok(())
        }
        _ => {
            println!("{:?}", m);
            Ok(())
        }
    }
}
