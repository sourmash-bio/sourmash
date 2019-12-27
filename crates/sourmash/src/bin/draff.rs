use clap::{load_yaml, App};
use exitfailure::ExitFailure;
//use human_panic::setup_panic;
use sourmash::cmd::{draff_compare, draff_index, draff_search, draff_signature};

fn main() -> Result<(), ExitFailure> {
    //setup_panic!();

    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let yml = load_yaml!("draff.yml");
    let m = App::from_yaml(yml).get_matches();

    match m.subcommand_name() {
        Some("compute") => {
            let cmd = m.subcommand_matches("compute").unwrap();
            let inputs = cmd
                .values_of("inputs")
                .map(|vals| vals.collect::<Vec<_>>())
                .unwrap();

            let ksize: usize = cmd.value_of("ksize").unwrap().parse().unwrap();
            let wsize: usize = cmd.value_of("wsize").unwrap().parse().unwrap();

            draff_signature(inputs, ksize, wsize)?;
        }
        Some("search") => {
            let cmd = m.subcommand_matches("search").unwrap();

            let index: &str = cmd.value_of("index").unwrap();
            let query: &str = cmd.value_of("query").unwrap();

            draff_search(index, query)?;
        }
        Some("compare") => {
            let cmd = m.subcommand_matches("compare").unwrap();
            let inputs = cmd
                .values_of("inputs")
                .map(|vals| vals.collect::<Vec<_>>())
                .unwrap();

            draff_compare(inputs)?;
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
        _ => {
            println!("{:?}", m);
        }
    }
    Ok(())
}
