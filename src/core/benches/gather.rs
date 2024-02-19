use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use sourmash::collection::Collection;
use sourmash::signature::Signature;
use sourmash::sketch::Sketch;
use sourmash::{index::calculate_gather_stats, storage::SigStore};

use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn gather_stats_benchmarks(c: &mut Criterion) {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/track_abund/47.fa.sig");
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");
    let orig_query = if let Sketch::MinHash(mh) = &sigs.swap_remove(0).sketches()[0] {
        mh.clone()
    } else {
        unimplemented!()
    };
    let query = orig_query.clone();

    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/track_abund/63.fa.sig");
    // load collection to get sig in sigstore
    let signatures = Signature::from_path(filename).expect("cant find file");
    let collection = Collection::from_sigs(signatures).expect("cant make collection");
    let match_sig: SigStore = collection.sig_for_dataset(0).expect("cant load sig");
    let test_cases = vec![(false, false), (true, false), (false, true), (true, true)];

    let mut group = c.benchmark_group("gather_stats");
    for (calc_abund_stats, calc_ani_ci) in test_cases {
        let test_name = format!(
            "abund{}_ani_ci{}",
            calc_abund_stats as u8, calc_ani_ci as u8
        );
        group.bench_function(&test_name, |b| {
            b.iter(|| {
                calculate_gather_stats(
                    black_box(&orig_query),
                    black_box(query.clone()),
                    black_box(match_sig.clone()),
                    black_box(42), // Example match_size
                    black_box(1),  // Example gather_result_rank
                    black_box(200),
                    black_box(calc_abund_stats),
                    black_box(calc_ani_ci),
                    black_box(None), // don't set custom confidence intervals
                )
                .expect("error calculating gather stats");
            });
        });
    }
    group.finish();
}

criterion_group!(gather, gather_stats_benchmarks);
criterion_main!(gather);
