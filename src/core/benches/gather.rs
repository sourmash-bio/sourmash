#[macro_use]
extern crate criterion;

use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::sketch::Sketch;
use sourmash::index::{FastGatherResult, calculate_gather_stats};
use sourmash::encodings::HashFunctions;
use sourmash::cmd::ComputeParameters;

use criterion::Criterion;

fn gather_stats_benchmarks(c: &mut Criterion) {
    let scaled = 10;
    let params = ComputeParameters::builder()
        .ksizes(vec![3])
        .scaled(scaled)
        .build();

    let mut match_sig = Signature::from_params(&params);
    // create two minhash
    let mut match_mh = KmerMinHash::new(scaled, 3, HashFunctions::Murmur64Dna, 42, true, 0);
    match_mh.add_hash(1);
    match_mh.add_hash(2);
    match_mh.add_hash(3);
    match_mh.add_hash(5);

    match_sig.reset_sketches();
    match_sig.push(Sketch::MinHash(match_mh.clone()));
    match_sig.set_filename("match-filename");
    match_sig.set_name("match-name");

    eprintln!("num_sketches: {:?}", match_sig.size());
    eprintln!("match_md5: {:?}", match_sig.md5sum());

    // Setup orig_query minhash with abundances and non-matching hash
    let mut orig_query = KmerMinHash::new(scaled, 3, HashFunctions::Murmur64Dna, 42, true, 0);
    orig_query.add_hash_with_abundance(1, 1);
    orig_query.add_hash_with_abundance(2, 3);
    orig_query.add_hash_with_abundance(4, 6); // Non-matching hash
    orig_query.add_hash_with_abundance(5, 9);
    orig_query.add_hash_with_abundance(6, 1);

    let mut query = orig_query.clone();
    let rm_hashes = vec![1];
    query.remove_many(rm_hashes.as_slice()).unwrap(); // remove hash 1
    // here define all of the options we want to test!
    // in particular, calc_ani and calc_ani_ci

    // set up single FastGatherResult
    let fgres = FastGatherResult::builder()
            .orig_query(orig_query)
            .query(query)
            .match_(match_sig.into())
            .match_mh(match_mh.clone())
            .match_size(2) // 2  -- only 2 hashes match, one was previously consumed
            .remaining_hashes(30) // arbitrary
            .gather_result_rank(5) // arbitrary
            // .total_orig_query_abund(20) // sum of orig_query abundances
            .build();
    // set up all the options we want to pass in
    let inputs = [(fgres, false, false)];//, fgres.clone(), true, false];

    let mut group = c.benchmark_group("gather_stats");
    
    for (i, input) in inputs.iter().enumerate() {
        group.bench_function(format!("calculate_gather_stats_{}", i), |b| {
            b.iter(|| {
                // to do -- need to get around this clone or impl clone for FastGatherResult
                calculate_gather_stats(input.0.clone(), input.1, input.2);
            });
        });
    }
    group.finish();
}

criterion_group!(gather, gather_stats_benchmarks);
criterion_main!(gather);