#[macro_use]
extern crate criterion;

use std::path::PathBuf;

use criterion::{Bencher, Criterion, Fun};
use sourmash::index::linear::LinearIndexBuilder;
use sourmash::index::nodegraph::Nodegraph;
use sourmash::index::sbt::{Node, MHBT, SBT};
use sourmash::index::search::search_minhashes;
use sourmash::index::{Index, Leaf};
use sourmash::Signature;

fn find_small_bench(c: &mut Criterion) {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/test-data/v5.sbt.json");

    let sbt: MHBT = SBT::from_path(filename).expect("Loading error");

    let leaf: Leaf<Signature> = (*sbt.leaves().first().unwrap()).clone();

    let mut linear = LinearIndexBuilder::default()
        .storage(sbt.storage())
        .build()
        .unwrap();
    for l in &sbt.leaves() {
        linear.insert(l);
    }

    let sbt_find = Fun::new(
        "sbt_find",
        move |b: &mut Bencher, leaf: &Leaf<Signature>| {
            b.iter(|| sbt.find(search_minhashes, leaf, 0.1))
        },
    );

    let linear_find = Fun::new(
        "linear_find",
        move |b: &mut Bencher, leaf: &Leaf<Signature>| {
            b.iter(|| linear.find(search_minhashes, leaf, 0.1))
        },
    );

    let functions = vec![sbt_find, linear_find];
    c.bench_functions("find_small", functions, leaf);
}

fn find_subset_bench(c: &mut Criterion) {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/test-data/subset.sbt.json");

    let sbt: MHBT = SBT::from_path(filename).expect("Loading error");

    let leaf: Leaf<Signature> = (*sbt.leaves().first().unwrap()).clone();

    let mut linear = LinearIndexBuilder::default()
        .storage(sbt.storage())
        .build()
        .unwrap();
    for l in &sbt.leaves() {
        linear.insert(l);
    }

    let sbt_find = Fun::new(
        "sbt_find",
        move |b: &mut Bencher, leaf: &Leaf<Signature>| {
            b.iter(|| sbt.find(search_minhashes, leaf, 0.1))
        },
    );

    let linear_find = Fun::new(
        "linear_find",
        move |b: &mut Bencher, leaf: &Leaf<Signature>| {
            b.iter(|| linear.find(search_minhashes, leaf, 0.1))
        },
    );

    let functions = vec![sbt_find, linear_find];
    c.bench_functions("find_subset", functions, leaf);
}

criterion_group!(benches, find_small_bench, find_subset_bench);
criterion_main!(benches);
