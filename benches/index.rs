#[macro_use]
extern crate criterion;

use std::path::PathBuf;

use criterion::{Bencher, Criterion, Fun};
use sourmash::index::bigsi::BIGSI;
use sourmash::index::linear::LinearIndex;
use sourmash::index::storage::ReadData;
use sourmash::index::MHBT;
use sourmash::index::{Dataset, Index};
use sourmash::signature::Signature;

fn find_small_bench(c: &mut Criterion) {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/test-data/v5.sbt.json");

    let sbt = MHBT::from_path(filename).expect("Loading error");

    let leaf: Dataset<Signature> = (*sbt.datasets().first().unwrap()).clone();

    let mut linear = LinearIndex::builder().storage(sbt.storage()).build();

    for l in &sbt.datasets() {
        linear.insert(l);
    }

    let mut bigsi = BIGSI::new(10000, 10);
    for l in &sbt.datasets() {
        let data = l.data().unwrap();
        bigsi.insert(data);
    }

    let sbt_find = Fun::new(
        "sbt_search",
        move |b: &mut Bencher, leaf: &Dataset<Signature>| b.iter(|| sbt.search(leaf, 0.1, false)),
    );

    let linear_find = Fun::new(
        "linear_search",
        move |b: &mut Bencher, leaf: &Dataset<Signature>| {
            b.iter(|| linear.search(leaf, 0.1, false))
        },
    );

    let bigsi_find = Fun::new(
        "bigsi_search",
        move |b: &mut Bencher, leaf: &Dataset<Signature>| {
            let data = leaf.data().unwrap();
            b.iter(|| bigsi.search(data, 0.1, false))
        },
    );

    let functions = vec![sbt_find, linear_find, bigsi_find];
    c.bench_functions("search_small", functions, leaf);
}

fn find_subset_bench(c: &mut Criterion) {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/test-data/subset.sbt.json");

    let sbt = MHBT::from_path(filename).expect("Loading error");

    let leaf: Dataset<Signature> = (*sbt.datasets().first().unwrap()).clone();

    let mut linear = LinearIndex::builder().storage(sbt.storage()).build();
    for l in &sbt.datasets() {
        linear.insert(l);
    }

    let mut bigsi = BIGSI::new(10000, 10);
    for l in &sbt.datasets() {
        let data = l.data().unwrap();
        bigsi.insert(data);
    }

    let sbt_find = Fun::new(
        "sbt_search",
        move |b: &mut Bencher, leaf: &Dataset<Signature>| b.iter(|| sbt.search(leaf, 0.1, false)),
    );

    let linear_find = Fun::new(
        "linear_search",
        move |b: &mut Bencher, leaf: &Dataset<Signature>| {
            b.iter(|| linear.search(leaf, 0.1, false))
        },
    );

    let bigsi_find = Fun::new(
        "bigsi_search",
        move |b: &mut Bencher, leaf: &Dataset<Signature>| {
            let data = leaf.data().unwrap();
            b.iter(|| bigsi.search(data, 0.1, false))
        },
    );

    let functions = vec![sbt_find, linear_find, bigsi_find];
    c.bench_functions("search_subset", functions, leaf);
}

criterion_group!(benches, find_small_bench, find_subset_bench);
criterion_main!(benches);
