#[macro_use]
extern crate criterion;

use std::fs::File;
use std::io::BufReader;
use std::path::PathBuf;

use sourmash::signature::{Signature, SigsTrait};
use sourmash::sketch::minhash::{KmerMinHash, KmerMinHashBTree};
use sourmash::sketch::Sketch;

use criterion::Criterion;

fn intersection(c: &mut Criterion) {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/gather-abund/genome-s10.fa.gz.sig");
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");
    let mh = if let Sketch::MinHash(mh) = &sigs.swap_remove(0).sketches()[0] {
        mh.clone()
    } else {
        unimplemented!()
    };

    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("../../tests/test-data/gather-abund/genome-s11.fa.gz.sig");
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut sigs: Vec<Signature> = serde_json::from_reader(reader).expect("Loading error");
    let mh2 = if let Sketch::MinHash(mh) = &sigs.swap_remove(0).sketches()[0] {
        mh.clone()
    } else {
        unimplemented!()
    };

    let mut group = c.benchmark_group("minhash");
    group.sample_size(10);

    group.bench_function("intersection", |b| {
        b.iter(|| {
            mh.intersection(&mh2).unwrap();
        });
    });

    group.bench_function("intersection_size", |b| {
        b.iter(|| {
            mh.intersection_size(&mh2).unwrap();
        });
    });

    let mut mh1 = KmerMinHash::builder()
        .num(0)
        .max_hash(1_000_000)
        .ksize(21)
        .build();
    let mut mh2 = KmerMinHash::builder()
        .num(0)
        .max_hash(1_000_000)
        .ksize(21)
        .build();

    let mut mh1_btree = KmerMinHashBTree::builder()
        .num(0)
        .max_hash(1_000_000)
        .ksize(21)
        .build();
    let mut mh2_btree = KmerMinHashBTree::builder()
        .num(0)
        .max_hash(1_000_000)
        .ksize(21)
        .build();

    for i in 0..=1_000_000 {
        if i % 2 == 0 {
            mh1.add_hash(i);
            mh1_btree.add_hash(i);
        }
        if i % 45 == 0 {
            mh2.add_hash(i);
            mh2_btree.add_hash(i);
        }
    }

    group.bench_function("large intersection", |b| {
        b.iter(|| {
            mh1.intersection(&mh2).unwrap();
        });
    });

    group.bench_function("large intersection_size", |b| {
        b.iter(|| {
            mh1.intersection_size(&mh2).unwrap();
        });
    });

    group.bench_function("large intersection btree", |b| {
        b.iter(|| {
            mh1_btree.intersection(&mh2_btree).unwrap();
        });
    });

    group.bench_function("large intersection_size btree", |b| {
        b.iter(|| {
            mh1_btree.intersection_size(&mh2_btree).unwrap();
        });
    });
}

criterion_group!(minhash, intersection);
criterion_main!(minhash);
