#[macro_use]
extern crate criterion;

use std::fs::File;
use std::io::{BufReader, Cursor, Read};
use std::path::PathBuf;

use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
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
}

criterion_group!(minhash, intersection);
criterion_main!(minhash);
