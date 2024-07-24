use std::fs::File;
use std::io::{BufWriter, Cursor, Read};

use sourmash::sketch::nodegraph::Nodegraph;

use codspeed_criterion_compat::{black_box, criterion_group, criterion_main, Criterion};

fn save_load(c: &mut Criterion) {
    let mut data: Vec<u8> = vec![];
    let mut f = File::open("../../tests/test-data/.sbt.v3/internal.0").unwrap();
    let _ = f.read_to_end(&mut data);

    let mut group = c.benchmark_group("nodegraph");
    group.sample_size(10);

    let mut reader = Cursor::new(data.clone());
    let ng = Nodegraph::from_reader(&mut reader).unwrap();

    group.bench_function("load nodegraph", |b| {
        b.iter(|| {
            let mut reader = Cursor::new(data.clone());
            let _ng = Nodegraph::from_reader(&mut reader).unwrap();
        });
    });

    group.bench_function("save nodegraph", |b| {
        b.iter(|| {
            let mut buf = Vec::new();
            let mut writer = BufWriter::new(&mut buf);
            ng.save_to_writer(&mut writer).unwrap();
        });
    });

    group.bench_function("save compressed nodegraph", |b| {
        b.iter(|| {
            let mut buf = Vec::new();
            let mut writer = niffler::get_writer(
                Box::new(&mut buf),
                niffler::compression::Format::Gzip,
                niffler::compression::Level::One,
            )
            .unwrap();

            ng.save_to_writer(&mut writer).unwrap();
        });
    });
}

criterion_group!(nodegraph, save_load);
criterion_main!(nodegraph);
