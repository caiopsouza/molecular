use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("1ppt.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("1ppt.nmr"));
    }));

    c.bench_function("2erl.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("2erl.nmr"));
    }));

    c.bench_function("1fs3.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("1fs3.nmr"));
    }));

    c.bench_function("1mqq.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("1mqq.nmr"));
    }));

    c.bench_function("1m40.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("1m40.nmr"));
    }));

    c.bench_function("1rwh.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("1rwh.nmr"));
    }));

    c.bench_function("3b34.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("3b34.nmr"));
    }));

    c.bench_function("2e7z.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("2e7z.nmr"));
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);