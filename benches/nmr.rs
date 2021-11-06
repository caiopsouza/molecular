use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("test.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("test.nmr"));
    }));

    c.bench_function("10-3.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("10-3.nmr"));
    }));

    c.bench_function("1ppt.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("1ppt.nmr"));
    }));

    c.bench_function("2erl.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("2erl.nmr"));
    }));

    c.bench_function("1fs3.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("1fs3.nmr"));
    }));

    c.bench_function("2e7z.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("2e7z.nmr"));
    }));

    c.bench_function("7nyz.nmr", |b| b.iter(|| {
        molecular::load_solve_and_format(black_box("7nyz.nmr"));
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);