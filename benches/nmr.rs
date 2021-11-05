use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("test", |b| b.iter(|| {
        let problem = molecular::load_problem("test.nmr");
        molecular::solve(black_box(&problem));
    }));

    c.bench_function("10-3.nmr", |b| b.iter(|| {
        let problem = molecular::load_problem("10-3.nmr");
        molecular::solve(black_box(&problem));
    }));

    c.bench_function("1ppt.nmr", |b| b.iter(|| {
        molecular::run(molecular::solve_1ppt);
    }));

    c.bench_function("2erl.nmr", |b| b.iter(|| {
        molecular::run(molecular::solve_2erl);
    }));

    c.bench_function("1fs3.nmr", |b| b.iter(|| {
        molecular::run(molecular::solve_1fs3);
    }));

    c.bench_function("2e7z.nmr", |b| b.iter(|| {
        molecular::run(molecular::solve_2e7z);
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);