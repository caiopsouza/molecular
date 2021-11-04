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
        let problem = molecular::load_problem("1ppt.nmr");
        let positions = molecular::solve(black_box(&problem));
        molecular::assert_solution(&problem, &positions);
    }));

    c.bench_function("2erl.nmr", |b| b.iter(|| {
        let problem = molecular::load_problem("2erl.nmr");
        molecular::solve(black_box(&problem));
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);