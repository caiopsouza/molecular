use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn criterion_benchmark(c: &mut Criterion) {
    let (problem, positions, cumulative_torsion_3) = molecular::load_problem("test.nmr");
    c.bench_function("test", |b| b.iter(|| {
        let res = molecular::solve(black_box(&problem), black_box(&positions), black_box(&cumulative_torsion_3));
    }));

    let (problem, positions, cumulative_torsion_3) = molecular::load_problem("10-3.nmr");
    c.bench_function("10-3.nmr", |b| b.iter(|| {
        let res = molecular::solve(black_box(&problem), black_box(&positions), black_box(&cumulative_torsion_3));
    }));

    let (problem, positions, cumulative_torsion_3) = molecular::load_problem("1ppt.nmr");
    c.bench_function("1ppt.nmr", |b| b.iter(|| {
        let res = molecular::solve(black_box(&problem), black_box(&positions), black_box(&cumulative_torsion_3));
        molecular::assert_solution(&problem, &res);
    }));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);