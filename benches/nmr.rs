use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

fn criterion_benchmark(c: &mut Criterion) {
    let instances = ["1ppt", "2erl", "1crn", "1jk2", "1pht", "1a70", "1fs3", "1hoe", "1poa", "1mbn", "1ptq", "1m40", "1n4w", "1mqq", "1bpm", "3b34", "2e7z", "1rwh", "1rgs"];

    let mut group = c.benchmark_group("molecular");
    group.sample_size(10);

    // Exact
    /*for file in instances.iter() {
        group.bench_with_input(BenchmarkId::new("exact", file), file,
                               |b, file| b.iter(|| molecular::solve_exact(file)));
    }*/

    // Heuristics: Greedy
    /*for file in instances.iter() {
        group.bench_with_input(BenchmarkId::new("greedy_one", file), file,
                               |b, file| b.iter(|| molecular::solve_greedy_one(file)));
    }

    for file in instances.iter() {
        group.bench_with_input(BenchmarkId::new("greedy_two", file), file,
                               |b, file| b.iter(|| molecular::solve_greedy_two(file)));
    }

    for file in instances.iter() {
        group.bench_with_input(BenchmarkId::new("greedy_worst_one", file), file,
                               |b, file| b.iter(|| molecular::solve_greedy_worst_one(file)));
    }*/

    // Heuristics: Local Search
    /*for file in instances.iter() {
        group.bench_with_input(BenchmarkId::new("local_forward", file), file,
                               |b, file| b.iter(|| molecular::solve_local_search_forward(file)));
    }

    for file in instances.iter() {
        group.bench_with_input(BenchmarkId::new("local_backward", file), file,
                               |b, file| b.iter(|| molecular::solve_local_search_backward(file)));
    }*/

    // Heuristics: Iterated Local Search
    for file in instances.iter() {
        group.bench_with_input(BenchmarkId::new("ils", file), file,
                               |b, file| b.iter(|| molecular::solve_iterated_local_search(file)));
    }

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);