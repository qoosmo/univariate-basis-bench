use criterion::{
    black_box,
    criterion_group,
    criterion_main,
    BenchmarkId,
    Criterion,
    Throughput,
};

use univariate_basis_bench::{
    field::F,
    kernel,
    monomial,
};

fn bench_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("single_point_eval_final");

    group.sample_size(50);
    group.warm_up_time(std::time::Duration::from_secs(5));
    group.measurement_time(std::time::Duration::from_secs(10));

    for log_n in [10usize, 12, 14, 16, 18, 20] {
        let n = 1usize << log_n;

        group.throughput(Throughput::Elements(n as u64));

        let monomial_coeffs: Vec<F> = (0..n)
            .map(|i| {
                F::from_canonical_u64(
                    (i as u64)
                        .wrapping_mul(17)
                        .wrapping_add(3),
                )
            })
            .collect();

        let kernel_coeffs = kernel::from_monomial(&monomial_coeffs);
        let x = F::from_canonical_u64(19);

        group.bench_with_input(
            BenchmarkId::new("monomial_horner", n),
            &n,
            |b, _| {
                b.iter(|| {
                    black_box(monomial::eval_horner(
                        black_box(&monomial_coeffs),
                        black_box(x),
                    ))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("monomial_parallel_blocks_4096", n),
            &n,
            |b, _| {
                b.iter(|| {
                    black_box(monomial::eval_parallel_blocks(
                        black_box(&monomial_coeffs),
                        black_box(x),
                        4096,
                    ))
                })
            },
        );

        group.bench_with_input(
            BenchmarkId::new("kernel_serial", n),
            &n,
            |b, _| {
                b.iter(|| {
                    black_box(kernel::eval_serial(
                        black_box(&kernel_coeffs),
                        black_box(x),
                    ))
                })
            },
        );

        // Only benchmark subtree parallelism when it actually performs
        // a decomposition into multiple subtrees.
        if n > 4096 {
            group.bench_with_input(
                BenchmarkId::new("kernel_subtree_4096", n),
                &n,
                |b, _| {
                    b.iter(|| {
                        black_box(kernel::eval_parallel_subtrees(
                            black_box(&kernel_coeffs),
                            black_box(x),
                            4096,
                        ))
                    })
                },
            );
        }
    }

    group.finish();
}

criterion_group!(benches, bench_eval);
criterion_main!(benches);
