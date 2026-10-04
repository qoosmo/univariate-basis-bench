use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rand::{RngCore, SeedableRng};
use rand::rngs::StdRng;
use univariate_basis_bench::{field::F, kernel, lagrange, monomial};

fn random_vec(n: usize, seed: u64) -> Vec<F> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n).map(|_| F::new(rng.next_u64() as u128)).collect()
}

fn bench_eval(c: &mut Criterion) {
    let mut group = c.benchmark_group("single_point_eval");
    group.sample_size(20);

    // Keep initial CI/runtime manageable. Extend to 2^24 on the target machine.
    for log_n in [10u32, 12, 14, 16, 18, 20] {
        let n = 1usize << log_n;
        group.throughput(Throughput::Elements(n as u64));

        let coeffs = random_vec(n, 0xA11CE + n as u64);
        let kernel_coeffs = random_vec(n, 0xBEEF + n as u64);
        let lagrange_values = random_vec(n, 0xCAFE + n as u64);
        let domain = lagrange::roots_of_unity_domain(n);
        let x = F::from_canonical_u64(123456789);

        group.bench_with_input(BenchmarkId::new("monomial_horner", n), &n, |b, _| {
            b.iter(|| black_box(monomial::eval_horner(black_box(&coeffs), black_box(x))))
        });

        group.bench_with_input(BenchmarkId::new("monomial_parallel_blocks_4096", n), &n, |b, _| {
            b.iter(|| black_box(monomial::eval_parallel_blocks(black_box(&coeffs), black_box(x), 4096)))
        });

        group.bench_with_input(BenchmarkId::new("lagrange_barycentric", n), &n, |b, _| {
            b.iter(|| black_box(lagrange::eval_roots_of_unity(black_box(&lagrange_values), black_box(&domain), black_box(x))))
        });

        group.bench_with_input(BenchmarkId::new("lagrange_barycentric_parallel", n), &n, |b, _| {
            b.iter(|| black_box(lagrange::eval_roots_of_unity_parallel(black_box(&lagrange_values), black_box(&domain), black_box(x))))
        });

        group.bench_with_input(BenchmarkId::new("kernel_serial", n), &n, |b, _| {
            b.iter(|| black_box(kernel::eval_serial(black_box(&kernel_coeffs), black_box(x))))
        });

        group.bench_with_input(BenchmarkId::new("kernel_parallel", n), &n, |b, _| {
            b.iter(|| black_box(kernel::eval_parallel(black_box(&kernel_coeffs), black_box(x), 2048)))
        });
    }

    group.finish();
}

criterion_group!(benches, bench_eval);
criterion_main!(benches);
