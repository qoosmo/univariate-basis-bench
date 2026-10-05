use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};

use univariate_basis_bench::{field::F, kernel, lagrange, monomial};

fn ntt_evaluate_monomial(coeffs: &[F]) -> Vec<F> {
    assert!(coeffs.len().is_power_of_two());

    let n = coeffs.len();
    let mut a = coeffs.to_vec();

    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            a.swap(i, j);
        }
    }

    let mut len = 2usize;
    while len <= n {
        let w_len = F::primitive_root_of_unity(len);

        for start in (0..n).step_by(len) {
            let mut w = F::ONE;
            let half = len / 2;

            for offset in 0..half {
                let u = a[start + offset];
                let v = a[start + offset + half] * w;

                a[start + offset] = u + v;
                a[start + offset + half] = u - v;

                w *= w_len;
            }
        }

        len <<= 1;
    }

    a
}

fn bench_final_bases(c: &mut Criterion) {
    let mut group = c.benchmark_group("final_basis_comparison");

    group.sample_size(50);
    group.warm_up_time(std::time::Duration::from_secs(5));
    group.measurement_time(std::time::Duration::from_secs(10));

    for log_n in [10usize, 12, 14, 16, 18, 20] {
        let n = 1usize << log_n;

        group.throughput(Throughput::Elements(n as u64));

        let monomial_coeffs: Vec<F> = (0..n)
            .map(|i| F::from_canonical_u64((i as u64).wrapping_mul(17).wrapping_add(3)))
            .collect();

        let kernel_coeffs = kernel::from_monomial(&monomial_coeffs);
        let domain = lagrange::roots_of_unity_domain(n);
        let lagrange_values = ntt_evaluate_monomial(&monomial_coeffs);

        let x = F::from_canonical_u64(19);

        let expected = monomial::eval_horner(&monomial_coeffs, x);

        assert_eq!(
            kernel::eval_parallel_subtrees_worker_scratch(&kernel_coeffs, x, 4096,),
            expected,
            "kernel mismatch at N={n}",
        );

        assert_eq!(
            lagrange::eval_roots_of_unity(&lagrange_values, &domain, x),
            expected,
            "serial Lagrange mismatch at N={n}",
        );

        assert_eq!(
            lagrange::eval_roots_of_unity_parallel(&lagrange_values, &domain, x,),
            expected,
            "parallel Lagrange mismatch at N={n}",
        );

        group.bench_with_input(BenchmarkId::new("kernel_optimized", n), &n, |b, _| {
            b.iter(|| {
                black_box(kernel::eval_parallel_subtrees_worker_scratch(
                    black_box(&kernel_coeffs),
                    black_box(x),
                    4096,
                ))
            })
        });

        group.bench_with_input(BenchmarkId::new("monomial_parallel_4096", n), &n, |b, _| {
            b.iter(|| {
                black_box(monomial::eval_parallel_blocks(
                    black_box(&monomial_coeffs),
                    black_box(x),
                    4096,
                ))
            })
        });

        group.bench_with_input(BenchmarkId::new("monomial_horner", n), &n, |b, _| {
            b.iter(|| {
                black_box(monomial::eval_horner(
                    black_box(&monomial_coeffs),
                    black_box(x),
                ))
            })
        });

        group.bench_with_input(BenchmarkId::new("lagrange_parallel", n), &n, |b, _| {
            b.iter(|| {
                black_box(lagrange::eval_roots_of_unity_parallel(
                    black_box(&lagrange_values),
                    black_box(&domain),
                    black_box(x),
                ))
            })
        });

        group.bench_with_input(BenchmarkId::new("lagrange_serial", n), &n, |b, _| {
            b.iter(|| {
                black_box(lagrange::eval_roots_of_unity(
                    black_box(&lagrange_values),
                    black_box(&domain),
                    black_box(x),
                ))
            })
        });
    }

    group.finish();
}

criterion_group!(benches, bench_final_bases);
criterion_main!(benches);
