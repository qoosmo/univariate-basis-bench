use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use rayon::prelude::*;

use univariate_basis_bench::{field::F, kernel, monomial};

const P: u64 = 0xffff_ffff_0000_0001;
const MONT_NINV: u64 = 0xffff_fffe_ffff_ffff;
const MONT_R2: u64 = 0xffff_fffe_0000_0001;

#[derive(Copy, Clone, Debug, Default, Eq, PartialEq)]
#[repr(transparent)]
struct Mont(u64);

impl Mont {
    const ZERO: Self = Self(0);
    const ONE: Self = Self(0x0000_0000_ffff_ffff);

    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        let sum = self.0 as u128 + rhs.0 as u128;
        let p = P as u128;
        if sum >= p {
            Self((sum - p) as u64)
        } else {
            Self(sum as u64)
        }
    }

    #[inline(always)]
    fn mul(self, rhs: Self) -> Self {
        Self(mont_redc(self.0 as u128 * rhs.0 as u128))
    }

    #[inline(always)]
    fn square(self) -> Self {
        self.mul(self)
    }

    #[inline(always)]
    fn from_canonical(x: u64) -> Self {
        Self(mont_redc((x as u128) * (MONT_R2 as u128)))
    }

    #[inline(always)]
    fn to_canonical(self) -> u64 {
        mont_redc(self.0 as u128)
    }

    #[inline]
    fn pow(mut self, mut e: u64) -> Self {
        let mut acc = Self::ONE;
        while e != 0 {
            if e & 1 == 1 {
                acc = acc.mul(self);
            }
            self = self.square();
            e >>= 1;
        }
        acc
    }
}

#[inline(always)]
fn mont_redc(t: u128) -> u64 {
    let m = (t as u64).wrapping_mul(MONT_NINV);
    let mp = (m as u128) * (P as u128);
    let (sum, carry) = t.overflowing_add(mp);

    let mut u = (sum >> 64) + ((carry as u128) << 64);
    if u >= P as u128 {
        u -= P as u128;
    }
    u as u64
}

#[inline(always)]
fn mont_fold(t: Mont, a: Mont, b: Mont) -> Mont {
    t.mul(a.add(b)).add(b)
}

fn mont_horner(coeffs: &[Mont], x: Mont) -> Mont {
    coeffs
        .iter()
        .rev()
        .copied()
        .fold(Mont::ZERO, |acc, a| acc.mul(x).add(a))
}

fn mont_parallel_blocks(coeffs: &[Mont], x: Mont, block_size: usize) -> Mont {
    assert!(block_size > 0);
    if coeffs.len() <= block_size {
        return mont_horner(coeffs, x);
    }

    let x_b = x.pow(block_size as u64);

    coeffs
        .par_chunks(block_size)
        .enumerate()
        .map(|(j, chunk)| {
            let local = mont_horner(chunk, x);
            local.mul(x_b.pow(j as u64))
        })
        .reduce(|| Mont::ZERO, |a, b| a.add(b))
}

fn mont_kernel_worker_scratch(coeffs: &[Mont], x: Mont, block_size: usize) -> Mont {
    assert!(coeffs.len().is_power_of_two());
    assert!(block_size.is_power_of_two());

    let block_size = block_size.min(coeffs.len());
    if coeffs.len() == 1 {
        return coeffs[0];
    }

    let total_levels = coeffs.len().trailing_zeros() as usize;
    let local_levels = block_size.trailing_zeros() as usize;

    let mut ts = Vec::with_capacity(total_levels);
    let mut t = x;
    for _ in 0..total_levels {
        ts.push(t);
        t = t.square();
    }

    if coeffs.len() <= block_size {
        let mut buf = coeffs.to_vec();
        let mut len = buf.len();
        for &level_t in &ts {
            let pairs = len / 2;
            for j in 0..pairs {
                let a = buf[2 * j];
                let b = buf[2 * j + 1];
                buf[j] = mont_fold(level_t, a, b);
            }
            len = pairs;
        }
        return buf[0];
    }

    let mut roots: Vec<Mont> = coeffs
        .par_chunks(block_size)
        .map_init(
            || Vec::<Mont>::with_capacity(block_size / 2),
            |buf, chunk| {
                buf.clear();

                for pair in chunk.chunks_exact(2) {
                    buf.push(mont_fold(ts[0], pair[0], pair[1]));
                }

                let mut len = block_size / 2;

                for &level_t in &ts[1..local_levels] {
                    let pairs = len / 2;
                    for j in 0..pairs {
                        let a = buf[2 * j];
                        let b = buf[2 * j + 1];
                        buf[j] = mont_fold(level_t, a, b);
                    }
                    len = pairs;
                }

                buf[0]
            },
        )
        .collect();

    for &level_t in &ts[local_levels..] {
        let pairs = roots.len() / 2;
        for j in 0..pairs {
            let a = roots[2 * j];
            let b = roots[2 * j + 1];
            roots[j] = mont_fold(level_t, a, b);
        }
        roots.truncate(pairs);
    }

    roots[0]
}

fn sanity() {
    let vals = [
        0u64,
        1,
        2,
        7,
        19,
        0xffff_ffff,
        P - 1,
        P - 2,
        0x1234_5678_9abc_def0,
    ];

    for &a in &vals {
        if a >= P {
            continue;
        }
        let am = Mont::from_canonical(a);
        assert_eq!(am.to_canonical(), a);

        for &b in &vals {
            if b >= P {
                continue;
            }
            let bm = Mont::from_canonical(b);
            let got = am.mul(bm).to_canonical();
            let want = ((a as u128 * b as u128) % P as u128) as u64;
            assert_eq!(got, want, "Montgomery mismatch: a={a}, b={b}");
        }
    }

    let mut a = 0x1234_5678_9abc_def0u64 % P;
    let mut b = 0x0fed_cba9_8765_4321u64 % P;
    for _ in 0..100_000 {
        a = a.wrapping_mul(6364136223846793005).wrapping_add(1) % P;
        b = b.wrapping_mul(1442695040888963407).wrapping_add(33) % P;

        let got = Mont::from_canonical(a)
            .mul(Mont::from_canonical(b))
            .to_canonical();
        let want = ((a as u128 * b as u128) % P as u128) as u64;
        assert_eq!(got, want);
    }
}

fn bench_mul_backends(c: &mut Criterion) {
    sanity();

    let a = F::from_canonical_u64(0x1234_5678_9abc_def0 % P);
    let b = F::from_canonical_u64(0x0fed_cba9_8765_4321 % P);

    let am = Mont::from_canonical(a.0);
    let bm = Mont::from_canonical(b.0);

    let mut g = c.benchmark_group("mul_backend");

    g.sample_size(100);
    g.warm_up_time(std::time::Duration::from_secs(5));
    g.measurement_time(std::time::Duration::from_secs(15));

    g.bench_function("goldilocks_special_mul", |bencher| {
        bencher.iter(|| black_box(black_box(a) * black_box(b)))
    });

    g.bench_function("montgomery_mul", |bencher| {
        bencher.iter(|| black_box(black_box(am).mul(black_box(bm))))
    });

    let t = F::from_canonical_u64(19);
    let aa = F::from_canonical_u64(1234567);
    let bb = F::from_canonical_u64(7654321);

    let tm = Mont::from_canonical(19);
    let aam = Mont::from_canonical(1234567);
    let bbm = Mont::from_canonical(7654321);

    g.bench_function("goldilocks_kernel_fold", |bencher| {
        bencher.iter(|| black_box(F::kernel_fold(black_box(t), black_box(aa), black_box(bb))))
    });

    g.bench_function("montgomery_kernel_fold", |bencher| {
        bencher.iter(|| black_box(mont_fold(black_box(tm), black_box(aam), black_box(bbm))))
    });

    g.finish();
}

fn bench_full_evaluators(c: &mut Criterion) {
    let n = 1usize << 20;
    let block = 4096usize;

    let monomial_coeffs: Vec<F> = (0..n)
        .map(|i| F::from_canonical_u64((i as u64).wrapping_mul(17).wrapping_add(3)))
        .collect();

    let kernel_coeffs = kernel::from_monomial(&monomial_coeffs);
    let x = F::from_canonical_u64(19);

    let monomial_m: Vec<Mont> = monomial_coeffs
        .iter()
        .map(|v| Mont::from_canonical(v.0))
        .collect();

    let kernel_m: Vec<Mont> = kernel_coeffs
        .iter()
        .map(|v| Mont::from_canonical(v.0))
        .collect();

    let xm = Mont::from_canonical(19);

    let expected = monomial::eval_horner(&monomial_coeffs, x).0;

    assert_eq!(
        kernel::eval_parallel_subtrees_worker_scratch(&kernel_coeffs, x, block).0,
        expected
    );
    assert_eq!(
        mont_parallel_blocks(&monomial_m, xm, block).to_canonical(),
        expected
    );
    assert_eq!(
        mont_kernel_worker_scratch(&kernel_m, xm, block).to_canonical(),
        expected
    );

    let mut g = c.benchmark_group("backend_full_eval_n1048576");
    g.sample_size(50);
    g.warm_up_time(std::time::Duration::from_secs(5));
    g.measurement_time(std::time::Duration::from_secs(12));
    g.throughput(Throughput::Elements(n as u64));

    g.bench_with_input(
        BenchmarkId::new("kernel_goldilocks_special", n),
        &n,
        |b, _| {
            b.iter(|| {
                black_box(kernel::eval_parallel_subtrees_worker_scratch(
                    black_box(&kernel_coeffs),
                    black_box(x),
                    block,
                ))
            })
        },
    );

    g.bench_with_input(BenchmarkId::new("kernel_montgomery", n), &n, |b, _| {
        b.iter(|| {
            black_box(mont_kernel_worker_scratch(
                black_box(&kernel_m),
                black_box(xm),
                block,
            ))
        })
    });

    g.bench_with_input(
        BenchmarkId::new("monomial_goldilocks_special", n),
        &n,
        |b, _| {
            b.iter(|| {
                black_box(monomial::eval_parallel_blocks(
                    black_box(&monomial_coeffs),
                    black_box(x),
                    block,
                ))
            })
        },
    );

    g.bench_with_input(BenchmarkId::new("monomial_montgomery", n), &n, |b, _| {
        b.iter(|| {
            black_box(mont_parallel_blocks(
                black_box(&monomial_m),
                black_box(xm),
                block,
            ))
        })
    });

    g.finish();
}

criterion_group!(benches, bench_mul_backends, bench_full_evaluators);
criterion_main!(benches);
