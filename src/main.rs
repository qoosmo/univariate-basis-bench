use rand::{RngCore, SeedableRng};
use rand::rngs::StdRng;
use std::time::Instant;
use univariate_basis_bench::{field::F, kernel, lagrange, monomial};

fn random_vec(n: usize, seed: u64) -> Vec<F> {
    let mut rng = StdRng::seed_from_u64(seed);
    (0..n).map(|_| F::new(rng.next_u64() as u128)).collect()
}

fn measure<FN: FnMut() -> F, F>(mut f: FN, reps: usize) -> f64 {
    let start = Instant::now();
    let mut sink = None;
    for _ in 0..reps { sink = Some(f()); }
    std::hint::black_box(sink);
    start.elapsed().as_secs_f64() * 1e3 / reps as f64
}

fn main() {
    println!("basis,n,ms,ns_per_element");
    for log_n in [10u32, 12, 14, 16, 18, 20] {
        let n = 1usize << log_n;
        let coeffs = random_vec(n, 1 + n as u64);
        let kcoeffs = random_vec(n, 2 + n as u64);
        let lvals = random_vec(n, 3 + n as u64);
        let domain = lagrange::roots_of_unity_domain(n);
        let x = F::from_canonical_u64(123456789);
        let reps = if n <= (1<<14) { 100 } else if n <= (1<<18) { 20 } else { 5 };

        let rows = [
            ("monomial_horner", measure(|| monomial::eval_horner(&coeffs, x), reps)),
            ("monomial_parallel_blocks", measure(|| monomial::eval_parallel_blocks(&coeffs, x, 4096), reps)),
            ("lagrange_barycentric", measure(|| lagrange::eval_roots_of_unity(&lvals, &domain, x), reps)),
            ("lagrange_barycentric_parallel", measure(|| lagrange::eval_roots_of_unity_parallel(&lvals, &domain, x), reps)),
            ("kernel_serial", measure(|| kernel::eval_serial(&kcoeffs, x), reps)),
            ("kernel_parallel", measure(|| kernel::eval_parallel(&kcoeffs, x, 2048), reps)),
        ];
        for (name, ms) in rows {
            println!("{name},{n},{ms:.6},{:.3}", ms * 1e6 / n as f64);
        }
    }
}
