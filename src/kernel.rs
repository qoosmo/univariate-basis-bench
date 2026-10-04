use rayon::prelude::*;
use crate::field::F;

/// Evaluate P(X)=sum_y c_y K_y(X),
/// K_y(X)=prod_i (X^(2^i)+y_i), y_i in {0,1}.
///
/// Coefficient ordering: integer y, bit i = y_i.
/// At level i, each adjacent pair (y_i=0, y_i=1) is folded as
///     c' = t*(a+b)+b,  t = x^(2^i).
pub fn eval_serial(coeffs: &[F], x: F) -> F {
    assert!(coeffs.len().is_power_of_two());
    let mut buf = coeffs.to_vec();
    let mut len = buf.len();
    let mut t = x;
    while len > 1 {
        for j in 0..(len / 2) {
            let a = buf[2*j];
            let b = buf[2*j + 1];
            buf[j] = t * (a + b) + b;
        }
        len /= 2;
        t = t.square();
    }
    buf[0]
}

/// Same arithmetic, but independent pair folds at each tree level use Rayon.
/// The threshold avoids thread-pool overhead once levels become small.
pub fn eval_parallel(coeffs: &[F], x: F, parallel_threshold: usize) -> F {
    assert!(coeffs.len().is_power_of_two());
    let mut current = coeffs.to_vec();
    let mut t = x;
    while current.len() > 1 {
        let pairs = current.len() / 2;
        let next: Vec<F> = if pairs >= parallel_threshold {
            current
                .par_chunks_exact(2)
                .map(|ab| t * (ab[0] + ab[1]) + ab[1])
                .collect()
        } else {
            current
                .chunks_exact(2)
                .map(|ab| t * (ab[0] + ab[1]) + ab[1])
                .collect()
        };
        current = next;
        t = t.square();
    }
    current[0]
}

/// Allocation-free parallel version using one work buffer.
/// Each level writes into the lower half, after computing results in parallel
/// into a temporary level buffer. Kept separate so we can benchmark algorithmic
/// parallelism independently of allocator effects.
pub fn eval_parallel_reuse(coeffs: &[F], x: F, parallel_threshold: usize) -> F {
    assert!(coeffs.len().is_power_of_two());
    let mut buf = coeffs.to_vec();
    let mut len = buf.len();
    let mut t = x;
    while len > 1 {
        let pairs = len / 2;
        if pairs >= parallel_threshold {
            let next: Vec<F> = buf[..len]
                .par_chunks_exact(2)
                .map(|ab| t * (ab[0] + ab[1]) + ab[1])
                .collect();
            buf[..pairs].copy_from_slice(&next);
        } else {
            for j in 0..pairs {
                let a = buf[2*j];
                let b = buf[2*j + 1];
                buf[j] = t * (a + b) + b;
            }
        }
        len = pairs;
        t = t.square();
    }
    buf[0]
}


/// Convert monomial coefficients a_e to kernel coefficients c_y.
///
/// For N=2^m and M=N-1,
///   a_e = sum_{y superset (M xor e)} c_y.
/// Thus after bit-complementing the monomial index, this is an inverse
/// superset-zeta (Boolean Möbius) transform.
pub fn from_monomial(monomial: &[F]) -> Vec<F> {
    assert!(monomial.len().is_power_of_two());
    let n = monomial.len();
    let mask = n - 1;
    let mut c: Vec<F> = (0..n).map(|z| monomial[mask ^ z]).collect();
    let m = n.trailing_zeros() as usize;
    for i in 0..m {
        let bit = 1usize << i;
        for z in 0..n {
            if z & bit == 0 {
                c[z] -= c[z | bit];
            }
        }
    }
    c
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::monomial;

    fn direct_kernel_eval(c: &[F], x: F) -> F {
        let n = c.len();
        let m = n.trailing_zeros() as usize;
        let mut out = F::ZERO;
        for y in 0..n {
            let mut k = F::ONE;
            let mut t = x;
            for i in 0..m {
                k *= if (y >> i) & 1 == 1 { t + F::ONE } else { t };
                t = t.square();
            }
            out += c[y] * k;
        }
        out
    }

    #[test]
    fn tree_matches_direct_definition() {
        let c: Vec<F> = (0..16).map(|i| F::from_canonical_u64((i * 17 + 3) as u64)).collect();
        let x = F::from_canonical_u64(19);
        assert_eq!(eval_serial(&c, x), direct_kernel_eval(&c, x));
    }

    #[test]
    fn monomial_conversion_preserves_polynomial() {
        let a: Vec<F> = (0..16).map(|i| F::from_canonical_u64((i * 11 + 5) as u64)).collect();
        let c = from_monomial(&a);
        for x0 in [2u64, 7, 19, 123] {
            let x = F::from_canonical_u64(x0);
            assert_eq!(monomial::eval_horner(&a, x), eval_serial(&c, x));
        }
    }
}
