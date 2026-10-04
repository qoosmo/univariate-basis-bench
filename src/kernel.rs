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


/// In-place parallel kernel tree.
///
/// Unlike `eval_parallel`, this performs no per-level allocation.
/// At level i, each independent chunk of size 2^(i+1) stores its
/// folded result in chunk[0].
pub fn eval_parallel_in_place(
    coeffs: &[F],
    x: F,
    parallel_threshold: usize,
) -> F {
    assert!(coeffs.len().is_power_of_two());

    if coeffs.len() == 1 {
        return coeffs[0];
    }

    let mut buf = coeffs.to_vec();
    let n = buf.len();

    let mut span = 2usize;
    let mut t = x;

    loop {
        let half = span / 2;
        let independent_chunks = n / span;

        if independent_chunks >= parallel_threshold {
            buf.par_chunks_mut(span).for_each(|chunk| {
                let a = chunk[0];
                let b = chunk[half];
                chunk[0] = t * (a + b) + b;
            });
        } else {
            for chunk in buf.chunks_mut(span) {
                let a = chunk[0];
                let b = chunk[half];
                chunk[0] = t * (a + b) + b;
            }
        }

        if span == n {
            break;
        }

        span *= 2;
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

/// Coarse-grained parallel kernel evaluation.
///
/// The coefficient array is split into contiguous power-of-two subtrees.
/// Each Rayon worker reduces one complete subtree locally, avoiding a
/// global synchronization barrier at every kernel-tree level.
pub fn eval_parallel_subtrees(
    coeffs: &[F],
    x: F,
    block_size: usize,
) -> F {
    assert!(coeffs.len().is_power_of_two());
    assert!(block_size.is_power_of_two());

    let block_size = block_size.min(coeffs.len());

    if coeffs.len() <= block_size {
        return eval_serial(coeffs, x);
    }

    let total_levels = coeffs.len().trailing_zeros() as usize;
    let local_levels = block_size.trailing_zeros() as usize;

    // t_i = x^(2^i)
    let mut ts = Vec::with_capacity(total_levels);
    let mut t = x;
    for _ in 0..total_levels {
        ts.push(t);
        t = t.square();
    }

    // Reduce each contiguous subtree independently.
    let mut roots: Vec<F> = coeffs
        .par_chunks(block_size)
        .map(|chunk| {
            let mut buf = chunk.to_vec();
            let mut len = block_size;

            for &level_t in &ts[..local_levels] {
                let pairs = len / 2;

                for j in 0..pairs {
                    let a = buf[2 * j];
                    let b = buf[2 * j + 1];
                    buf[j] = level_t * (a + b) + b;
                }

                len = pairs;
            }

            buf[0]
        })
        .collect();

    // Finish the small top tree serially.
    for &level_t in &ts[local_levels..] {
        let pairs = roots.len() / 2;

        for j in 0..pairs {
            let a = roots[2 * j];
            let b = roots[2 * j + 1];
            roots[j] = level_t * (a + b) + b;
        }

        roots.truncate(pairs);
    }

    roots[0]
}

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
                let rhs = c[z | bit];
                c[z] -= rhs;
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
    fn parallel_in_place_matches_serial() {
        let c: Vec<F> = (0..65536)
            .map(|i| F::from_canonical_u64((i * 17 + 3) as u64))
            .collect();

        let x = F::from_canonical_u64(19);

        assert_eq!(
            eval_serial(&c, x),
            eval_parallel_in_place(&c, x, 2048)
        );
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

    #[test]
    fn parallel_subtrees_matches_serial() {
        let coeffs: Vec<F> = (0..65536)
            .map(|i| F::from_canonical_u64((i as u64).wrapping_mul(17).wrapping_add(3)))
            .collect();

        let x = F::from_canonical_u64(19);

        let expected = eval_serial(&coeffs, x);

        for block_size in [512usize, 1024, 2048, 4096, 8192, 16384] {
            assert_eq!(
                expected,
                eval_parallel_subtrees(&coeffs, x, block_size),
                "block_size={block_size}"
            );
        }
    }


}
