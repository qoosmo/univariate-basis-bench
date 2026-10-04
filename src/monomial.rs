use rayon::prelude::*;
use crate::field::F;

/// Coefficients are [a_0, a_1, ..., a_{n-1}].
#[inline]
pub fn eval_horner(coeffs: &[F], x: F) -> F {
    coeffs.iter().rev().copied().fold(F::ZERO, |acc, a| acc * x + a)
}

/// Parallel single-polynomial evaluation by splitting coefficients into blocks.
///
/// P(x) = sum_j x^(jB) * P_j(x), where each block P_j is evaluated by Horner.
/// This is an important fair baseline: monomial evaluation is not forced to be
/// purely serial if extra cores are available.
pub fn eval_parallel_blocks(coeffs: &[F], x: F, block_size: usize) -> F {
    assert!(block_size > 0);
    if coeffs.len() <= block_size {
        return eval_horner(coeffs, x);
    }
    let x_b = x.pow(block_size as u64);
    coeffs
        .par_chunks(block_size)
        .enumerate()
        .map(|(j, chunk)| {
            let local = eval_horner(chunk, x);
            local * x_b.pow(j as u64)
        })
        .reduce(|| F::ZERO, |a, b| a + b)
}
