use rayon::prelude::*;
use crate::field::F;

/// Evaluation domain H = {1, omega, ..., omega^(n-1)}.
pub fn roots_of_unity_domain(n: usize) -> Vec<F> {
    let omega = F::primitive_root_of_unity(n);
    let mut out = Vec::with_capacity(n);
    let mut cur = F::ONE;
    for _ in 0..n {
        out.push(cur);
        cur *= omega;
    }
    out
}

/// Evaluate from Lagrange values v_i=P(h_i) on a roots-of-unity domain.
///
/// For h_i^n=1:
/// L_i(x) = (x^n-1) * h_i / (n * (x-h_i)).
/// We batch-invert all denominators with one field inversion.
pub fn eval_roots_of_unity(values: &[F], domain: &[F], x: F) -> F {
    assert_eq!(values.len(), domain.len());
    let n = values.len();
    assert!(n.is_power_of_two());

    if let Some((idx, _)) = domain.iter().enumerate().find(|(_, h)| **h == x) {
        return values[idx];
    }

    let mut denoms: Vec<F> = domain.iter().map(|&h| x - h).collect();
    batch_inverse_in_place(&mut denoms);

    let n_inv = F::from_canonical_u64(n as u64).inverse();
    let factor = (x.pow(n as u64) - F::ONE) * n_inv;

    values.iter().zip(domain).zip(&denoms)
        .fold(F::ZERO, |acc, ((&v, &h), &inv)| acc + v * h * inv) * factor
}

/// Parallel numerator accumulation; denominator inversion is still the usual
/// one-inversion batch algorithm.
pub fn eval_roots_of_unity_parallel(values: &[F], domain: &[F], x: F) -> F {
    assert_eq!(values.len(), domain.len());
    let n = values.len();
    assert!(n.is_power_of_two());

    if let Some((idx, _)) = domain.iter().enumerate().find(|(_, h)| **h == x) {
        return values[idx];
    }

    let mut inv_denoms: Vec<F> = domain.par_iter().map(|&h| x - h).collect();
    batch_inverse_in_place(&mut inv_denoms);

    let n_inv = F::from_canonical_u64(n as u64).inverse();
    let factor = (x.pow(n as u64) - F::ONE) * n_inv;

    values.par_iter().zip(domain.par_iter()).zip(inv_denoms.par_iter())
        .map(|((&v, &h), &inv)| v * h * inv)
        .reduce(|| F::ZERO, |a,b| a+b) * factor
}

pub fn batch_inverse_in_place(xs: &mut [F]) {
    let n = xs.len();
    if n == 0 { return; }
    let mut prefix = Vec::with_capacity(n);
    let mut acc = F::ONE;
    for &x in xs.iter() {
        prefix.push(acc);
        acc *= x;
    }
    let mut inv_acc = acc.inverse();
    for i in (0..n).rev() {
        let x = xs[i];
        xs[i] = inv_acc * prefix[i];
        inv_acc *= x;
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::monomial;

    #[test]
    fn lagrange_eval_recovers_polynomial() {
        let n = 16usize;
        let a: Vec<F> = (0..n).map(|i| F::from_canonical_u64((i * 13 + 9) as u64)).collect();
        let domain = roots_of_unity_domain(n);
        let values: Vec<F> = domain.iter().map(|&h| monomial::eval_horner(&a, h)).collect();
        let x = F::from_canonical_u64(123456789);
        assert_eq!(eval_roots_of_unity(&values, &domain, x), monomial::eval_horner(&a, x));
        assert_eq!(eval_roots_of_unity_parallel(&values, &domain, x), monomial::eval_horner(&a, x));
    }
}
