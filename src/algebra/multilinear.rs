use ark_ff::{AdditiveGroup, Field};

use crate::algebra::embedding::{Embedding, Identity};

/// Evaluate the multi-linear extension of `evals` in `point`.
pub fn multilinear_extend<F: Field>(evals: &[F], point: &[F]) -> F {
    mixed_multilinear_extend(&Identity::<F>::new(), evals, point)
}

/// Evaluate the multi-linear extension of `evals` in `point`.
///
/// Supports implicit zero-padding: when `evals.len() < 1 << point.len()`,
/// the missing tail entries are treated as zeros.
#[allow(clippy::too_many_lines)]
pub fn mixed_multilinear_extend<M: Embedding>(
    embedding: &M,
    evals: &[M::Source],
    point: &[M::Target],
) -> M::Target {
    #[inline]
    fn eval_exact<M: Embedding>(
        embedding: &M,
        evals: &[M::Source],
        point: &[M::Target],
    ) -> M::Target {
        debug_assert_eq!(evals.len(), 1 << point.len());

        // Helper to compute (a + (b - a) * c) efficiently with a, b in source field.
        let mixed = |a, b, c| embedding.mixed_add(embedding.mixed_mul(c, b - a), a);

        match point {
            [] => embedding.map(evals[0]),
            [x] => mixed(evals[0], evals[1], *x),
            [x0, x1] => {
                let a0 = mixed(evals[0], evals[1], *x1);
                let a1 = mixed(evals[2], evals[3], *x1);
                a0 + (a1 - a0) * *x0
            }
            [x0, x1, x2] => {
                let a00 = mixed(evals[0], evals[1], *x2);
                let a01 = mixed(evals[2], evals[3], *x2);
                let a10 = mixed(evals[4], evals[5], *x2);
                let a11 = mixed(evals[6], evals[7], *x2);
                let a0 = a00 + (a01 - a00) * *x1;
                let a1 = a10 + (a11 - a10) * *x1;
                a0 + (a1 - a0) * *x0
            }
            [x0, x1, x2, x3] => {
                let a000 = mixed(evals[0], evals[1], *x3);
                let a001 = mixed(evals[2], evals[3], *x3);
                let a010 = mixed(evals[4], evals[5], *x3);
                let a011 = mixed(evals[6], evals[7], *x3);
                let a100 = mixed(evals[8], evals[9], *x3);
                let a101 = mixed(evals[10], evals[11], *x3);
                let a110 = mixed(evals[12], evals[13], *x3);
                let a111 = mixed(evals[14], evals[15], *x3);
                let a00 = a000 + (a001 - a000) * *x2;
                let a01 = a010 + (a011 - a010) * *x2;
                let a10 = a100 + (a101 - a100) * *x2;
                let a11 = a110 + (a111 - a110) * *x2;
                let a0 = a00 + (a01 - a00) * *x1;
                let a1 = a10 + (a11 - a10) * *x1;
                a0 + (a1 - a0) * *x0
            }
            [x, tail @ ..] => {
                let (f0, f1) = evals.split_at(evals.len() / 2);
                #[cfg(not(feature = "parallel"))]
                let (f0, f1) = (
                    eval_exact(embedding, f0, tail),
                    eval_exact(embedding, f1, tail),
                );

                #[cfg(feature = "parallel")]
                let (f0, f1) = {
                    use crate::utils::workload_size;
                    if evals.len() > workload_size::<M::Source>() {
                        rayon::join(
                            || eval_exact(embedding, f0, tail),
                            || eval_exact(embedding, f1, tail),
                        )
                    } else {
                        (
                            eval_exact(embedding, f0, tail),
                            eval_exact(embedding, f1, tail),
                        )
                    }
                };

                f0 + (f1 - f0) * *x
            }
        }
    }

    #[inline]
    fn eval_partial<M: Embedding>(
        embedding: &M,
        evals: &[M::Source],
        point: &[M::Target],
    ) -> M::Target {
        let size = 1 << point.len();
        debug_assert!(evals.len() <= size);
        if evals.is_empty() {
            return M::Target::ZERO;
        }
        if evals.len() == size {
            return eval_exact(embedding, evals, point);
        }

        match point {
            [] => embedding.map(evals[0]),
            [x, tail @ ..] => {
                let half = size / 2;

                // Only low half has data; high half is all implicit zeros.
                if evals.len() <= half {
                    let f0 = eval_partial(embedding, evals, tail);
                    return f0 * (M::Target::ONE - *x);
                }

                // Low subtree is exact/full, high subtree is partial.
                let (low, high) = evals.split_at(half);

                #[cfg(not(feature = "parallel"))]
                let (f0, f1) = (
                    eval_exact(embedding, low, tail),
                    eval_partial(embedding, high, tail),
                );

                #[cfg(feature = "parallel")]
                let (f0, f1) = {
                    use crate::utils::workload_size;
                    if evals.len() > workload_size::<M::Source>() {
                        rayon::join(
                            || eval_exact(embedding, low, tail),
                            || eval_partial(embedding, high, tail),
                        )
                    } else {
                        (
                            eval_exact(embedding, low, tail),
                            eval_partial(embedding, high, tail),
                        )
                    }
                };

                f0 + (f1 - f0) * *x
            }
        }
    }

    eval_partial(embedding, evals, point)
}

/// Evaluate the MLE of `evals` at `point` using a mixed binary/ternary
/// folding schedule.
///
/// Clones `evals` once, then folds in-place using the mixed schedule.
/// For pow2 sizes where `evals.len() == 1 << point.len()`, delegates to
/// the optimized butterfly `multilinear_extend`.
///
/// `ternary_start`: index in `point` at which ternary rounds begin.
/// Challenges before this index always use binary folding.
pub fn smooth_multilinear_extend<F: Field>(evals: &[F], point: &[F], ternary_start: usize) -> F {
    if evals.is_empty() {
        return F::ZERO;
    }
    if evals.len() == 1 || point.is_empty() {
        return evals[0];
    }
    let mut w = evals.to_vec();
    smooth_multilinear_extend_in_place(&mut w, point, ternary_start)
}

/// In-place variant of `smooth_multilinear_extend`. Consumes the caller's
/// buffer and returns the final fold result — skips the `evals.to_vec()`
/// clone, saving an O(N) allocation when the caller already owns the vector.
pub fn smooth_multilinear_extend_in_place<F: Field>(
    w: &mut Vec<F>,
    point: &[F],
    ternary_start: usize,
) -> F {
    use crate::algebra::sumcheck::{fold, fold3};
    use crate::protocols::sumcheck::is_ternary_round;

    if w.is_empty() {
        return F::ZERO;
    }
    if w.len() == 1 || point.is_empty() {
        return w[0];
    }

    let mut size = w.len();
    for (i, &r) in point.iter().enumerate() {
        if size <= 1 {
            break;
        }
        if i >= ternary_start && is_ternary_round(size) {
            fold3(w, r);
            size /= 3;
        } else {
            fold(w, r);
            size = (size + 1) / 2;
        }
    }
    debug_assert_eq!(w.len(), 1);
    w[0]
}

fn smooth_mle_inner<F: Field>(
    evals: &[F],
    point: &[F],
    challenge_offset: usize,
    ternary_start: usize,
) -> F {
    use crate::protocols::sumcheck::is_ternary_round;

    if evals.is_empty() {
        return F::ZERO;
    }
    if evals.len() == 1 || point.is_empty() {
        return evals[0];
    }

    let size = evals.len();
    let r = point[0];
    let rest = &point[1..];
    let is_tern = challenge_offset >= ternary_start && is_ternary_round(size);

    if is_tern {
        // Ternary: split into exact thirds — no allocation.
        let third = size / 3;
        let e0 = &evals[..third];
        let e1 = &evals[third..2 * third];
        let e2 = &evals[2 * third..];

        let half_inv = F::from(2u64).inverse().expect("char != 2");
        let two = F::from(2u64);
        let l0 = (r - F::ONE) * (r - two) * half_inv;
        let l1 = r * (two - r);
        let l2 = r * (r - F::ONE) * half_inv;

        let f0 = smooth_mle_inner(e0, rest, challenge_offset + 1, ternary_start);
        let f1 = smooth_mle_inner(e1, rest, challenge_offset + 1, ternary_start);
        let f2 = smooth_mle_inner(e2, rest, challenge_offset + 1, ternary_start);

        f0 * l0 + f1 * l1 + f2 * l2
    } else if size % 2 == 0 {
        // Even binary: clean split — no allocation.
        let half = size / 2;
        let lo = &evals[..half];
        let hi = &evals[half..];

        let f0 = smooth_mle_inner(lo, rest, challenge_offset + 1, ternary_start);
        let f1 = smooth_mle_inner(hi, rest, challenge_offset + 1, ternary_start);

        f0 + (f1 - f0) * r
    } else {
        // Odd binary: small allocation for the folded vector.
        // Odd sizes are rare in practice (only from partial binary reduction
        // of 3^b residuals). The allocation is small (half+1 elements).
        let half = size / 2;
        let output_len = half + 1;
        let mut folded = Vec::with_capacity(output_len);
        for i in 0..half {
            folded.push(evals[i] + (evals[half + i] - evals[i]) * r);
        }
        folded.push(evals[half + half] * r);
        smooth_mle_inner(&folded, rest, challenge_offset + 1, ternary_start)
    }
}

/// Computes eq(points, p) on the hypercube for all p ∈ {0,1}^k.
pub fn eq_weights<F: Field>(point: &[F]) -> Vec<F> {
    let mut result = vec![F::ZERO; 1 << point.len()];
    eval_eq(&mut result, point, F::ONE);
    result
}

/// Accumulates a scaled evaluation of the equality function.
///
/// Given an evaluation point `point`, the function computes
/// the equality polynomial recursively using the formula:
///
/// ```text
/// eq(X) = ∏ (1 - X_i + 2X_i z_i)
/// ```
///
/// where `z_i` are the points.
pub fn eval_eq<F: Field>(accumulator: &mut [F], point: &[F], scalar: F) {
    assert_eq!(accumulator.len(), 1 << point.len());
    if let [x0, xs @ ..] = point {
        let (acc_0, acc_1) = accumulator.split_at_mut(1 << xs.len());
        let s1 = scalar * x0; // Contribution when `X_i = 1`
        let s0 = scalar - s1; // Contribution when `X_i = 0`

        #[cfg(feature = "parallel")]
        {
            use crate::utils::workload_size;
            if acc_0.len() > workload_size::<F>() {
                rayon::join(|| eval_eq(acc_0, xs, s0), || eval_eq(acc_1, xs, s1));
                return;
            }
        }
        eval_eq(acc_0, xs, s0);
        eval_eq(acc_1, xs, s1);
    } else {
        accumulator[0] += scalar;
    }
}

#[cfg(test)]
mod tests {
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    use proptest::proptest;

    use super::*;
    use crate::algebra::{random_vector, sumcheck::tests::zero_pad};

    pub type F = crate::algebra::fields::Field64;

    #[test]
    fn test_multilinear_zero_extend() {
        proptest!(|(seed:u64, length in 0_usize..(1 << 14))| {
            let mut rng = StdRng::seed_from_u64(seed);
            let vector: Vec<F> = random_vector(&mut rng, length);
            let extended_vector = zero_pad(&vector);
            let num_variables = length.next_power_of_two().trailing_zeros();
            let point = random_vector(&mut rng, num_variables as usize);
            assert_eq!(
                multilinear_extend(&vector, &point),
                multilinear_extend(&extended_vector, &point)
            );
        });
    }

    #[test]
    fn test_smooth_multilinear_extend() {
        use crate::algebra::sumcheck::{fold, fold3};
        use crate::protocols::sumcheck::is_ternary_round;
        use crate::smooth_domain::is_smooth;

        proptest!(|(seed: u64, length in 1_usize..300)| {
            let mut rng = StdRng::seed_from_u64(seed);
            let vector: Vec<F> = random_vector(&mut rng, length);

            // Generate enough challenges to fold to size 1
            let mut num_rounds = 0;
            let mut s = length;
            while s > 1 {
                if is_ternary_round(s) { s /= 3; } else { s = (s + 1) / 2; }
                num_rounds += 1;
            }
            let point: Vec<F> = random_vector(&mut rng, num_rounds);

            // Reference: fold in-place
            let mut v = vector.clone();
            let mut sz = v.len();
            for &r in &point {
                if is_ternary_round(sz) {
                    fold3(&mut v, r);
                    sz /= 3;
                } else {
                    fold(&mut v, r);
                    sz = (sz + 1) / 2;
                }
            }
            assert_eq!(v.len(), 1);
            let expected = v[0];

            // Test: smooth_multilinear_extend with ternary_start=0
            let actual = smooth_multilinear_extend(&vector, &point, 0);
            assert_eq!(actual, expected, "mismatch at length={length}");
        });
    }

    #[test]
    fn test_multilinear_extra_variables() {
        proptest!(|(seed:u64, length in 0_usize..(1 << 10), excess_variables in 0_usize..3)| {
            let mut rng = StdRng::seed_from_u64(seed);
            let vector: Vec<F> = random_vector(&mut rng, length);
            let num_variables = length.next_power_of_two().trailing_zeros() as usize + excess_variables;
            let point = random_vector(&mut rng, num_variables);
            let mut extended_vector = vector.clone();
            extended_vector.resize(1 << num_variables, F::ZERO);
            assert_eq!(
                multilinear_extend(&vector, &point),
                multilinear_extend(&extended_vector, &point)
            );
        });
    }
}
