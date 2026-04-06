use ark_ff::Field;
#[cfg(feature = "parallel")]
use rayon::join;

use crate::algebra::{dot, embedding::Embedding, scalar_mul};
#[cfg(feature = "parallel")]
use crate::utils::workload_size;

/// Compute the split point for the sumcheck fold.
///
/// Uses `len / 2` (clean halving) so that the fold boundary aligns with
/// interleaving block boundaries for smooth-{2,3} sizes like 12 = 4 × 3.
///
/// For even `len` this gives two equal halves.
/// For odd `len` the low half gets `floor(len/2)` and the high half gets
/// `ceil(len/2)` elements; the extra high element is paired with an implicit
/// zero on the low side.
#[inline]
pub fn fold_half(len: usize) -> usize {
    len / 2
}

/// Computes the constant and quadratic coefficient of the sumcheck polynomial.
///
/// The vector is split at `fold_half(len)`. When the halves are unequal,
/// the shorter side's missing elements are treated as implicit zeros.
pub fn compute_sumcheck_polynomial<F: Field>(a: &[F], b: &[F]) -> (F, F) {
    fn recurse<F: Field>(a0: &[F], a1: &[F], b0: &[F], b1: &[F]) -> (F, F) {
        debug_assert_eq!(a0.len(), b0.len());
        debug_assert_eq!(a1.len(), b1.len());
        debug_assert!(a0.len() == a1.len());

        #[cfg(feature = "parallel")]
        if a0.len() * 4 > workload_size::<F>() {
            let mid = a0.len() / 2;
            let (a0l, a0r) = a0.split_at(mid);
            let (b0l, b0r) = b0.split_at(mid);
            let (a1l, a1r) = a1.split_at(mid);
            let (b1l, b1r) = b1.split_at(mid);
            let (left, right) = join(
                || recurse(a0l, a1l, b0l, b1l),
                || recurse(a0r, a1r, b0r, b1r),
            );
            return (left.0 + right.0, left.1 + right.1);
        }
        let mut acc0 = F::ZERO;
        let mut acc2 = F::ZERO;
        for ((&a0, &a1), (&b0, &b1)) in a0.iter().zip(a1).zip(b0.iter().zip(b1)) {
            acc0 += a0 * b0;
            acc2 += (a1 - a0) * (b1 - b0);
        }
        (acc0, acc2)
    }

    let non_padded = a.len().min(b.len());
    let a = &a[..non_padded];
    let b = &b[..non_padded];
    if a.is_empty() {
        return (F::ZERO, F::ZERO);
    }
    if a.len() == 1 {
        return (a[0] * b[0], F::ZERO);
    }

    let half = fold_half(a.len());
    let (a0, a1) = a.split_at(half);
    let (b0, b1) = b.split_at(half);

    // Paired portion: min(low, high) elements
    let paired = a0.len().min(a1.len());
    let (acc0, acc2) = recurse(&a0[..paired], &a1[..paired], &b0[..paired], &b1[..paired]);

    // Low-side tail: extra elements in the low half paired with implicit zero
    // from the high side. Term: a0*(1-t) * b0*(1-t) → c(0) += a0*b0, c(2) += a0*b0
    let low_tail = dot(&a0[paired..], &b0[paired..]);

    // High-side tail: extra elements in the high half paired with implicit zero
    // from the low side. Term: a1*t * b1*t → c(0) += 0, c(2) += a1*b1
    let high_tail = dot(&a1[paired..], &b1[paired..]);

    (acc0 + low_tail, acc2 + low_tail + high_tail)
}

/// Folds evaluations by linear interpolation at the given weight, in place.
///
/// Splits at `fold_half(len)` so the fold boundary aligns with interleaving
/// block boundaries for smooth-{2,3} sizes.
///
/// For the paired portion: `low[i] += (high[i] - low[i]) * weight`.
/// Low-side tail (extra low elements): `low[i] *= (1 - weight)`.
/// High-side tail (extra high elements): `result[i] = high[i] * weight`.
pub fn fold<F: Field>(values: &mut Vec<F>, weight: F) {
    fn recurse_both<F: Field>(low: &mut [F], high: &[F], weight: F) {
        #[cfg(feature = "parallel")]
        if low.len() > workload_size::<F>() {
            let split = low.len() / 2;
            let (ll, lr) = low.split_at_mut(split);
            let (hl, hr) = high.split_at(split);
            rayon::join(
                || recurse_both(ll, hl, weight),
                || recurse_both(lr, hr, weight),
            );
            return;
        }

        for (low, high) in low.iter_mut().zip(high) {
            *low += (*high - *low) * weight;
        }
    }

    if values.len() <= 1 {
        return;
    }

    let len = values.len();
    let half = fold_half(len);
    let high_len = len - half;
    let output_len = half.max(high_len);

    // Save any extra high-side elements before the in-place mutation.
    let high_extras: Vec<F> = if high_len > half {
        values[half + half..].iter().map(|&v| v * weight).collect()
    } else {
        Vec::new()
    };

    {
        let (low, high) = values.split_at_mut(half);

        let paired = low.len().min(high.len());
        recurse_both(&mut low[..paired], &high[..paired], weight);

        // Low-side tail: low[i] *= (1 - weight)
        scalar_mul(&mut low[paired..], F::ONE - weight);
    }

    values.truncate(output_len);
    for (i, &val) in high_extras.iter().enumerate() {
        values[half + i] = val;
    }
    values.shrink_to_fit();
}

pub fn fold_and_compute_polynomial<F: Field>(a: &mut Vec<F>, b: &mut Vec<F>, weight: F) -> (F, F) {
    // TODO: Replace with a single pass implementation.
    fold(a, weight);
    fold(b, weight);
    compute_sumcheck_polynomial(a, b)
}

/// Evaluate a coefficient vector at a multilinear point in the target field.
pub fn mixed_eval<M: Embedding>(
    embedding: &M,
    coeff: &[M::Source],
    eval: &[M::Target],
    scalar: M::Target,
) -> M::Target {
    debug_assert_eq!(coeff.len(), 1 << eval.len());

    if let Some((&x, tail)) = eval.split_first() {
        let (low, high) = coeff.split_at(coeff.len() / 2);

        #[cfg(feature = "parallel")]
        if low.len() > workload_size::<M::Source>() {
            let (a, b) = join(
                || mixed_eval(embedding, low, tail, scalar),
                || mixed_eval(embedding, high, tail, scalar * x),
            );
            return a + b;
        }

        mixed_eval(embedding, low, tail, scalar) + mixed_eval(embedding, high, tail, scalar * x)
    } else {
        embedding.mixed_mul(scalar, coeff[0])
    }
}

#[cfg(test)]
pub(crate) mod tests {
    use ark_ff::AdditiveGroup;
    use ark_std::rand::{rngs::StdRng, Rng, SeedableRng};
    use proptest::proptest;

    use super::*;
    use crate::algebra::{fields::Field64, random_vector};

    type F = Field64;

    /// Zero-pad to the next power of two.
    pub fn zero_pad<F: Field>(values: &[F]) -> Vec<F> {
        if values.is_empty() {
            return Vec::new();
        }
        let mut vec = values.to_vec();
        vec.resize(vec.len().next_power_of_two(), F::ZERO);
        vec
    }

    /// For power-of-2 lengths, fold_half == next_power_of_two >> 1,
    /// so the result must match the zero-padded version.
    #[test]
    fn sumcheck_poly_power_of_two() {
        proptest!(|(seed:u64, log_length in 0_usize..14)| {
            let length = 1 << log_length;
            let mut rng = StdRng::seed_from_u64(seed);
            let vector: Vec<F> = random_vector(&mut rng, length);
            let covector: Vec<F> = random_vector(&mut rng, length);
            let expected = compute_sumcheck_polynomial(&vector, &covector);
            let extended_vector = zero_pad(&vector);
            let extended_covector = zero_pad(&covector);
            assert_eq!(compute_sumcheck_polynomial(&extended_vector, &extended_covector), expected);
        });
    }

    /// Fold output length must be ceil(len/2).
    #[test]
    fn fold_output_length() {
        proptest!(|(seed:u64, length in 2_usize..1000)| {
            let mut rng = StdRng::seed_from_u64(seed);
            let mut vector: Vec<F> = random_vector(&mut rng, length);
            let weight = rng.gen::<F>();
            fold(&mut vector, weight);
            let expected_len = (length + 1) / 2;
            assert_eq!(vector.len(), expected_len);
        });
    }

    /// For power-of-2 lengths, fold matches the zero-padded version.
    #[test]
    fn fold_power_of_two() {
        proptest!(|(seed:u64, log_length in 1_usize..14)| {
            let length = 1 << log_length;
            let mut rng = StdRng::seed_from_u64(seed);
            let mut vector: Vec<F> = random_vector(&mut rng, length);
            let mut extended_vector = zero_pad(&vector);
            let weight = rng.gen::<F>();
            fold(&mut vector, weight);
            fold(&mut extended_vector, weight);
            assert_eq!(vector, extended_vector);
        });
    }

    /// Sumcheck invariant: c(0) + c(1) == dot(a, b).
    #[test]
    fn sumcheck_poly_invariant() {
        proptest!(|(seed:u64, length in 0_usize..500)| {
            let mut rng = StdRng::seed_from_u64(seed);
            let vector: Vec<F> = random_vector(&mut rng, length);
            let covector: Vec<F> = random_vector(&mut rng, length);
            let sum = crate::algebra::dot(&vector, &covector);
            let (c0, c2) = compute_sumcheck_polynomial(&vector, &covector);
            let c1 = sum - c0.double() - c2;
            assert_eq!(c0 + (c0 + c1 + c2), sum);
        });
    }

    /// Fold+sumcheck round trip: after folding, dot(a', b') == c(r).
    #[test]
    fn fold_sumcheck_round_trip() {
        proptest!(|(seed:u64, length in 2_usize..500)| {
            let mut rng = StdRng::seed_from_u64(seed);
            let mut a: Vec<F> = random_vector(&mut rng, length);
            let mut b: Vec<F> = random_vector(&mut rng, length);
            let sum = crate::algebra::dot(&a, &b);
            let (c0, c2) = compute_sumcheck_polynomial(&a, &b);
            let c1 = sum - c0.double() - c2;
            let r = rng.gen::<F>();
            let expected_new_sum = (c2 * r + c1) * r + c0;
            fold(&mut a, r);
            fold(&mut b, r);
            let actual = crate::algebra::dot(&a, &b);
            assert_eq!(actual, expected_new_sum);
        });
    }
}
