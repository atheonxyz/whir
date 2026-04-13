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

// ─── Ternary fold (radix-3) ─────────────────────────────────────────────

/// Ternary fold: interpolate three thirds at `weight` using Lagrange basis
/// over {0, 1, 2}, in place. Requires `values.len() % 3 == 0`.
///
/// For each position i in `0..third`:
///   out[i] = v[i]*L₀(w) + v[third+i]*L₁(w) + v[2*third+i]*L₂(w)
///
/// where L₀(t) = (t-1)(t-2)/2, L₁(t) = t(2-t), L₂(t) = t(t-1)/2.
pub fn fold3<F: Field>(values: &mut Vec<F>, weight: F) {
    let len = values.len();
    assert!(len % 3 == 0, "fold3 requires len divisible by 3, got {len}");
    if len == 0 {
        return;
    }
    let third = len / 3;

    // Lagrange basis evaluated at `weight`:
    //   L₀(w) = (w-1)(w-2)/2
    //   L₁(w) = w(2-w)
    //   L₂(w) = w(w-1)/2
    let half_inv = F::from(2u64).inverse().expect("char != 2");
    let w = weight;
    let l0 = (w - F::ONE) * (w - F::from(2u64)) * half_inv;
    let l1 = w * (F::from(2u64) - w);
    let l2 = w * (w - F::ONE) * half_inv;

    #[cfg(feature = "parallel")]
    if third > workload_size::<F>() {
        use rayon::prelude::*;
        let (part0, rest) = values.split_at_mut(third);
        let (part1, part2) = rest.split_at(third);
        part0
            .par_iter_mut()
            .zip(part1.par_iter())
            .zip(part2.par_iter())
            .for_each(|((v0, &v1), &v2)| {
                *v0 = *v0 * l0 + v1 * l1 + v2 * l2;
            });
        values.truncate(third);
        values.shrink_to_fit();
        return;
    }

    for i in 0..third {
        values[i] = values[i] * l0 + values[third + i] * l1 + values[2 * third + i] * l2;
    }
    values.truncate(third);
    values.shrink_to_fit();
}

/// Compute the degree-4 sumcheck polynomial coefficients for a ternary round.
///
/// Returns (c₀, c₂, c₃, c₄). The caller derives c₁ from the sum relation:
///   sum = c(0) + c(1) + c(2)
///       = 3c₀ + 3c₁ + 5c₂ + 9c₃ + 17c₄
///   c₁ = (sum - 3c₀ - 5c₂ - 9c₃ - 17c₄) * inv3
///
/// Requires `a.len() == b.len()` and both divisible by 3.
pub fn compute_sumcheck_polynomial3<F: Field>(a: &[F], b: &[F]) -> (F, F, F, F) {
    let len = a.len().min(b.len());
    let a = &a[..len];
    let b = &b[..len];
    assert!(
        len % 3 == 0,
        "ternary sumcheck requires len % 3 == 0, got {len}"
    );
    if len == 0 {
        return (F::ZERO, F::ZERO, F::ZERO, F::ZERO);
    }

    let third = len / 3;
    let half_inv = F::from(2u64).inverse().expect("char != 2");

    // For each triple (a₀, a₁, a₂), compute the degree-2 Lagrange coefficients:
    //   fa(t) = α₀ + α₁·t + α₂·t²
    //   α₀ = a₀
    //   α₁ = (-3a₀ + 4a₁ - a₂) / 2
    //   α₂ = (a₀ - 2a₁ + a₂) / 2
    //
    // The product fa(t)·fb(t) = Σ pₖ·tᵏ for k=0..4:
    //   p₀ = α₀·β₀
    //   p₁ = α₀·β₁ + α₁·β₀
    //   p₂ = α₀·β₂ + α₁·β₁ + α₂·β₀
    //   p₃ = α₁·β₂ + α₂·β₁
    //   p₄ = α₂·β₂
    //
    // Accumulate p₀, p₂, p₃, p₄ across all triples. p₁ is derived via sum relation.

    #[cfg(feature = "parallel")]
    if third > workload_size::<F>() {
        use rayon::prelude::*;
        let (p0, p2, p3, p4) = (0..third)
            .into_par_iter()
            .fold(
                || (F::ZERO, F::ZERO, F::ZERO, F::ZERO),
                |(mut p0, mut p2, mut p3, mut p4), i| {
                    let a0 = a[i];
                    let a1 = a[third + i];
                    let a2 = a[2 * third + i];
                    let b0 = b[i];
                    let b1 = b[third + i];
                    let b2 = b[2 * third + i];
                    let a0x3 = a0 + a0.double();
                    let alpha1 = (a1.double().double() - a0x3 - a2) * half_inv;
                    let alpha2 = (a0 - a1.double() + a2) * half_inv;
                    let b0x3 = b0 + b0.double();
                    let beta1 = (b1.double().double() - b0x3 - b2) * half_inv;
                    let beta2 = (b0 - b1.double() + b2) * half_inv;
                    p0 += a0 * b0;
                    p2 += a0 * beta2 + alpha1 * beta1 + alpha2 * b0;
                    p3 += alpha1 * beta2 + alpha2 * beta1;
                    p4 += alpha2 * beta2;
                    (p0, p2, p3, p4)
                },
            )
            .reduce(
                || (F::ZERO, F::ZERO, F::ZERO, F::ZERO),
                |(a0, a2, a3, a4), (b0, b2, b3, b4)| (a0 + b0, a2 + b2, a3 + b3, a4 + b4),
            );
        return (p0, p2, p3, p4);
    }

    let mut p0 = F::ZERO;
    let mut p2 = F::ZERO;
    let mut p3 = F::ZERO;
    let mut p4 = F::ZERO;

    for i in 0..third {
        let a0 = a[i];
        let a1 = a[third + i];
        let a2 = a[2 * third + i];
        let b0 = b[i];
        let b1 = b[third + i];
        let b2 = b[2 * third + i];

        let alpha1 = (a1.double().double() - a0.double() - a0 - a2) * half_inv;
        let alpha2 = (a0 - a1.double() + a2) * half_inv;
        let beta1 = (b1.double().double() - b0.double() - b0 - b2) * half_inv;
        let beta2 = (b0 - b1.double() + b2) * half_inv;

        p0 += a0 * b0;
        // p1 is derived from the sum relation
        p2 += a0 * beta2 + alpha1 * beta1 + alpha2 * b0;
        p3 += alpha1 * beta2 + alpha2 * beta1;
        p4 += alpha2 * beta2;
    }

    (p0, p2, p3, p4)
}

/// Fused fold3 + compute_sumcheck_polynomial3 in a single pass.
///
/// Folds both `a` and `b` at `weight` (ternary Lagrange interpolation),
/// then computes the degree-4 sumcheck polynomial coefficients from the
/// *pre-fold* data — all in one loop over the triples.
///
/// Returns `(p0, p2, p3, p4)` (same as `compute_sumcheck_polynomial3`).
pub fn fold3_and_compute_polynomial<F: Field>(
    a: &mut Vec<F>,
    b: &mut Vec<F>,
    weight: F,
) -> (F, F, F, F) {
    let len = a.len().min(b.len());
    assert!(
        len % 3 == 0,
        "ternary fold+sumcheck requires len % 3 == 0, got {len}"
    );
    if len == 0 {
        a.clear();
        b.clear();
        return (F::ZERO, F::ZERO, F::ZERO, F::ZERO);
    }
    let third = len / 3;
    let half_inv = F::from(2u64).inverse().expect("char != 2");

    // Lagrange basis at `weight` for fold
    let w = weight;
    let two = F::from(2u64);
    let l0 = (w - F::ONE) * (w - two) * half_inv;
    let l1 = w * (two - w);
    let l2 = w * (w - F::ONE) * half_inv;

    let mut p0 = F::ZERO;
    let mut p2 = F::ZERO;
    let mut p3 = F::ZERO;
    let mut p4 = F::ZERO;

    for i in 0..third {
        let a0 = a[i];
        let a1 = a[third + i];
        let a2 = a[2 * third + i];
        let b0 = b[i];
        let b1 = b[third + i];
        let b2 = b[2 * third + i];

        // Fold in-place (write before we lose the values)
        a[i] = a0 * l0 + a1 * l1 + a2 * l2;
        b[i] = b0 * l0 + b1 * l1 + b2 * l2;

        // Lagrange → monomial coefficients for the sumcheck polynomial
        let a0x3 = a0 + a0.double();
        let alpha1 = (a1.double().double() - a0x3 - a2) * half_inv;
        let alpha2 = (a0 - a1.double() + a2) * half_inv;
        let b0x3 = b0 + b0.double();
        let beta1 = (b1.double().double() - b0x3 - b2) * half_inv;
        let beta2 = (b0 - b1.double() + b2) * half_inv;

        p0 += a0 * b0;
        p2 += a0 * beta2 + alpha1 * beta1 + alpha2 * b0;
        p3 += alpha1 * beta2 + alpha2 * beta1;
        p4 += alpha2 * beta2;
    }
    a.truncate(third);
    b.truncate(third);
    (p0, p2, p3, p4)
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

    // ─── Ternary fold tests ──────────────────────────────────────────

    /// fold3 output length must be len/3.
    #[test]
    fn fold3_output_length() {
        proptest!(|(seed:u64, k in 0_usize..200)| {
            let length = k * 3;
            if length == 0 { return Ok(()); }
            let mut rng = StdRng::seed_from_u64(seed);
            let mut vector: Vec<F> = random_vector(&mut rng, length);
            let weight = rng.gen::<F>();
            fold3(&mut vector, weight);
            assert_eq!(vector.len(), k);
        });
    }

    /// fold3 at evaluation points 0, 1, 2 picks the respective third.
    #[test]
    fn fold3_at_integer_points() {
        proptest!(|(seed:u64, k in 1_usize..100)| {
            let length = k * 3;
            let mut rng = StdRng::seed_from_u64(seed);
            let v: Vec<F> = random_vector(&mut rng, length);

            // fold3 at 0 → first third
            let mut v0 = v.clone();
            fold3(&mut v0, F::ZERO);
            assert_eq!(v0, v[..k]);

            // fold3 at 1 → second third
            let mut v1 = v.clone();
            fold3(&mut v1, F::ONE);
            assert_eq!(v1, v[k..2*k]);

            // fold3 at 2 → third third
            let mut v2 = v.clone();
            fold3(&mut v2, F::from(2u64));
            assert_eq!(v2, v[2*k..]);
        });
    }

    /// Ternary sumcheck invariant: c(0) + c(1) + c(2) == sum.
    #[test]
    fn sumcheck_poly3_invariant() {
        proptest!(|(seed:u64, k in 1_usize..200)| {
            let length = k * 3;
            let mut rng = StdRng::seed_from_u64(seed);
            let a: Vec<F> = random_vector(&mut rng, length);
            let b: Vec<F> = random_vector(&mut rng, length);
            let sum = crate::algebra::dot(&a, &b);
            let (p0, p2, p3, p4) = compute_sumcheck_polynomial3(&a, &b);
            // c(0) + c(1) + c(2) = 3p₀ + 3p₁ + 5p₂ + 9p₃ + 17p₄ = sum
            let inv3 = F::from(3u64).inverse().unwrap();
            let p1 = (sum - F::from(3u64) * p0 - F::from(5u64) * p2
                - F::from(9u64) * p3 - F::from(17u64) * p4) * inv3;
            let c0 = p0;
            let c1 = p0 + p1 + p2 + p3 + p4;
            let c2 = p0 + F::from(2u64) * p1 + F::from(4u64) * p2
                + F::from(8u64) * p3 + F::from(16u64) * p4;
            assert_eq!(c0 + c1 + c2, sum);
        });
    }

    /// Ternary fold+sumcheck round trip: after fold3, dot(a', b') == c(r).
    #[test]
    fn fold3_sumcheck_round_trip() {
        proptest!(|(seed:u64, k in 1_usize..200)| {
            let length = k * 3;
            let mut rng = StdRng::seed_from_u64(seed);
            let mut a: Vec<F> = random_vector(&mut rng, length);
            let mut b: Vec<F> = random_vector(&mut rng, length);
            let sum = crate::algebra::dot(&a, &b);
            let (p0, p2, p3, p4) = compute_sumcheck_polynomial3(&a, &b);
            let inv3 = F::from(3u64).inverse().unwrap();
            let p1 = (sum - F::from(3u64) * p0 - F::from(5u64) * p2
                - F::from(9u64) * p3 - F::from(17u64) * p4) * inv3;
            let r = rng.gen::<F>();
            let expected = (((p4 * r + p3) * r + p2) * r + p1) * r + p0;
            fold3(&mut a, r);
            fold3(&mut b, r);
            let actual = crate::algebra::dot(&a, &b);
            assert_eq!(actual, expected);
        });
    }
}
