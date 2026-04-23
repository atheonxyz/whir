//! Quadratic sumcheck protocol with mixed binary/ternary rounds.
//!
//! When the current size is divisible by 3, a ternary round is used
//! (degree-4 polynomial, size → size/3). Otherwise a binary round is
//! used (degree-2 polynomial, size → ceil(size/2)).

use std::fmt;

use ark_ff::Field;
use ark_std::rand::{CryptoRng, RngCore};
use serde::{Deserialize, Serialize};
#[cfg(feature = "tracing")]
use tracing::instrument;

use crate::{
    algebra::{
        dot,
        sumcheck::{
            compute_sumcheck_polynomial, compute_sumcheck_polynomial3, fold, fold3,
            fold_and_compute_polynomial,
        },
        univariate_evaluate,
    },
    protocols::proof_of_work,
    transcript::{
        codecs::U64, Codec, Decoding, DuplexSpongeInterface, ProverState, VerificationResult,
        VerifierMessage, VerifierState,
    },
    type_info::Type,
    utils::chunks_exact_or_empty,
};

/// Returns true if this round should be a ternary fold (size divisible by 3).
#[inline]
pub fn is_ternary_round(size: usize) -> bool {
    size >= 3 && size % 3 == 0
}

/// Compute the size after one fold step: ternary if divisible by 3, else binary.
#[inline]
pub fn fold_one_step(size: usize) -> usize {
    if size <= 1 {
        size
    } else if is_ternary_round(size) {
        size / 3
    } else {
        (size + 1) / 2
    }
}

/// Compute the size after `rounds` of binary-only folding.
///
/// Each fold: `size → ceil(size / 2)`.
pub fn fold_n_times_binary(initial_size: usize, rounds: usize) -> usize {
    let mut size = initial_size;
    for _ in 0..rounds {
        if size <= 1 {
            break;
        }
        size = (size + 1) / 2;
    }
    size
}

/// Compute the size of a vector after `rounds` of mixed folding.
///
/// Each round: if size % 3 == 0, ternary fold (size → size/3);
/// otherwise binary fold (size → ceil(size/2)).
pub fn fold_n_times(initial_size: usize, rounds: usize) -> usize {
    let mut size = initial_size;
    for _ in 0..rounds {
        if size <= 1 {
            break;
        }
        size = fold_one_step(size);
    }
    size
}

/// Fold a vector one step using the mixed ternary/binary schedule.
/// Returns the new size after folding.
#[inline]
pub fn fold_mixed<F: Field>(v: &mut Vec<F>, r: F, size: usize) -> usize {
    use crate::algebra::sumcheck::{fold, fold3};
    if is_ternary_round(size) {
        fold3(v, r);
        size / 3
    } else {
        fold(v, r);
        (size + 1) / 2
    }
}

/// Count the number of rounds needed to reduce `size` to 1 using mixed folding.
pub fn rounds_to_one(size: usize) -> usize {
    let mut s = size;
    let mut rounds = 0;
    while s > 1 {
        s = fold_one_step(s);
        rounds += 1;
    }
    rounds
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(bound = "")]
pub struct Config<F>
where
    F: Field,
{
    pub field: Type<F>,
    pub initial_size: usize,
    pub round_pow: proof_of_work::Config,
    pub num_rounds: usize,
    pub mask_length: usize,
    /// When true, use mixed ternary/binary folding (ternary when size % 3 == 0).
    /// When false, use binary-only folding (the default for non-final sumchecks).
    #[serde(default)]
    pub ternary: bool,
}

/// Field constants reused across sumcheck rounds.
struct RoundConsts<F: Field> {
    half: F,
    inv3: F,
    three: F,
    five: F,
    nine: F,
    seventeen: F,
}

impl<F: Field> RoundConsts<F> {
    fn new() -> Self {
        Self {
            half: F::from(2u64).inverse().unwrap(),
            inv3: F::from(3u64).inverse().expect("char != 3"),
            three: F::from(3u64),
            five: F::from(5u64),
            nine: F::from(9u64),
            seventeen: F::from(17u64),
        }
    }
}

/// Assemble a masked univariate U(t) = mask_rlc · c(t) + β_r · s(t) + K, where
/// K is chosen so that U summed over the fold points equals `mask_sum`.
///
/// * `out` is written to length `max(poly.len(), mask.len())`; higher slots zeroed.
/// * `mask_sum_over_fold_pts` is `eval_01(mask)` for binary rounds and
///   `eval_012(mask)` for ternary rounds.
/// * `inv_fold_count` is `1/2` for binary rounds, `1/3` for ternary.
fn apply_mask<F: Field>(
    out: &mut [F],
    mask: &[F],
    mask_sum: F,
    mask_rlc: F,
    beta_r: F,
    poly: &[F],
    mask_sum_over_fold_pts: F,
    inv_fold_count: F,
) {
    for (slot, m) in out.iter_mut().zip(mask.iter()) {
        *slot = beta_r * *m;
    }
    for slot in out.iter_mut().skip(mask.len()) {
        *slot = F::ZERO;
    }
    out[0] += (mask_sum - beta_r * mask_sum_over_fold_pts) * inv_fold_count;
    for (i, &c) in poly.iter().enumerate() {
        out[i] += mask_rlc * c;
    }
}

impl<F: Field> Config<F> {
    /// Compute the vector size after `num_rounds` of folding.
    ///
    /// When `ternary` is true, uses mixed binary/ternary schedule.
    /// Otherwise, each fold: `size → ceil(size / 2)`.
    pub fn final_size(&self) -> usize {
        if self.ternary {
            fold_n_times(self.initial_size, self.num_rounds)
        } else {
            fold_n_times_binary(self.initial_size, self.num_rounds)
        }
    }

    /// Runs the sumcheck protocol with mixed binary/ternary rounds.
    ///
    /// It reduces a claim of the form `dot(a, b) == sum` to an exponentially
    /// smaller claim `dot(a', b') == sum'` where `a'` is `a` folded in place
    /// and similarly for `b`.
    ///
    /// When the current vector size is divisible by 3, a ternary round is used
    /// (degree-4 polynomial). Otherwise a binary round (degree-2 polynomial).
    ///
    /// This function:
    /// - Samples random values to progressively reduce the polynomial.
    /// - Applies proof-of-work grinding if required.
    /// - Returns the sampled folding randomness values used in each reduction step.
    #[cfg_attr(feature = "tracing", instrument(skip_all))]
    pub fn prove<H, R>(
        &self,
        prover_state: &mut ProverState<H, R>,
        a: &mut Vec<F>,
        b: &mut Vec<F>,
        sum: &mut F,
        masks: &[F],
    ) -> (Vec<F>, F, F)
    where
        H: DuplexSpongeInterface,
        R: CryptoRng + RngCore,
        F: Codec<[H::U]>,
        [u8; 32]: Decoding<[H::U]>,
        U64: Codec<[H::U]>,
    {
        assert!(
            self.num_rounds == 0 || self.final_size() >= 1,
            "too many rounds ({}) for initial_size {}",
            self.num_rounds,
            self.initial_size,
        );
        let has_ternary_round = self.ternary
            && fold_schedule(self.initial_size, self.num_rounds, true)
                .iter()
                .any(|&f| f == 3);
        if has_ternary_round {
            assert!(
                self.mask_length == 0 || self.mask_length >= 5,
                "ternary sumcheck masking requires mask_length >= 5 (got {})",
                self.mask_length,
            );
        } else {
            assert!(self.mask_length == 0 || self.mask_length >= 3);
        }
        assert_eq!(a.len(), self.initial_size);
        assert_eq!(b.len(), self.initial_size);
        debug_assert_eq!(dot(a, b), *sum);
        assert_eq!(masks.len(), self.num_rounds * self.mask_length);
        let k = RoundConsts::<F>::new();

        // Precompute the fold schedule and per-round β scaling.
        let schedule = fold_schedule(self.initial_size, self.num_rounds, self.ternary);
        // β_r = ∏_{j>r} ff(j) as field element.
        let mut beta: Vec<F> = vec![F::ONE; self.num_rounds];
        for r in (0..self.num_rounds).rev() {
            let next = r + 1;
            if next < self.num_rounds {
                beta[r] = beta[next] * F::from(schedule[next] as u64);
            }
        }

        // Send mask sum and get combination randomness.
        let mut mask_sum = F::ZERO;
        let mut mask_rlc = F::ONE;
        if !masks.is_empty() {
            // Initial mask_sum = Σ_r prefix_r · β_r · eval_fold_r(s_r),
            // where prefix_r = ∏_{j<r} ff(j), β_r = ∏_{j>r} ff(j), and
            // eval_fold_r sums s_r over its round's fold points (eval_01 or eval_012).
            let mut prefix = F::ONE;
            for (r, m) in masks.chunks_exact(self.mask_length).enumerate() {
                let e = if schedule[r] == 3 { eval_012(m) } else { eval_01(m) };
                mask_sum += prefix * beta[r] * e;
                prefix *= F::from(schedule[r] as u64);
            }
            prover_state.prover_message(&mask_sum);
            mask_rlc = prover_state.verifier_message();
        }

        // Univariate buffer: ternary → degree 4 (5 coeffs); binary only → degree 2 (3 coeffs).
        let max_deg = if self.ternary { 4 } else { 2 };
        let buf_len = (max_deg + 1).max(self.mask_length);
        let mut univariate: Vec<F> = vec![F::ZERO; buf_len];

        let mut res = Vec::with_capacity(self.num_rounds);
        let mut current_size = self.initial_size;
        // Track: (is_ternary, folding_randomness) from previous round.
        let mut prev_fold: Option<(bool, F)> = None;
        for (round, mask) in
            chunks_exact_or_empty(masks, self.mask_length, self.num_rounds).enumerate()
        {
            let ternary = self.ternary && is_ternary_round(current_size);
            let r = if ternary {
                // Apply any pending fold from the previous round before computing
                // the ternary polynomial (ternary has no fused fold-and-compute).
                if let Some((prev_ternary, w)) = prev_fold {
                    if prev_ternary { fold3(a, w); fold3(b, w); }
                    else { fold(a, w); fold(b, w); }
                }
                let r = self.prove_ternary_round(
                    prover_state, a, b, sum, mask,
                    &mut mask_sum, mask_rlc, beta[round], &mut univariate, &k,
                );
                current_size /= 3;
                r
            } else {
                // Binary rounds can fuse the previous binary fold into the polynomial compute.
                let (c0, c2) = match prev_fold {
                    Some((true, w)) => {
                        fold3(a, w); fold3(b, w);
                        compute_sumcheck_polynomial(a, b)
                    }
                    Some((false, w)) => fold_and_compute_polynomial(a, b, w),
                    None => compute_sumcheck_polynomial(a, b),
                };
                let r = self.prove_binary_round(
                    prover_state, (c0, c2), sum, mask,
                    &mut mask_sum, mask_rlc, beta[round], &mut univariate, &k,
                );
                current_size = (current_size + 1) / 2;
                r
            };
            res.push(r);
            prev_fold = Some((ternary, r));
        }
        if let Some((prev_ternary, w)) = prev_fold {
            if prev_ternary { fold3(a, w); fold3(b, w); }
            else { fold(a, w); fold(b, w); }
        }

        *sum = mask_sum + mask_rlc * *sum;
        (res, mask_sum, mask_rlc)
    }

    /// One ternary round: compute degree-4 polynomial, send it (optionally
    /// combined with a mask polynomial), consume PoW, sample folding
    /// randomness, and update `*sum` (and `*mask_sum` when masking is on).
    #[allow(clippy::too_many_arguments)]
    fn prove_ternary_round<H, R>(
        &self,
        prover_state: &mut ProverState<H, R>,
        a: &[F],
        b: &[F],
        sum: &mut F,
        mask: &[F],
        mask_sum: &mut F,
        mask_rlc: F,
        beta_r: F,
        univariate: &mut [F],
        k: &RoundConsts<F>,
    ) -> F
    where
        H: DuplexSpongeInterface,
        R: CryptoRng + RngCore,
        F: Codec<[H::U]>,
        [u8; 32]: Decoding<[H::U]>,
        U64: Codec<[H::U]>,
    {
        let (p0, p2, p3, p4) = compute_sumcheck_polynomial3(a, b);
        // c(0)+c(1)+c(2) = sum  ⇒  3p₀+3p₁+5p₂+9p₃+17p₄ = sum
        let p1 = (*sum - k.three * p0 - k.five * p2 - k.nine * p3 - k.seventeen * p4) * k.inv3;

        if mask.is_empty() {
            prover_state.prover_messages(&[p0, p2, p3, p4]);
        } else {
            apply_mask(
                univariate, mask, *mask_sum, mask_rlc, beta_r,
                &[p0, p1, p2, p3, p4], eval_012(mask), k.inv3,
            );
            prover_state.prover_message(&univariate[0]);
            prover_state.prover_messages(&univariate[2..]);
        }

        self.round_pow.prove(prover_state);
        let r = prover_state.verifier_message::<F>();
        *sum = (((p4 * r + p3) * r + p2) * r + p1) * r + p0;
        if !mask.is_empty() {
            *mask_sum = univariate_evaluate(univariate, r) - mask_rlc * *sum;
        }
        r
    }

    /// One binary round: take the precomputed `(c0, c2)` degree-2 polynomial
    /// (possibly fused with the previous fold by the caller), send it
    /// (optionally combined with a mask polynomial), consume PoW, sample
    /// folding randomness, and update `*sum` (and `*mask_sum` when masking is on).
    #[allow(clippy::too_many_arguments)]
    fn prove_binary_round<H, R>(
        &self,
        prover_state: &mut ProverState<H, R>,
        (c0, c2): (F, F),
        sum: &mut F,
        mask: &[F],
        mask_sum: &mut F,
        mask_rlc: F,
        beta_r: F,
        univariate: &mut [F],
        k: &RoundConsts<F>,
    ) -> F
    where
        H: DuplexSpongeInterface,
        R: CryptoRng + RngCore,
        F: Codec<[H::U]>,
        [u8; 32]: Decoding<[H::U]>,
        U64: Codec<[H::U]>,
    {
        let c1 = *sum - c0.double() - c2;

        if mask.is_empty() {
            prover_state.prover_messages(&[c0, c2]);
        } else {
            apply_mask(
                univariate, mask, *mask_sum, mask_rlc, beta_r,
                &[c0, c1, c2], eval_01(mask), k.half,
            );
            prover_state.prover_message(&univariate[0]);
            prover_state.prover_messages(&univariate[2..]);
        }

        self.round_pow.prove(prover_state);
        let r = prover_state.verifier_message::<F>();
        *sum = (c2 * r + c1) * r + c0;
        if !mask.is_empty() {
            *mask_sum = univariate_evaluate(univariate, r) - mask_rlc * *sum;
        }
        r
    }

    #[cfg_attr(feature = "tracing", instrument(skip_all))]
    pub fn verify<H>(
        &self,
        verifier_state: &mut VerifierState<H>,
        sum: &mut F,
    ) -> VerificationResult<(Vec<F>, F)>
    where
        H: DuplexSpongeInterface,
        F: Codec<[H::U]>,
        [u8; 32]: Decoding<[H::U]>,
        U64: Codec<[H::U]>,
    {
        assert!(
            self.num_rounds == 0 || self.final_size() >= 1,
            "too many rounds ({}) for initial_size {}",
            self.num_rounds,
            self.initial_size,
        );
        let has_ternary_round = self.ternary
            && fold_schedule(self.initial_size, self.num_rounds, true)
                .iter()
                .any(|&f| f == 3);
        if has_ternary_round {
            assert!(
                self.mask_length == 0 || self.mask_length >= 5,
                "ternary sumcheck masking requires mask_length >= 5 (got {})",
                self.mask_length,
            );
        } else {
            assert!(self.mask_length == 0 || self.mask_length >= 3);
        }

        let mut mask_rlc = F::ONE;
        if self.mask_length > 0 && self.num_rounds > 0 {
            let mask_sum: F = verifier_state.prover_message()?;
            mask_rlc = verifier_state.verifier_message();
            *sum = mask_sum + mask_rlc * *sum;
        }

        let k = RoundConsts::<F>::new();
        let max_deg = if self.ternary { 4 } else { 2 };
        let buf_len = (max_deg + 1).max(self.mask_length);
        let mut univariate = vec![F::ZERO; buf_len];
        let has_mask = self.mask_length > 0;
        let mut res = Vec::with_capacity(self.num_rounds);
        let mut current_size = self.initial_size;
        for _ in 0..self.num_rounds {
            let ternary = self.ternary && is_ternary_round(current_size);
            let r = if ternary {
                let r = self.verify_ternary_round(verifier_state, sum, has_mask, &mut univariate, &k)?;
                current_size /= 3;
                r
            } else {
                let r = self.verify_binary_round(verifier_state, sum, has_mask, &mut univariate)?;
                current_size = (current_size + 1) / 2;
                r
            };
            res.push(r);
        }
        Ok((res, mask_rlc))
    }

    /// One ternary verify round: receive the round polynomial (4 coeffs unmasked,
    /// or U₀ + U₂..U_{buf-1} masked), derive the missing coefficient from the
    /// sum-over-fold-points check, sample randomness, and update `*sum`.
    fn verify_ternary_round<H>(
        &self,
        verifier_state: &mut VerifierState<H>,
        sum: &mut F,
        has_mask: bool,
        univariate: &mut [F],
        k: &RoundConsts<F>,
    ) -> VerificationResult<F>
    where
        H: DuplexSpongeInterface,
        F: Codec<[H::U]>,
        [u8; 32]: Decoding<[H::U]>,
        U64: Codec<[H::U]>,
    {
        if !has_mask {
            let p0: F = verifier_state.prover_message()?;
            let p2: F = verifier_state.prover_message()?;
            let p3: F = verifier_state.prover_message()?;
            let p4: F = verifier_state.prover_message()?;
            let p1 = (*sum - k.three * p0 - k.five * p2 - k.nine * p3 - k.seventeen * p4) * k.inv3;
            self.round_pow.verify(verifier_state)?;
            let r = verifier_state.verifier_message::<F>();
            *sum = (((p4 * r + p3) * r + p2) * r + p1) * r + p0;
            Ok(r)
        } else {
            univariate[0] = verifier_state.prover_message()?;
            for c in &mut univariate[2..] {
                *c = verifier_state.prover_message()?;
            }
            // Derive U₁ from U(0)+U(1)+U(2) = sum.
            // Coefficient of U_i in U(0)+U(1)+U(2): 3 if i∈{0,1}; else 1+2^i.
            let two = F::from(2u64);
            let mut tail = k.three * univariate[0];
            // t2_pow tracks 2^offset. Skip starts at offset=2, init to 2,
            // bump before use → 2^2=4 first iter, 2^3=8 next, etc.
            let mut t2_pow = two;
            for &c in univariate.iter().skip(2) {
                t2_pow *= two;
                tail += (F::ONE + t2_pow) * c;
            }
            univariate[1] = (*sum - tail) * k.inv3;
            self.round_pow.verify(verifier_state)?;
            let r = verifier_state.verifier_message::<F>();
            *sum = univariate_evaluate(univariate, r);
            Ok(r)
        }
    }

    /// One binary verify round: receive the round polynomial (2 coeffs unmasked,
    /// or U₀ + U₂..U_{buf-1} masked), derive the linear coefficient from the
    /// sum-over-fold-points check, sample randomness, and update `*sum`.
    fn verify_binary_round<H>(
        &self,
        verifier_state: &mut VerifierState<H>,
        sum: &mut F,
        has_mask: bool,
        univariate: &mut [F],
    ) -> VerificationResult<F>
    where
        H: DuplexSpongeInterface,
        F: Codec<[H::U]>,
        [u8; 32]: Decoding<[H::U]>,
        U64: Codec<[H::U]>,
    {
        if !has_mask {
            let c0: F = verifier_state.prover_message()?;
            let c2: F = verifier_state.prover_message()?;
            let c1 = *sum - c0.double() - c2;
            self.round_pow.verify(verifier_state)?;
            let r = verifier_state.verifier_message::<F>();
            *sum = (c2 * r + c1) * r + c0;
            Ok(r)
        } else {
            univariate[0] = verifier_state.prover_message()?;
            for c in &mut univariate[2..] {
                *c = verifier_state.prover_message()?;
            }
            univariate[1] = *sum - univariate[0].double() - univariate[2..].iter().sum::<F>();
            self.round_pow.verify(verifier_state)?;
            let r = verifier_state.verifier_message::<F>();
            *sum = univariate_evaluate(univariate, r);
            Ok(r)
        }
    }
}

impl<F: Field> fmt::Display for Config<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "size {} rounds {} pow {:.2} ℓ_zk {}",
            self.initial_size,
            self.num_rounds,
            self.round_pow.difficulty(),
            self.mask_length
        )
    }
}

// Evaluated a univariate as p(0) + p(1)
fn eval_01<F: Field>(coefficients: &[F]) -> F {
    if coefficients.is_empty() {
        return F::ZERO;
    }
    coefficients[0] + coefficients.iter().sum::<F>()
}

// Evaluate a univariate as p(0) + p(1) + p(2).
// For p(t) = Σ_{i≤4} p_i t^i this is 3p₀ + 3p₁ + 5p₂ + 9p₃ + 17p₄.
fn eval_012<F: Field>(coefficients: &[F]) -> F {
    if coefficients.is_empty() {
        return F::ZERO;
    }
    let p_at_0 = coefficients[0];
    let p_at_1 = coefficients.iter().copied().sum::<F>();
    let mut p_at_2 = F::ZERO;
    let mut two_pow = F::ONE;
    let two = F::from(2u64);
    for &c in coefficients {
        p_at_2 += c * two_pow;
        two_pow *= two;
    }
    p_at_0 + p_at_1 + p_at_2
}

/// Compute per-round fold factors for a mixed-radix sumcheck schedule.
///
/// For each of `num_rounds`, returns 2 if the round is binary or 3 if ternary.
/// Starts from `initial_size` and follows `is_ternary_round` dispatch (gated
/// by `ternary_enabled`).
pub fn fold_schedule(initial_size: usize, num_rounds: usize, ternary_enabled: bool) -> Vec<usize> {
    let mut out = Vec::with_capacity(num_rounds);
    let mut s = initial_size;
    for _ in 0..num_rounds {
        if ternary_enabled && is_ternary_round(s) {
            out.push(3);
            s /= 3;
        } else {
            out.push(2);
            s = (s + 1) / 2;
        }
    }
    out
}

#[cfg(test)]
mod tests {
    // TODO: Proptest based tests checking invariants and post conditions.
    use ark_std::rand::{
        distributions::{Distribution, Standard},
        rngs::StdRng,
        SeedableRng,
    };
    use proptest::{prelude::Just, prop_oneof, proptest, strategy::Strategy};
    #[cfg(feature = "tracing")]
    use tracing::instrument;

    use super::*;
    use crate::{
        algebra::{
            fields::{self, Field64},
            random_vector,
        },
        transcript::DomainSeparator,
        utils::zip_strict,
    };

    impl<F: Field> Config<F> {
        pub fn arbitrary() -> impl Strategy<Value = Self> {
            let mask_length = prop_oneof![
                3 => Just(0_usize),
                7 => 3_usize..100,
            ];
            (0_usize..(1 << 12), 0_usize..16, mask_length).prop_map(
                |(initial_size, num_rounds, mask_length)| {
                    let max_rounds = rounds_to_one(initial_size);
                    let num_rounds = num_rounds.min(max_rounds);
                    // Ternary sumcheck masking needs mask_length >= 5. If any round
                    // in this config would be ternary, bump mask_length to at least 5.
                    let has_ternary = {
                        let mut s = initial_size;
                        let mut found = false;
                        for _ in 0..num_rounds {
                            if is_ternary_round(s) {
                                found = true;
                                break;
                            }
                            s = fold_one_step(s);
                        }
                        found
                    };
                    let mask_length = if has_ternary && mask_length != 0 {
                        mask_length.max(5)
                    } else {
                        mask_length
                    };
                    Self {
                        field: Type::new(),
                        initial_size,
                        num_rounds,
                        round_pow: proof_of_work::Config::none(),
                        mask_length,
                        ternary: true,
                    }
                },
            )
        }
    }

    #[cfg_attr(feature = "tracing", instrument)]
    fn test_config<F>(seed: u64, config: &Config<F>)
    where
        F: Field + Codec,
        Standard: Distribution<F>,
    {
        // Pseudo-random Instance
        let instance = U64(seed);
        let ds = DomainSeparator::protocol(config)
            .session(&format!("Test at {}:{}", file!(), line!()))
            .instance(&instance);
        let mut rng = StdRng::seed_from_u64(seed);
        let initial_vector = random_vector(&mut rng, config.initial_size);
        let initial_covector = random_vector(&mut rng, config.initial_size);
        let initial_sum = dot(&initial_vector, &initial_covector);
        let masks = random_vector(&mut rng, config.mask_length * config.num_rounds);

        // Prover
        let mut vector = initial_vector.clone();
        let mut covector = initial_covector.clone();
        let mut sum = initial_sum;
        let mut prover_state = ProverState::new_std(&ds);
        let (point, mask_sum, mask_rlc) = config.prove(
            &mut prover_state,
            &mut vector,
            &mut covector,
            &mut sum,
            &masks,
        );
        let expected_mask_sum = zip_strict(
            chunks_exact_or_empty(&masks, config.mask_length, config.num_rounds),
            &point,
        )
        .map(|(m, x)| univariate_evaluate(m, *x))
        .sum::<F>();
        assert_eq!(
            vector.len(),
            config.final_size(),
            "vector.len()={} != final_size()={} for initial_size={} num_rounds={}",
            vector.len(),
            config.final_size(),
            config.initial_size,
            config.num_rounds
        );
        assert_eq!(covector.len(), config.final_size());
        assert_eq!(mask_sum, expected_mask_sum);
        assert_eq!(mask_sum + mask_rlc * dot(&vector, &covector), sum);
        if config.final_size() == 1 {
            // Verify by replaying the fold on the initial vector.
            let mut v_check = initial_vector.clone();
            let mut sz = config.initial_size;
            for &r in &point {
                sz = fold_mixed(&mut v_check, r, sz);
            }
            assert_eq!(v_check.len(), 1);
            assert_eq!(v_check[0], vector[0]);

            let mut c_check = initial_covector.clone();
            sz = config.initial_size;
            for &r in &point {
                sz = fold_mixed(&mut c_check, r, sz);
            }
            assert_eq!(c_check[0], covector[0]);
        } else {
            // TODO: Check correct folding.
        }
        let proof = prover_state.proof();

        // Verifier
        let mut verifier_sum = initial_sum;
        let mut verifier_state = VerifierState::new_std(&ds, &proof);
        let (verifier_point, verifier_mask_rlc) = config
            .verify(&mut verifier_state, &mut verifier_sum)
            .unwrap();
        assert_eq!(verifier_point, point);
        assert_eq!(verifier_mask_rlc, mask_rlc);
        assert_eq!(verifier_sum, sum);
        verifier_state.check_eof().unwrap();
    }

    fn test_sumcheck<F>()
    where
        F: Field + Codec,
        Standard: Distribution<F>,
    {
        crate::tests::init();
        proptest!(|(seed: u64, config in Config::arbitrary())| {
            test_config(seed, &config);
        });
    }

    #[test]
    fn test_single_round() {
        test_config(
            0,
            &Config::<Field64> {
                field: Type::new(),
                initial_size: 2,
                round_pow: proof_of_work::Config::none(),
                num_rounds: 1,
                mask_length: 3,
                ternary: true,
            },
        );
    }

    #[test]
    fn test_two_rounds() {
        // initial_size=3 triggers a ternary round (3→1 in 1 round).
        test_config(
            0,
            &Config::<Field64> {
                field: Type::new(),
                initial_size: 3,
                round_pow: proof_of_work::Config::none(),
                num_rounds: 1, // rounds_to_one(3) = 1 with ternary fold
                mask_length: 0,
                ternary: true,
            },
        );
    }

    #[test]
    fn test_three_rounds() {
        // initial_size=5 → 3 → 1: round 1 binary (5→3), round 2 ternary (3→1).
        test_config(
            0,
            &Config::<Field64> {
                field: Type::new(),
                initial_size: 5,
                round_pow: proof_of_work::Config::none(),
                num_rounds: 2, // rounds_to_one(5) = 2 with mixed folding
                mask_length: 0,
                ternary: true,
            },
        );
    }

    #[test]
    fn test_field64_1() {
        test_sumcheck::<fields::Field64>();
    }

    #[test]
    #[ignore = "Somewhat expensive and redundant"]
    fn test_field64_2() {
        test_sumcheck::<fields::Field64_2>();
    }

    #[test]
    #[ignore = "Somewhat expensive and redundant"]
    fn test_field64_3() {
        test_sumcheck::<fields::Field64_3>();
    }

    #[test]
    #[ignore = "Somewhat expensive and redundant"]
    fn test_field128() {
        test_sumcheck::<fields::Field128>();
    }

    #[test]
    #[ignore = "Somewhat expensive and redundant"]
    fn test_field192() {
        test_sumcheck::<fields::Field192>();
    }

    #[test]
    #[ignore = "Somewhat expensive and redundant"]
    fn test_field256() {
        test_sumcheck::<fields::Field256>();
    }
}
