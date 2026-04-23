# Smooth-domain support with mixed binary/ternary sumcheck folding

## Background: what WHIR is doing

WHIR is a polynomial commitment scheme. The prover commits to a vector of N field elements (treated as the evaluations of a multilinear polynomial on the boolean hypercube), and later proves statements of the form "the polynomial evaluates to v at point r" without revealing the vector.

Internally, WHIR works by interpreting the vector as a Reed-Solomon codeword over an NTT-friendly domain, running a sumcheck to reduce the claim, folding the codeword in half (or thirds), and repeating. Each round halves the problem size. After enough rounds the codeword is small enough to send in full, and the verifier spot-checks consistency against Merkle commitments.

Two mathematical tricks underpin this:

- **The sumcheck protocol** reduces a claim of the form "the sum of f(x) over all x ∈ {0,1}^k equals S" to a claim about f at a single point. Each round sends a small univariate polynomial and receives a random challenge; the sum gets restricted to one fewer variable. After k rounds, the verifier has a single-point evaluation claim.

- **The NTT (Number Theoretic Transform)** is the finite-field analog of the FFT. It lets us move between coefficient and evaluation representations of polynomials in O(n log n). NTTs work cleanly only on domains whose size divides `p-1` for the field characteristic p — for BN254 this means sizes like 2^a, or more generally 2^a × 3^b × 5^c × ...

## Problem

WHIR currently requires the committed vector size to be a power of two. When a witness has N elements where N is not a power of two, the prover must zero-pad to the next power of two before committing. For sizes like 3 × 2^a, this wastes up to 33% of the domain — the prover runs the full NTT, Merkle tree, and sumcheck over those padding zeros.

Concretely: a witness of size 6,291,456 (= 3 × 2^21) gets padded to 8,388,608. The prover does 8.4M elements of work when only 6.3M carry real data.

## Solution overview

This PR enables WHIR to operate on smooth-{2,3} sizes (N = 2^a × 3^b) directly, without padding. The commitment domain matches the witness size exactly.

Making this work required three layers of changes, because the power-of-two assumption was baked in at multiple levels: the NTT, the sumcheck folding schedule, and the verifier's constraint evaluation.

### Layer 1: smooth NTT domain support (commit `0654a4d`)

**Macro view:** teach all the data-processing primitives to handle "weird" sizes like 3 × 2^a instead of just powers of two.

The NTT at smooth sizes was already available (Cooley-Tukey with coset support). What needed changing was everything that *consumed* the NTT output: the sumcheck fold, the polynomial reconstruction, the Merkle tree indexing.

Specific changes:
- `Config::new` now accepts smooth-{2,3} sizes instead of rejecting them.
- `fold()` in `algebra/sumcheck.rs` generalized to split at `len/2` (integer division) instead of the next power-of-two boundary. For odd sizes this creates an asymmetric split with a "paired" middle, a "low tail", and a "high tail" that get handled separately.
- `compute_sumcheck_polynomial()` reworked to produce correct coefficients across this asymmetric split.
- Added `smooth_domain.rs` with helpers: `is_smooth()`, `odd_part()`, `extra_rounds()`.
- Fixed `challenge_indices.rs` so random query indices aren't biased when `num_leaves` isn't a power of two (we added 8 extra bytes of entropy for modular reduction).

### Layer 2: verifier optimization (commit `4cb6af4`)

**Macro view:** the verifier was doing O(n) work to evaluate constraints on the smooth domain; we replaced this with O(log² n).

Inside WHIR, each round accumulates a "constraint" — essentially, an evaluation of some structured polynomial (like `f(x) = x^i`) at the current folding randomness. For power-of-two sizes, these structured evaluations have a closed-form tensor product identity that the verifier can compute in O(log n) field operations.

For smooth sizes, the tensor identity breaks because the fold boundaries don't align to bit positions. A naive verifier would materialize the full weight vector (O(n) memory) and fold it down (O(n) work). That defeats the purpose of succinct verification.

Specific additions:
- `sum_x_fold_eq()` recursively computes `Σ x^i · fold_eq(i)` where `fold_eq(i)` is the coefficient position i receives after all fold rounds. The recursion does one `x^half` exponentiation per level; depth is O(log n), per-level cost is O(log n), total O(log² n).
- `fold_coeff()` traces a single position through all binary fold rounds in O(k) where k is the number of remaining challenges.

### Layer 3: mixed binary/ternary sumcheck (commit `e4dfe9d`)

**Macro view:** when the remaining vector size is divisible by 3, fold by thirds instead of halves. Reduces round count by ~37% for the 3^b residual portion.

After the NTT-based folding reduces the committed domain to its smooth residual (e.g., a factor of 3^b remains), the final sumcheck still has to fold this down to size 1. With only binary folding, this takes ceil(log₂(3^b)) rounds. Ternary folding reduces 3^b in exactly b rounds.

**How a sumcheck round works (the big picture):**

In each round, the prover and verifier have a current "sum" S claimed to equal `Σ f(x) · g(x)` over some domain of size n. The prover sends a polynomial c(t) that encodes what the sum would be if we fixed the next variable to t:

```
c(0) = sum over the "0 half" of the domain
c(1) = sum over the "1 half" of the domain
c(0) + c(1) = S   (consistency check the verifier does)
```

The verifier samples a random challenge r, and the new claim becomes `S' = c(r)`, with the domain reduced by half. Repeat until the domain has size 1.

In the **binary** version, c(t) has degree 2 (product of two linear functions of t). The prover sends the coefficients (c₀, c₂) — 2 field elements — and the verifier derives c₁ from the sum consistency check.

In the **ternary** version, we fold by 3 instead of 2. Evaluation at three points {0, 1, 2} instead of {0, 1}. The round polynomial c(t) now has degree 4 (product of two degree-2 Lagrange interpolants). The prover sends (p₀, p₂, p₃, p₄) — 4 field elements. The verifier derives p₁ from the consistency relation:

```
c(0) + c(1) + c(2) = S
```

Expanding in the monomial basis:
```
c(0) = p₀
c(1) = p₀ + p₁ + p₂ + p₃ + p₄
c(2) = p₀ + 2p₁ + 4p₂ + 8p₃ + 16p₄
Sum  = 3p₀ + 3p₁ + 5p₂ + 9p₃ + 17p₄ = S
```

So `p₁ = (S - 3p₀ - 5p₂ - 9p₃ - 17p₄) / 3`. The verifier plugs in the challenge, evaluates c(r) using all 5 coefficients, and gets the new sum.

**Prover implementation:**
- `fold3()` — ternary Lagrange interpolation at weight w, using basis polynomials L₀(w) = (w-1)(w-2)/2, L₁(w) = w(2-w), L₂(w) = w(w-1)/2. These satisfy Lᵢ(j) = δᵢⱼ on the evaluation points {0, 1, 2}.
- `compute_sumcheck_polynomial3()` — takes the current a, b vectors (split into thirds a₀|a₁|a₂ and b₀|b₁|b₂) and computes (p₀, p₂, p₃, p₄). Parallel via rayon for large inputs.
- `fold3_and_compute_polynomial()` — fused version that folds and computes the next round's polynomial in a single pass over the data (saves one memory pass).
- The sumcheck prover loop dynamically selects ternary vs binary per round: ternary when `size >= 3 && size % 3 == 0`, binary otherwise.

**Verifier implementation:**
- Receives (p₀, p₂, p₃, p₄), reconstructs p₁ using the sum relation above.
- Evaluates c(r) at the challenge r and sets `new_sum = c(r)`.
- `sum_x_fold_eq` extended to handle ternary rounds: at ternary levels, the factor becomes `L₀(r) + L₁(r)·x^third + L₂(r)·x^{2·third}` instead of the binary `(1-r) + r·x^half`.

**Proof size tradeoff:**
- Per round: ternary = 4 field elements, binary = 2.
- Rounds to reduce n to 1: ternary = log₃(n), binary = log₂(n). So ~37% fewer ternary rounds.
- Net for a size-9 residual: ternary sends 2 × 4 = 8 elements; binary sends 4 × 2 = 8 elements. Wash.
- For 3 × 2^a shaped residuals, ternary saves the final 2-3 binary rounds, which dominates at small residual sizes.

### Additional optimizations

Beyond the three main layers, this PR includes several performance and code quality improvements:

- **`fold_mixed` helper** — one function replaces the duplicated `if is_ternary_round { fold3 } else { fold }` pattern across 7 call sites.
- **Precomputed field constants** — `F::from(3/5/9/17)` and `inv3` hoisted once per sumcheck instead of recomputed per round.
- **Precomputed `half_inv` in `sum_x_fold_eq`** — one field inversion per verify call instead of per recursion level.
- **Debug assertion on `fold_coeff`** — catches misuse if someone passes a ternary-eligible size to the binary-only fold coefficient function.
- **`smooth_multilinear_extend`** — evaluates the MLE of an arbitrary smooth-size vector at a point using the mixed fold schedule. Used by `fold_based_mle_evaluate` for `Covector` linear forms. The implementation downcasts the `LinearForm` trait object via `std::any::Any` and, for Covectors, skips the zero-alloc + accumulate dance and folds the stored vector directly.
- **Removed `final_vector.clone()` in verifier** — the smooth verify path now moves the final vector into the fold instead of cloning it.
- **Removed dead code** — `fold_eq_weights` and `fold_based_dot_evaluate` from earlier iterations are superseded by `sum_x_fold_eq`.

## Files changed

| File | Lines | What |
|------|-------|------|
| `src/algebra/sumcheck.rs` | +451 | `fold3`, `compute_sumcheck_polynomial3`, fused variant, parallel ternary |
| `src/algebra/multilinear.rs` | +141 | `smooth_multilinear_extend` + property test |
| `src/algebra/mod.rs` | +4 | Export `smooth_multilinear_extend` |
| `src/protocols/sumcheck.rs` | +388 | Mixed ternary/binary prove/verify, `is_ternary_round`, `fold_mixed`, constant hoisting |
| `src/protocols/whir/config.rs` | +53 | Smooth size acceptance, `num_binary_variables`, `current_size` tracking, `ternary` flag |
| `src/protocols/whir/mod.rs` | +370 | `sum_x_fold_eq` (with ternary), `fold_coeff`, `fold_based_mle_evaluate`, `FinalClaim` metadata |
| `src/protocols/whir/prover.rs` | +7 | Thread ternary metadata |
| `src/protocols/whir/verifier.rs` | +59 | Mixed fold in `poly_eval`, `smooth_multilinear_extend` for final vector |
| `src/protocols/basecase.rs` | +29 | `fold_mixed` / `smooth_multilinear_extend` usage |
| `src/protocols/challenge_indices.rs` | +18 | Non-power-of-two bias correction |
| `src/smooth_domain.rs` | +41 | `is_smooth`, `odd_part`, `extra_rounds` |
| `src/lib.rs` | +1 | `pub mod smooth_domain` |

**Total: +1993 lines, -151 lines across 15 files.**

## Benchmarks

Parameters match what ProveKit uses in production: BN254 field, folding factor 3, 128-bit security, 10-bit proof-of-work, rate 1/4. Sizes from 2^17 to 2^24, including non-power-of-two sizes (3 × 2^a and 9 × 2^a).

Red = pow2 path (zero-pad to next power of two, current behavior). Blue = smooth + ternary (commit at exact size, this PR).

### Prove time

Blue tracks or beats red at every data point. The win comes from committing to a smaller vector — the NTT, Merkle tree, and sumcheck all scale with the commit size. At the pow2 stepping points (e.g., N = 9,437,184 pads up to 16,777,216), red pays for nearly 2× the real data.

| N | Red (pow2) | Blue (smooth) | Speedup |
|---|---|---|---|
| 786,432 | 215ms (pads to 1M) | 101ms | 2.1× |
| 1,572,864 | 423ms (pads to 2M) | 185ms | 2.3× |
| 6,291,456 | 2,469ms (pads to 8M) | 833ms | 3.0× |

At pure power-of-two sizes both lines coincide (same commit size = same work).

### Prove memory

Same staircase pattern: red uses the padded commit size, blue uses the actual smooth size. At N = 6.3M: red uses 3,136MB, blue uses 2,800MB.

### Verify time

Smooth verify is slower at non-pow2 sizes due to two structural costs:

1. **O(log² n) constraint evaluation.** `sum_x_fold_eq` does a `x^half` exponentiation per recursion level. Each exponentiation is O(log n); recursion depth is O(log n). The power-of-two path uses the tensor product identity in O(log n) total.

2. **O(n) MLE evaluation for Covectors.** `smooth_multilinear_extend` clones the covector vector and folds in-place. The power-of-two path uses `multilinear_extend`'s zero-allocation butterfly.

At N > 4M, smooth verify is ~1.5–2× slower than pow2. Below N = 1M the gap is negligible.

### Verify memory

Smooth uses more verify memory at non-pow2 sizes due to the covector clone. At pow2 sizes both paths are identical.

## Tests

121 tests pass. One new property-based test for `smooth_multilinear_extend` (verifies it agrees with fold-in-place for all smooth sizes 1..300 over Field64).

## Known limitations

1. **Verify cost at smooth sizes.** The O(log² n) tensor evaluation and O(n) MLE evaluation are structural — they come from the mixed fold schedule not admitting a clean tensor product factorization. Closing this fully is an open research question; we'd need a smooth-domain analog of the tensor identity.

2. **ProveKit integration.** ProveKit pins an older WHIR rev (`790bdf0`) and uses `whir_zk::Config` (the ZK wrapper). These changes are in the inner WHIR layer, so they apply transparently once ProveKit updates its dependency. ProveKit's custom `PrefixCovector` and `OffsetCovector` linear forms would either need `Any`-downcast support added or would use the fallback materialization path.

3. **Only final sumcheck uses ternary.** Initial and mid-round sumchecks remain binary-only because the NTT folding factor is pow2. Ternary folding could extend to intermediate rounds if the coset NTT supports smooth-size domains natively.
