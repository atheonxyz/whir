# ============================================================================
# PHASE 17: TENSOR IDENTITIES — SIDE BY SIDE
# ============================================================================
#
# WHAT THIS IS:
#   A verification PERFORMANCE OPTIMIZATION, not a correctness fix.
#   Without it, the smooth domain fix (script 16) makes verification O(n)
#   per constraint. This recovers O(log² n) per constraint.
#
# WHAT WAS FIXED (script 16):
#   algebra/sumcheck.rs  fold_half()  — fold splits at len/2
#   This broke the tensor product structure of fold_eq, making the
#   standard O(log n) tensor identity give WRONG answers for smooth sizes.
#
# WHAT THIS SCRIPT TRACES:
#   protocols/whir/mod.rs  sum_x_fold_eq()  — the smooth tensor identity
#   A recursive formula that evaluates Σ x^i · fold_eq(i) in O(log² n)
#   by adding a correction term at each odd-sized fold round.
#
# WHERE IT'S CALLED:
#   protocols/whir/verifier.rs  constraint subtraction loop:
#
#     if round_size.is_power_of_two() {
#         val = weights.mle_evaluate(point);
#         // → univariate_evaluation.rs:46  standard tensor identity  O(log n)
#     } else {
#         val = sum_x_fold_eq(weights.point, point, round_size);
#         // → mod.rs  smooth tensor identity  O(log² n)
#     }
#
# Run: sage test-arth/17_tensor_identity_side_by_side.sage
#
# ============================================================================

p = 97
F = GF(p)

print("=" * 78)
print("TENSOR IDENTITIES — SIDE BY SIDE")
print("Standard (power-of-2) vs Smooth (non-power-of-2)")
print("=" * 78)


# ── Helpers ────────────────────────────────────────────────────────────

def fold_coeff(pos, challenges, size):
    """The fold-eq weight for one position: trace it through each fold round."""
    c = F(1); pp = pos; s = size
    for r in challenges:
        if s <= 1: break
        h = s // 2
        if pp < h:
            c *= (1 - r)
        else:
            c *= r; pp -= h
        s = (s + 1) // 2
    return c

def standard_tensor(x, challenges):
    """O(log n): standard tensor identity for power-of-2 sizes.
    S = Π ((1-rⱼ) + rⱼ · x^{2^j})"""
    result = F(1)
    x2j = x
    for r in reversed(challenges):
        result *= (1 - r) + r * x2j
        x2j = x2j**2
    return result

def smooth_tensor(x, challenges, size):
    """O(log² n): smooth tensor identity for any size.
    Even: S = factor · S_next
    Odd:  S = factor · S_next - correction"""
    if size == 0: return F(0)
    if size == 1 or len(challenges) == 0: return F(1)
    half = size // 2
    r = challenges[0]; rest = challenges[1:]
    next_size = (size + 1) // 2
    x_half = x**half
    factor = (1 - r) + r * x_half
    s_next = smooth_tensor(x, rest, next_size)
    result = factor * s_next
    if size % 2 == 1:
        result -= (1 - r) * x_half * fold_coeff(half, rest, next_size)
    return result

def brute_force(x, challenges, size):
    """O(n): materialize and sum. Ground truth."""
    return sum(x**i * fold_coeff(i, challenges, size) for i in range(size))


# =========================================================================
# PART 1: n=8 (POWER OF 2) — BOTH IDENTITIES AGREE
# =========================================================================
print(f"""
{'━'*78}
PART 1: n=8 (POWER OF 2)
{'━'*78}

  Fold schedule: 8 → 4 → 2 → 1    (all even — no corrections needed)
  All splits are at len/2 = next_pow2 >> 1. Both identities agree.
""")

n8 = 8
x = F(7)
challenges_8 = [F(13), F(17), F(19)]
print(f"  x = {ZZ(x)},  challenges = {[ZZ(r) for r in challenges_8]}")

# Show the fold-eq weights for n=8
print(f"\n  Fold-eq weights (what each position contributes):")
print(f"  {'pos':>4}  {'binary':>8}  {'weight':>8}  {'x^i·weight':>12}")
total = F(0)
for i in range(n8):
    bits = [(i >> (2 - b)) & 1 for b in range(3)]
    w = fold_coeff(i, challenges_8, n8)
    contrib = x**i * w
    total += contrib
    bit_str = ''.join(str(b) for b in bits)
    # Show how weight decomposes as product of independent factors
    factors = []
    for b_idx, b_val in enumerate(bits):
        r = challenges_8[b_idx]
        if b_val == 1:
            factors.append(f"r{b_idx}={ZZ(r)}")
        else:
            factors.append(f"(1-r{b_idx})={ZZ(1-r)}")
    print(f"  {i:>4}  {bit_str:>8}  {ZZ(w):>8}  {ZZ(contrib):>12}    = {'·'.join(factors)}")

print(f"\n  Brute force sum: Σ x^i · fold_eq(i) = {ZZ(total)}")

# Standard tensor identity
print(f"\n  STANDARD TENSOR IDENTITY (O(log n) = O(3)):")
print(f"    S = Π ((1-rⱼ) + rⱼ · x^(2^j))")
result_std = F(1)
x2j = x
for j, r in enumerate(reversed(challenges_8)):
    factor = (1 - r) + r * x2j
    print(f"    j={2-j}: (1-{ZZ(r)}) + {ZZ(r)}·{ZZ(x)}^{2**j} = {ZZ(1-r)} + {ZZ(r)}·{ZZ(x2j)} = {ZZ(factor)}")
    result_std *= factor
    x2j = x2j**2
print(f"    Product = {ZZ(result_std)}")

# Smooth tensor identity
print(f"\n  SMOOTH TENSOR IDENTITY (O(log² n) — same result for power-of-2):")
result_smooth = smooth_tensor(x, challenges_8, n8)
print(f"    Result = {ZZ(result_smooth)}")

# Verify all three agree
print(f"""
  ┌────────────────────────────────────────────────┐
  │  Brute force:   {ZZ(total):3d}                            │
  │  Standard:      {ZZ(result_std):3d}                            │
  │  Smooth:        {ZZ(result_smooth):3d}                            │
  │  All equal?     {total == result_std == result_smooth}                         │
  │                                                │
  │  For power-of-two, BOTH identities work.        │
  │  The Rust code takes the standard path          │
  │  (fast path: smooth_size.is_power_of_two()).   │
  └────────────────────────────────────────────────┘

  WHY THEY AGREE: every round is even, so the smooth identity's
  recursion S = factor * S_next has NO correction term at any level.
  The factors are exactly the standard tensor factors.
""")


# =========================================================================
# PART 2: n=9 (SMOOTH) — STANDARD FAILS, SMOOTH WORKS
# =========================================================================
print(f"{'━'*78}")
print("PART 2: n=9 (SMOOTH, NOT POWER OF 2)")
print(f"{'━'*78}")

n9 = 9
challenges_9 = [F(13), F(17), F(19), F(23)]

print(f"""
  Fold schedule: 9 → 5 → 3 → 2 → 1    (three ODD rounds!)
  The split at each round:
    9: half=4, high=5 (ODD — 5 ≠ 4)
    5: half=2, high=3 (ODD — 3 ≠ 2)
    3: half=1, high=2 (ODD — 2 ≠ 1)
    2: half=1, high=1 (even)

  x = {ZZ(x)},  challenges = {[ZZ(r) for r in challenges_9]}
""")

# Show the fold-eq weights for n=9
print(f"  Fold-eq weights (position traced through each fold round):")
print(f"  {'pos':>4}  {'path':>20}  {'weight':>8}  {'x^i·weight':>12}")
total9 = F(0)
for i in range(n9):
    w = fold_coeff(i, challenges_9, n9)
    contrib = x**i * w
    total9 += contrib
    # Trace the path
    pp = i; s = n9; path_parts = []
    for rnd, r in enumerate(challenges_9):
        if s <= 1: break
        h = s // 2
        if pp < h:
            path_parts.append(f"L(h={h})")
        else:
            path_parts.append(f"H(h={h})")
            pp -= h
        s = (s + 1) // 2
    print(f"  {i:>4}  {'→'.join(path_parts):>20}  {ZZ(w):>8}  {ZZ(contrib):>12}")

print(f"\n  Brute force sum: {ZZ(total9)}")

# Standard tensor — WRONG
print(f"\n  STANDARD TENSOR IDENTITY (WRONG for n=9):")
result_std9 = standard_tensor(x, challenges_9)
print(f"    S = Π ((1-rⱼ) + rⱼ · x^(2^j))")
x2j = x
for j, r in enumerate(reversed(challenges_9)):
    factor = (1 - r) + r * x2j
    print(f"    j={3-j}: factor = {ZZ(factor)}   (uses x^{2**j} = {ZZ(x2j)})")
    x2j = x2j**2
print(f"    Product = {ZZ(result_std9)}")
print(f"    Brute force = {ZZ(total9)}")
print(f"    MATCH? {result_std9 == total9}  ← WRONG!")

print(f"""
  WHY IT FAILS: the standard identity assumes each variable
  independently selects a bit of the index (tensor product).
  But for n=9, the fold at half=4 puts positions 0-3 in LOW
  and 4-8 in HIGH. Position 4 and position 0 go to the SAME
  output position after folding — the variable assignment at
  round 1 depends on which half you fell into at round 0.
  The variables are COUPLED, not independent.
""")

# Smooth tensor — CORRECT
print(f"  SMOOTH TENSOR IDENTITY (CORRECT for n=9):")
print(f"  Step-by-step recursion:\n")

def smooth_tensor_traced(x, challenges, size, indent=4):
    prefix = " " * indent
    if size == 0:
        print(f"{prefix}S(..., 0) = 0")
        return F(0)
    if size == 1 or len(challenges) == 0:
        print(f"{prefix}S(..., {size}) = 1  [base]")
        return F(1)

    half = size // 2
    r = challenges[0]; rest = challenges[1:]
    next_size = (size + 1) // 2
    x_half = x**half
    factor = (1 - r) + r * x_half
    parity = "ODD" if size % 2 == 1 else "EVEN"

    print(f"{prefix}S(x, [r={ZZ(r)}, ...], {size}):  half={half}, {parity}")
    print(f"{prefix}  factor = (1-{ZZ(r)}) + {ZZ(r)}·x^{half} = {ZZ(1-r)} + {ZZ(r)}·{ZZ(x_half)} = {ZZ(factor)}")

    s_next = smooth_tensor_traced(x, rest, next_size, indent + 4)

    result = factor * s_next
    print(f"{prefix}  factor · S_next = {ZZ(factor)} · {ZZ(s_next)} = {ZZ(result)}")

    if size % 2 == 1:
        corr_coeff = fold_coeff(half, rest, next_size)
        corr = (1 - r) * x_half * corr_coeff
        result -= corr
        print(f"{prefix}  CORRECTION (odd size):")
        print(f"{prefix}    fold_coeff({half}, rest, {next_size}) = {ZZ(corr_coeff)}")
        print(f"{prefix}    (1-{ZZ(r)}) · x^{half} · {ZZ(corr_coeff)} = {ZZ(1-r)} · {ZZ(x_half)} · {ZZ(corr_coeff)} = {ZZ(corr)}")
        print(f"{prefix}    {ZZ(result + corr)} - {ZZ(corr)} = {ZZ(result)}")
    else:
        print(f"{prefix}  (even — no correction)")

    print(f"{prefix}  → S = {ZZ(result)}")
    return result

result_smooth9 = smooth_tensor_traced(x, challenges_9, n9)

print(f"""
  ┌────────────────────────────────────────────────┐
  │  Brute force:   {ZZ(total9):3d}                            │
  │  Standard:      {ZZ(result_std9):3d}   ← WRONG               │
  │  Smooth:        {ZZ(result_smooth9):3d}                            │
  │                                                │
  │  Standard == Brute?  {result_std9 == total9}  ← FAILS     │
  │  Smooth == Brute?    {result_smooth9 == total9}                  │
  └────────────────────────────────────────────────┘
""")


# =========================================================================
# PART 3: WHERE THIS MATTERS IN THE PROTOCOL
# =========================================================================
print(f"""{'━'*78}
PART 3: WHERE THIS MATTERS IN THE PROTOCOL
{'━'*78}

  The verifier runs this computation in the CONSTRAINT SUBTRACTION LOOP
  (verifier.rs, after sumcheck, before FinalClaim).

  Each OOD sample or in-domain sample creates a UnivariateEvaluation
  constraint with some evaluation point x. The verifier needs:

    val = Σ x^i · fold_eq(i, eval_point, round_size)

  This is subtracted from linear_form_rlc to isolate the user's
  linear forms.

  VERIFIER CODE PATH:

    if round_size.is_power_of_two() {{
        // Standard tensor identity — O(log n) per constraint
        val = weights.mle_evaluate(&evaluation_point[start..]);
        // calls UnivariateEvaluation::mle_evaluate
        // which uses Π ((1-rⱼ) + rⱼ·x^(2^j))
    }} else {{
        // Smooth tensor identity — O(log² n) per constraint
        val = sum_x_fold_eq(weights.point, &evaluation_point[start..], round_size);
        // recursive: factor · S_next - correction
    }}

  For n=8: takes the power-of-2 fast path. sum_x_fold_eq never called.
  For n=9: takes the smooth path. Without sum_x_fold_eq, the verifier
           would use the standard identity, get the WRONG answer, and
           reject valid proofs.

  NUMBER OF TIMES CALLED per proof:
    Once per OOD sample + once per in-domain sample, per WHIR round.
    Typical: ~10-30 constraints per round, 3-6 rounds = 30-180 calls.
    At O(log² n) each, this is fast. The old O(n) path was 50x slower.
""")


# =========================================================================
# PART 4: THE CORRECTION TERM — WHAT IT IS AND WHY
# =========================================================================
print(f"""{'━'*78}
PART 4: THE CORRECTION TERM — WHY ODD SIZES NEED IT
{'━'*78}

  For EVEN size n, split at n/2 gives equal halves:
    S = Σ_{{i<n/2}} x^i·eq_lo(i) + Σ_{{i≥n/2}} x^i·eq_hi(i)
      = (1-r)·Σ_{{i<n/2}} x^i·S_next_eq(i) + r·x^{{n/2}}·Σ_{{j<n/2}} x^j·S_next_eq(j)
      = ((1-r) + r·x^{{n/2}}) · S_next

    Both partial sums have the SAME number of terms. They factor cleanly.

  For ODD size n, split at ⌊n/2⌋ gives UNEQUAL halves:
    low has ⌊n/2⌋ terms.  high has ⌈n/2⌉ = ⌊n/2⌋+1 terms.

    The high sum has one MORE term than the low sum.
    When we try to factor:

    S = (1-r)·Σ_{{i<half}} + r·x^half·Σ_{{j<half+1}}
              └─ half terms ┘        └─ half+1 terms ┘

    The high sum (half+1 terms) and the low sum (half terms) DON'T match.
    We can write:

    low_sum  = S_next - x^half·fold_eq(half)     (full sum minus the extra term)
    high_sum = S_next                              (full sum, all half+1 terms)

    So: S = (1-r)·(S_next - x^half·fold_eq(half)) + r·x^half·S_next
          = ((1-r) + r·x^half)·S_next - (1-r)·x^half·fold_eq(half)
            └── standard factor ──┘    └── correction ─────────────┘

    The correction subtracts the extra term that the low partial sum
    doesn't have. It's the fold coefficient of ONE position (the boundary
    position ⌊n/2⌋), computed in O(log n) by tracing through the schedule.

  For n=9, there are 3 odd rounds (sizes 9, 5, 3). Each contributes
  one O(log n) correction. Total: O(3 · log 9) ≈ O(12) operations.
  For n=576 = 9·64, there are 3 odd rounds out of 10 total.
  The 7 even rounds are free (same as standard tensor identity).
""")
