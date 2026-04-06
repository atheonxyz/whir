# ============================================================================
# PHASE 16: THE BUG AND THE FIX — SIDE BY SIDE
# ============================================================================
#
# One script. Same witness, same randomness. Two fold strategies.
# Shows exactly where the old fold breaks and the new fold works.
#
# Left column:  OLD fold (next_power_of_two >> 1)  — the bug
# Right column: NEW fold (len / 2)                  — the fix
#
# Witness: n=12, d=4, k=3, ff=2
# Circuit: w = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
#
# Run: sage test-arth/16_bug_and_fix_side_by_side.sage
#
# ============================================================================

p = 97
F = GF(p)
g = F.multiplicative_generator()
R.<X> = PolynomialRing(F)
R2.<t> = PolynomialRing(F)

n = 12; d = 4; k = 3; ff = 2

w = [F(i+1) for i in range(n)]
blocks = [[w[j*k+i] for i in range(k)] for j in range(d)]
polys = [sum(F(blocks[j][i])*X**i for i in range(k)) for j in range(d)]

# Simple covector for clarity
covector = [F(1)]*n
S = sum(w[i]*covector[i] for i in range(n))

# Codeword (committed before any fold — same for both strategies)
cw_len = 6
w_cw = g**((p-1)//cw_len)
codeword = [[poly(w_cw**i) for poly in polys] for i in range(cw_len)]

# Fold challenges (same for both)
r1, r2 = F(7), F(11)

print("=" * 78)
print("THE BUG AND THE FIX — SIDE BY SIDE")
print("Same witness, same randomness, two fold strategies")
print("=" * 78)

# =========================================================================
# PART 1: SETUP — WHAT BOTH SIDES SHARE
# =========================================================================
print(f"""
{'━'*78}
PART 1: SETUP (identical for both strategies)
{'━'*78}

  Witness: w = {[ZZ(x) for x in w]}   (n={n})

  Interleaving: d={d} blocks of k={k} elements
    Block 0: w[0:3]  = {[ZZ(x) for x in blocks[0]]}  →  f₀(X) = {polys[0]}
    Block 1: w[3:6]  = {[ZZ(x) for x in blocks[1]]}  →  f₁(X) = {polys[1]}
    Block 2: w[6:9]  = {[ZZ(x) for x in blocks[2]]}  →  f₂(X) = {polys[2]}
    Block 3: w[9:12] = {[ZZ(x) for x in blocks[3]]}  →  f₃(X) = {polys[3]}

  Codeword (RS encode, committed via Merkle — SAME for both):
    row[0] = {[ZZ(v) for v in codeword[0]]}   at ω⁰ = {ZZ(w_cw**0)}
    row[1] = {[ZZ(v) for v in codeword[1]]}   at ω¹ = {ZZ(w_cw**1)}
    row[2] = {[ZZ(v) for v in codeword[2]]}   at ω² = {ZZ(w_cw**2)}
    row[3] = {[ZZ(v) for v in codeword[3]]}   at ω³ = {ZZ(w_cw**3)}
    row[4] = {[ZZ(v) for v in codeword[4]]}   at ω⁴ = {ZZ(w_cw**4)}
    row[5] = {[ZZ(v) for v in codeword[5]]}   at ω⁵ = {ZZ(w_cw**5)}

  Claim: dot(covector, w) = {ZZ(S)}
  Fold challenges: r₁={ZZ(r1)}, r₂={ZZ(r2)}
""")

# =========================================================================
# PART 2: FOLD ROUND 1 — WHERE THE STRATEGIES DIVERGE
# =========================================================================
print(f"{'━'*78}")
print("PART 2: FOLD ROUND 1 — THE CRITICAL DIFFERENCE")
print(f"{'━'*78}")

old_half = 8   # next_power_of_two(12) >> 1
new_half = 6   # 12 // 2

print(f"""
  ┌─────────────────────────────────┬─────────────────────────────────┐
  │  OLD: half = next_pow2(12)>>1   │  NEW: half = 12 / 2            │
  │       half = 16 >> 1 = 8        │       half = 6                  │
  ├─────────────────────────────────┼─────────────────────────────────┤
  │                                 │                                 │
  │  Block 0 [0-2]  ──┐            │  Block 0 [0-2]  ──┐            │
  │  Block 1 [3-5]  ──┤            │  Block 1 [3-5]  ──┤            │
  │  Block 2 [6-7]  ──┤ LOW (8)    │                    │ LOW (6)    │
  │  ─ ─ ─ ─ ─ ─ ─ ─ ┤            │  ── boundary ──────┤            │
  │  Block 2 [8]    ──┤ HIGH (4)   │                    │            │
  │  Block 3 [9-11] ──┘            │  Block 2 [6-8]  ──┤ HIGH (6)   │
  │                                 │  Block 3 [9-11] ──┘            │
  │  ✗ Block 2 is SPLIT!           │  ✓ Every block intact           │
  │    Low has 2/3 of Block 2      │    Blocks 0,1 in LOW            │
  │    High has 1/3 of Block 2     │    Blocks 2,3 in HIGH           │
  └─────────────────────────────────┴─────────────────────────────────┘
""")

# ── OLD fold round 1 ──
old_w = list(w)
old_cov = list(covector)
old_lo = old_w[:old_half]
old_hi = old_w[old_half:]  # only 4 elements
old_clo = old_cov[:old_half]
old_chi = old_cov[old_half:]

print(f"  OLD Round 1: split at {old_half}")
print(f"    w_lo[0:{old_half}] = {[ZZ(x) for x in old_lo]}")
print(f"    w_hi[{old_half}:{n}] = {[ZZ(x) for x in old_hi]}  (only {len(old_hi)} elements!)")
print(f"    Positions {len(old_hi)}..{old_half-1} in high are IMPLICIT ZEROS")

# Fold old
old_folded = []
for i in range(old_half):
    hi_val = F(old_hi[i]) if i < len(old_hi) else F(0)
    old_folded.append((1-r1)*F(old_lo[i]) + r1*hi_val)
print(f"    Folded: {[ZZ(x) for x in old_folded]}  (len={len(old_folded)})")

# ── NEW fold round 1 ──
new_w = list(w)
new_cov = list(covector)
new_lo = new_w[:new_half]
new_hi = new_w[new_half:]
new_clo = new_cov[:new_half]
new_chi = new_cov[new_half:]

print(f"\n  NEW Round 1: split at {new_half}")
print(f"    w_lo[0:{new_half}] = {[ZZ(x) for x in new_lo]}")
print(f"    w_hi[{new_half}:{n}] = {[ZZ(x) for x in new_hi]}  (exactly {len(new_hi)} elements)")
print(f"    Both halves are full — no implicit zeros")

new_folded = [(1-r1)*F(new_lo[i]) + r1*F(new_hi[i]) for i in range(new_half)]
print(f"    Folded: {[ZZ(x) for x in new_folded]}  (len={len(new_folded)})")

# ── Show what happened to each block ──
print(f"""
  WHAT HAPPENED TO BLOCK 2 (positions 6, 7, 8):

  OLD fold at 8:
    w[6] → low[6],  paired with high[6] (implicit zero)
      result[6] = (1-{ZZ(r1)})·{ZZ(w[6])} + {ZZ(r1)}·0 = {ZZ((1-r1)*w[6])}
    w[7] → low[7],  paired with high[7] (implicit zero)
      result[7] = (1-{ZZ(r1)})·{ZZ(w[7])} + {ZZ(r1)}·0 = {ZZ((1-r1)*w[7])}
    w[8] → high[0], paired with low[0]
      result[0] += {ZZ(r1)}·{ZZ(w[8])} (mixed into Block 0's position!)

    Block 2's data is scattered: two parts scaled by (1-r1) at positions 6,7,
    one part mixed into position 0. The block is destroyed.

  NEW fold at 6:
    w[6] → high[0], paired with low[0] = w[0]
      result[0] = (1-{ZZ(r1)})·{ZZ(w[0])} + {ZZ(r1)}·{ZZ(w[6])} = {ZZ(new_folded[0])}
    w[7] → high[1], paired with low[1] = w[1]
      result[1] = (1-{ZZ(r1)})·{ZZ(w[1])} + {ZZ(r1)}·{ZZ(w[7])} = {ZZ(new_folded[1])}
    w[8] → high[2], paired with low[2] = w[2]
      result[2] = (1-{ZZ(r1)})·{ZZ(w[2])} + {ZZ(r1)}·{ZZ(w[8])} = {ZZ(new_folded[2])}

    Block 2 pairs cleanly with Block 0. Each element goes to the correct
    coefficient position. The block structure is preserved.
""")

# ── Sumcheck check — both pass! ──
# OLD
old_cp = sum(((1-t)*F(old_clo[i])+t*(F(old_chi[i]) if i < len(old_chi) else F(0)))*
             ((1-t)*F(old_lo[i])+t*(F(old_hi[i]) if i < len(old_hi) else F(0)))
             for i in range(old_half))
# NEW
new_cp = sum(((1-t)*F(new_clo[i])+t*F(new_chi[i]))*
             ((1-t)*F(new_lo[i])+t*F(new_hi[i]))
             for i in range(new_half))

print(f"  SUMCHECK INVARIANT — BOTH PASS (that's the trap!):")
print(f"    OLD: c(0)+c(1) = {ZZ(old_cp(0))}+{ZZ(old_cp(1))} = {ZZ(old_cp(0)+old_cp(1))} = S = {ZZ(S)} ✓")
print(f"    NEW: c(0)+c(1) = {ZZ(new_cp(0))}+{ZZ(new_cp(1))} = {ZZ(new_cp(0)+new_cp(1))} = S = {ZZ(S)} ✓")
print(f"""
    The sumcheck invariant holds for ANY split point — it just partitions
    the dot product into two halves. This is why the prover never catches
    the bug. The bug only surfaces at the in-domain check.
""")

# =========================================================================
# PART 3: FOLD ROUND 2 — COMPLETE THE FOLD
# =========================================================================
print(f"{'━'*78}")
print("PART 3: FOLD ROUND 2")
print(f"{'━'*78}")

# OLD round 2 (8 elements → 4)
old_half2 = 4
old_lo2 = old_folded[:old_half2]
old_hi2 = old_folded[old_half2:]
old_folded2 = [(1-r2)*F(old_lo2[i]) + r2*F(old_hi2[i]) for i in range(old_half2)]

# NEW round 2 (6 elements → 3)
new_half2 = 3
new_lo2 = new_folded[:new_half2]
new_hi2 = new_folded[new_half2:]
new_folded2 = [(1-r2)*F(new_lo2[i]) + r2*F(new_hi2[i]) for i in range(new_half2)]

print(f"""
  ┌─────────────────────────────────┬─────────────────────────────────┐
  │  OLD: {len(old_folded)} → {old_half2} (split at {old_half2})       │  NEW: {len(new_folded)} → {new_half2} (split at {new_half2})        │
  ├─────────────────────────────────┼─────────────────────────────────┤
  │  input:  {str([ZZ(x) for x in old_folded[:4]]):24s}│  input:  {str([ZZ(x) for x in new_folded[:3]]):24s}│
  │          {str([ZZ(x) for x in old_folded[4:]]):24s}│          {str([ZZ(x) for x in new_folded[3:]]):24s}│
  │  output: {str([ZZ(x) for x in old_folded2]):24s}│  output: {str([ZZ(x) for x in new_folded2]):24s}│
  │  length: {len(old_folded2)} elements                │  length: {len(new_folded2)} elements                │
  │  needed: k = {k}                      │  needed: k = {k}                      │
  │  {len(old_folded2)} ≠ {k}  WRONG SIZE              │  {len(new_folded2)} = {k}  CORRECT ✓                │
  └─────────────────────────────────┴─────────────────────────────────┘
""")

# =========================================================================
# PART 4: THE IN-DOMAIN CHECK — WHERE THE BUG KILLS THE PROOF
# =========================================================================
print(f"{'━'*78}")
print("PART 4: THE IN-DOMAIN CHECK — OLD FAILS, NEW PASSES")
print(f"{'━'*78}")

# eq weights (same for both — they depend on randomness, not fold strategy)
eq_w = [(1-r1)*(1-r2), (1-r1)*r2, r1*(1-r2), r1*r2]

# The correct algebraic folded polynomial
f_folded = sum(eq_w[j]*polys[j] for j in range(d))
correct_coeffs = [f_folded[i] for i in range(k)]

print(f"""
  eq_weights = {[ZZ(x) for x in eq_w]}
  (from challenges r₁={ZZ(r1)}, r₂={ZZ(r2)})

  Correct folded polynomial (what the codeword represents):
    f_folded = Σ eq_w[j]·fⱼ(X) = {f_folded}
    Coefficients: {[ZZ(c) for c in correct_coeffs]}

  OLD folded vector: {[ZZ(x) for x in old_folded2]}  (length {len(old_folded2)})
  NEW folded vector: {[ZZ(x) for x in new_folded2]}  (length {len(new_folded2)})

  Do they match the correct coefficients?
""")

old_match = len(old_folded2) == k and all(F(old_folded2[i]) == correct_coeffs[i] for i in range(k))
new_match = len(new_folded2) == k and all(F(new_folded2[i]) == correct_coeffs[i] for i in range(k))

print(f"    OLD: {[ZZ(x) for x in old_folded2]} vs {[ZZ(c) for c in correct_coeffs]} → ", end="")
if len(old_folded2) != k:
    print(f"WRONG LENGTH ({len(old_folded2)} ≠ {k})")
elif not old_match:
    print("WRONG VALUES")
else:
    print("MATCH ✓")

print(f"    NEW: {[ZZ(x) for x in new_folded2]} vs {[ZZ(c) for c in correct_coeffs]} → ", end="")
if new_match:
    print("MATCH ✓")
else:
    print("MISMATCH")

# ── Check at every codeword position ──
print(f"""
  In-domain check at ALL {cw_len} codeword positions:

  ┌─────┬─────┬────────────────────────────┬────────────────────────────┐
  │  ωⁱ │  pt │  OLD: LEFT vs RIGHT        │  NEW: LEFT vs RIGHT        │
  ├─────┼─────┼────────────────────────────┼────────────────────────────┤""")

old_pass_count = 0
new_pass_count = 0
for idx in range(cw_len):
    pt = w_cw**idx
    left = sum(eq_w[j]*codeword[idx][j] for j in range(d))
    old_right = sum(F(old_folded2[i])*pt**i for i in range(len(old_folded2)))
    new_right = sum(F(new_folded2[i])*pt**i for i in range(len(new_folded2)))
    old_ok = left == old_right
    new_ok = left == new_right
    if old_ok: old_pass_count += 1
    if new_ok: new_pass_count += 1
    old_status = "✓ PASS" if old_ok else "✗ FAIL"
    new_status = "✓ PASS" if new_ok else "✗ FAIL"
    print(f"  │  ω^{idx} │ {ZZ(pt):3d} │  {ZZ(left):2d} vs {ZZ(old_right):2d}  {old_status:8s}      │  {ZZ(left):2d} vs {ZZ(new_right):2d}  {new_status:8s}      │")

print(f"  └─────┴─────┴────────────────────────────┴────────────────────────────┘")
print(f"         OLD: {old_pass_count}/{cw_len} pass                    NEW: {new_pass_count}/{cw_len} pass")

# =========================================================================
# PART 5: WHY THE OLD FOLD CORRUPTS ALL COEFFICIENTS
# =========================================================================
print(f"""
{'━'*78}
PART 5: WHY THE OLD FOLD CORRUPTS ALL COEFFICIENTS (not just one extra)
{'━'*78}

  The fold at position 8 doesn't just produce an extra coefficient.
  It corrupts ALL of them. Here's why:

  The correct combination is:
    f_folded = eq_w[0]·f₀ + eq_w[1]·f₁ + eq_w[2]·f₂ + eq_w[3]·f₃

  This requires Block j to be fully in one fold-half so that eq_w[j]
  applies uniformly to all k coefficients of that block.

  With the old fold (boundary at 8):
    Block 0 [0-2]:  fully in LOW  → gets factor (1-r₁)  ✓
    Block 1 [3-5]:  fully in LOW  → gets factor (1-r₁)  ✓
    Block 2 [6-8]:  SPLIT!
      positions 6,7: in LOW  → get factor (1-r₁)
      position 8:    in HIGH → gets factor r₁
    Block 3 [9-11]: fully in HIGH → gets factor r₁      ✓

  Block 2's coefficient at position i in the folded vector:
    coeff[i] should be eq_w[2] · block2[i]
    but actually:
      for i=0 (position 6): (1-r₁)·w[6]  — Block 2 data with WRONG weight
      for i=1 (position 7): (1-r₁)·w[7]  — Block 2 data with WRONG weight
      for i=2 (position 8): mixed into position 0 of the high half
                             — Block 2 data at WRONG position

  This means Block 2's contribution to the folded polynomial is WRONG
  at every coefficient, not just the extra one. And since the fold
  combines blocks additively, the corruption propagates to ALL output
  coefficients.

  OLD output = {[ZZ(x) for x in old_folded2]}
  Correct    = {[ZZ(c) for c in correct_coeffs]}
  Difference = {[ZZ(F(old_folded2[i]) - correct_coeffs[i]) for i in range(min(len(old_folded2), k))]}
  Every coefficient is wrong.
""")

# =========================================================================
# PART 6: WHY THE NEW FOLD WORKS
# =========================================================================
print(f"""{'━'*78}
PART 6: WHY THE NEW FOLD PRESERVES BLOCK STRUCTURE
{'━'*78}

  With the new fold (boundary at 6):

  Round 1 (split at 6):
    LOW  = [w[0], w[1], w[2], w[3], w[4], w[5]]  = Block 0 + Block 1
    HIGH = [w[6], w[7], w[8], w[9], w[10], w[11]] = Block 2 + Block 3

    fold(i) = (1-r₁)·low[i] + r₁·high[i]

    Position 0: (1-r₁)·w[0] + r₁·w[6]   = (1-r₁)·block0[0] + r₁·block2[0]
    Position 1: (1-r₁)·w[1] + r₁·w[7]   = (1-r₁)·block0[1] + r₁·block2[1]
    Position 2: (1-r₁)·w[2] + r₁·w[8]   = (1-r₁)·block0[2] + r₁·block2[2]
    Position 3: (1-r₁)·w[3] + r₁·w[9]   = (1-r₁)·block1[0] + r₁·block3[0]
    Position 4: (1-r₁)·w[4] + r₁·w[10]  = (1-r₁)·block1[1] + r₁·block3[1]
    Position 5: (1-r₁)·w[5] + r₁·w[11]  = (1-r₁)·block1[2] + r₁·block3[2]

    Block 0 pairs with Block 2.  Block 1 pairs with Block 3.
    Within each pair, coefficient i pairs with coefficient i.

  Round 2 (split at 3):
    LOW  = [result[0], result[1], result[2]]   = Block(0+2) folded
    HIGH = [result[3], result[4], result[5]]   = Block(1+3) folded

    Final[i] = (1-r₂)·[(1-r₁)·block0[i]+r₁·block2[i]]
             + r₂·[(1-r₁)·block1[i]+r₁·block3[i]]
             = (1-r₁)(1-r₂)·block0[i] + (1-r₁)r₂·block1[i]
             + r₁(1-r₂)·block2[i]     + r₁r₂·block3[i]
             = eq_w[0]·block0[i] + eq_w[1]·block1[i]
             + eq_w[2]·block2[i] + eq_w[3]·block3[i]

    This IS the correct algebraic combination f_folded(X) evaluated
    at coefficient i. The fold output equals the polynomial coefficients.""")

for i in range(k):
    terms = " + ".join(f"{ZZ(eq_w[j])}·{ZZ(blocks[j][i])}" for j in range(d))
    print(f"    coeff[{i}] = {terms} = {ZZ(correct_coeffs[i])} = fold[{i}] = {ZZ(new_folded2[i])} ✓")

# =========================================================================
# PART 7: IMPACT ON SUBSEQUENT STAGES
# =========================================================================
print(f"""
{'━'*78}
PART 7: WHAT ELSE CHANGED (the cascade)
{'━'*78}

  The fold fix (one line: len/2 instead of next_pow2>>1) breaks an
  assumption: that fold() and multilinear_extend() use the same splits.

  For power-of-2: len/2 = next_pow2>>1. Nothing changes. Backward compatible.

  For smooth sizes: fold splits at len/2, multilinear_extend splits at
  the power-of-2 boundary. They produce different variable orderings.
  This means fold_eq is NOT a tensor product — variables are coupled.

  ┌──────────────────────────────────────────────────────────────────────┐
  │  EVERY PLACE THE VERIFIER CALLED multilinear_extend ON A FOLD      │
  │  RESULT MUST NOW USE fold-based EVALUATION INSTEAD.                │
  │                                                                      │
  │  1. verifier.rs: poly_eval for non-pow2 final_vector               │
  │     → fold-replay instead of MultilinearExtension                   │
  │                                                                      │
  │  2. verifier.rs: constraint subtraction loop                        │
  │     → sum_x_fold_eq() O(log²n) smooth tensor identity              │
  │     (script 15 derives this)                                        │
  │                                                                      │
  │  3. mod.rs: FinalClaim::verify                                      │
  │     → fold_based_mle_evaluate for non-pow2 linear forms             │
  │                                                                      │
  │  4. config.rs: final_size computation                               │
  │     → iterative ceil(n/2) instead of next_pow2 >> rounds            │
  │                                                                      │
  │  5. challenge_indices.rs: remove power-of-2 assertion               │
  │     → smooth codeword lengths are valid num_leaves                  │
  └──────────────────────────────────────────────────────────────────────┘

  And for odd-length folds (k=9: 9→5→3→2→1):
    Each fold of an odd-length vector has unequal halves.
    The extra high element gets result[j] = r·high[j].
    The sumcheck invariant still holds: the extra element contributes
    only t²·a·b to the polynomial (zero at t=0 and t=1).

    See script 14 for full trace, script 15 for the O(log²n) identity.
""")

# =========================================================================
# PART 8: FINAL SUMMARY
# =========================================================================
print(f"""{'━'*78}
SUMMARY
{'━'*78}

  ┌───────────────────────┬───────────────────────┬───────────────────────┐
  │                       │ OLD (next_pow2 >> 1)  │ NEW (len / 2)         │
  ├───────────────────────┼───────────────────────┼───────────────────────┤
  │ Split point (n=12)    │ 8                     │ 6                     │
  │ Block 2 intact?       │ NO — split at pos 8   │ YES — fully in HIGH   │
  │ Folded length         │ 4 (wrong)             │ 3 = k (correct)       │
  │ Coefficients correct? │ ALL wrong             │ ALL correct           │
  │ Sumcheck c(0)+c(1)=S? │ YES (always)          │ YES (always)          │
  │ In-domain check       │ FAILS at all points   │ PASSES at all points  │
  │ Power-of-2 behavior   │ (baseline)            │ IDENTICAL             │
  │ Backward compatible?  │ —                     │ YES                   │
  └───────────────────────┴───────────────────────┴───────────────────────┘

  The bug: fold preserves the DOT PRODUCT (sumcheck passes) but
  corrupts the POLYNOMIAL STRUCTURE (in-domain check fails).

  The prover never notices because it only needs the dot product.
  The verifier catches it because it checks the polynomial structure
  against the committed codeword.

  Script 13 diagnosed it. Script 14 traces the fix. Script 15 recovers
  fast verification via the O(log²n) smooth tensor identity.
""")
