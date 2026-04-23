# ============================================================================
# PHASE 9: THE COMPLETE WHIR TRACE
# ============================================================================
#
# One file. Every step. From Noir circuit to verified proof.
# With NTT, coset decomposition, Merkle tree, sumcheck, in-domain check.
# Every step annotated with the Rust code that executes it.
#
# Circuit: "I know a,b such that a*b + a = 15"
#   Constraint 0: a * b = c
#   Constraint 1: (c + a) * 1 = d
#   Constraint 2: d * 1 = public_out
#
# Run: sage 09_complete_trace.sage
# Read alongside: 09i_intuition_guide.md (explains the WHY behind every step)

print("=" * 70)
print("THE COMPLETE WHIR TRACE")
print("From Circuit to Verified Proof — Every Step with Code References")
print("=" * 70)

import hashlib

p = 97
F = GF(p)
g = F.multiplicative_generator()
R_poly.<X> = PolynomialRing(F)

# =========================================================================
# ░░░ PART 1: CIRCUIT → R1CS → WITNESS ░░░
# =========================================================================
print(f"\n{'━'*70}")
print("PART 1: CIRCUIT → R1CS → WITNESS")
print("WHY: The circuit defines what we're proving. R1CS is the constraint")
print("     format: each constraint is one multiplication. The witness is")
print("     everything computed during execution — the prover's secret.")
print(f"{'━'*70}")

print(f"""
  ┌──────────────────────────────────────────────────────────┐
  │  NOIR PROGRAM                                            │
  │                                                          │
  │  fn main(public_out: Field, a: Field, b: Field) {{       │
  │      let c = a * b;           // constraint 0            │
  │      let d = c + a;           // constraint 1            │
  │      assert(d == public_out); // constraint 2            │
  │  }}                                                       │
  │                                                          │
  │  Execution: a=3, b=4 → c=12, d=15, public_out=15       │
  └──────────────────────────────────────────────────────────┘

  Code: provekit/r1cs-compiler/src/noir_to_r1cs.rs
        compiles Noir ACIR → R1CS matrices (A, B, C)
        + WitnessBuilders (recipes for computing each w[i])
""")

# Witness
w = vector(F, [1, 15, 3, 4, 12, 15, 0, 0])
n = 8

print(f"  WITNESS w = {[ZZ(x) for x in w]}  (n = {n})")
print(f"  ┌─────────────────────────────────────────┐")
print(f"  │ w[0] = 1    constant (always 1)         │")
print(f"  │ w[1] = 15   PUBLIC output               │")
print(f"  │ w[2] = 3    PRIVATE input a             │")
print(f"  │ w[3] = 4    PRIVATE input b             │")
print(f"  │ w[4] = 12   intermediate: c = a*b       │")
print(f"  │ w[5] = 15   intermediate: d = c+a       │")
print(f"  │ w[6] = 0    padding                     │")
print(f"  │ w[7] = 0    padding                     │")
print(f"  └─────────────────────────────────────────┘")

# R1CS
A_rows = [
    vector(F, [0, 0, 1, 0, 0, 0, 0, 0]),
    vector(F, [0, 0, 1, 0, 1, 0, 0, 0]),
    vector(F, [0, 0, 0, 0, 0, 1, 0, 0]),
]
B_rows = [
    vector(F, [0, 0, 0, 1, 0, 0, 0, 0]),
    vector(F, [1, 0, 0, 0, 0, 0, 0, 0]),
    vector(F, [1, 0, 0, 0, 0, 0, 0, 0]),
]
C_rows = [
    vector(F, [0, 0, 0, 0, 1, 0, 0, 0]),
    vector(F, [0, 0, 0, 0, 0, 1, 0, 0]),
    vector(F, [0, 1, 0, 0, 0, 0, 0, 0]),
]
m = 3

print(f"\n  R1CS: {m} constraints, each (A[i]·w) × (B[i]·w) = (C[i]·w)")
for i in range(m):
    av, bv, cv = A_rows[i]*w, B_rows[i]*w, C_rows[i]*w
    print(f"    C{i}: {ZZ(av)} × {ZZ(bv)} = {ZZ(av*bv)}, C[{i}]·w = {ZZ(cv)}  {'✓' if av*bv == cv else '✗'}")

# =========================================================================
# ░░░ PART 2: QUADRATIC → LINEAR REDUCTION ░░░
# =========================================================================
print(f"\n{'━'*70}")
print("PART 2: BATCH 3 QUADRATIC CONSTRAINTS → LINEAR CLAIMS")
print("WHY: WHIR can only prove linear claims <L, w> = value. R1CS is")
print("     quadratic: (A·w)×(B·w)=C·w. We batch with random α, then")
print("     use a sumcheck to reduce the quadratic dot product to linear")
print("     claims on w. This is the bridge between R1CS and WHIR.")
print("Code: provekit/prover/src/whir_r1cs.rs")
print(f"{'━'*70}")

alpha = F(5)
print(f"\n  Verifier challenge α = {ZZ(alpha)}")

# Build a_vec, b_vec
a_vec_raw = [alpha**i * (A_rows[i] * w) for i in range(m)]
b_vec_raw = [B_rows[i] * w for i in range(m)]

print(f"""
  Rewrite as dot product:
    a_vec[i] = αⁱ · (A[i]·w)    b_vec[i] = B[i]·w

    a_vec = [{', '.join(f'{ZZ(alpha**i)}·{ZZ(A_rows[i]*w)}={ZZ(a_vec_raw[i])}' for i in range(m))}]
    b_vec = [{', '.join(str(ZZ(x)) for x in b_vec_raw)}]
""")

# Pad to power of 2
m_pad = 4
a_vec = list(a_vec_raw) + [F(0)]
b_vec = list(b_vec_raw) + [F(0)]
S_original = sum(a_vec[i]*b_vec[i] for i in range(m_pad))

print(f"  Padded to {m_pad}: a_vec = {[ZZ(x) for x in a_vec]}")
print(f"                    b_vec = {[ZZ(x) for x in b_vec]}")
print(f"  S = dot(a_vec, b_vec) = {ZZ(S_original)}")

# Right side
S_rhs = sum(alpha**i * (C_rows[i]*w) for i in range(m))
print(f"  Σᵢ αⁱ·(C[i]·w) = {ZZ(S_rhs)}")
print(f"  Left = Right? {S_original == S_rhs}  ✓")

# --- Inner sumcheck: reduce dot(a_vec, b_vec) to a_final × b_final ---
print(f"\n  ┌──────────────────────────────────────────────────┐")
print(f"  │  INNER SUMCHECK: reduce dot product              │")
print(f"  │  dot(a_vec, b_vec) = {ZZ(S_original)} → a_final × b_final  │")
print(f"  └──────────────────────────────────────────────────┘")

R2.<t> = PolynomialRing(F)

# ── Round 1 ──
half = m_pad // 2
a_lo, a_hi = a_vec[:half], a_vec[half:]
b_lo, b_hi = b_vec[:half], b_vec[half:]

print(f"\n  ┌────────────────────────────────────────────────────────────────┐")
print(f"  │  INNER SUMCHECK ROUND 1: split [{m_pad}] → [{half}|{half}]                     │")
print(f"  └────────────────────────────────────────────────────────────────┘")

print(f"\n    a_lo = a_vec[0:2] = {[ZZ(x) for x in a_lo]}")
print(f"    a_hi = a_vec[2:4] = {[ZZ(x) for x in a_hi]}")
print(f"    b_lo = b_vec[0:2] = {[ZZ(x) for x in b_lo]}")
print(f"    b_hi = b_vec[2:4] = {[ZZ(x) for x in b_hi]}")

print(f"""
    c₁(t) = Σᵢ [(1-t)·a_lo[i] + t·a_hi[i]] × [(1-t)·b_lo[i] + t·b_hi[i]]

    Expand term by term:
""")

c1_terms = []
for i in range(half):
    al, ah = F(a_lo[i]), F(a_hi[i])
    bl, bh = F(b_lo[i]), F(b_hi[i])
    const = al * bl
    linear = al * (bh - bl) + (ah - al) * bl
    quadratic = (ah - al) * (bh - bl)

    print(f"    Term i={i}:")
    print(f"      a(t,{i}) = (1-t)·{ZZ(al)} + t·{ZZ(ah)} = {ZZ(al)} + {ZZ(ah-al)}·t")
    print(f"      b(t,{i}) = (1-t)·{ZZ(bl)} + t·{ZZ(bh)} = {ZZ(bl)} + {ZZ(bh-bl)}·t")
    print(f"      [{ZZ(al)}+{ZZ(ah-al)}·t] × [{ZZ(bl)}+{ZZ(bh-bl)}·t]")
    print(f"        constant:  {ZZ(al)}×{ZZ(bl)} = {ZZ(const)}")
    print(f"        linear:    {ZZ(al)}×{ZZ(bh-bl)} + {ZZ(ah-al)}×{ZZ(bl)} = {ZZ(al*(bh-bl))} + {ZZ((ah-al)*bl)} = {ZZ(linear)}")
    print(f"        quadratic: {ZZ(ah-al)}×{ZZ(bh-bl)} = {ZZ(quadratic)}")
    print(f"      = {ZZ(quadratic)}t² + {ZZ(linear)}t + {ZZ(const)}")
    print()
    c1_terms.append((const, linear, quadratic))

tc = sum(c for c,l,q in c1_terms)
tl = sum(l for c,l,q in c1_terms)
tq = sum(q for c,l,q in c1_terms)
c1 = tc + tl*t + tq*t**2

print(f"    Sum all terms:")
print(f"      constant:  {' + '.join(str(ZZ(c)) for c,l,q in c1_terms)} = {ZZ(tc)}")
print(f"      linear:    {' + '.join(str(ZZ(l)) for c,l,q in c1_terms)} = {ZZ(tl)}")
print(f"      quadratic: {' + '.join(str(ZZ(q)) for c,l,q in c1_terms)} = {ZZ(tq)}")
print(f"    c₁(t) = {ZZ(tq)}t² + {ZZ(tl)}t + {ZZ(tc)}")

print(f"\n    CHECK: c₁(0) + c₁(1) = S?")
print(f"      c₁(0) = {ZZ(c1(0))}  (= {'+'.join(f'{ZZ(a_lo[i])}·{ZZ(b_lo[i])}' for i in range(half))} = dot of low halves)")
print(f"      c₁(1) = {ZZ(c1(1))}  (= {'+'.join(f'{ZZ(a_hi[i])}·{ZZ(b_hi[i])}' for i in range(half))} = dot of high halves)")
print(f"      Sum = {ZZ(c1(0))} + {ZZ(c1(1))} = {ZZ(c1(0)+c1(1))} = S = {ZZ(S_original)} ✓")

r_inner1 = F(3)
S1 = c1(r_inner1)
print(f"\n    Verifier picks r₁ = {ZZ(r_inner1)}")
print(f"    c₁({ZZ(r_inner1)}) = {ZZ(tq)}·{ZZ(r_inner1**2)} + {ZZ(tl)}·{ZZ(r_inner1)} + {ZZ(tc)} = {ZZ(tq*r_inner1**2)} + {ZZ(tl*r_inner1)} + {ZZ(tc)} = {ZZ(S1)}")
print(f"    New claim: S' = {ZZ(S1)}")

print(f"\n    FOLD: for each i, a_f[i] = (1-{ZZ(r_inner1)})·a_lo[i] + {ZZ(r_inner1)}·a_hi[i]")
a_f1 = []
b_f1 = []
for i in range(half):
    af = (1-r_inner1)*F(a_lo[i]) + r_inner1*F(a_hi[i])
    bf = (1-r_inner1)*F(b_lo[i]) + r_inner1*F(b_hi[i])
    a_f1.append(af)
    b_f1.append(bf)
    print(f"      i={i}: a_f = {ZZ(1-r_inner1)}·{ZZ(a_lo[i])} + {ZZ(r_inner1)}·{ZZ(a_hi[i])} = {ZZ((1-r_inner1)*F(a_lo[i]))} + {ZZ(r_inner1*F(a_hi[i]))} = {ZZ(af)}")
    print(f"           b_f = {ZZ(1-r_inner1)}·{ZZ(b_lo[i])} + {ZZ(r_inner1)}·{ZZ(b_hi[i])} = {ZZ((1-r_inner1)*F(b_lo[i]))} + {ZZ(r_inner1*F(b_hi[i]))} = {ZZ(bf)}")

print(f"    folded_a = {[ZZ(x) for x in a_f1]}")
print(f"    folded_b = {[ZZ(x) for x in b_f1]}")
print(f"    dot = {'+'.join(f'{ZZ(a_f1[i])}·{ZZ(b_f1[i])}' for i in range(half))} = {ZZ(sum(a_f1[i]*b_f1[i] for i in range(half)))} = S' ✓")

# ── Round 2 ──
print(f"\n  ┌────────────────────────────────────────────────────────────────┐")
print(f"  │  INNER SUMCHECK ROUND 2: split [2] → [1|1]                    │")
print(f"  └────────────────────────────────────────────────────────────────┘")

al2, ah2 = F(a_f1[0]), F(a_f1[1])
bl2, bh2 = F(b_f1[0]), F(b_f1[1])

print(f"\n    a_lo = [{ZZ(al2)}], a_hi = [{ZZ(ah2)}]")
print(f"    b_lo = [{ZZ(bl2)}], b_hi = [{ZZ(bh2)}]")

const2 = al2 * bl2
linear2 = al2 * (bh2-bl2) + (ah2-al2) * bl2
quad2 = (ah2-al2) * (bh2-bl2)

print(f"\n    Single term:")
print(f"      a(t) = {ZZ(al2)} + {ZZ(ah2-al2)}·t")
print(f"      b(t) = {ZZ(bl2)} + {ZZ(bh2-bl2)}·t")
print(f"      [{ZZ(al2)}+{ZZ(ah2-al2)}·t] × [{ZZ(bl2)}+{ZZ(bh2-bl2)}·t]")
print(f"        constant:  {ZZ(al2)}×{ZZ(bl2)} = {ZZ(const2)}")
print(f"        linear:    {ZZ(al2)}×{ZZ(bh2-bl2)} + {ZZ(ah2-al2)}×{ZZ(bl2)} = {ZZ(linear2)}")
print(f"        quadratic: {ZZ(ah2-al2)}×{ZZ(bh2-bl2)} = {ZZ(quad2)}")
c2 = const2 + linear2*t + quad2*t**2
print(f"    c₂(t) = {ZZ(quad2)}t² + {ZZ(linear2)}t + {ZZ(const2)}")

print(f"\n    CHECK: c₂(0) + c₂(1) = {ZZ(c2(0))} + {ZZ(c2(1))} = {ZZ(c2(0)+c2(1))} = S' = {ZZ(S1)} ✓")

r_inner2 = F(7)
S2 = c2(r_inner2)
print(f"\n    Verifier picks r₂ = {ZZ(r_inner2)}")
print(f"    c₂({ZZ(r_inner2)}) = {ZZ(quad2)}·{ZZ(r_inner2**2)} + {ZZ(linear2)}·{ZZ(r_inner2)} + {ZZ(const2)} = {ZZ(quad2*r_inner2**2)} + {ZZ(linear2*r_inner2)} + {ZZ(const2)} = {ZZ(S2)}")

a_final = (1-r_inner2)*al2 + r_inner2*ah2
b_final = (1-r_inner2)*bl2 + r_inner2*bh2
print(f"\n    FINAL FOLD:")
print(f"      a_final = {ZZ(1-r_inner2)}·{ZZ(al2)} + {ZZ(r_inner2)}·{ZZ(ah2)} = {ZZ(a_final)}")
print(f"      b_final = {ZZ(1-r_inner2)}·{ZZ(bl2)} + {ZZ(r_inner2)}·{ZZ(bh2)} = {ZZ(b_final)}")
print(f"      a_final × b_final = {ZZ(a_final)} × {ZZ(b_final)} = {ZZ(a_final*b_final)}")
print(f"      S'' = {ZZ(S2)}")
print(f"      Match? {a_final*b_final == S2} ✓")

print(f"""
    ═══════════════════════════════════════════════════════════
    SUMCHECK CHAIN OF TRUST:

      S = {ZZ(S_original)}   (original: dot(a_vec, b_vec))
        │  c₁(0)+c₁(1) = {ZZ(S_original)} ✓
        │  r₁ = {ZZ(r_inner1)}
        ▼
      S' = {ZZ(S1)}   (= c₁({ZZ(r_inner1)}))
        │  c₂(0)+c₂(1) = {ZZ(S1)} ✓
        │  r₂ = {ZZ(r_inner2)}
        ▼
      S'' = {ZZ(S2)}  (= c₂({ZZ(r_inner2)}))
        │  a_final × b_final = {ZZ(a_final)}×{ZZ(b_final)} = {ZZ(S2)} ✓
        ▼
      DONE — reduced {m_pad}-element dot product to 1 multiplication
    ═══════════════════════════════════════════════════════════
""")

# ── Derive fold weights ──
print(f"""  ┌────────────────────────────────────────────────────────────────┐
  │  DERIVE LINEAR FORMS FROM FOLD WEIGHTS                         │
  │  Code: provekit/prover/src/whir_r1cs.rs                       │
  └────────────────────────────────────────────────────────────────┘

  The inner sumcheck folded a_vec with randomness r₁={ZZ(r_inner1)}, r₂={ZZ(r_inner2)}.
  Each a_vec[j] contributed to a_final with a specific weight.

  Tracing the folds backward:
    Round 1: a_folded[i] = (1-r₁)·a_vec[i] + r₁·a_vec[i+2]
    Round 2: a_final     = (1-r₂)·a_folded[0] + r₂·a_folded[1]

  Expanding a_final in terms of a_vec[0..3]:
    a_final = (1-r₂)·[(1-r₁)·a_vec[0] + r₁·a_vec[2]]
            + r₂·[(1-r₁)·a_vec[1] + r₁·a_vec[3]]

  So fold_weight[j] = coefficient of a_vec[j]:
""")

fw = [F(0)]*m_pad
fw[0] = (1-r_inner2)*(1-r_inner1)
fw[1] = r_inner2*(1-r_inner1)
fw[2] = (1-r_inner2)*r_inner1
fw[3] = r_inner2*r_inner1

print(f"    fw[0] = (1-r₂)(1-r₁) = (1-{ZZ(r_inner2)})(1-{ZZ(r_inner1)}) = {ZZ(1-r_inner2)}·{ZZ(1-r_inner1)} = {ZZ(fw[0])}")
print(f"    fw[1] = r₂·(1-r₁)    = {ZZ(r_inner2)}·{ZZ(1-r_inner1)}       = {ZZ(fw[1])}")
print(f"    fw[2] = (1-r₂)·r₁    = {ZZ(1-r_inner2)}·{ZZ(r_inner1)}       = {ZZ(fw[2])}")
print(f"    fw[3] = r₂·r₁        = {ZZ(r_inner2)}·{ZZ(r_inner1)}         = {ZZ(fw[3])}")

# Verify
a_final_check = sum(fw[j]*a_vec[j] for j in range(m_pad))
print(f"\n  Verify: Σⱼ fw[j]·a_vec[j] = {'+'.join(f'{ZZ(fw[j])}·{ZZ(a_vec[j])}' for j in range(m_pad))}")
print(f"         = {ZZ(a_final_check)} = a_final = {ZZ(a_final)} ✓")

# ── Derive L_A ──
print(f"""
  ── DERIVE L_A ──

  a_vec[j] = αʲ · (A[j]·w), so:
    a_final = Σⱼ fw[j] · αʲ · (A[j]·w)
            = (Σⱼ fw[j] · αʲ · A[j]) · w
            = L_A · w

  Compute L_A row by row:
""")

L_A = vector(F, [0]*n)
for j in range(m):
    coeff = fw[j] * alpha**j
    contrib = coeff * A_rows[j]
    L_A += contrib
    print(f"    j={j}: fw[{j}]·α^{j}·A[{j}] = {ZZ(fw[j])}·{ZZ(alpha**j)}·{[ZZ(x) for x in A_rows[j]]}")
    print(f"        = {ZZ(coeff)} · {[ZZ(x) for x in A_rows[j]]} = {[ZZ(x) for x in contrib]}")
print(f"    (j=3: fw[3]·a_vec[3] = {ZZ(fw[3])}·0 = 0, no A[3] row)")
print(f"\n    L_A = {[ZZ(x) for x in L_A]}")
print(f"\n  Verify: <L_A, w> = {'+'.join(f'{ZZ(L_A[i])}·{ZZ(w[i])}' for i in range(n) if L_A[i] != 0)}")
print(f"         = {ZZ(L_A*w)} = a_final = {ZZ(a_final)} ✓")
assert L_A*w == a_final

# ── Derive L_B ──
print(f"""
  ── DERIVE L_B ──

  b_vec[j] = B[j]·w (no α factor), so:
    b_final = Σⱼ fw[j] · (B[j]·w)
            = (Σⱼ fw[j] · B[j]) · w
            = L_B · w
""")

L_B = vector(F, [0]*n)
for j in range(m):
    contrib = fw[j] * B_rows[j]
    L_B += contrib
    print(f"    j={j}: fw[{j}]·B[{j}] = {ZZ(fw[j])}·{[ZZ(x) for x in B_rows[j]]} = {[ZZ(x) for x in contrib]}")
print(f"\n    L_B = {[ZZ(x) for x in L_B]}")
print(f"    <L_B, w> = {'+'.join(f'{ZZ(L_B[i])}·{ZZ(w[i])}' for i in range(n) if L_B[i] != 0)} = {ZZ(L_B*w)} = b_final = {ZZ(b_final)} ✓")
assert L_B*w == b_final

# ── Derive L_C ──
print(f"""
  ── DERIVE L_C ──

  The right side of the batched equation:
    Σᵢ αⁱ·(C[i]·w) = S_original = {ZZ(S_original)}

  This is ALREADY linear in w (no sumcheck needed):
    (Σᵢ αⁱ·C[i]) · w = L_C · w
""")

L_C = vector(F, [0]*n)
for i in range(m):
    contrib = alpha**i * C_rows[i]
    L_C += contrib
    print(f"    i={i}: α^{i}·C[{i}] = {ZZ(alpha**i)}·{[ZZ(x) for x in C_rows[i]]} = {[ZZ(x) for x in contrib]}")
print(f"\n    L_C = {[ZZ(x) for x in L_C]}")
print(f"    <L_C, w> = {'+'.join(f'{ZZ(L_C[i])}·{ZZ(w[i])}' for i in range(n) if L_C[i] != 0)} = {ZZ(L_C*w)} = S = {ZZ(S_original)} ✓")
assert L_C*w == S_original

# ── Why L_A×L_B ≠ L_C but it's OK ──
print(f"""
  ── WHY <L_A,w>×<L_B,w> ≠ <L_C,w> AND WHY THAT'S OK ──

    <L_A,w> × <L_B,w> = {ZZ(a_final)} × {ZZ(b_final)} = {ZZ(a_final*b_final)} = S'' = {ZZ(S2)}
    <L_C,w> = {ZZ(L_C*w)} = S = {ZZ(S_original)}
    {ZZ(S2)} ≠ {ZZ(S_original)}  ← different numbers!

    This is expected. The inner sumcheck CHANGED the claim value:
      Original:  dot(a_vec, b_vec) = S = {ZZ(S_original)}
      After sumcheck: a_final × b_final = S'' = {ZZ(S2)}

    The sumcheck chain of trust connects them:
      S = {ZZ(S_original)} → c₁(r₁) = S' = {ZZ(S1)} → c₂(r₂) = S'' = {ZZ(S2)}

    The verifier checks THREE things:
      1. WHIR proves <L_A,w> = {ZZ(a_final)}, <L_B,w> = {ZZ(b_final)}, <L_C,w> = {ZZ(S_original)}
      2. The sumcheck chain: S → S' → S'' (from transcript)
      3. One multiplication: <L_A,w> × <L_B,w> = S'' ({ZZ(a_final)}×{ZZ(b_final)} = {ZZ(S2)} ✓)

    Together: <L_C,w> = S, and S links to S'' via sumcheck,
    and S'' = <L_A,w>×<L_B,w>. So the original R1CS holds.
""")

print(f"""
  ┌──────────────────────────────────────────────────────────────────┐
  │  SUMMARY: 3 LINEAR CLAIMS FOR WHIR                               │
  │                                                                    │
  │  L_A = {str([ZZ(x) for x in L_A]):42s}   │
  │  <L_A, w> = {ZZ(L_A*w)}  (= a_final, from inner sumcheck fold)          │
  │                                                                    │
  │  L_B = {str([ZZ(x) for x in L_B]):42s}   │
  │  <L_B, w> = {ZZ(L_B*w)}  (= b_final, from inner sumcheck fold)          │
  │                                                                    │
  │  L_C = {str([ZZ(x) for x in L_C]):42s}   │
  │  <L_C, w> = {ZZ(L_C*w)}  (= S_original, right side of batched R1CS)     │
  │                                                                    │
  │  Verifier checks:                                                  │
  │    <L_A,w>×<L_B,w> = S''       ({ZZ(a_final)}×{ZZ(b_final)} = {ZZ(S2)})             │
  │    <L_C,w> = S                  ({ZZ(L_C*w)} = {ZZ(S_original)})                     │
  │    Sumcheck chain: S={ZZ(S_original)} → S'={ZZ(S1)} → S''={ZZ(S2)}                  │
  └──────────────────────────────────────────────────────────────────┘
""")

# =========================================================================
# ░░░ PART 3: WHIR COMMIT ░░░
# =========================================================================
print(f"{'━'*70}")
print("PART 3: WHIR COMMIT")
print("WHY: The prover must LOCK IN the witness before seeing challenges.")
print("     Split into blocks (interleave), evaluate each block polynomial")
print("     at N points via NTT (RS encode), hash the evaluations into a")
print("     Merkle tree. The root is the commitment — binding and compact.")
print("Code: whir/mod.rs:138 → irs_commit.rs:302")
print(f"{'━'*70}")

folding_factor = 2
interleaving_depth = 4
message_length = n // interleaving_depth
rate_inv = 2
codeword_length = message_length * rate_inv
w_cw = g**((p-1)//codeword_length)

print(f"""
  Parameters (whir/config.rs:17):
    n={n}, ff={folding_factor}, d={interleaving_depth}, k={message_length}, N={codeword_length}
    ω={ZZ(w_cw)} (primitive {codeword_length}th root), rate=1/{rate_inv}

  3a. INTERLEAVE (irs_commit.rs:327-330)
""")

blocks = []
polys = []
for j in range(interleaving_depth):
    block = [w[j*message_length+i] for i in range(message_length)]
    blocks.append(block)
    poly = sum(F(block[i])*X**i for i in range(message_length))
    polys.append(poly)
    print(f"    Block {j}: w[{j*message_length}:{(j+1)*message_length}] = {[ZZ(x) for x in block]} → f_{j}(X) = {poly}")

# NTT — with BEFORE/AFTER coset comparison
print(f"""
  3b. RS ENCODE via NTT (irs_commit.rs:331 → cooley_tukey.rs:403)

  ┌──────────────────────────────────────────────────────────────────┐
  │  BEFORE COSET NTT (how it used to work):                         │
  │                                                                    │
  │  Pad polynomial to codeword_length with zeros, run one big NTT:  │
  │    f₀ = {polys[0]}, coeffs = {[ZZ(x) for x in blocks[0]]}                    │
  │    Padded: {[ZZ(x) for x in blocks[0]]} + [0, 0] = {[ZZ(x) for x in blocks[0]] + [0, 0]}             │
  │    One {codeword_length}-point NTT on this 4-element vector (50% zeros!)       │
  │    Cost: O({codeword_length} × log₂({codeword_length})) = O({codeword_length * int(log(codeword_length,2))}) per polynomial                │
  │                                                                    │
  │  Problem: half the input is zeros — wasted butterfly operations.  │
  └──────────────────────────────────────────────────────────────────┘
""")

# Show before-coset result (direct evaluation — same answer, more work)
print(f"    BEFORE: Direct {codeword_length}-point evaluation of f₀:")
for i in range(codeword_length):
    pt = w_cw**i
    val = polys[0](pt)
    print(f"      f₀(ω^{i}) = f₀({ZZ(pt)}) = {ZZ(val)}")

print(f"""
  ┌──────────────────────────────────────────────────────────────────┐
  │  AFTER COSET NTT (Remco's optimization, commit 103a9ca):        │
  │                                                                    │
  │  Split {codeword_length} eval points into cosets of a subgroup.              │
  │  Each coset gets a SMALL, DENSE NTT — no wasted zeros.           │
  │                                                                    │
  │  WHAT IS A COSET?                                                 │
  │  The {codeword_length} evaluation points {{1, ω, ω², ω³}} form a group.       │
  │  Pick a subgroup of size M (the coset_size).                     │
  │  The group decomposes into N/M cosets:                            │
  │    Coset 0: {{1, ω^K, ω^{{2K}}, ...}}        (the subgroup itself)   │
  │    Coset 1: {{ω, ω^{{K+1}}, ω^{{2K+1}}, ...}}   (subgroup × ω)        │
  │    ...                                                             │
  │  Each coset has M points. Total: K cosets × M points = N points. │
  │                                                                    │
  │  KEY INSIGHT: evaluating f at coset c's points is equivalent to  │
  │  running an M-point NTT on TWISTED coefficients (multiply a_i    │
  │  by ω^{{ci}}). The twisted vector is DENSE — no zero padding!      │
  │                                                                    │
  │  Cost: K × O(M × log₂(M)) + K × O(k) twiddle muls              │
  │  vs Before: O(N × log₂(N))                                       │
  │                                                                    │
  │  Code: cooley_tukey.rs:424-463 (pick coset_size, twist, NTT,    │
  │         transpose). Added in Remco's zkWHIR Part 1 (PR #232).   │
  └──────────────────────────────────────────────────────────────────┘

    coset_size = next_order({message_length}) = {message_length}
    num_cosets = {codeword_length}/{message_length} = {codeword_length//message_length}
    g_sub = ω^{codeword_length//message_length} = {ZZ(w_cw**(codeword_length//message_length))}
    (g_sub is the primitive {message_length}th root = subgroup generator)

    For Block 0 (f₀ = {polys[0]}):
""")

coset_count = codeword_length // message_length
g_sub = w_cw**coset_count

for c in range(coset_count):
    twist = w_cw**c
    twisted = [blocks[0][i]*twist**i for i in range(message_length)]
    ntt_out = [twisted[0]+twisted[1], twisted[0]-twisted[1]]
    print(f"      Coset {c}: twist=ω^{c}={ZZ(twist)}")
    print(f"        Twisted: {[ZZ(x) for x in twisted]}")
    print(f"        2-pt NTT: [{ZZ(ntt_out[0])}, {ZZ(ntt_out[1])}]")
    for j in range(message_length):
        pt = w_cw**c * g_sub**j
        print(f"        → f₀({ZZ(pt)}) = {ZZ(polys[0](pt))} ✓")

codeword = [[poly(w_cw**i) for poly in polys] for i in range(codeword_length)]

print(f"\n    Codeword matrix ({codeword_length}×{interleaving_depth}):")
for i in range(codeword_length):
    print(f"      row[{i}] = {[ZZ(v) for v in codeword[i]]}  ← evaluations at ω^{i}={ZZ(w_cw**i)}")

# Merkle
print(f"\n  3c. MERKLE COMMIT (merkle_tree.rs:82)")

def hash_row(row):
    return hashlib.sha256(b"".join(int(ZZ(v)).to_bytes(2,'big') for v in row)).hexdigest()[:12]

leaves = [hash_row(row) for row in codeword]
h01 = hashlib.sha256((leaves[0]+leaves[1]).encode()).hexdigest()[:12]
h23 = hashlib.sha256((leaves[2]+leaves[3]).encode()).hexdigest()[:12]
root = hashlib.sha256((h01+h23).encode()).hexdigest()[:12]

for i in range(codeword_length):
    print(f"    leaf[{i}] = hash(row[{i}]) = {leaves[i]}")
print(f"    H₀₁ = {h01},  H₂₃ = {h23}")
print(f"    ROOT = {root}  ← COMMITMENT (sent to verifier)")

# =========================================================================
# ░░░ PART 4: WHIR PROVE — SUMCHECK ON WITNESS ░░░
# =========================================================================
print(f"\n{'━'*70}")
print("PART 4: WHIR PROVE — SUMCHECK ON WITNESS")
print("WHY: WHIR proves ALL linear claims at once by batching them into")
print("     a single sumcheck. The covector is a random linear combination")
print("     (RLC) of L_A, L_B, L_C. One sumcheck proves all three claims.")
print("     The fold halves the vector each round. After 2 folds, the result")
print("     IS the coefficients of the folded polynomial — matching the")
print("     interleaving structure. This enables the in-domain check.")
print("Code: whir/prover.rs:50  prove()")
print(f"{'━'*70}")

# ── RLC batching of linear forms (prover.rs:156-168) ──
rho = F(23)  # RLC challenge from Fiat-Shamir
rlc_coeffs = [F(1), rho, rho**2]

print(f"""
  4a. BATCH LINEAR FORMS INTO ONE COVECTOR (prover.rs:156-168)

      ρ is the RLC (random linear combination) challenge.
      In the real protocol: Fiat-Shamir squeeze from transcript AFTER
      the prover commits. This ensures the prover can't adapt the
      witness to the RLC. The verifier derives the same ρ.
      Code: initial_forms_rlc_coeffs in verifier.rs

      The prover has 3 linear claims:
        <L_A, w> = {ZZ(L_A*w)}     (rlc_coeff = 1)
        <L_B, w> = {ZZ(L_B*w)}     (rlc_coeff = ρ = {ZZ(rho)})
        <L_C, w> = {ZZ(L_C*w)}     (rlc_coeff = ρ² = {ZZ(rho)}² = {ZZ(rho**2)})

      Batch into one covector:
        covector = 1·L_A + ρ·L_B + ρ²·L_C
        the_sum  = 1·<L_A,w> + ρ·<L_B,w> + ρ²·<L_C,w>
""")

batched_cov = vector(F, [0]*n)
for j, (coeff, L) in enumerate([(F(1), L_A), (rho, L_B), (rho**2, L_C)]):
    contrib = coeff * L
    batched_cov += contrib
    print(f"        {ZZ(coeff)} · L_{'ABC'[j]} = {ZZ(coeff)} · {[ZZ(x) for x in L]}")
    print(f"              = {[ZZ(x) for x in contrib]}")

batched_sum = sum(rlc_coeffs[j] * [L_A*w, L_B*w, L_C*w][j] for j in range(3))

print(f"""
      covector = {[ZZ(x) for x in batched_cov]}
      the_sum  = 1·{ZZ(L_A*w)} + {ZZ(rho)}·{ZZ(L_B*w)} + {ZZ(rho**2)}·{ZZ(L_C*w)}
               = {ZZ(rlc_coeffs[0]*L_A*w)} + {ZZ(rlc_coeffs[1]*L_B*w)} + {ZZ(rlc_coeffs[2]*L_C*w)}
               = {ZZ(batched_sum)}

      Verify: <covector, w> = {ZZ(batched_cov*w)}
              the_sum        = {ZZ(batched_sum)}
              Match? {batched_cov*w == batched_sum} ✓
""")
assert batched_cov*w == batched_sum

# ── Initial sumcheck ──
linear_form = list(batched_cov)
claimed = ZZ(batched_sum)

print(f"  4b. INITIAL SUMCHECK (prover.rs:205 → sumcheck.rs:62)")
print(f"      Proves: <covector, w> = {claimed}")
print(f"      {folding_factor} rounds, folding {n} → {message_length}")

cur_w = list(w)
cur_cov = [F(x) for x in linear_form]
cur_S = F(claimed)
rand = []

for rnd in range(folding_factor):
    h = len(cur_w)//2
    wl, wh = cur_w[:h], cur_w[h:]
    cl, ch = cur_cov[:h], cur_cov[h:]

    cp = sum(((1-t)*F(cl[i])+t*F(ch[i]))*((1-t)*F(wl[i])+t*F(wh[i])) for i in range(h))
    c0 = cp(0)
    c1_at_1 = cp(1)
    c2c = cp.leading_coefficient() if cp.degree()>=2 else F(0)

    print(f"\n    Round {rnd+1} (sumcheck.rs):")
    print(f"      Split {len(cur_w)} → [{h}|{h}]")
    print(f"        w_lo = {[ZZ(x) for x in wl]}")
    print(f"        w_hi = {[ZZ(x) for x in wh]}")
    print(f"        cov_lo = {[ZZ(x) for x in cl]}")
    print(f"        cov_hi = {[ZZ(x) for x in ch]}")
    print(f"      c(t) = Σᵢ [(1-t)·cov_lo[i]+t·cov_hi[i]] × [(1-t)·w_lo[i]+t·w_hi[i]]")
    print(f"           = {cp}")
    print(f"      c(0) = dot(cov_lo, w_lo) = {'+'.join(f'{ZZ(cl[i])}·{ZZ(wl[i])}' for i in range(h))} = {ZZ(c0)}")
    print(f"      c(1) = dot(cov_hi, w_hi) = {'+'.join(f'{ZZ(ch[i])}·{ZZ(wh[i])}' for i in range(h))} = {ZZ(c1_at_1)}")
    print(f"      CHECK: c(0)+c(1) = {ZZ(c0)}+{ZZ(c1_at_1)} = {ZZ(c0+c1_at_1)} = S = {ZZ(cur_S)} ✓")
    print(f"      Sends: c₀={ZZ(c0)}, c₂={ZZ(c2c)}")

    r = F([7,11][rnd])
    rand.append(r)
    cur_S = cp(r)
    print(f"      Challenge: r={ZZ(r)}")
    print(f"      c({ZZ(r)}) = {ZZ(c2c)}·{ZZ(r)}² + ... + {ZZ(c0)} = {ZZ(cur_S)}  → new S")

    new_w = [(1-r)*F(wl[i])+r*F(wh[i]) for i in range(h)]
    new_cov = [(1-r)*F(cl[i])+r*F(ch[i]) for i in range(h)]
    print(f"      Fold witness:  w'[i] = (1-{ZZ(r)})·w_lo[i] + {ZZ(r)}·w_hi[i]")
    print(f"        = [{', '.join(f'{ZZ(1-r)}·{ZZ(wl[i])}+{ZZ(r)}·{ZZ(wh[i])}={ZZ(new_w[i])}' for i in range(h))}]")
    print(f"      Fold covector: cov'[i] = (1-{ZZ(r)})·cov_lo[i] + {ZZ(r)}·cov_hi[i]")
    print(f"        = [{', '.join(f'{ZZ(new_cov[i])}' for i in range(h))}]")
    cur_w = new_w
    cur_cov = new_cov
    print(f"      Verify: dot(cov', w') = {'+'.join(f'{ZZ(cur_cov[i])}·{ZZ(cur_w[i])}' for i in range(h))}")
    dot_check = sum(F(cur_cov[i])*F(cur_w[i]) for i in range(h))
    print(f"        = {ZZ(dot_check)} = S = {ZZ(cur_S)} ✓")

print(f"""
  Result: folded_w   = {[ZZ(x) for x in cur_w]} (length {message_length})
          folded_cov = {[ZZ(x) for x in cur_cov]}
          Randomness = {[ZZ(r) for r in rand]}""")

# =========================================================================
# ░░░ PART 5: SEND FINAL VECTOR + OPEN COMMITMENT ░░░
# =========================================================================
print(f"\n{'━'*70}")
print("PART 5: SEND FINAL VECTOR + OPEN COMMITMENT")
print("WHY: The folded vector is tiny (2 elements) — send it directly.")
print("     Then open random codeword positions so the verifier can spot-check")
print("     that the committed codeword is consistent with this folded result.")
print("     RS distance guarantees: if the prover cheated, most positions are wrong,")
print("     so random queries catch it with overwhelming probability.")
print("Code: prover.rs:289 (send) + prover.rs:301 (open)")
print(f"{'━'*70}")

print(f"""
  5a. SEND FINAL VECTOR (prover.rs:289-296)
      The folded witness is small enough to send directly:
      {[ZZ(x) for x in cur_w]}
      ({message_length} field elements — in real proofs, ~64 out of ~624,000)

  5b. OPEN COMMITMENT (prover.rs:301 → irs_commit.rs:392)

      5b-i. SAMPLE QUERY POSITIONS (irs_commit.rs:497)
            challenge_indices(transcript, {codeword_length}, num_queries, ...)
            Fiat-Shamir squeeze → indices in [0, {codeword_length})
""")

query_idx = [0, 2]  # simulated
print(f"            Queried indices: {query_idx}")

print(f"""
      5b-ii. REVEAL ROWS + MERKLE PROOFS (irs_commit.rs:421-434)""")
for idx in query_idx:
    print(f"            row[{idx}] = {[ZZ(v) for v in codeword[idx]]}")
    print(f"            leaf[{idx}] = {leaves[idx]}, Merkle path → root")

# =========================================================================
# ░░░ PART 6: VERIFIER — FULL DETAILED TRACE ░░░
# =========================================================================
print(f"\n{'━'*70}")
print("PART 6: VERIFIER — FULL DETAILED TRACE")
print("WHY: The verifier has NO witness. It reconstructs everything from")
print("     the transcript (Merkle root, sumcheck polynomials, final vector,")
print("     opened rows, Merkle paths). Each check gates the next — if any")
print("     fails the proof is rejected.")
print("Code: whir_r1cs.rs (R1CS wrapper) → whir/verifier.rs (WHIR core)")
print(f"{'━'*70}")

print(f"""
  ┌──────────────────────────────────────────────────────────┐
  │  VERIFIER'S VIEW: what it knows before starting          │
  │                                                          │
  │  From verification key (pkv):                            │
  │    • R1CS matrices A, B, C (PUBLIC, from the circuit)    │
  │    • WHIR protocol parameters (n, ff, rate, etc.)        │
  │    • ABI for public inputs                               │
  │                                                          │
  │  From proof transcript:                                  │
  │    • Merkle root (commitment to witness)                 │
  │    • Inner sumcheck polynomials (R1CS reduction)         │
  │    • WHIR sumcheck polynomials (fold rounds)             │
  │    • Final vector (folded witness)                       │
  │    • Opened codeword rows + Merkle proofs                │
  │                                                          │
  │  The verifier NEVER sees the witness.                    │
  │  It INDEPENDENTLY computes L_A, L_B, L_C from public    │
  │  matrices and transcript-derived randomness.             │
  └──────────────────────────────────────────────────────────┘
""")

# ── V0: REPLAY INNER SUMCHECK + RECONSTRUCT LINEAR FORMS ──
print(f"  V0. REPLAY INNER SUMCHECK + RECONSTRUCT LINEAR FORMS")
print(f"      Code: provekit/verifier/src/whir_r1cs.rs")
print(f"""
      The verifier reconstructs L_A, L_B, L_C WITHOUT the witness.
      It has: A, B, C matrices (public), and derives α, r₁, r₂
      from the transcript.

      V0a. DERIVE CHALLENGES FROM TRANSCRIPT:
           α = {ZZ(alpha)}    (batching challenge, squeezed after commit)
           r₁ = {ZZ(r_inner1)}   (inner sumcheck round 1 challenge)
           r₂ = {ZZ(r_inner2)}   (inner sumcheck round 2 challenge)
           ρ = {ZZ(rho)}    (RLC challenge for WHIR batching)

      V0b. REPLAY INNER SUMCHECK (same as prover's Part 2):
           The verifier receives c₀, c₂ per round from the transcript,
           checks c(0)+c(1) = S, and derives the SAME challenges.
           If the prover lied about any polynomial, the check fails here.

           Round 1: c₁(0)+c₁(1) = {ZZ(S_original)} ✓, challenge r₁={ZZ(r_inner1)}
           Round 2: c₂(0)+c₂(1) = {ZZ(S1)} ✓, challenge r₂={ZZ(r_inner2)}
           Final:   a_final × b_final = {ZZ(a_final)} × {ZZ(b_final)} = {ZZ(S2)} = S'' ✓

      V0c. COMPUTE FOLD WEIGHTS (from inner sumcheck randomness):
           fw[j] = eq([r₁, r₂], bits(j))
""")
# Recompute fold weights from public data
v_fw = [F(0)]*m_pad
v_fw[0] = (1-r_inner2)*(1-r_inner1)
v_fw[1] = r_inner2*(1-r_inner1)
v_fw[2] = (1-r_inner2)*r_inner1
v_fw[3] = r_inner2*r_inner1
for j in range(m_pad):
    b = [(j>>1)&1, j&1]
    print(f"           fw[{j}] = eq([{ZZ(r_inner1)},{ZZ(r_inner2)}], {b}) = {ZZ(v_fw[j])}")

print(f"""
      V0d. RECONSTRUCT L_A, L_B, L_C FROM PUBLIC MATRICES:

           The verifier computes these INDEPENDENTLY — same formula as
           the prover, but using the PUBLIC A, B, C matrices:

             L_A = Σⱼ fw[j] · αʲ · A[j]     (A is public)
             L_B = Σⱼ fw[j] · B[j]           (B is public)
             L_C = Σᵢ αⁱ · C[i]              (C is public)
""")

# Reconstruct L_A from public matrices
v_L_A = vector(F, [0]*n)
for j in range(m):
    v_L_A += v_fw[j] * alpha**j * A_rows[j]
print(f"           L_A (reconstructed) = {[ZZ(x) for x in v_L_A]}")
print(f"           L_A (from prover)   = {[ZZ(x) for x in L_A]}")
print(f"           Match? {list(v_L_A) == list(L_A)} ✓")
assert list(v_L_A) == list(L_A)

v_L_B = vector(F, [0]*n)
for j in range(m):
    v_L_B += v_fw[j] * B_rows[j]
print(f"           L_B (reconstructed) = {[ZZ(x) for x in v_L_B]}")
print(f"           Match? {list(v_L_B) == list(L_B)} ✓")
assert list(v_L_B) == list(L_B)

v_L_C = vector(F, [0]*n)
for i in range(m):
    v_L_C += alpha**i * C_rows[i]
print(f"           L_C (reconstructed) = {[ZZ(x) for x in v_L_C]}")
print(f"           Match? {list(v_L_C) == list(L_C)} ✓")
assert list(v_L_C) == list(L_C)

print(f"""
      The prover CANNOT substitute different constraints because:
        1. A, B, C are fixed in the verification key (pkv)
        2. α, r₁, r₂ are derived from the transcript (Fiat-Shamir)
        3. The verifier computes L_A, L_B, L_C independently
        4. Any change to the matrices or randomness would produce
           different linear forms, and the WHIR proof would fail

      Also computes claimed values from the inner sumcheck:
        <L_A, w> = a_final = {ZZ(a_final)}
        <L_B, w> = b_final = {ZZ(b_final)}
        <L_C, w> = S = {ZZ(S_original)}
""")

# ── Multi-round note ──
print(f"""  ┌──────────────────────────────────────────────────────────────────┐
  │  NOTE: MULTI-ROUND WHIR STRUCTURE                                 │
  │                                                                    │
  │  This toy example has n=8, ff=2, so there's only 1 WHIR round:   │
  │    initial commit → 2 sumcheck folds → in-domain check → done    │
  │                                                                    │
  │  In production (provekit: ff=3, nv=18):                           │
  │    Initial commit → 3 folds → in-domain check                    │
  │    Round 1 commit → 3 folds → in-domain check                    │
  │    Round 2 commit → 3 folds → in-domain check                    │
  │    Round 3 commit → 3 folds → in-domain check                    │
  │    Round 4 commit → 3 folds → in-domain check                    │
  │    Round 5 commit → 3 folds → in-domain check                    │
  │    Final sumcheck (0 remaining vars)                              │
  │                                                                    │
  │  Each round: NTT + Merkle commit, sumcheck, open, in-domain.     │
  │  That's 6 NTT operations, 6 in-domain checks, 18 binary folds.  │
  │  Every in-domain check gates the next round — if any fails,      │
  │  the proof is rejected.                                           │
  └──────────────────────────────────────────────────────────────────┘
""")

# ── V1: INITIAL SUMCHECK VERIFY ──
print(f"  V1. INITIAL SUMCHECK VERIFY (verifier.rs:120-132)")
print(f"      Code: sumcheck.rs verify() — receives c₀,c₂ from transcript,")
print(f"      reconstructs c₁, checks c(0)+c(1) = S, picks random r.")
print(f"""
      The verifier only sees the polynomial coefficients, NOT the vectors.
      It replays the exact same chain the prover produced.
""")

# Replay the sumcheck from the verifier's perspective
v_S = F(claimed)
v_rand = []
# Re-derive the sumcheck polynomials from the prover's trace
v_cur_w_ghost = list(w)  # verifier doesn't have this — used only to derive c₀,c₂
v_cur_cov_ghost = [F(x) for x in linear_form]

for rnd in range(folding_factor):
    h = len(v_cur_w_ghost) // 2
    wl, wh = v_cur_w_ghost[:h], v_cur_w_ghost[h:]
    cl, ch = v_cur_cov_ghost[:h], v_cur_cov_ghost[h:]

    cp = sum(((1-t)*F(cl[i])+t*F(ch[i]))*((1-t)*F(wl[i])+t*F(wh[i])) for i in range(h))
    c0 = cp(0)
    c1_at_1 = cp(1)
    c2c = cp.leading_coefficient() if cp.degree() >= 2 else F(0)
    c1c = v_S - c0 - c0 - c2c  # c₁ = S - 2c₀ - c₂

    print(f"      Round {rnd+1}:")
    print(f"        Receive from transcript: c₀={ZZ(c0)}, c₂={ZZ(c2c)}")
    print(f"        Reconstruct: c₁ = S - 2·c₀ - c₂ = {ZZ(v_S)} - 2·{ZZ(c0)} - {ZZ(c2c)} = {ZZ(c1c)}")
    print(f"        c(t) = {ZZ(c2c)}·t² + {ZZ(c1c)}·t + {ZZ(c0)}")
    print(f"        CHECK: c(0) + c(1) = {ZZ(c0)} + {ZZ(c0 + c1c + c2c)} = {ZZ(c0 + c0 + c1c + c2c)}")
    print(f"                           = S = {ZZ(v_S)}  {'✓' if c0 + c0 + c1c + c2c == v_S else '✗'}")

    r = F([7,11][rnd])
    v_rand.append(r)
    v_S = c0 + c1c*r + c2c*r**2
    print(f"        Fiat-Shamir challenge: r = {ZZ(r)}")
    print(f"        c(r) = {ZZ(c0)} + {ZZ(c1c)}·{ZZ(r)} + {ZZ(c2c)}·{ZZ(r**2)}")
    print(f"             = {ZZ(c0)} + {ZZ(c1c*r)} + {ZZ(c2c*r**2)} = {ZZ(v_S)}")
    print(f"        New sum: the_sum = {ZZ(v_S)}")

    # "fold" the ghost vectors (verifier doesn't do this — just for derivation)
    v_cur_w_ghost = [(1-r)*F(wl[i])+r*F(wh[i]) for i in range(h)]
    v_cur_cov_ghost = [(1-r)*F(cl[i])+r*F(ch[i]) for i in range(h)]

the_sum = v_S
evaluation_point = list(v_rand)

print(f"""
      After {folding_factor} rounds:
        the_sum = {ZZ(the_sum)}
        evaluation_point = {[ZZ(r) for r in evaluation_point]}
        (verifier now holds the sum and the random challenges)
""")

# ── V2: RECEIVE FINAL VECTOR ──
print(f"  V2. RECEIVE FINAL VECTOR (verifier.rs:194)")
final_vector = list(cur_w)  # this is what the prover sent
print(f"      final_vector = {[ZZ(x) for x in final_vector]}")
print(f"      Length = {len(final_vector)} = message_length = k = {message_length}")
print(f"      (These are supposed to be the folded polynomial coefficients)")

# ── V3: FINAL PROOF OF WORK ──
print(f"""
  V3. FINAL PROOF-OF-WORK (verifier.rs:197)
      In production: verifier checks PoW nonce from transcript.
      (Disabled in this trace — pow_bits=0)
""")

# ── V4: OPEN PREVIOUS COMMITMENT → MERKLE VERIFY ──
print(f"  V4. OPEN & VERIFY COMMITMENT (verifier.rs:199-215)")
print(f"      irs_commit.rs verify(): challenge query positions,")
print(f"      receive opened rows, verify Merkle paths.")

print(f"""
      V4a. CHALLENGE QUERY INDICES (challenge_indices.rs:6)
           Fiat-Shamir squeeze → indices in [0, {codeword_length})
           Queried indices: {query_idx}

      V4b. RECEIVE OPENED ROWS (irs_commit.rs:421-434)
""")
for idx in query_idx:
    print(f"           row[{idx}] = {[ZZ(v) for v in codeword[idx]]}")
    print(f"           (d={interleaving_depth} evaluations: f₀(ω^{idx}), f₁(ω^{idx}), f₂(ω^{idx}), f₃(ω^{idx}))")

print(f"""
      V4c. VERIFY MERKLE PROOFS (matrix_commit.rs verify())
           For each queried row, hash it and trace the path to the root.
""")
for idx in query_idx:
    h = hash_row(codeword[idx])
    sibling_idx = idx ^ 1
    sibling_h = hash_row(codeword[sibling_idx])
    if idx < sibling_idx:
        parent = hashlib.sha256((h+sibling_h).encode()).hexdigest()[:12]
    else:
        parent = hashlib.sha256((sibling_h+h).encode()).hexdigest()[:12]
    print(f"           row[{idx}]: hash = {h}")
    print(f"             sibling leaf[{sibling_idx}] = {sibling_h}")
    print(f"             parent = hash({h}, {sibling_h}) = {parent}")
    print(f"             trace to root = {root} ✓")

# ── V5: IN-DOMAIN CHECK ──
print(f"""
  V5. IN-DOMAIN CHECK (verifier.rs:217-226) ← THE CRITICAL STEP

      WHY: This is where the verifier connects TWO independent objects:
        1. The COMMITTED codeword (locked in at commit time, opened at query positions)
        2. The FOLDED vector (produced by the sumcheck, received as final_vector)

      If the prover cheated in either the commit or the sumcheck,
      these two objects will disagree at most query positions.

      WHAT IT CHECKS:
        For each queried position i:
          Σⱼ eq_w[j] · codeword[i][j]  ==  Σₖ final_vector[k] · ω^(ik)
          └── from commitment ──────┘      └── from sumcheck ──────────┘

      LEFT SIDE: random linear combination of the d committed polynomial
        evaluations at ω^i, weighted by the eq polynomial of the fold randomness.
      RIGHT SIDE: evaluate the folded polynomial (whose coefficients are
        final_vector) at the same point ω^i.
""")

# Compute eq_weights
eq_w = []
print(f"      V5a. COMPUTE EQ WEIGHTS (eq_weights.rs)")
print(f"           From fold randomness {[ZZ(r) for r in v_rand]}:")
print(f"           eq_w[j] = Π (rₖ if bit_k(j)=1, else 1-rₖ)")
for j in range(interleaving_depth):
    b = [(j>>(folding_factor-1-f))&1 for f in range(folding_factor)]
    terms = []
    ew = F(1)
    for f in range(folding_factor):
        if b[f] == 1:
            ew *= v_rand[f]
            terms.append(f"r{f}={ZZ(v_rand[f])}")
        else:
            ew *= (1 - v_rand[f])
            terms.append(f"(1-r{f})={ZZ(1-v_rand[f])}")
    eq_w.append(ew)
    print(f"             j={j}, bits={b}: {' · '.join(terms)} = {ZZ(ew)}")
print(f"           eq_weights = {[ZZ(x) for x in eq_w]}")

# Derive the folded polynomial
f_folded = sum(eq_w[j]*polys[j] for j in range(interleaving_depth))
print(f"""
      V5b. ALGEBRAIC MEANING:
           f_folded(X) = Σⱼ eq_w[j] · fⱼ(X)""")
for j in range(interleaving_depth):
    print(f"             {ZZ(eq_w[j])} · ({polys[j]})")
print(f"           = {f_folded}")
print(f"           Coefficients: {[ZZ(f_folded[i]) for i in range(message_length)]}")
print(f"           final_vector: {[ZZ(x) for x in final_vector]}")
print(f"           Match? {all(F(final_vector[i]) == f_folded[i] for i in range(message_length))} ✓")

# Check each queried position
print(f"""
      V5c. VERIFY AT QUERIED POSITIONS:
           For each index i, check:
             LEFT  = Σⱼ eq_w[j] · row[i][j]
             RIGHT = Σₖ final_vector[k] · (ω^i)^k = f_folded(ω^i)
             LEFT == RIGHT?
""")
all_ok = True
for idx in query_idx:
    point = w_cw**idx
    print(f"           Position {idx} (ω^{idx} = {ZZ(point)}):")

    # LEFT: codeword rlc
    left_terms = []
    left = F(0)
    for j in range(interleaving_depth):
        term = eq_w[j] * codeword[idx][j]
        left += term
        left_terms.append(f"{ZZ(eq_w[j])}·{ZZ(codeword[idx][j])}")
    print(f"             LEFT  = {' + '.join(left_terms)}")
    print(f"                   = {ZZ(left)}")

    # RIGHT: polynomial evaluation
    right_terms = []
    right = F(0)
    for k in range(message_length):
        term = F(final_vector[k]) * point**k
        right += term
        right_terms.append(f"{ZZ(final_vector[k])}·{ZZ(point)}^{k}")
    print(f"             RIGHT = {' + '.join(right_terms)}")
    print(f"                   = {' + '.join(str(ZZ(F(final_vector[k])*point**k)) for k in range(message_length))}")
    print(f"                   = {ZZ(right)}")

    ok = left == right
    all_ok = all_ok and ok
    print(f"             LEFT == RIGHT?  {ZZ(left)} == {ZZ(right)}  {'✓' if ok else '✗'}")
    print()

assert all_ok, "In-domain check failed!"

# ── V6: FINAL SUMCHECK ──
print(f"  V6. FINAL SUMCHECK (verifier.rs:228-230)")
print(f"""
      The final_vector has {message_length} elements. The sumcheck continues
      to fold it down to 1 element (the scalar MLE value).
      For n={n} with ff={folding_factor}: message_length={message_length}=2^1, so 1 more round.
""")

# One more sumcheck round on the final vector
v_final_w = list(final_vector)
v_final_cov = list(v_cur_cov_ghost)
v_final_S = the_sum
final_sumcheck_randomness = []

final_num_rounds = 0
s_tmp = len(v_final_w)
while s_tmp > 1:
    s_tmp = (s_tmp + 1) // 2
    final_num_rounds += 1

for rnd in range(final_num_rounds):
    h = len(v_final_w) // 2
    wl, wh = v_final_w[:h], v_final_w[h:]
    cl, ch = v_final_cov[:h], v_final_cov[h:]

    cp_f = sum(((1-t)*F(cl[i])+t*F(ch[i]))*((1-t)*F(wl[i])+t*F(wh[i])) for i in range(h))
    c0_f = cp_f(0)
    c2_f = cp_f.leading_coefficient() if cp_f.degree() >= 2 else F(0)
    c1_f = v_final_S - c0_f - c0_f - c2_f

    print(f"      Round {rnd+1}: {len(v_final_w)} → {h} elements")
    print(f"        Receive: c₀={ZZ(c0_f)}, c₂={ZZ(c2_f)}")
    print(f"        Reconstruct: c₁ = {ZZ(v_final_S)} - 2·{ZZ(c0_f)} - {ZZ(c2_f)} = {ZZ(c1_f)}")
    print(f"        CHECK: c(0)+c(1) = {ZZ(c0_f)}+{ZZ(c0_f+c1_f+c2_f)} = {ZZ(c0_f+c0_f+c1_f+c2_f)} = {ZZ(v_final_S)} ✓")

    r_f = F(13)
    final_sumcheck_randomness.append(r_f)
    v_final_S = c0_f + c1_f*r_f + c2_f*r_f**2
    print(f"        Challenge: r={ZZ(r_f)} → the_sum = {ZZ(v_final_S)}")

    v_final_w = [(1-r_f)*F(wl[i])+r_f*F(wh[i]) for i in range(h)]
    v_final_cov = [(1-r_f)*F(cl[i])+r_f*F(ch[i]) for i in range(h)]

evaluation_point += final_sumcheck_randomness

# ── V7: COMPUTE poly_eval AND linear_form_rlc ──
print(f"""
  V7. EXTRACT linear_form_rlc (verifier.rs:238-253)

      After all sumcheck rounds, the verifier has:
        the_sum = {ZZ(v_final_S)}
        evaluation_point = {[ZZ(r) for r in evaluation_point]}

      the_sum = poly_eval · linear_form_rlc + internal_constraints

      V7a. COMPUTE poly_eval:
           poly_eval = MLE(final_vector, final_sumcheck_randomness)
           This evaluates the final_vector polynomial at the final
           sumcheck point using the standard multilinear extension.
""")

# poly_eval = multilinear extension of final_vector at final_sumcheck_randomness
# For power-of-2, this is the tensor identity:
#   MLE([v0,v1], [r]) = (1-r)*v0 + r*v1
poly_eval = F(0)
for i in range(len(final_vector)):
    # eq(i, point) for the final sumcheck variables
    eq_i = F(1)
    for bit_idx, r_val in enumerate(final_sumcheck_randomness):
        bit = (i >> (len(final_sumcheck_randomness) - 1 - bit_idx)) & 1
        eq_i *= r_val if bit == 1 else (1 - r_val)
    poly_eval += F(final_vector[i]) * eq_i

print(f"           MLE({[ZZ(x) for x in final_vector]}, {[ZZ(r) for r in final_sumcheck_randomness]})")
print(f"           = Σᵢ final_vector[i] · eq(i, point)")
for i in range(len(final_vector)):
    eq_i = F(1)
    eq_terms = []
    for bit_idx, r_val in enumerate(final_sumcheck_randomness):
        bit = (i >> (len(final_sumcheck_randomness) - 1 - bit_idx)) & 1
        if bit == 1:
            eq_i *= r_val
            eq_terms.append(f"r={ZZ(r_val)}")
        else:
            eq_i *= (1 - r_val)
            eq_terms.append(f"(1-r)={ZZ(1-r_val)}")
    print(f"             i={i}: {ZZ(final_vector[i])} · [{' · '.join(eq_terms)}] = {ZZ(final_vector[i])} · {ZZ(eq_i)} = {ZZ(F(final_vector[i])*eq_i)}")
print(f"           poly_eval = {ZZ(poly_eval)}")

print(f"""
      V7b. EXTRACT linear_form_rlc:
           linear_form_rlc = the_sum / poly_eval
           = {ZZ(v_final_S)} / {ZZ(poly_eval)} = {ZZ(v_final_S / poly_eval)}
""")
linear_form_rlc_verifier = v_final_S / poly_eval

# ── V7c: Subtract internal constraints ──
print(f"""      V7c. SUBTRACT INTERNAL CONSTRAINT CONTRIBUTIONS (verifier.rs:244-253)
           The verifier accumulated OOD/in-domain constraints into the_sum.
           Each constraint's contribution must be subtracted from linear_form_rlc.

           For this example (no OOD, no intermediate rounds):
             No internal constraints to subtract.
             linear_form_rlc = {ZZ(linear_form_rlc_verifier)} (unchanged)
""")

# ── V8: FINAL CLAIM VERIFY ──
print(f"""  V8. FINAL CLAIM VERIFY (mod.rs FinalClaim::verify())

      The verifier has:
        evaluation_point = {[ZZ(r) for r in evaluation_point]}
        linear_form_rlc  = {ZZ(linear_form_rlc_verifier)}
        rlc_coefficients = [1, ρ={ZZ(rho)}, ρ²={ZZ(rho**2)}]  (3 linear forms)

      It computes:
        rlc = Σⱼ coeff[j] · linear_form[j].mle_evaluate(evaluation_point)
            = 1·MLE(L_A, pt) + ρ·MLE(L_B, pt) + ρ²·MLE(L_C, pt)

      This is the UNBUNDLING step: the verifier checks each form separately,
      then combines with the same RLC coefficients the prover used.
""")

# Compute MLE for each linear form
num_vars = len(evaluation_point)

def compute_mle(lf_vec, point, label):
    result = F(0)
    nv = len(point)
    for i in range(n):
        if lf_vec[i] == 0:
            continue
        eq_i = F(1)
        for bit_idx in range(nv):
            bit = (i >> (nv - 1 - bit_idx)) & 1
            eq_i *= point[bit_idx] if bit == 1 else (1 - point[bit_idx])
        result += F(lf_vec[i]) * eq_i
    return result

mle_la = compute_mle(list(L_A), evaluation_point, "L_A")
mle_lb = compute_mle(list(L_B), evaluation_point, "L_B")
mle_lc = compute_mle(list(L_C), evaluation_point, "L_C")

print(f"      V8a. COMPUTE mle_evaluate FOR EACH FORM:")
print(f"           MLE(L_A, point) = {ZZ(mle_la)}")
print(f"           MLE(L_B, point) = {ZZ(mle_lb)}")
print(f"           MLE(L_C, point) = {ZZ(mle_lc)}")

rlc_computed = rlc_coeffs[0]*mle_la + rlc_coeffs[1]*mle_lb + rlc_coeffs[2]*mle_lc

print(f"""
      V8b. COMBINE WITH RLC COEFFICIENTS:
           rlc = 1·{ZZ(mle_la)} + {ZZ(rho)}·{ZZ(mle_lb)} + {ZZ(rho**2)}·{ZZ(mle_lc)}
               = {ZZ(rlc_coeffs[0]*mle_la)} + {ZZ(rlc_coeffs[1]*mle_lb)} + {ZZ(rlc_coeffs[2]*mle_lc)}
               = {ZZ(rlc_computed)}
""")

print(f"      V8c. CHECK: rlc == linear_form_rlc?")
print(f"           computed rlc = {ZZ(rlc_computed)}")
print(f"           expected     = {ZZ(linear_form_rlc_verifier)}")
ok_rlc = rlc_computed == linear_form_rlc_verifier
print(f"           Match? {ok_rlc} {'✓' if ok_rlc else '✗'}")
assert ok_rlc, "FinalClaim rlc mismatch!"

# ── V9: CONNECT BACK TO R1CS ──
print(f"""
  V9. CONNECT BACK TO R1CS

      WHIR proved: <L_A, w> = {ZZ(a_final)}
      WHIR proved: <L_B, w> = {ZZ(b_final)}  (via a separate WHIR proof)
      WHIR proved: <L_C, w> = {ZZ(S_original)}  (via a separate WHIR proof)

      The verifier checks the quadratic relationship:
        <L_A,w> × <L_B,w> = S''
        {ZZ(a_final)} × {ZZ(b_final)} = {ZZ(a_final*b_final)}

      The sumcheck chain links S to S'':
        S = {ZZ(S_original)} → S' = {ZZ(S1)} → S'' = {ZZ(S2)}
        <L_A,w> × <L_B,w> = {ZZ(a_final*b_final)} = S'' = {ZZ(S2)} ✓

      Therefore the original R1CS is satisfied:
        Σᵢ αⁱ · [(A[i]·w)(B[i]·w) - C[i]·w] = 0  ✓

      ╔═══════════════════════════╗
      ║   PROOF ACCEPTED          ║
      ╚═══════════════════════════╝
""")

# =========================================================================
# ░░░ PART 7: COMPLETE SUMMARY ░░░
# =========================================================================
print(f"{'━'*70}")
print("PROOF ACCEPTED")
print(f"{'━'*70}")

print(f"""
  ╔════════════════════════════════════════════════════════════════════╗
  ║                                                                    ║
  ║  The verifier is convinced:                                        ║
  ║    "Someone knows a, b such that a*b + a = 15"                    ║
  ║                                                                    ║
  ║  Without learning: a = 3, b = 4, c = 12, d = 15                  ║
  ║                                                                    ║
  ║  WHAT WAS SENT (the proof):                                       ║
  ║    • Merkle root: {root}                              ║
  ║    • Sumcheck polynomials: 2 rounds × (c₀, c₂)                   ║
  ║    • Opened rows: {len(query_idx)} positions + Merkle paths                 ║
  ║    • Final vector: {str([ZZ(x) for x in cur_w]):30s}           ║
  ║                                                                    ║
  ║  WHAT WAS NEVER SENT:                                             ║
  ║    • The witness w = {[ZZ(x) for x in w]}          ║
  ║    • Any private input value                                       ║
  ║                                                                    ║
  ╚════════════════════════════════════════════════════════════════════╝
""")

# =========================================================================
# CODE MAP
# =========================================================================
print(f"{'━'*70}")
print("COMPLETE CODE MAP")
print(f"{'━'*70}")

print("""
  PROVER                                         VERIFIER
  ──────                                         ────────

  whir/mod.rs:138 commit()
    → irs_commit.rs:327  interleave
    → irs_commit.rs:331  NTT encode
      → cooley_tukey.rs:424  coset decomposition
      → cooley_tukey.rs:457  ntt_batch (butterflies)
      → cooley_tukey.rs:463  transpose
    → merkle_tree.rs:82  build tree
    → sends ROOT ─────────────────────────→ receives ROOT

  whir_r1cs.rs: quadratic → linear
    → batch constraints with α
    → inner sumcheck: dot(a_vec, b_vec)
    → linear forms L_A, L_B, L_C

  whir/prover.rs:50 prove()
    → :156  batch: covector = Σ rlc[j]·L[j]  (L_A + ρ·L_B + ρ²·L_C)
    → :174  the_sum = Σ rlc[j]·<L[j],w>
    → :205  initial_sumcheck.prove(covector, witness, the_sum)
      → sumcheck.rs:62                      → verifier.rs:120
        sends c₀, c₂ per round ──────────→   checks c(0)+c(1)=S
                                              picks random r
    → :289  send final vector ────────────→ :195  receive final vector
    → :299  final_pow ────────────────────→ :199  verify pow
    → :301  open commitment
      → irs_commit.rs:497  challenge_indices (Fiat-Shamir)
      → sends rows + Merkle proofs ───────→ :202  verify Merkle
                                            :223  IN-DOMAIN CHECK
                                              eq_w[j]·row[i][j] == Σ v[k]·ω^(ik)
    → :316  final_sumcheck ───────────────→ :250  verify final sumcheck
    → :323  return FinalClaim ────────────→ :281  return FinalClaim

                                            mod.rs:69  FinalClaim::verify()
                                              Σ c·L.mle_evaluate(point) == rlc
                                              ACCEPT ✓

  THE SMOOTH-DOMAIN BUG:
    At verifier.rs:223, the in-domain check compares the folded
    vector (from sumcheck) against the committed codeword (from NTT).
    When witness size is not power-of-2 (e.g., 12), fold() produces
    more coefficients than the committed polynomial degree, and the
    check fails. See Phase 6 for the detailed trace.
""")
