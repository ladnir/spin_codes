# Smaller quarter-rate outer: target operating points

Status, 2026-09-11: **both operating points are now outward-certified**, with
the existing inner unchanged. The [certificate and replay](SMALLER_OUTWARD_CERTIFICATE.md)
supersede the earlier binary64-only evidence.

| Relative distance | Certified margin target | Outward margin diagnostic |
|---|---:|---:|
| 16.5% | 40 bits | 41.08348 bits |
| 19% | 30 bits | 30.05251 bits |

Displayed margin diagnostics are truncated; exact dyadic bounds establish
setup failure below 2^-40 and 2^-30. `SMALLER_MARGIN_CERTIFICATE.json` is the
authoritative upper-bound receipt, with a 512-bit replay. The earlier
`SMALLER_MARGIN_CLOSED.json` remains the final binary64 exploration receipt.

This work tests relative distance 16.5% at 40-bit margin and relative distance
19% at 30-bit margin. Fix K=2^20 message bits and the existing optimized
RM2Sub inner with t=128 and s=19. The proposed outer is a fixed BCH-derived
[128,32,32] code, giving N=4K output bits.

The certificates evaluate analytic first-moment bounds outward over
the same setup randomness as `FIXED_INNER_EVALUATION.md`. All reported
failure events include output weight at most floor(delta N). No inner
parameters, permutations, state transitions, or end conditions are changed.
Outward certification is complete; quarter-rate encoder benchmarks remain.

## Fixed construction and exact spectrum

Let P and Q be the primitive binary BCH codes of length 127 with designed
distances 31 and 43. Their dimensions are 36 and 29. `smaller_outer.py`
constructs their generator polynomials p and q over GF(128), using modulus
0x83 and primitive element 2. It checks q=p h and the dimensions directly.

Multiplication by X modulo h cycles through all 127 nonzero residues.
Consequently, cyclic coordinate shifts act transitively on the nonzero
cosets of Q in P. Their weight enumerators are equal. Parity extension
preserves this action and adds coordinate 127.

Define C32 as the span of the extended rows q, Xq, ..., X^28 q, together
with p, Xp, X^2 p. The three added residue classes are independent. The
checker verifies rank 32, parent containment, even parity, and the all-one
word. If C29 and C36 denote the two extended BCH codes, then

    W_C32 = W_C29 + 7 (W_C36 - W_C29) / 127.

The [published C29 spectrum](https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/EBCH/desaki/EBCH128_29.wd)
and [published C36 spectrum](https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/EBCH/desaki/EBCH128_36.wd)
are retained in `sources/`. Their numerical coefficients were compared
against the author-hosted files on 2026-09-11. The
[authors' index](https://isec.ec.okayama-u.ac.jp/home/kusaka/wd/)
attributes these tables to Desaki, Fujiwara, and Kasami's 1997 work.

The reconstructed spectrum has mass 2^32, minimum weight 32, and 588 words
of weight 32. Every coefficient is an integer. Complement symmetry, even
parity, and exact MacWilliams checks pass. `SMALLER_OUTER_AUDIT.json` retains
the fixed generator and complete spectrum. This is an intermediate subcode,
not a claim that the narrow-sense BCH family has dimension exactly 32.

## Earlier numerical calculation and checks

The message has L=32768 outer rows. The evaluator uses the existing exact
selected inner map, including its independently audited image and kernel
spectra. Its occupation Q counts the nonzero outer rows.

- Q=1 uses exact outer-weight support averaging.
- Q=2 uses the existing native pair-support recurrence. This averages two
  independent row supports exactly and uses distinct positions inside each
  region. It replaces the earlier, much looser spectrum-band approximation.
- Q=3 through 64 use the kernel-aware regional transfer and deterministic
  spectrum-band bounds. Composition enumeration separately refines Q=3,4.
- The initial dense calculation uses Q=65 through L and a disjoint cover of zero, ordinary, and all-one row
  counts. The all-one word is an exact atom, not part of the ordinary band.

The sparse log-surprisal grid has spacing 1/8. Its earlier spacing of 1/2
left artificial dips at individual occupations. Refining the numerical
grid changes the bound, not the code or inner map.

At 19%, the remaining dense difficulty was confined to small occupations.
`refine_smaller_dense.py` improves continuous numerical witnesses and splits
boxes without changing the counting measures. `extend_smaller_sparse.py`
then computes Q=65 through 128 directly and intersects the dense cover with
Q>=129. The final dense contribution gives 241.78 bits; the sparse
contribution gives 30.0529 bits. Both ranges meet without a gap or overlap.
At 16.5%, the split remains Q<=64 versus Q>=65; the dense contribution
gives 307.92 bits and the sparse contribution gives 41.0835 bits.

`evaluate_smaller_margins.py` checks every selected dense witness and the
exact integer cover before combining the sparse and dense contributions.
`verify_smaller_margins.py` independently reruns that retained-witness replay.
Sparse bounds are reproduced by rerunning the producer; the standalone
verifier only checks their retained union arithmetic and source identities.
Those exploration calculations use nearest binary64. The subsequent
[outward calculation](SMALLER_OUTWARD_CERTIFICATE.md) recomputes every bound
in Arb, using simpler Q2 composition bounds instead of the native pair kernel.

## Reproduction

Use Python with NumPy and SciPy. The existing inner's generated implementation
manifest is needed for the map-identity comparison; if absent, run
`python workstreams/bare_bch_rm2sub/generate.py` first.

Build the unchanged native pair-support helper on Windows:

```powershell
& workstreams/finite_asymptotic_theory/landscape_db/build_activation_q2.ps1
```

The native helper's existing tests include exact support enumeration and
agreement with the Q1 recurrence on the axes. It is a diagnostic arithmetic
kernel, not an encoder benchmark. Then run:

```sh
python workstreams/rate_quarter_bch/smaller_outer.py
python workstreams/rate_quarter_bch/evaluate_smaller_margins.py
python workstreams/rate_quarter_bch/refine_smaller_dense.py
python workstreams/rate_quarter_bch/extend_smaller_sparse.py
python workstreams/rate_quarter_bch/verify_smaller_margins.py
python -m unittest discover -s workstreams/rate_quarter_bch -p 'test_*.py'
python -m unittest discover -s workstreams/finite_asymptotic_theory/landscape_db -p test_activation_q2.py
```

`SMALLER_MARGIN_SCREEN.json` retains the initial arguments, numerical search grids,
per-occupation margins, dense cover and witnesses, and source hashes. Its
`spectrum` field lists nonzero outer words only; the separate outer audit
contains the zero word as well. No generated binaries are part of the
retained research data. `SMALLER_MARGIN_REFINED.json` retains the continuous
refinement, and `SMALLER_MARGIN_CLOSED.json` retains the final sparse extension
and clipped dense cover. Source hashes bind that chain. The final receipt's
per-row metadata records the additional sparse grid and changed boundary.

## Next step

Implement and benchmark the quarter-rate encoder using the certified
t128_s19 baseline. The selected t64_s16 alternative remains numerical-only.
The t256 investigation is parked; its [deterministic ceiling](inner_calibration/T256_DEDICATED_RUN.md)
remains useful but does not affect the completed t128_s19 certificates.
