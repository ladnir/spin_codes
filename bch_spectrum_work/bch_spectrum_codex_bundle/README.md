# BCH(256) spectrum project

Mathematical arguments and research notes on the intermediate `[256,128,38]`
extended BCH code. This Git snapshot contains notes only; computational
evidence and experiments remain local. See [REPOSITORY_SCOPE.md](REPOSITORY_SCOPE.md)
before using the artifact paths or replay commands below.

Latest: [RANDOM_INNER_M22_CLOSURE.md](RANDOM_INNER_M22_CLOSURE.md).
The original ideal RandomStepConv-M22 goal is closed: k=20, rate 1/2,
relative distance > 0.1, and exact failure upper bound < 0.61 * 2^-40.
The diagnostic margin is 40.72497637925699 bits. No statistical BCH shell
assumption or variance estimate is used. The result does not transfer to
RM2Sub, a PRG, or the draft's different convolution without another argument.

The proof combines [independent Fourier rank certificates](BCH_SHIFT_RANK_CERTIFICATES.md)
with a joint weighted-shell LP and all-occupation transfer bounds.
Run `python -B code/verify_bch_m22_closure.py`; see the final theorem note
for the complete sequential replay list. Next: independent review and paper
integration of this finite theorem.

## Historical milestones

The entries below describe earlier results and open obligations at those times.
The final M22 theorem above supersedes their statements that the goal is open.

Newest lead: [BCH_AFFINE_CANONICAL_AND_SCHAUBPLUS_LEAD.md](BCH_AFFINE_CANONICAL_AND_SCHAUBPLUS_LEAD.md).
Faster exact orbit verification gives a partial weight20 bound. A stronger
dual-distance table suggests a meaningful improvement, but its primary-source
audit is pending. The resulting38.64383-bit calculation is PROVISIONAL;
the certified M22 margin remains37.57062bits and the40bit goal remains open.

Latest: [BCH_HULL_DUAL_WEIGHT18.md](BCH_HULL_DUAL_WEIGHT18.md).
Complete affine enumeration proves A18(HQdual)=0. The deterministic M22
margin is now37.57062bits, with A38(C)<=226659579388193. Exact feasible
hypothetical spectra show that even the full intersection model plus these
exact low shells remains insufficient. The original40bit goal stays open.

Latest: [BCH_HULL_DUAL_LOW_SHELLS.md](BCH_HULL_DUAL_LOW_SHELLS.md).
Complete geometric enumerations prove A14(HQdual)=0 and A16(HQdual)=16592.
The certified M22 margin improves to37.56826bits; A38(C)<=227085449947170.
The original40bit target remains open. All new counts are deterministic,
using a published exact spectrum to certify enumeration completeness.

Latest method: [BCH_HULL_INTERSECTION_CUTS.md](BCH_HULL_INTERSECTION_CUTS.md).
An exact Farkas witness turns intersection constraints into a small valid cut.
The first cut improves the A38 cap only negligibly, to 228870861644645;
the certified M22 margin remains 37.55842 bits. The goal is still open.

Latest: [BCH_HULL_SPECTRUM_CONSTRAINTS.md](BCH_HULL_SPECTRUM_CONSTRAINTS.md).
Coupled hull spectra and a derived exact dimension-69 anchor improve the
deterministic M22 bound to 37.55842 bits. The new A38 cap is 228870861653604.
The original 40-bit target remains open. Exact replay:
`python -B code/audit_bch_hull_anchor_cap.py --verify`.

Latest deterministic refinement: [BCH_LOCAL_SLACK_AND_MOD4.md](BCH_LOCAL_SLACK_AND_MOD4.md).
Exact signed weight identities improve the M22 bound to 36.12681 bits.
The remaining sufficient A38<=10^13 condition is still unproved.

Newest M22 reduction: [M22_ONE_SHELL_REDUCTION.md](M22_ONE_SHELL_REDUCTION.md).
Sharper certified one-row tails reduce the remaining assumption to A38<=10^13;
weights 40 and 42 now use deterministic caps. The full bound is about 40.14372
bits conditional on that single inequality, or 36.11873 bits with only current
deterministic information. The remaining A38 inequality is not proved.

Latest result: [RANDOM_INNER_M25_CLOSURE.md](RANDOM_INNER_M25_CLOSURE.md)
gives a full deterministic-spectrum certificate for ideal RandomStepConv-M25
at the same 2^20 message bits and 10% distance target. Its exact failure bound
is below 2^-41. It uses no statistical shell assumptions, but changes memory
from 22 to 25. Run `python -B code/verify_bch_m25_closure.py` to check it.
The original M22 goal is not declared closed by this parameter change.

The M22 result: [RANDOM_MODEL_FULL_DISTANCE_BOUND.md](RANDOM_MODEL_FULL_DISTANCE_BOUND.md)
gives an all-occupation RandomStepConv-M22 distance bound below 2^-40 for
2^20 message bits, conditional on the three statistically tested shell caps.
The q>=4 tail is now certified below 2^-83. This does not establish deterministic
shell caps or transfer the result to RM2Sub. Verify the saved dependencies and
exact aggregate with `python -B code/verify_bch_full_random_inner.py`.

The direct endpoint test and its exact sample budgets are documented in
`ENDPOINT_SAMPLING_CERTIFICATE.md`. Run `python code/analyze_endpoint_sampling.py`
to check the counting identities and prospective IID zero-hit test budgets.
The fixed 150,361,035-sample run and full independent replay completed with zero
endpoint hits. `generated/endpoint_certificate_20260904/certificate.json` records
conditional statistical acceptance at error 2^-40 under ideal independent
bytes; actual execution used retained Windows CNG pseudorandom bytes.
See `ENDPOINT_APPLICATION_BUDGET_CORRECTION.md`: this shell result alone does
not close the full occupation-one budget.

## First commands

```bash
python code/check_all.py
python code/export_exact_constraints.py
```

`check_all.py` uses only the Python standard library and verifies:

- the complete `EBCH(256,71)` anchor sums to `2^71`;
- its exact MacWilliams transform has dual minimum distance 16;
- the complete `EBCH(256,187)^perp` anchor sums to `2^69`;
- MacWilliams reconstruction reproduces the published low-weight
  `EBCH(256,187)` coefficients and the derived weights 52--60;
- the two Wambach minimum words have punctured weights 37 and 39;
- under the field/coordinate conventions described in the handoff, all four
  natural octal bit-order conventions give trivial affine stabilizer, hence
  orbit size 65,280 for the extended seeds.
- the two BCH generators have dimensions 131 and 123, differ by exactly the
  degree-8 factor for the binary cyclotomic coset of 37, and cyclic shift has
  one orbit on all 255 nonzero elements of the resulting quotient;
- the exact locator-reduction constants, the derivative-zero family bound,
  the saved five-million-sample incidence receipt, and the exhaustive binary
  and GF(4) coefficient-family receipts.

`export_exact_constraints.py` writes `generated/coupled_lp_exact.json`.
It is an integer-coefficient representation of the coupled `q_w, h_w`
BCH-sandwich LP.  This deliberately avoids silently converting the core
constraints to double precision.  It also exports the Wambach seed lower bounds,
the quotient/dual integrality conditions, and the affine 1- and 2-design shell
divisibilities as explicit constraints.

`solve_float_diagnostic.py` reads that model and applies binomial variable scaling
and per-row normalization before calling HiGHS.  By default it replaces the exact
falling-factorial OA equations with their equivalent low-degree Krawtchouk-zero
basis for better numerical conditioning.  Its output is deliberately marked
non-rigorous; use it only to locate candidate optima for subsequent exact checking.

`export_lp.py` converts the continuous part of the integer JSON model to exact
CPLEX LP syntax for QSopt_ex or rational SoPlex.  The separate divisibility
constraints are intentionally not weakened into floating-point rows.  Its default
OA basis contains only degrees 0,2,...,14: complement symmetry makes the odd
Krawtchouk rows identically zero, so this removes exact linear dependencies.

`solve_z3_exact.py` checks feasibility or optimizes the same independent continuous
model directly in exact linear-rational arithmetic.  This avoids accepting a
floating-point infeasibility status or rounded objective as a certificate.

## Mathematical objects

```text
Q = EBCH(256,123,40)
P = EBCH(256,131,38)
Q <= C <= P, dim(C)=128
```

For the common nonzero coset distribution `h_w` of `Q` in `P`:

```text
P_w = Q_w + 255 h_w
C_w = Q_w + 31 h_w
```

Read `handoff/bch_spectrum_codex_handoff.md` before changing the model.

## Highest-value next target

Compute or tightly certify

```text
A_37(BCH(255,131,37))
```

which has the form

```text
A_37 = 4845 m,  m >= 2.
```

Then

```text
A_38(P) = 32640 m
h_38     = 128 m
A_38(C) = 3968 m.
```

The LP should then be rerun with that shell fixed.

Equivalently, let \(L_{37}\) be the normalized locator-polynomial list size.
Then \(L_{37}=19m\), and \(L_{37}<2^{34.093}\) suffices for the application.
Every candidate has at most 37 agreements. A directed-rounding transfer
certificate gives the exact integer sufficient cap
\(A_{38}(C)\leq3827351840403\). Equivalently, a sixth-incidence factor
\(R_6\leq6\) suffices with 0.02347 bits of certified margin. A
five-million-sample importance sampler estimates \(R_6=1.032\), with a normal
95-percent upper endpoint of 1.136. The estimate is evidence, not a certified
point bound. Exact exhaustion of the nested GF(4)-coefficient family gives
only \(1.134\cdot10^{-8}\) of the global \(R_6\) budget and certifies ten new
affine weight-38 orbits. A ten-million-sample invariant-base experiment gives
GF(16)-stratum ratio \(1.00624\), with normal 95-percent upper endpoint
1.01432. Its nine new exact affine orbits improve
\(A_{38}(C)\geq1{,}095{,}168\).

## Dependencies

The verification code is standard-library only.

For an exact production LP solve, use an exact rational solver such as
QSopt_ex or rational SoPlex.  A floating-point solve is useful only as a
diagnostic / candidate optimizer because the low-shell scales are badly
conditioned.

## Files

- `handoff/bch_spectrum_codex_handoff.md` — full research handoff.
- `paper/e104-a_9_1321.pdf` — Fujiwara--Kusaka IEICE paper.
- `data/ebch256_71_half.csv` — exact half-spectrum of the k=71 anchor.
- `data/ebch256_187_dual_half.csv` — exact half-spectrum of the 69-dim dual.
- `data/ebch256_187_low_38_60.csv` — reconstructed k=187 coefficients.
- `data/wambach_min_words.txt` — explicit minimum-word seeds.
- `code/spectrum_exact.py` — exact Krawtchouk/MacWilliams utilities.
- `code/anchors.py` — exact anchor tables.
- `code/affine_wambach.py` — GF(256) and affine-stabilizer checks.
- `code/bch_quotient.py` — exact BCH generator and quotient-action checks.
- `code/export_exact_constraints.py` — exact coupled-LP constraint exporter.
- `code/export_lp.py` — exact continuous-model CPLEX LP exporter.
- `code/solve_float_diagnostic.py` — scaled floating relaxation driver (diagnostic only).
- `code/solve_z3_exact.py` — exact rational feasibility/optimization driver.
- `code/check_all.py` — reproducibility checks.
- `code/locator_incidence_bounds.py` — exact thresholds and degenerate-locus bound.
- `code/sample_locator_extensions.cpp` — targeted 24-point incidence sampler.
- `code/certify_random_inner_threshold_outward.py` — directed application-threshold certificate.
- `code/analyze_low_degree_146_slices.py` — exact low-degree perturbation and Walsh-ridge check.
- `code/analyze_rs_coset_moment_lp.py` — exact primal-dual universal moment bound.
- `code/verify_symmetric_incidence.py` — symmetric-function incidence identity checks.
- `code/analyze_binary_146_family.py` — exhaustive Frobenius-fixed coefficient family.
- `code/analyze_frobenius_invariant_incidence.py` — exhaustive invariant-set and new-orbit receipt.
- `code/sample_frobenius4_incidence.cpp` — AVX2 exhaustive GF(4) invariant-set kernel.
- `code/analyze_frobenius4_invariant_incidence.py` — exact GF(4)-family and new-orbit receipt.
- `code/sample_frobenius16_locator_extensions.cpp` — invariant-base GF(16) importance sampler.
- `code/analyze_frobenius16_locator_sample.py` — sampled GF(16) receipt and exact endpoint closures.
