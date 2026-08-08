# Proof Status And High-Risk TODOs

Updated: 2026-08-07

Purpose: keep the post-restructure proof state explicit.  This is not a new
theorem source; it is the working checklist for what is in good narrative
shape, what is checkable, and what remains too risky to hide in prose.

## Current Narrative State

The broad paper restructure is done.  The compiled manuscript now has one
clean spine:

1. introduction and preliminaries;
2. modular first-moment framework;
3. accumulator warmup;
4. random sliding dense construction;
5. structured local-code construction.

The original exploratory branches remain preserved but uncompiled:
`innerDense.tex`, `integration.tex`, `outerExpandAcc.tex`, and
`innerSparse.tex`.  Topic notes under `explorations/` retain the dense finite
diagnostics, block-recursive BCH inner attempts, BCH spectra work, outer
spectrum comparisons, sparse inner notes, and expander/accumulator outer notes.

## Theorem-Facing Or Checkable Claims

### Accumulator Warmup

Status: theorem-facing warmup.

The accumulator exact enumerator and contraction are in `innerAcc.tex`.  The
combined warmup theorem says a log-memory random sliding dense outer plus the
accumulator inner gives linear distance for sufficiently small constants.  The
role is explanatory; it is not the intended concrete construction.

### Random Sliding Dense Construction

Status: theorem-facing asymptotic line.

The random sliding dense parent section is `randomDenseConstruction.tex`; it
inputs `outerDense.tex`, `innerDenseScalar.tex`, and `integrationDense.tex`.
The key claims are:

- `thm:dense-dense-integration`: linear distance under explicit parameter
  inequalities.
- `thm:dense-dense-explicit-constant`: concrete `0.109` asymptotic checkpoint.

Current verifier:

```powershell
python scripts\verify_dense_claims.py --delta 0.109
```

Latest checked status: sampled-grid verification passes with worst gap
`-0.003107731647`.

### Structured Local-Code Construction

Status: finite checkable certificate line.

The structured parent section is `localCodeStructured.tex`; it inputs the
local-code outer interface, full-split inner interface, RM/EBCH certificate,
and BCH projection note.  The proved finite instantiation is:

- outer: direct sum of `4096` copies of `RM(4,9) [512,256,32]`;
- inner: full-split EBCH `[128,64,22]` with `b=64`;
- length: `N=2^21`;
- target: `delta=.09`, `d=188743`;
- theorem-safe rational/outward first moment: `E[Z_d] <= 2^-37.27`;
- conclusion: some realization is a binary
  `[2^21,2^20,d_min >= 188744]` code.

Current verifier:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
```

Latest checked status: verifier passes and reports
`current_checked_ledger_total_log2,-37.278528`; the underlying stored value is
`-37.278527626...`.  This remains a sharper checked-log diagnostic.  The
theorem-facing rational/outward certificate instead reports
`fullsplit_complete_rational_37_27_bits_status,PASS` and uses `2^-37.27`.

### BCH And Better Local-Code Rows

Status: heuristic projection only.

The BCH256 and BCH512 rows are spectrum-model projections, not theorem claims.
They should not be promoted until an exact spectrum or rigorous low-weight
envelope replaces the modeled spectrum and the affected tails are recomputed.

## High-Risk Items Already Handled

- The broad section structure is now real in the compiled paper, not only in
  planning notes.
- The structured local-code line is now one parent section with local-code
  independent outer and inner interfaces before the RM/EBCH instantiation.
- The random sliding dense line is now one parent section containing the outer,
  inner, and asymptotic integration.
- The dense `0.109` tiny-window proof now uses an explicit
  `xi_tiny=1/4` handoff instead of hiding the linear-window `xi=8`.
- The dense `0.109` low-weight window now uses the polynomial-free global
  outer envelope `A_h <= 3^h`, avoiding the loose `n^{O(1)}` prefactor.
- The scalar dense inner section now includes a dependency map identifying the
  ON/OFF lemmas as local machinery and pointing to the packaged envelopes
  consumed by the integration theorem.
- The dense integration handoff now uses an explicit checked endpoint
  `eta_1=0.99` in the linear-window criterion and treats the top endpoint
  `h=n` separately instead of applying the outer linear exponent at `eta=1`.
- The full-split finite certificate now states the block-time convention
  (`T=B-i+1` from the first active block), the `late_blocks=5949` split, and
  the unconstrained terminal-state convention in both paper prose and manifest
  metadata.
- The generic first-moment theorem and the block-outer interface now bound the
  joint event that the inner restriction is noninjective or the image has
  minimum distance at most the target.  This repairs the former implicit
  assumption that every square, length-preserving inner map is invertible.
- The full-split one-step turnoff law now states its probability space:
  conditional on block occupancies and preceding transition weights, the
  input and state supports are independent uniform subsets of their specified
  sizes.
- The finite theorem is factored through named complete `h<=500` and `h>=501`
  certificate lemmas; the theorem itself only combines their exact dyadic
  thresholds and concludes the certified dimension and distance.
- The finite certificate states its audit layers explicitly.  The legacy
  manifest total still re-sums floating logarithmic rows within stored
  tolerances, while independent exact/rational and outward-rounded artifacts
  now cover both the complete `h<=500` family and every `h>=501` post-prefix
  family.
- The finite manifest now gives explicit `PASS` status fields for the exact
  RM support check, the prefix interior audit, and the prefix placement-ratio
  certificate, rather than leaving those checks implicit in successful
  execution.
- The finite certificate prose now distinguishes the default manifest verifier
  from full regeneration: the default pass re-sums committed row artifacts and
  manifest interval constants, while the high/complement/early/late interval
  families have opt-in recomputation flags.
- A full opt-in interval recomputation audit is now tracked at
  `scripts/fullsplit_interval_recompute_audit.md`.  It recomputes all 45
  high, complement-high, early-postprefix, early-accelerated, and late
  postprefix rows and records `PASS` under upper-bound semantics; no row
  recomputed above its stored manifest bound, and several high-interval rows
  recomputed to slightly safer values.
- The dominant RM/EBCH prefix mechanism is now isolated in the manifest:
  after exact RM reweighting the \(h=32,r=1,\mathrm{gap}=1..4000\) row is the
  peak, the full \(32\le h\le500,e\le8\) prefix family is only about
  `0.002423` bits above that peak, and the companion ridge-shape check records
  the first-ridge ratio thresholds.
- The same peak-to-family inflation is now a thresholded verifier check:
  total family inflation, total non-peak remainder, \(h=32\)-slice inflation,
  \(h=32\)-slice remainder, and \(h>32\) separation all have explicit gates.
- The dominant peak row arithmetic is now recomputed as its own manifest check:
  exact \(A_{32}=4096\,A^{RM}_{32}\), exact first-gap placement
  numerator/denominator, and the checked `T=5949,H=31` inner knot table value.
- The dominant peak now also has an independent rational upper-bound audit.
  It constructs the EBCH split law with exact fractions, proves the finite
  termination atom is below the rational cap `2^-63` at the monotone endpoint,
  uses the rational Chernoff pole `2333/2373`, rounds the survival probability
  upward to an explicit 96-bit dyadic, bounds the `e=1..8` termination terms by
  `sum_e binom(H+1,e)2^(-63e)`, and checks the complete RM/placement/inner peak against
  `49/2^43` by exact cross multiplication.  The resulting displayed upper
  exponent is about `-37.385767`.  This remains a useful independent peak
  cross-check; the complete low-weight family is now covered by the aggregate
  rational certificate below.
- The dominant `T=5949,H=31` inner knot is now regenerated inside the verifier
  from the EBCH spectrum and the `e<=1` full-split episode formula under the
  stored effective-turnoff convention, then compared against the knot table.
- The ultra-late prefix `T<5949` is now a thresholded verifier check using
  exact RM direct-sum outer coefficients; the total, peak, split, above-split
  mass, and above-split gap all have explicit gates.
- The same ultra-late prefix now has a fully rational aggregate check: the
  verifier sums all 114 supported terms
  `A_h binom(64*5949,h)/binom(2^21,h)` as exact fractions and proves the total
  is below `2^-41` by integer cross multiplication.  Its displayed exact-sum
  exponent is about `-41.113442`; logarithms are used only for reporting.
- The complete `h<=500` first-active family now has an independent
  exact/outward-rounded certificate.  It covers all six gap buckets, keeps
  `e=0..8` in the three prefix buckets and `e=0..16` in the three early
  buckets, and covers the remaining `e>=9`/`e>=17` tails.  Exact RM counts,
  exact placement sums, a 1024-bit outward dyadic inner envelope, and exact
  final cross multiplication give a rational upper bound whose diagnostic
  logarithm is about `-37.276548513006`.  The theorem-safe threshold is the
  exact rational `6793/2^50`; the integer inequality
  `6793^100 <= 2^1273` proves it is at most `2^-37.27` without evaluating a
  transcendental logarithm.
- The complete `h>=501` post-prefix aggregate now has an independent
  rational/outward-rounded certificate.  Exact rational arithmetic handles the
  EBCH split law, MGFs, turnoff envelope, adjacent-ratio peak locations, and 56
  endpoint-in-`T` monotonicity checks.  Decimal logarithms are evaluated with
  directed outward rounding at precision 100.  The diagnostic aggregate upper
  logarithm is `-183.800029543464`, and the theorem-facing gate is `2^-180`.
  The default verifier validates the SHA-256-protected artifact; full
  regeneration is available with
  `python scripts\certify_fullsplit_postprefix_rational.py --recompute-artifact`.
- Both rational certificate paths structurally validate the committed RM and
  EBCH weight tables and record their SHA-256 fingerprints.  The post-prefix
  artifact is rejected if either source table changes.  These checks bind and
  internally validate the inputs; they do not independently derive either
  spectrum from a code definition.  In particular, the EBCH table's
  mathematical provenance remains a declared input.
- The two independent bounds are now combined exactly:
  `(6793*2^130+1)/2^180`.  The default verifier checks
  `(6793*2^130+1)^100 <= 2^14273`, proving that the complete `1<=h<=N`
  rational/outward first moment is at most `2^-37.27` with integer arithmetic.
- The RM/EBCH finite result is stated with the combined rational/outward first
  moment as its theorem-facing bound; the sharper checked-log ledger remains a
  diagnostic, and BCH rows remain separated as projections.
- The current LaTeX log is clean after two passes: no warnings, overfulls,
  undefined references, multiply-defined labels, or rerun requests were found.

## Explicit TODOs Before Stronger Claims

### P0: Must Be Resolved Before Theorem Upgrades

- Replace BCH/random-like projection spectra with exact spectra or rigorous
  low-weight envelopes.
- Recompute the `h>2000` tail under any new BCH-like outer model.
- Keep the BCH projection rows visibly heuristic until the two items above are
  done.

### P1: Proof-Hardening And Audit Items

- For a fully formal audit packet, store raw regenerated interval artifacts or
  make the canonical verifier require the full opt-in recomputation pass under
  outward-rounded interval arithmetic.
- Replace the sampled-grid dense `0.109` exponent-gap check by an interval
  certificate if a fully formal computer-assisted proof is needed.
- Further compress the scalar dense inner ON/OFF proof chain if the main paper
  remains too long; the current narrative now marks which lemmas are local
  scaffolding versus consumed interfaces.
- Keep the bare `verify_dense_claims.py` default from being confused with the
  manuscript theorem: the current theorem check is explicitly
  `--delta 0.109`.

### P2: Presentation And Future-Construction Items

- Shorten dense outer proof machinery only if the main paper remains too long.
- Build a cleaner outer-spectrum comparison note for sliding banded, random
  block, RM, BCH-like, and possible sliding-structured outers.
- Preserve `outerExpandAcc.tex` and `innerSparse.tex` as separate construction
  lines until they are either revived or intentionally retired.

## Verification Gates

Use these as the basic post-edit checks:

```powershell
python scripts\verify_dense_claims.py --delta 0.109
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
python scripts\compare_outer_modes_fullsplit.py --delta 0.09
pdflatex -interaction=nonstopmode main_permConv.tex
pdflatex -interaction=nonstopmode main_permConv.tex
```

Do not run multiple long verification or benchmark commands at the same time.
