# Proof Status And High-Risk TODOs

Generated: 2026-06-13

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
- checked first moment: `E[Z_d] <= 2^-37.278528`.

Current verifier:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
```

Latest checked status: verifier passes and reports
`current_checked_ledger_total_log2,-37.278528`.

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
- The RM/EBCH finite result is stated as a manifest-backed checked first
  moment, with BCH rows separated as projections.
- The current LaTeX log is clean after two passes: no warnings, overfulls,
  undefined references, multiply-defined labels, or rerun requests were found.

## Explicit TODOs Before Stronger Claims

### P0: Must Be Resolved Before Theorem Upgrades

- Replace BCH/random-like projection spectra with exact spectra or rigorous
  low-weight envelopes.
- Recompute the `h>2000` tail under any new BCH-like outer model.
- Keep the BCH projection rows visibly heuristic until the two items above are
  done.
- Clarify the full-split construction and boundary convention enough that the
  paper, manifest, and scripts visibly evaluate the same object.
- Decide whether the finite exact-support prefix rows remain accepted numerical
  certificates or are replaced by analytic monotonicity/support lemmas.

### P1: Proof-Hardening And Audit Items

- Harden the most important finite certificate rows with interval, rational,
  or exact-integer safeguards.
- Audit the dense integration high-weight and linear-window handoff one more
  time after the recent tiny-window and low-window fixes.
- Compress the scalar dense inner ON/OFF proof chain so the narrative exposes
  only the lemmas consumed by the integration theorem.
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
