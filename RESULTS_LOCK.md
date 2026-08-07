# Results Lock

Generated: 2026-06-12

Purpose: lock the current dense+dense research state before paper cleanup or
prose movement. This file is a status snapshot, not a new proof.

## Main Result Lines

### Accumulator warmup line

Status: proved framework warmup, useful for exposition rather than concrete
performance.

Primary locations:

- `innerAcc.tex`: exact accumulator input-output enumerator and contraction.
- `outerDense.tex`: random sliding dense outer generating-function envelope.
- `framework.tex`: generic inner-interface first-moment consumption.

Locked claims:

- The basic rate-1 accumulator has exact enumerator
  `A^{Acc}_{w,h} = binom(h-1, ceil(w/2)-1) binom(n-h, w-ceil(w/2))`
  over the feasible range.
- For `delta < 1/4` and `w <= n/2`,
  `p_w(delta) <= (4e delta)^{ceil(w/2)}`.
- Combined with the random sliding dense outer generating-function envelope,
  this gives a clean linear-distance warmup theorem for sufficiently small
  `delta`, but with a weak constant compared with the dense recursive inner.

### Analytic dense-band outer x dense recursive inner

Status: proved asymptotic theorem in the current manuscript.

Primary locations:

- `framework.tex`: serial-concatenation first-moment framework.
- `randomDenseConstruction.tex`: compiled parent section for the random
  sliding dense construction.
- `outerDense.tex`: systematic dense/banded outer spectrum interface.
- `innerDenseScalar.tex`: scalar dense recursive inner interface in the current
  compiled spine, nested under `randomDenseConstruction.tex`.
- `integrationDense.tex`: dense+dense integration theorem in the current
  compiled spine, nested under `randomDenseConstruction.tex`.
- `innerDense.tex` and `integration.tex`: preserved original sources.

Locked claims:

- The systematic dense/banded outer plus dense recursive inner has linear
  minimum distance with high probability under the parameter conditions in
  `thm:dense-dense-integration`.
- The concrete sampled-grid verification supports the manuscript's current
  `0.109` asymptotic constant theorem, `thm:dense-dense-explicit-constant`.
- The command `python scripts\verify_dense_claims.py --delta 0.109` reports
  `STATUS: VERIFIED ON THE SAMPLED GRID` with worst gap
  `-0.003107731647` near `eta=0.499993218813`.

Important command note:

- The requested bare command `python scripts\verify_dense_claims.py` currently
  uses default `delta=0.12` and fails with worst gap `0.029360865155`.
  Treat this as a stale/harsher default, not as a failure of the current
  `0.109` manuscript claim.

### Finite RM outer x full-split EBCH inner

Status: current main finite construction and checked certificate.

Primary locations:

- `localCodeStructured.tex`: compiled parent section for the structured
  local-code construction.
- `localCodeOuter.tex`: direct-sum outer first-moment interface.
- `localCodeInner.tex`: full-split local-code recursive inner definition and
  one-step split law.
- `localCodeCertificate.tex`, especially:
  - the manifest-backed checked ledger summary in the RM/EBCH certificate
    section;
  - `thm:fullsplit-rm-finite-certificate-009`.
- `explorations/fullsplit_certificate_development.tex`: detailed development
  notes copied from the original dense-inner workspace, including the former
  row-level wrappers:
  - `lem:rm-prefix-support-split`
  - `cor:fullsplit-small-prefix-all-first-active`
  - `cor:fullsplit-postprefix-checked-ledger`
  - `cor:fullsplit-current-finite-checkpoint`
- `innerDense.tex`: preserved original source containing the pre-cleanup
  full-split ledger material.
- `AUDIT.md`: full audit packet and command log.
- `scripts/fullsplit_finite_ledger_manifest.json`: machine-readable checked
  ledger snapshot.

Locked parameters:

- Outer: direct sum of `4096` copies of `RM(4,9)`.
- Local outer code: binary `[512,256,32]`.
- Outer length: `N=2^21=2097152`.
- Message length: `K=2^20=1048576`.
- Inner: full-codeword random-split dense recursive inner.
- Inner local spectrum: EBCH `[128,64,22]` table from
  `scripts/EBCH128_64.wd`.
- Inner block/state/output half-size: `b=64`.
- Target relative distance: `delta=.09`.
- Target integer distance: `d=floor(.09*N)=188743`.

Locked certificate:

```text
E[Z_d] <= 2^-37.2785
Pr[d_min <= d] <= 2^-37.2785
d_min >= 188744 with positive probability
188744 / 2^21 > .09
```

Verification status:

- `python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json`
  passed.
- The verifier reported:
  - `current_checked_ledger_total_log2,-37.278528`
  - unrounded stored total `-37.278527626...`
  - theorem-safe outward-rounded exponent `-37.2785`
  - `current_checked_ledger_margin_bits,37.278528`
  - `h501_plus_checked_rows_log2,-182.259739`
  - `prefix_32_500_e_le8_dominant_term_log2,-37.385767`
  - `dominant_peak_rational_upper_bound_status,PASS`
  - exact rational dominant-peak threshold `49/2^43`
  - `late_prefix_exact_rational_status,PASS`
  - exact rational ultra-late-prefix threshold `2^-41`
  - `fullsplit_h500_complete_rational_status,PASS`
  - complete rational/outward `h<=500` diagnostic log `-37.276548513006`
  - exact threshold `6793/2^50 <= 2^-37.27`
  - dominant row `outer_weight=32, first_r=1, gap=1..4000`

Canonical artifacts:

- `scripts/fullsplit_finite_ledger_manifest.json`
- `scripts/rm512_256_spectrum.csv`
- `scripts/EBCH128_64.wd`
- `scripts/certify_fullsplit_h500_rational.py`
- `scripts/fullsplit_h500_gap_sums_exact.json`
- `scripts/fullsplit_h500_inner_bounds_dyadic.json`
- `scripts/fullsplit_piecewise_h32_500_csv.csv`
- `scripts/fullsplit_piecewise_early_h32_500_e16_uniformsurv.csv`
- `scripts/fullsplit_piecewise_h501_2000_eall_hsummary.csv`

## Scalar Dense Finite Diagnostics

Status: diagnostic and partially superseded by the RM/full-split finite
certificate, but still important for understanding the older dense-band x
dense-conv line.

Locked facts:

- The old fixed-tap random-banded outer x fixed-tap dense recursive inner scan
  at `k=2^20` showed that `sigma=80` crosses about `40` bits of margin.
- From `integration.tex` and saved scan files:
  - at `delta=.106`, `sigma=80`, `-log2 mu_{<=30}=40.993`, peak `h=16`;
  - at `delta=.09`, `sigma=80`, `-log2 mu_{<=30}=44.651`, peak `h=14`;
  - the larger prefix check at `delta=.106`, `sigma=80` gave
    `log2 mu_{<=500} <= -40.991450`, peak `h=16`.

Interpretation:

- These are not the current finite main theorem.
- They remain useful as a comparison point for sliding/banded outer low-weight
  enumerators.
- Do not confuse `sigma=80` in this scalar construction with the local
  dimension `256` in block outers.

## Outer-Mode Projections

Status: design diagnostics only. BCH rows are heuristic until an exact local
spectrum or rigorous low-weight envelope replaces the model.

Command:

```powershell
python scripts\compare_outer_modes_fullsplit.py --delta 0.09
```

Artifacts:

- `scripts/outer_mode_comparison_delta009.csv`
- `scripts/outer_mode_comparison_delta009.md`

Locked projection rows at `N=2^21`, `delta=.09`, `d=188743`:

| mode | status | local outer | blocks | log2 mu | margin bits | dominant h |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `rm512_exact` | proved | RM `[512,256,32]` | 4096 | -37.278528 | 37.278528 | 32 |
| `bch256_heuristic_plus0` | heuristic | BCH-like `[256,128,38]` | 8192 | -37.089134 | 37.089134 | 59 |
| `bch512_heuristic_plus0` | heuristic | extended BCH-like `[512,256,>=62]` | 4096 | -68.171304 | 68.171304 | 119 |

Sensitivity warning:

- BCH256 is fragile under low-weight inflation: `+10` bits leaves
  `27.167072` margin, while `+20` bits makes the `501..2000` projection
  positive.
- BCH512 is more robust in this projection: `+10,+20,+40` give margins
  `65.802609`, `56.111333`, and `36.111659`.
- The `>2000` tail is not recomputed for BCH projections in this driver.

## Conceptual Conclusions

Current working conclusions to preserve during cleanup:

- The first moment depends on the low-weight enumerator shape, not just the
  distance floor.
- Direct-sum block outers have thick low-weight spectra because many input
  bits share the same local output range and local low-weight words are reused
  across blocks.
- Sliding/banded outers have a different run/span geometry; their expected
  extreme low-weight enumerators can be much thinner even at comparable
  locality.
- The current finite RM/full-split bottleneck is low/intermediate outer weight
  plus late placement, dominated by `h=32`, `r=1`, `gap=1..4000`.
- The current BCH512 projection suggests a real constant-factor upgrade could
  exist, but only after the local spectrum or an envelope is certified.

## Proof Debt

Do not promote these to proved claims until resolved:

- Keep `PROOF_STATUS.md` current as the explicit high-risk TODO ledger.
- Replace BCH/random-like projection spectra with exact spectra or rigorous
  low-weight envelopes.
- Recompute the `h>2000` tail under any new BCH-like outer model.
- Continue rational/interval hardening of the `h>=501` interval families.
  The complete `h<=500` prefix-family inflation and all first-active positions
  now have an independent 37.27-bit rational/outward certificate.
- Clarify the full-split construction definition and boundary convention so
  all scripts visibly evaluate the same object.
- Decide whether finite exact-support prefix lemmas should remain finite
  certificates or be replaced by analytic monotonicity/support lemmas.
- Build an explicit outer-spectrum comparison audit for banded, random-block,
  RM, BCH-like, and possible sliding-structured outers.
- Preserve alternate construction ideas, especially `outerExpandAcc.tex` and
  `innerSparse.tex`, as separate notes or appendices rather than deleting them
  from the project.

## Cleanup Flow To Preserve

The intended main-paper spine after inventory review is:

1. Intro, preliminaries, and modular first-moment framework.
2. Accumulator warmup: log-memory outer plus accumulator inner gives linear
   distance with a weak constant.
3. Random sliding dense construction: random/log-memory banded outer plus dense
   recursive inner, including the `0.109` asymptotic theorem.
4. Structured local-code construction: local-code independent outer interface,
   local-code independent full-split recursive inner interface, finite
   certificate template, then RM/EBCH as the first proved instantiation.
5. BCH and better local-code instantiations: projections and spectrum
   requirements, explicitly not theorem-facing yet.
6. Preserved explorations: scalar dense finite diagnostics, block-recursive
   BCH inner, BCH spectra, outer spectrum comparison, expander/accumulator
   outer, and sparse refresh inner.

## Verification Commands Run In This Lock Pass

```powershell
python scripts\verify_dense_claims.py
python scripts\verify_dense_claims.py --delta 0.109
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
python scripts\compare_outer_modes_fullsplit.py --delta 0.09
pdflatex -interaction=nonstopmode main_permConv.tex
pdflatex -interaction=nonstopmode main_permConv.tex
```

Results:

- Bare `verify_dense_claims.py`: failed at default `delta=0.12`; see note
  above.
- `verify_dense_claims.py --delta 0.109`: passed sampled-grid check.
- Full-split finite ledger verifier: passed and wrote the manifest.
- Outer-mode comparison: refreshed reports.
- LaTeX: compiled twice and produced `main_permConv.pdf`; existing undefined
  references and multiply-defined labels were present in the initial lock pass.

## Current Restructure Verification Addendum

After the physical section-wrapper cleanup, the compiled paper spine is:

1. `intro.tex`
2. `prelim.tex`
3. `framework.tex`
4. `innerAcc.tex`
5. `randomDenseConstruction.tex`, which inputs `outerDense.tex`,
   `innerDenseScalar.tex`, and `integrationDense.tex`
6. `localCodeStructured.tex`, which inputs `localCodeOuter.tex`,
   `localCodeInner.tex`, `localCodeCertificate.tex`, and
   `localCodeProjections.tex`

The latest cleanup verification ran:

```powershell
pdflatex -interaction=nonstopmode main_permConv.tex
pdflatex -interaction=nonstopmode main_permConv.tex
python scripts\verify_dense_claims.py --delta 0.109
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
```

The LaTeX log scan for `Warning`, `Overfull`, `Underfull`, `undefined`,
`multiply defined`, and `Rerun` is now clean after the second pass.  The dense
`0.109` verifier still reports worst gap `-0.003107731647`, and the full-split
finite ledger still reports `current_checked_ledger_total_log2,-37.278528`.

## Exclusions

Generated scratch CSV/PNG files, `tmp_*` files, and `__pycache__` files are not
part of this result lock unless explicitly listed above or named in
`AUDIT.md` as a canonical artifact.
