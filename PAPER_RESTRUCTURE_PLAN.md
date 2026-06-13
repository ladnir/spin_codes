# Paper Restructure Plan

Generated: 2026-06-12

Purpose: one working source of truth for the paper cleanup. This consolidates
the result lock, paper inventory, revised paper flow, preservation policy, and
the current implementation status of the split.

Implementation checkpoint: the first real split and compact structured-section
rewrite have now been applied.  The compiled main spine is
`intro.tex`, `prelim.tex`, `framework.tex`, `innerAcc.tex`, `outerDense.tex`,
`innerDenseScalar.tex`, `integrationDense.tex`, `localCodeOuter.tex`,
`localCodeInner.tex`, `localCodeCertificate.tex`, and
`localCodeProjections.tex`.  The original large source files are preserved, and
the demoted working material has been copied under `explorations/`.  In
particular, the long full-split finite-ledger development was preserved in
`explorations/fullsplit_certificate_development.tex`, while the structured
construction now has separate compiled files for the local-code outer
interface, full-split inner interface, RM/EBCH certificate, and BCH/local-code
projection path.
The follow-up cleanup pass removed visible research-log residue from the
compiled spine: the accumulator TODOs now live in
`explorations/accumulator_warmup_notes.tex`, the dense integration section is
named around the random sliding dense construction, and the local-code section
explicitly separates the checked RM/EBCH certificate from spectrum-model
projections.
The scalar dense inner has also been tightened: positional refinements that are
not consumed by the analytic theorem now live in
`explorations/scalar_dense_inner_refinements.tex`, while
`innerDenseScalar.tex` keeps the ON/OFF lemmas, run-tail exponent,
forced-termination contraction, and packaged inner envelopes.
The random sliding outer has likewise been tightened: fixed-tap tiny-weight
facts now live in `explorations/random_sliding_outer_low_weight_notes.tex`, and
`outerDense.tex` keeps only the theorem-facing span law and spectrum envelopes.
The structured local-code section now treats the finite RM/EBCH result as a
certificate theorem: the main spine keeps the certificate inputs, prefix and
post-prefix totals, and final first-moment theorem, while detailed post-prefix
row tables remain in the full-split exploration/audit notes.
The structured certificate has now been tightened further: the row-level
RM-prefix/post-prefix wrapper statements were collapsed into a manifest-backed
checked ledger summary in `localCodeCertificate.tex`, while the former wrappers
remain in `explorations/fullsplit_certificate_development.tex`.
The compiled main spine has also had a language-audit cleanup: cleanup-era
phrases such as "preserved in notes", "former wrappers", and visible TODO
configuration were removed or rewritten as paper-facing references, while the
BCH spectrum-model subsection still explicitly marks heuristic projections as
non-theorem material.
A final structural pass moved the accumulator warmup immediately after the
framework in the compiled paper and added the missing combined theorem:
random sliding dense outer plus accumulator inner gives linear distance for
sufficiently small constants.  Thus the warmup is now a real proof example, not
just an isolated accumulator enumerator.
The scalar dense inner proof has also started its proof-polish pass: the
run-tail exponent proof in the main spine was compressed, its derivative sign
was corrected, and the expanded Stirling derivation was preserved in
`explorations/scalar_dense_inner_refinements.tex`.
The scalar dense inner was tightened again: the unused polynomial-prefactor
geometric envelope, the ON-to-OFF explanatory remark, and the parameter-use
summary were demoted to `explorations/scalar_dense_inner_refinements.tex`;
`innerDenseScalar.tex` keeps the constant-factor geometric envelope that is
actually consumed by the dense+dense theorem.
The same scalar dense pass now demotes unused fixed-tap finite-diagnostic
bookkeeping from the compiled spine: exact final-survivor averaging, episode
budgeting, multi-candidate termination, and isolated-pair candidate counts live
in the scalar dense-inner refinement note instead of `innerDenseScalar.tex`.
First-start-aware no-OFF tails, zero-gap survival, and the single-run base
factor have also been demoted to the scalar dense-inner refinement note because
they are finite-diagnostic/calibration tools rather than dependencies of the
current dense+dense theorem.
The final unused scalar dense declaration in the compiled path,
`OffAfter(T)`, has been removed, and the original monolithic `innerDense.tex`
and `integration.tex` files now carry explicit archival/provenance headers.
The dense integration spine has also been tightened: the explicit admissible
baseline region, one concrete baseline choice, and memory-scaling interpretation
were moved to `explorations/scalar_dense_finite_diagnostics.tex`, while
`integrationDense.tex` keeps the qualitative dense+dense theorem, the
linear-weight criterion, and the concrete \(0.109\) theorem.
The integration spine has been tightened again: the duplicate local
first-moment theorem and explanatory role/criterion paragraphs were moved to
`explorations/scalar_dense_finite_diagnostics.tex`; `integrationDense.tex`
now cites the framework theorem directly.
The structured projection subsection has now been tightened and physically
split as well: `localCodeProjections.tex` keeps only the theorem/projection
status boundary and a compact qualitative summary of BCH-like outer
projections, while the detailed BCH256/BCH512 table and low-weight inflation
stress tests live in `explorations/outer_spectrum_comparison.tex`.
The random sliding outer section has been renamed and tightened:
`outerDense.tex` is now headed as the random sliding dense outer, and duplicated
orientation/summary prose was moved to
`explorations/random_sliding_outer_low_weight_notes.tex`.
The scalar dense inner was tightened again: remaining compiled commentary about
non-main refinements, linear OFF-budget directions, and envelope interpretation
was moved out of `innerDenseScalar.tex`; the useful interpretation note now
lives in `explorations/scalar_dense_inner_refinements.tex`.
The framework section was also tightened: optional IOWE viewpoint and the
piecewise/single-envelope bookkeeping corollaries now live in
`explorations/framework_flexible_corollaries.tex`, while the compiled
`framework.tex` keeps the first-moment theorem and basic envelope corollary.
The random sliding outer section now keeps the concise fixed-tap band-cluster
low-weight law because `integrationDense.tex` consumes it in the tiny-window
late-placement bound; the exploration note keeps only auxiliary diagnostics.

## Summary

The cleaned paper should become a theorem spine plus preserved research notes.
The main paper keeps the active proof lines:

1. A modular serial-concatenation framework.
2. A simple accumulator warmup.
3. The random sliding dense construction with the asymptotic theorem.
4. A local-code independent structured construction.
5. The proved RM/EBCH finite instantiation.
6. BCH and better local codes as heuristic upgrade paths until spectra are
   certified.

Everything else should be preserved, but moved out of the main paper flow once
the split is reviewed.

## Main Paper Flow

### 1. Intro / Preliminaries / Framework

Keep the general setup and first-moment framework as the front of the paper.

Main source files:

- `intro.tex`
- `prelim.tex`
- `framework.tex`

Target content:

- serial concatenation setup;
- outer spectrum interface;
- inner slice-to-tail interface;
- first-moment minimum-distance reduction;
- short explanation that later sections instantiate the interfaces.

### 2. Accumulator Warmup

Promote the accumulator inner as a clean framework example.

Main source file:

- `innerAcc.tex`

Locked result:

- The binary accumulator has exact input-output enumerator
  ```text
  A^{Acc}_{w,h}
  = binom(h-1, ceil(w/2)-1) binom(n-h, w-ceil(w/2))
  ```
  over the feasible range.
- For `delta < 1/4` and `w <= n/2`,
  ```text
  p_w(delta) <= (4e delta)^{ceil(w/2)}.
  ```
- With a log-memory outer spectrum envelope, this gives a simple
  linear-distance theorem for sufficiently small `delta`; in the compiled
  paper this is stated as the random sliding dense outer plus accumulator
  inner warmup theorem.

Main-paper role:

- Demonstrates the framework transparently.
- Shows why recursive inners help.
- Explicitly note that constants are weak and this is not the performance
  construction.

### 3. Random Sliding Dense Construction

This is the original dense-band outer x dense recursive inner theory line.

Main source files:

- `outerDense.tex`
- scalar/theorem-facing parts of `innerDense.tex`
- theorem-facing parts of `integration.tex`

Target content:

- random/log-memory banded outer;
- dense recursive/convolution inner;
- dense+dense first-moment integration;
- asymptotic linear-distance theorem;
- concrete `0.109` theorem.

Locked result:

- `thm:dense-dense-integration`: systematic dense/banded outer plus dense
  recursive inner has linear minimum distance with high probability.
- `thm:dense-dense-explicit-constant`: concrete `0.109` asymptotic theorem.
- Verification command matching the current theorem:
  ```powershell
  python scripts\verify_dense_claims.py --delta 0.109
  ```
  reports sampled-grid verification with worst gap `-0.003107731647`.

Notes:

- The bare verifier default currently checks `delta=0.12` and fails. Treat
  that as a stale/harsher default, not as the manuscript theorem.
- Scalar finite scans belong in exploration notes, not the main theorem flow.

### 4. Structured Local-Code Construction

This section should be local-code independent first. RM and EBCH are
instantiations, not definitions of the construction.

Target content:

- local-code outer interface:
  ```text
  C_out = C_loc^{oplus B},     W_out(z) = W_loc(z)^B
  ```
- local-code full-split recursive inner interface:
  - choose a rate-1/2 local inner code of length `2b`, dimension `b`;
  - scramble a nonzero state/input block to a uniform nonzero local codeword;
  - randomly split `2b` coordinates into `b` output coordinates and `b` next
    state coordinates;
  - analyze through the local weight spectrum and split law.
- finite certificate template:
  ```text
  sum_h A_h^out p_h^in(delta)
  ```
  with outer supplied by a spectrum/envelope and inner supplied by an episode
  ledger/certificate.

Main source files:

- local-code interface material from `integration.tex`;
- full-split definition and certificate-interface material from `innerDense.tex`.
- implemented split spine: `localCodeOuter.tex`, `localCodeInner.tex`,
  `localCodeCertificate.tex`, and `localCodeProjections.tex`.
- implemented demotion: detailed modeled BCH projection rows live in
  `explorations/outer_spectrum_comparison.tex`, with only a short status note
  remaining in the compiled paper.

Main-paper role:

- Explain what data a local code must provide.
- Make clear that distance floor alone is insufficient; low-weight spectrum
  shape matters.
- Prepare the RM/EBCH instantiation.

### 5. Proved Instantiation: RM Outer + EBCH Full-Split Inner

This is the current main finite construction.

Main source files:

- certificate-facing full-split material in `innerDense.tex`;
- `AUDIT.md` for audit provenance;
- `scripts/fullsplit_finite_ledger_manifest.json` for machine-readable ledger
  state.

Locked parameters:

- Outer local code: `RM(4,9) [512,256,32]`.
- Number of outer blocks: `4096`.
- Inner local code: EBCH `[128,64,22]`.
- Full-split inner block size: `b=64`.
- Global length: `N=2^21=2097152`.
- Target: `delta=.09`, `d=188743`.

Locked certificate:

```text
E[Z_d] <= 2^-37.278528
Pr[d_min <= d] <= 2^-37.278528
d_min >= 188744 with positive probability
188744 / 2^21 > .09
```

Canonical verification:

```powershell
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
```

The verifier reports:

- `current_checked_ledger_total_log2,-37.278528`
- `current_checked_ledger_margin_bits,37.278528`
- dominant row `outer_weight=32, first_r=1, gap=1..4000`.

Main-paper role:

- State the finite theorem cleanly.
- Keep proof details enough to show every row family is covered.
- Move long command logs and row-family development history to audit/exploration
  notes.

### 6. BCH And Better Local-Code Instantiations

This is an upgrade/projection section, not theorem text yet.

Main source artifacts:

- `scripts/compare_outer_modes_fullsplit.py`
- `scripts/outer_mode_comparison_delta009.csv`
- `scripts/outer_mode_comparison_delta009.md`

Locked projection rows at `N=2^21`, `delta=.09`, `d=188743`:

| mode | status | local outer | blocks | log2 mu | margin bits | dominant h |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `rm512_exact` | proved | RM `[512,256,32]` | 4096 | -37.278528 | 37.278528 | 32 |
| `bch256_heuristic_plus0` | heuristic | BCH-like `[256,128,38]` | 8192 | -37.089134 | 37.089134 | 59 |
| `bch512_heuristic_plus0` | heuristic | extended BCH-like `[512,256,>=62]` | 4096 | -68.171304 | 68.171304 | 119 |

Main-paper role:

- Explain the likely upgrade path.
- Emphasize that BCH rows require exact spectra or rigorous low-weight
  envelopes before becoming theorem-facing.
- Explain the structural lesson: direct-sum block outers and sliding banded
  outers can have very different low-weight enumerator shapes.

### 7. Preserved Explorations

These ideas should not be deleted. They should be moved out of the main paper
flow after review.

Proposed destinations:

- `explorations/scalar_dense_finite_diagnostics.tex`
  - sigma scans;
  - isolated-slice and early-pair material;
  - old scalar finite proof targets.
- `explorations/fullsplit_certificate_development.tex`
  - row-family development history;
  - detailed finite ledger derivations not needed in main flow.
- `explorations/block_recursive_bch_inner.tex`
  - block-recursive BCH inner attempts;
  - column-profile, turnoff, and permutation diagnostics.
- `explorations/bch_outer_spectra.tex`
  - BCH spectrum search;
  - MacWilliams, affine-coset, public-table, and small-ladder material.
- `explorations/outer_spectrum_comparison.tex`
  - random banded versus random block versus RM/BCH local-code comparison.
- `explorations/outer_expand_accumulator.tex`
  - preserved material from `outerExpandAcc.tex`.
- `explorations/inner_sparse.tex`
  - preserved material from `innerSparse.tex`.

## Locked Results Summary

### Accumulator Warmup

Status: proved framework warmup.

- Exact accumulator enumerator is in `innerAcc.tex`.
- Low-output contraction:
  ```text
  p_w(delta) <= (4e delta)^{ceil(w/2)}
  ```
- Use as a clean example, not as the current best construction.

### Dense Sliding Asymptotic Theorem

Status: theorem-facing.

- Random/log-memory banded outer plus dense recursive inner gives linear
  distance.
- Concrete sampled theorem target: `0.109`.
- The `0.109` verifier passes on the sampled grid.

### Finite RM/EBCH Certificate

Status: current main finite certificate.

- `N=2^21`, `delta=.09`, `d=188743`.
- `E[Z_d] <= 2^-37.278528`.
- Dominant mode: outer `h=32`, first active occupancy `r=1`, near late window
  `gap=1..4000`.

### Scalar Dense Sigma-80 Diagnostic

Status: diagnostic/superseded for the main finite certificate.

- Old fixed-tap random-banded outer x fixed-tap dense recursive inner scan:
  - `delta=.106`, `sigma=80`: about `40.993` bits through `h<=30`, peak
    `h=16`;
  - `delta=.09`, `sigma=80`: about `44.651` bits through `h<=30`, peak
    `h=14`;
  - larger prefix at `delta=.106`, `sigma=80`:
    `log2 mu_{<=500} <= -40.991450`.

### BCH/RM Projections

Status: heuristic except the RM row.

- BCH256 no-inflation projection is near RM but sensitive to low-weight
  inflation.
- BCH512 no-inflation projection is stronger and more robust.
- No BCH projection becomes theorem text until the local spectrum/envelope and
  matching tail recomputation are certified.

## Restructure Map

| source | action |
| --- | --- |
| `intro.tex` | Keep main; rewrite later after section split stabilizes. |
| `prelim.tex` | Keep main. |
| `framework.tex` | Keep main as first-moment framework. |
| `innerAcc.tex` | Promote as accumulator warmup. |
| `outerDense.tex` | Keep dense outer theorem/interface. |
| `innerDense.tex` | Preserve original source. Scalar dense theorem extracted to `innerDenseScalar.tex`; full-split development copied to `explorations/fullsplit_certificate_development.tex`; block-recursive BCH copied to `explorations/block_recursive_bch_inner.tex`. |
| `integration.tex` | Preserve original source. Analytic dense+dense theorem extracted to `integrationDense.tex`; scalar finite scans copied to `explorations/scalar_dense_finite_diagnostics.tex`; BCH spectra/history copied to `explorations/bch_outer_spectra.tex`. |
| `outerExpandAcc.tex` | Preserve outside main paper via `explorations/outer_expand_accumulator.tex`. |
| `innerSparse.tex` | Preserve outside main paper via `explorations/inner_sparse.tex`. |
| `AUDIT.md` | Keep as audit packet and source of truth for verification details, not main prose. |
| `RESULTS_LOCK.md` | Keep as current result snapshot until restructure is complete. |
| `PAPER_INVENTORY.md` | Keep as provenance for this plan. |

## Phase 3 Guardrails

- Do not delete exploratory material on the first restructure pass.
- Do not move labels gratuitously; keep stable labels where possible.
- If labels move or are renamed, add a reference map.
- Do not present BCH/random-like rows as theorem text.
- Do not present old scalar finite diagnostics as the current finite theorem.
- Preserve `outerExpandAcc.tex` and `innerSparse.tex` as separate construction
  lines.
- Keep generated scratch files and `__pycache__` out of the paper flow unless
  explicitly named as canonical artifacts.

## Implemented Checkpoint

The first physical split is complete:

1. `explorations/` exists and contains topic notes for the demoted branches.
2. Scalar dense theorem material now compiles from `innerDenseScalar.tex`.
3. Analytic dense+dense integration now compiles from `integrationDense.tex`.
4. The local-code interface, RM/EBCH certificate, and BCH projection note now
   compile from `localCodeOuter.tex`, `localCodeInner.tex`,
   `localCodeCertificate.tex`, and `localCodeProjections.tex`.
5. `main_permConv.tex` uses the cleaned spine and no longer inputs the original
   monolithic `innerDense.tex`, `integration.tex`, `outerExpandAcc.tex`, or
   `innerSparse.tex`.
6. The accumulator warmup no longer exposes its optional TODO list in the main
   paper; those notes live in `explorations/accumulator_warmup_notes.tex`.
7. Scalar dense-inner refinement lemmas not needed by the theorem-facing
   integration now live in `explorations/scalar_dense_inner_refinements.tex`.
8. Random sliding outer tiny-weight diagnostics now live in
   `explorations/random_sliding_outer_low_weight_notes.tex`.
9. Structured local-code post-prefix row details are delegated to the
   full-split exploration/audit notes; the main spine keeps compact audited
   totals and row-cover categories.
10. A main-spine language audit removed compiled `\stan{...}` comments, stale
    working-save-point wording, and the final framework overfull. Remaining
    diagnostic/projection/heuristic language marks explicit non-theorem status
    boundaries.
11. The scalar dense run-tail exponent proof was shortened in the main spine,
    the monotonicity sign was corrected, and the expanded derivation was moved
    to the scalar dense-inner refinement note.
12. Fixed-tap finite-diagnostic bookkeeping lemmas not consumed by the
    dense+dense theorem were moved from `innerDenseScalar.tex` to
    `explorations/scalar_dense_inner_refinements.tex`.
13. First-start-aware and single-run scalar dense diagnostics were also moved
    out of the compiled theorem spine; `innerDenseScalar.tex` now proceeds from
    the three-term tail directly to forced-termination contraction.
14. The stale `OffAfter(T)` declaration was removed from the compiled scalar
    dense spine, and the original monolithic `innerDense.tex` and
    `integration.tex` sources were marked archival.
15. Baseline dense+dense parameter-bookkeeping corollaries were demoted from
    `integrationDense.tex` to `explorations/scalar_dense_finite_diagnostics.tex`.
16. Scalar dense-inner commentary paragraphs not used by the proof were removed
    from `innerDenseScalar.tex` and preserved as exploration context where useful.
17. Framework optional corollaries were demoted from `framework.tex` to
    `explorations/framework_flexible_corollaries.tex`, leaving only the core
    first-moment reduction and reusable envelope corollary in the compiled spine.
18. The fixed-tap band-cluster low-weight law was promoted back into
    `outerDense.tex`, and the dense integration proof now references it directly
    instead of relying on cleanup/provenance prose.

Remaining cleanup is proof-facing polish, not the initial physical split:

- continue compressing the scalar dense inner proof presentation, especially
  the long ON/OFF lemma chain;
- shorten dense outer proof machinery further only if the main paper remains
  too long;
- harden finite-certificate manifests for external audit;
- replace BCH projections with exact spectra or rigorous envelopes before any
  theorem upgrade.

## Verification Plan

For full restructure verification, run sequentially:

```powershell
python scripts\verify_dense_claims.py --delta 0.109
python scripts\verify_fullsplit_finite_ledger.py --write-manifest-json scripts\fullsplit_finite_ledger_manifest.json
python scripts\compare_outer_modes_fullsplit.py --delta 0.09
pdflatex -interaction=nonstopmode main_permConv.tex
pdflatex -interaction=nonstopmode main_permConv.tex
```

Expected caveats:

- The bare `verify_dense_claims.py` default checks `delta=0.12` and is not the
  current theorem check.
- Do not run multiple long verification/benchmark commands at the same time.

## Assumptions

- This plan remains the working cleanup source of truth.
- `RESULTS_LOCK.md` and `PAPER_INVENTORY.md` remain present for provenance.
- The current main finite construction is RM outer plus EBCH full-split inner.
- BCH and other stronger local-code rows remain projections until certified.
