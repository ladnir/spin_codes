# Claude Audit Report

Command used:

```powershell
claude -p --model claude-fable-5 --no-session-persistence --permission-mode bypassPermissions --tools "Read,Grep,Glob" --max-budget-usd 1 --output-format text
```

Scope:

- Read-only audit.
- No shell commands from Claude.
- Files requested: `AUDIT.md`, the relevant full-split sections of `innerDense.tex`, and the core verification scripts.

## Verdict

The certificate is structurally coherent: the script formulas for the e=1 Cauchy envelope, the e>=2 same-lambda multiplier, the selected-gap closed forms, and the sufficient endpoint monotonicity reduction all match the manuscript lemmas, and the four-piece ledger arithmetic in `verify_fullsplit_finite_ledger.py` is a correct log-sum of its inputs. However, the headline `log2 mu_finite <= -34.767174` is a numerical certificate with at least three load-bearing links that are currently assumed rather than checked: the bucket-wide reuse of left-endpoint-T exact grids in the dominant piece, the provenance of the termination atom `p_term`, and the hard-coded interval rows that the verifier merely re-sums. It is correctly described as not yet a polished global theorem.

## Findings

### F1 High: Exact-grid T monotonicity is not certified

The dominant piece uses `csv` mode on `32 <= h <= 500`. The wrapper looks up the inner bound by `H` only and ignores `T`, pairing grids computed at bucket left endpoints `T = 5949, 9949, 13949` with placements throughout each bucket.

This is sound only if the exact `e <= 8` episode sum is nonincreasing in `T`. The current monotonicity script certifies the fixed-pole envelope for the `h >= 2001` interval certificate, not the exact tiny-prefix grids.

Claude classified this as a soundness assumption, not just presentation debt.

### F2 High: Termination atom provenance is not audit-grade

The constant

```text
log2 p_term <= -63.8926492
```

is load-bearing. The lemma says `p_term` should cover both the `q'=0` split event and exact input/state cancellation if charged. The audited scripts expose an optional `--extra-turnoff-log2`, but the CSV provenance does not make clear whether/how cancellation was included.

Resolution note: `scripts/verify_fullsplit_turnoff_atom.py` now derives this atom for the finite-prefix domain
`0 <= H <= 499`, `T >= 5949`. It computes `log2 p0 = -64.00457432724919`,
`log2 eta = -67.63641076470455`, `log2(p0+eta) = -63.89264922047382`, and uses the rounded-up certificate atom
`-63.8926492`. Remaining issue: the same exact-cancellation add-on is not yet certified for larger `H`.

### F3 High: High interval rows are hard-coded in the ledger verifier

`verify_fullsplit_finite_ledger.py` re-sums literal constants for the interval certificate. It does not recompute them or document exact regeneration commands.

This is acceptable for a ledger manifest, but not yet an audit-grade proof artifact.

### F4 Medium-High: Infeasible left endpoint can silently drop mass

In the wrapper, `e01_cauchy_log2` / `eall_ratio_log2` can return `-inf` when the left endpoint `T_min` cannot hold `H`, even though later `T` in the same bucket can be feasible. The current important high-tail runs use `feasible-min`, but the general tool is unsafe in `min` mode for high `H`.

### F5 Medium: Current number is a late-window slice

The `-34.767174` claim is for the audited late-window split, not the entire dense+dense construction unless the other placement regimes are separately covered.

The paper and audit packet should keep this qualifier visible.

### F6 Medium: Product form for per-gap terminations needs a lemma

The exact grids and crude tail use products like `p_term^e`. This needs a stated independence/conditioning lemma for per-block dense scramblers and gap terminations.

### F7 Medium: Floating-point rigor

No interval/rational arithmetic is present. Claude highlighted two small examples:

- near-one handling in `log2_arith_geom_sum`
- the early break in `log2_ht_tail_envelope`

These are probably numerically tiny, but should be bounded explicitly or moved into interval arithmetic.

### F8 Low: Presentation and hygiene

- `innerDense.tex` had a stale `-34.771796` combined-prefix number while the current ledger is `-34.767174`.
- The older `e <= 1` diagnostic numbers coexist with the theorem-facing all-episode row and should be clearly labeled.
- `__pycache__` status noise should be cleaned/ignored.
- The monotonicity intervals and ledger intervals are consistent, but the relationship is not obvious.

## Author Questions

1. Outside the finite-prefix range, how is occupied-block exact cancellation handled in the all-episode lemma?
2. Which command generated the fifteen interval rows, and was `feasible-min` used?
3. What justifies applying the exact `e <= 8` grids at bucket left endpoints across entire buckets?
4. Where is the product bound for multiple terminations proved?
5. Which section covers first-active placements outside the audited late-window split?

## Suggested Next Fixes

1. Certify exact-grid `T` monotonicity for `32 <= h <= 500`, or keep the current full interior finite audit as an explicit finite lemma.
2. Extend the `p_term` exact-cancellation derivation beyond the finite prefix, or split larger-`H` cancellations into a separate bound.
3. Make the interval rows reproducible from documented commands, or add a spot-recompute mode to the ledger verifier.
4. Make `sum_fullsplit_piecewise_certificate.py` assert or use `feasible-min` when endpoint placement is paired with high `H`.
5. Correct stale prose numbers and separate old diagnostics from theorem-facing rows.
6. Add interval/rational arithmetic for the final ledger.
