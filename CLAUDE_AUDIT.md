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

The certificate is structurally coherent: the script formulas for the e=1 Cauchy envelope, the e>=2 same-lambda multiplier, the selected-gap closed forms, and the sufficient endpoint monotonicity reduction all match the manuscript lemmas, and the four-piece ledger arithmetic in `verify_fullsplit_finite_ledger.py` is a correct log-sum of its inputs. However, the headline `log2 mu_finite <= -34.767174` is still a numerical certificate, not a polished global theorem: the tiny-prefix exact-grid control is now checked by a finite interior artifact, and the high-interval manifest rows are command-reproducible, but both still need final theorem packaging or fully regenerated proof artifacts.

## Findings

### F1 Medium: Exact-grid T monotonicity is a finite artifact, not theorem text

The dominant piece uses `csv` mode on `32 <= h <= 500`. The wrapper looks up the inner bound by `H` only and ignores `T`, pairing grids computed at bucket left endpoints `T = 5949, 9949, 13949` with placements throughout each bucket.

This is sound only if the exact `e <= 8` episode sum is nonincreasing in `T`, or if every interior `T` is directly audited. The current fixed-pole monotonicity script certifies the `h >= 2001` interval certificate, not the exact tiny-prefix grids.

Resolution note: `scripts/check_fullsplit_exact_grid_interior.py` now audits every `H=0..499` and every `T` in the three placement buckets. The tracked artifact has `1500` rows and worst positive difference `2.19205276153e-09`, at the left endpoint itself. `verify_fullsplit_finite_ledger.py` checks this artifact by default. The manuscript now records this as Lemma `lem:fullsplit-tiny-prefix-interiorT`. Remaining debt: decide whether to keep this finite lemma or replace it with analytic monotonicity.

### F2 High: Termination atom provenance is not audit-grade

The constant

```text
log2 p_term <= -63.8926492
```

is load-bearing. The lemma says `p_term` should cover both the `q'=0` split event and exact input/state cancellation if charged. The audited scripts expose an optional `--extra-turnoff-log2`, but the CSV provenance does not make clear whether/how cancellation was included.

Resolution note: `scripts/verify_fullsplit_turnoff_atom.py` now derives this atom for the finite-prefix domain
`0 <= H <= 499`, `T >= 5949`. It computes `log2 p0 = -64.00457432724919`,
`log2 eta = -67.63641076470455`, `log2(p0+eta) = -63.89264922047382`, and uses the rounded-up certificate atom
`-63.8926492`. For larger `H`, `scripts/verify_fullsplit_global_turnoff_envelope.py` now proves the Bernstein
envelope `log2 p_term <= -62.4078758`. The high-interval rows carry a checked `+1.484774` bit adjustment for this
wider atom; `scripts/verify_fullsplit_high_interval_adjustment.py` computes the required adjustment as
`1.484773404861` bits, including the tail multiplier, leaving about `5.95e-7` bits of slack.

### F3 Medium: High interval rows are a command-reproducible manifest

`verify_fullsplit_finite_ledger.py` re-sums literal constants for the interval certificate. It now records the exact theorem-facing regeneration command for each high interval and exposes opt-in selected-row recomputation through `--check-high-intervals`, but it still does not recompute every interval by default.

This is acceptable for a ledger manifest, but not yet a fully regenerated audit-grade proof artifact.

### F4 Medium-High: Infeasible left endpoint can silently drop mass

In the wrapper, `e01_cauchy_log2` / `eall_ratio_log2` can return `-inf` when the left endpoint `T_min` cannot hold `H`, even though later `T` in the same bucket can be feasible. The current important high-tail runs use `feasible-min`, but the general tool is unsafe in `min` mode for high `H`.

### F5 Medium: First-active coverage is now explicit, but not global

The `-34.767174` claim is for the audited window `5949 <= T <= 17948`. The verifier also checks the placement-only ultra-late prefix `T < 5949`, `h <= 500`, and now checks the early small-weight region `T > 17948`, `h <= 500` by combining a bucket-uniform survival ledger through `e <= 16` with a crude `e >= 17` tail. The combined all-first-active-position value for `32 <= h <= 500` remains `-34.732187`.

Remaining coverage debt: the early region for `h > 500` and the post-prefix ultra-late cap beyond `h > 500`.

Update: the early high-density branch has a clear bookkeeping hazard. A separable endpoint-placement scan can look
positive because it pairs placement at `T_max` with survival at `T_min`. The new probe
`scripts/probe_fullsplit_early_paired_survival.py` keeps `T` paired; with a deliberately pessimistic `20000`-bit
episode allowance, sparse probes through `h=N` are still below `-111645` bits. This is evidence for the next theorem
target, not yet a full interval certificate.

Further update: the forced near-full-density endpoint also exposed an occupancy-coefficient pole issue. At
`T=32767`, `H=bT=2097088`, `lambda=2.3`, the eall-ratio inner bound improves from `-893425.892127` with `rho=10` to
`-1181735.473865` with `rho=1e6`, equal to the no-turnoff survival value. This supports treating the endpoint by a
large-`rho`/deficit-side coefficient envelope rather than by adding heuristic episode slack.

Outer-side update: the RM512 local spectrum is complement-symmetric, so the direct-sum outer can use
`A_h=A_{N-h}`. The paired probe and piecewise certificate now expose this via `--outer-complement-symmetry`; an
`h=N`, `r=64` smoke gives `-1181735.473865` in the paired probe and `-1181723.239347` in the endpoint piecewise
wrapper, the latter differing by the expected bucket-size overcount.

Paired all-episode update: `probe_fullsplit_early_paired_survival.py --inner-mode eallratio` now attaches the fixed
all-episode wrapper to each exact `T` placement term. Sparse high-density checks with complement symmetry have the
current sampled worst near `h=1200501`, local `r=32..43`, at `total_log2=-132647.058774`. This is not an interval
certificate yet, but it strongly suggests the remaining high-density task is packaging, not a new obstruction.

### F6 Medium: Product form for per-gap terminations now has a manuscript lemma

Resolved in the current manuscript draft by Lemma `lem:fullsplit-selected-gap-product`, which uses sequential exposure and the chain rule rather than unconditional independence. The numeric local atoms are now packaged in Lemma `lem:fullsplit-certified-termination-atoms`.

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

1. Should the high-interval manifest remain command-reproducible, or should every interval be regenerated into checked artifacts?
2. Should the opt-in high-interval recomputation be expanded into a scheduled/full certificate run?
3. Should Lemma `lem:fullsplit-tiny-prefix-interiorT` remain as a finite certificate, or should it be replaced by analytic monotonicity?
4. What theorem will cover the remaining first-active regimes `T > 17948, h > 500` and ultra-late `h > 500`?

## Suggested Next Fixes

1. Prove or certify the remaining first-active regimes `T > 17948, h > 500` and ultra-late `h > 500`.
2. Decide whether to keep the interval rows as a manifest with spot recomputation or regenerate all rows as explicit artifacts.
3. Make `sum_fullsplit_piecewise_certificate.py` assert or use `feasible-min` when endpoint placement is paired with high `H`.
4. Add interval/rational arithmetic for the final ledger.
