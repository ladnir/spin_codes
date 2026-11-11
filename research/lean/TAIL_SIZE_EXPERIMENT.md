# Tail certificate size experiment (2026-11-10)

Question: can the generated Lean tail certificates for the Turbo/RS8 s22 ladder
be made much smaller? This was a narrow experiment. Nothing was productionized, and no existing file
was modified. New files: `scripts/size_exp/{diag_marker,qstar_corrected,marker_probe}.py`,
`scripts/size_exp/qstar_out.txt`, `scripts/size_exp/_bcache/` (float caches),
`SpinCodes/Turbo/SizeExp/G22MarkerProbe.lean` (generated probe).

## Headline

The "75 bits per region" gap between the marker and exact bounds mostly comes from a
**planner bug**. The rate function is fine. Both untrusted planners score the marker term with
`C(L,q)^-62`. The Lean theorem (`pairedRegionDP22_marker_pow64_le`,
final check `tailMarkerFracFinalCheck22`: `... <= a * C^63 * tn0^(64q)`) and the
producer (`stage_marker`: `- 63 * log2(Cq)`) use `C(L,q)^-63`. So the planners
under-credit every marker term by log2 C(L,q) bits:

| q (G22) | reported gap marker-exact | log2 C(32768,q) | remainder | remainder / 64 |
|---|---|---|---|---|
| 352 | 3593 | 2802 | 791 | 12.4 |
| 495 | 4527 | 3697 | 830 | 13.0 |
| 546 | 4921 | 4001 | 920 | 14.4 |
| 760 | 5835 | 5204 | 631 | 9.9 |

Locations: `scripts/ladder_tail22.py::marker_term` uses `lch(L,q) + ... - 63*lch(L,q)`,
for a net coefficient of -62. `scripts/ladder_marker_plan22.py::best` uses `a0 = lch(L,q0) * -62`.
`exact_term` correctly nets -63. This is a planner-only issue: the Lean checks are
correct, and every certified bound stands. The bug only made the plans far more expensive
than necessary.

## (A) Marker vs exact: diagnosis

`scripts/size_exp/diag_marker.py` uses G22, tilt 26214425/2^30, alpha 3/5, n = 1024. It computes the
tilted distribution of the degree under `P~ = t*^g P[g] / rho(t*)` at the Chernoff
saddle t*, in exact fixed point (flint).

| q | t* | -log2 rho(T_q) (per-region marker loss) | max-entry loss | tilted mean g |
|---|---|---|---|---|
| 64 | 2^-6.62 | 4.34 bits | 4.3 | 64.7 |
| 352 | 2^-4.81 | 7.55 bits | 1.9 | 369.2 |
| 760 | 2^-4.76 | 8.37 bits | 2.2 | 769.7 |

The row-0 mass is unimodal and peaks at q. There are no competing phases and no non-convexity.
The marker's true loss is the expected Chernoff width loss of about 4-8 bits per region (x64 = 280-540
bits total). The leftover 10-14 bits per region in the table above come from comparing at a fixed
tilt/t rather than at each method's optimum. **A phase-split marker is not needed. NO-GO on that idea.**

## (A) Corrected crossover q*

`scripts/size_exp/qstar_corrected.py` re-screens over all 481 cached families x 13 alphas x
1088 rational t. It uses the same float model as the fine planner and the same targets
`-(q1_bits+20)`, and runs once with each normalizer. Output: `scripts/size_exp/qstar_out.txt`.

| rung | target | q* planner (C^-62) | q* in use | **q* Lean (C^-63)** | marker ranges from q* |
|---|---|---|---|---|---|
| G22 (L=2^15) | -76.37 | 405 pointwise (546 with ranges) | 546 | **53** | 10 |
| G24 (L=2^17) | -74.37 | 1174 | 1300 | **72** | 11 |
| G26 (L=2^19) | -72.35 | 4229 | 5000 | **105** | 13 |

Sample values (C^-63, log2 of the term): G26 q=128 -179, q=256 -885, q=1024 -6489. Below q*
the marker fails because log2 C(L,q) is small and the Chernoff loss of 64 x 4-8 bits dominates.
Those q stay exact. They are the cheap pieces.

### Producer / Lean check

`scripts/size_exp/marker_probe.py` reuses the built G22 family and pair modules read-only. It runs the
real producer (`stage_marker`, exact big-integer arithmetic, cache in scratchpad) on a
new marker range piece over q values the current plan covers with exact pieces:

| range | screen (float) | producer (exact ints) | verdict |
|---|---|---|---|
| 352..495 (replaces exact piece P12, 258 MB src / 1.65 GB olean) | -2185.7 | -2180.2 (a = 1, i.e. <= 2^-256) | passes by about 2100 bits |
| 496..545 (P13) | -3279.4 | -3274 | passes |
| 64..70 | -89.1 | -83.4 | passes (target -76.37) |
| 64..123, 124..247, 248..545 (one t) | +587 / +732 / +1996 | fail | the range is too wide for one t |

The producer is about 5-6 bits **less** favourable than the float screen, because Lean bounds the row sum of
`M_t^(64n)` and not `rho^(64n)`, and adds rounding. So ranges that the screen accepts
with less than ~8 bits of margin need splitting: G22 53..66 (-76.9), G24 72..87, 88..124, 125..207,
G26 105..116, 230..387, 7835..19295. Expect q* to rise slightly (about 10-20%) and a few more
ranges. Nothing else changes.

Lean compile of `SpinCodes/Turbo/SizeExp/G22MarkerProbe.lean` (q = 352..495), done locally with
`compile.ps1` at its default `-M 6000`: **kernel "excessive memory consumption" at `baseChecked`**
after 375 s. That module is structurally identical to the production `TailP14MBase` modules,
which are green on Peach under the ladder's standard `-M 14000`. Only the range and t differ. I did not retry at
14000 locally because only about 7.4 GB was free next to the planners. **To do:** compile the probe at the
standard limit on Peach or locally once memory is free. This is the one unverified step.
Size of a marker range piece (measured, G22 P14-P22): 0.4 MB source, 7.9 MB olean each.

## Size estimate with the corrected planner

Measured G22 today (`Generated/G22Tail`): exact pieces P1-P13 = 1037 MB source / 7.5 GB olean.
Family+pair modules = 136 MB / 1.7 GB (21 families). Markers 9 x 0.4 MB.
Exact pieces up to q = 46 (P1-P5) cost only 44 MB source in total.

| rung | exact pieces today (source) | exact pieces with corrected q* (source) |
|---|---|---|
| G22 | 3..545: 1037 MB (7.5 GB olean) | 3..52: about 60 MB (about 0.6 GB olean) |
| G24 | 3..1299: about 1 GB est. | 3..71: about 0.1 GB est. |
| G26 | 3..4999: about 4-5 GB est. | 3..104: about 0.15-0.2 GB est. |

The estimates use the measured cost of about 6 MB + 0.33 MB per degree per piece, with one extra squaring per rung.
Family/pair overhead (about 6.5 MB source / 80 MB olean per family) then dominates. Reusing tilts across ranges
cuts it further. Whole tails: about 0.2-0.4 GB source per rung, against 1.2 GB (G22) and several GB (G24/G26) today.
The sum of the per-q bounds does not change in kind: marker ranges already contribute
`max * count` near the target, as the current G22 markers do.

## (B) Kernel-computed chains: NO-GO for now (not prototyped)

With q* <= about 105 at every rung, the largest exact chain has degree about 104 and the stored
stages total about 0.1-0.2 GB per rung. Kernel-computed Kronecker-packed chains would remove at most that.
They would also need new soundness lemmas (packing, carry-free bounds, rounding-up division) and
kernel memory profiling. That is not worth it until the corrected plan has been tried. I would only revisit it if
family overhead is reduced and the exact pieces become the bottleneck again.

## (C) Cheaper stage encoding: NO-GO (unnecessary)

## Recommendation

**GO on (A) as a planner fix.** Correct the -62 -> -63 normalizer in
`ladder_tail22.marker_term` and `ladder_marker_plan22.best`. Require about 8 bits of screen margin, or
check range endpoints with the integer producer, as `marker_probe.py` does. Then re-plan.
No Lean or theorem changes are needed: all pieces use the existing `tail_term_fracAB_*` and
`TailMarkerRange22` interfaces.

Next steps:
1. Compile `SpinCodes/Turbo/SizeExp/G22MarkerProbe.lean` at the standard `-M 14000`. Expect green
   with axioms [propext, Classical.choice, Quot.sound].
2. Fix both planners, re-plan G24 with exact pieces 3..~80 and markers above, and build it.
   Then do G26 the same way.
3. Optional: re-plan G22 to shrink its regenerated size from 1.2 GB to about 0.3 GB. Its endpoint
   2^-56.3691 is already certified, so this matters only for reproducibility cost.
