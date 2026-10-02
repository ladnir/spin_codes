# Shared GF16 proof progress

The lower-distance two-update follow-on is now separate and successful:
at K=2^20, its fresh whole-code replays prove >7% relative distance with
more than 44.61 bits of setup-failure margin, and >5% with >73.38 bits, for the same
measured 5.863 ms encoder.
The [current record](SHARED_SHUFFLE.md#lower-distance-search-2026-09-29)
tracks those results and the still-open 7.5%/8% frontier. The earlier
>6% / >40.68-bit proof is preserved but dominated by the 7% result.
The final 7% replay on 2026-09-30 checked all 303 dense intervals and all
sparse occupancies q=1..32 at 384-bit precision. Its whole-code margin is
44.6143782060373124 bits, with a separate scope/hash/sum audit. All 405
proof-tool tests pass. The original four-update 10%
goal and its plateau are preserved unchanged below.

The goal is a complete certificate for the shared-shuffle four-update
encoder at K=2^20: distance greater than 10% and bad-setup probability
below 2^-40. The distribution and IMT(128,19) parameters are fixed during
this goal. A different distance or encoder is a separately labeled result,
not completion of this target. Old constructions and their certificates
remain indexed in [the proof index](../PROOF_INDEX.md).

**Outcome, 2026-09-29: stopped at the agreed plateau. The whole-code 10%
proof remains open.** The verified prefix advanced from q=1–18 to q=1–22.
The final search still leaves q=23 unresolved, and there is no complete
dense cover. This stops the bounded investigation, not the construction's
proof obligation. It does not show that the construction has low distance.

## Baseline

At the start of this goal, all supports at occupancies q=1,...,18 are
covered. Their aggregate margin exceeds 44.20 bits. Occupancy 19 is the
first gap; its best checked search produced a floating 12.7567-bit
proposal, not an outward certificate. The corresponding proposals at
q=20,24,32 have log2 uppers +42.53, +268.01, and +748.76.
No occupancy beyond eighteen has a complete certificate for this route.
The measured encoder time remains 6.532 ms; this goal does not retime it.

The most recent useful refinement kept output weight and return to zero
coupled. It raised q=16 from 40.16 to 113.71 bits and closed q=17,18.
A separate, exactly checked positive mixture of shared-support shells
was too loose in the basic scalar dense analysis. Its correctness is
established, but it did not improve occupancy coverage.

## Current checkpoint

The complete prefix is now **q=1–22**, with aggregate restricted-class
margin **44.20717305700 bits**. Occupancy 19 closes at 64.7197006266 bits
with the coupled counts. Its fresh 256-bit cover has 39 leaves and two
splits. The receipt is `tmp/shared-goal/coupled-counts-q19-32.json`, SHA256
`3ba099c055fb4b4bf4038eea84369fa65fbf2cd1a02df65847b741e0d28ad499`.
An independent 384-bit regeneration reproduces that result in
`tmp/shared-goal/q19-p384.json`, SHA256
`bbe2fc6350fc8bb60ba4d1a5d005e9849e37a3f9550ca619ba5630ad65497c27`.
The exact dyadic sum for q=1–19 is less than 2^-44.2; this combines the
earlier q=1–18 receipts with the new q=19 bound, not a dense certificate.
The earlier coupled-count run tested every occupancy through 32 but did
not close 20–32. Joint CDF constraints now close **q=20 at 79.71067 bits**
and **q=21 at 58.19630 bits**, using fresh 256-bit arithmetic. Their
receipt is `tmp/shared-goal/joint-counts-q20-24-bounded.json`, SHA256
`0140a157aa80609c695a331fe37bc05b7ecdc364ffe1c04e1bd43cf9f0a6de3d`.
An independent 384-bit regeneration confirms both occupancies and improves
their bounds to **81.84346 and 62.87364 bits**, respectively. The gain is
from additional exact count witnesses, not a claim that higher arithmetic
precision itself improves the construction. That receipt is
`tmp/shared-goal/joint-counts-q20-22-p384.json`, SHA256
`303198c78ac5b6c9d785edd33a9506179fe502bb7a1e52fedad1cef3891b9a12`.

The wider joint model closes **q=22 at 54.12979782879 bits**, with fresh
256-bit arithmetic over all supports. Its receipt is
`tmp/shared-goal/wide-counts-q22-32.json`, SHA256
`e319cfd4034313f6ff0e7d1401ff51122653e1e547ddcc8d8f89c6b0fab119b4`.
An independent 384-bit regeneration confirms closure, at **48.292552 bits**.
It found weaker count witnesses under the bounded numerical search; this
is not an arithmetic-precision regression. The independently regenerated
q=22 contribution gives aggregate **44.12599953929 bits** when combined
with the retained q=1–21 receipts. The stronger 256-bit contribution remains
valid. The new receipt is `tmp/shared-goal/wide-counts-q22-23-p384.json`,
SHA256 `c3e9184aa222d7a7a19c6f84545f23e63fce5c0f73cbfbdfbf92448879ea6baa`.
The first gap is **q=23**, with a final floating complete-cover proposal
of **23.20064 bits** after 512 splits, versus 21.83767 after 48 splits.
No outward certificate was produced at q=23. The final receipt is
`tmp/shared-goal/q23-retained-control.json`, SHA256
`ccf492344aeefe2160cb09d141ce4b5579fdd1ba65c59950a27c441953fe79b4`.
The earlier proposals at q=24,25,26,28,32
have log2 uppers +20.4055, +64.1411, +106.1074, +202.9473, and +390.8440.
No outward result is claimed for these unresolved occupancies.
The exact dyadic aggregate through q=22 is below 2^-44.2. No class beyond
q=22 has been incorporated into that guarantee.

## Progress and stopping rule

After each bounded experiment, record its mathematical change, receipt,
quantitative result, verification status, and next decision. A run is not
progress merely because it consumes more work or creates another file.
Meaningful progress is one of the following:

- Newly covered occupancies or a smaller uncovered range, with outward
  checks and correct setup scope.
- A material improvement to an unresolved bound, distinguished from
  floating search evidence until checked outward.
- A proved counting or transfer inequality that improves the bottleneck
  on the actual instance, with exact or outward numerical checks.

Minor tilt refinements are controls for search artifacts, not distinct
mathematical approaches. Coverage and complete-cover bounds take priority
over favorable singleton probes. Aggregate bounds must include every
covered occupancy; selected points cannot fill gaps.

Consider a plateau only after three consecutive substantive experiments
spanning distinct mathematical approaches give no meaningful improvement.
Before stopping, check for inadequate witness search, precision, and
implementation errors, and review untried high-value approaches. Any useful
improvement resets that count. A plateau report must state that the proof
is open, identify the first gap, and explain the evidence for stopping.

## Experiment ledger

| Experiment | Question or change | Result and verification | Decision |
|---|---|---|---|
| Baseline, 2026-09-29 | Shared counts plus coupled GF16 returns through local occupancy three | Complete prefix q=1–18, >44.20 bits; q=17 independently replayed at 384 bits | Continue from q=19; no plateau |
| Coupled counting, 2026-09-29 | Exponential basis-weight bound, then alternate primal/dual identities with fresh containment bounds after each exchange | Exact caps improve by 10.4366 bits at support 120, 16.4639 at 128, 15.0909 at 136, and 13.5010 at 144. Eight-iteration run took 112 seconds. A fourth iteration is within 1e-9 bits of the eighth at these points; no exact fixed point is claimed | Material counting progress; regenerate all-support covers q=19–32. The basis bound alone gave no gain at 120 or 128 |
| Coupled-count covers, 2026-09-29 | Use four coupled-count iterations with the existing three-packet return kernel | q=19 fully verified at 64.7197 bits; q=20–32 unresolved. No setup change | Contiguous prefix advances to 19; nonprogress streak resets to zero |
| Four-packet returns, 2026-09-29 | Extend the exact joint feedback/output census through local occupancy four; keep the baseline outer bounds as a control | All 1,820,475,000 local inputs checked against the independent feedback census. q=19 proposal improves from 12.7567 to 12.7654 bits, about 0.009 bits; no new coverage | Stop extending this census in the current envelope; the return tail is not the main bottleneck |
| Empty-component diagnosis, 2026-09-29 | Inspect the positive comparison measure itself before refining its transfer kernel | Its active-component empty-input mass is about 2^106.69 per group. The comparison therefore cannot certify the target, even with an exact inner | Identifies a mathematical obstruction in the comparison, not the encoder |
| Empty-budget mixture, 2026-09-29 | Require each active component to contribute at most 2^-64 to its empty-input atom, with all shell inequalities still checked exactly | Aggregate empty-active mass falls below 2^-61.22. Outward singleton diagnostics improve substantially but remain positive across the interior; no new covered occupancy | Keep the constraint; optimize the tilted comparison mass next. Selected means are not a complete-domain improvement claim |
| Tilted mixture, 2026-09-29 | Optimize component cost at input tilt 1/4, retaining the exact empty-input and shell checks | Basic outward singleton log2 bounds include +1906.10 at mean .0005 and +67336.67 at .032; no new coverage | Use as a controlled input to the stronger regional analysis, not as a distance certificate |
| Longer shortening LP, 2026-09-29 | Extend exact-verified primal and dual shortening LPs from length 104 through 160, then repeat coupled containment | Gain at support 96 is 1.307 bits beyond the new baseline. Across the current bottleneck 107–120, the relative count decrease is at most 4.18e-26 | No material gain where needed; extending these same LPs is low priority |
| q=20 witness-search control, 2026-09-29 | Finer output tilts and four times the split budget | Floating complete-cover proposal changes from 20.9053 to 20.9643 bits, only 0.059 bits. 785 retained leaves | Not a precision or coarse-grid explanation for the roughly 19-bit gap to the target |
| Regional dense comparison, 2026-09-29 | Retain conditional regional counts, variance, feedback classes, and lazy-density bounds for the tilted majorant | At the same mean .0005, the outward log2 bound falls from +1906.10 to +615.6485; at .032, from +67336.67 to +66953.35. Still no certified point or interval | Material conditional-bound improvement, but no dense coverage. Do not discard this approach merely because it has not closed |
| GF16 extension moments, 2026-09-29 | View BCH quadruples as GF16 codewords; use dual distance 30 and a squared degree-14 reproducing kernel | Exact counting bound at support 112 is 2^408.615 versus the existing 2^303.616 cap. No improvement at any of 96,104,108,112,116,120,128 | Dual-distance moments alone are weaker than the coupled binary-subspace bounds; retain as a checked negative result |
| Joint CDF constraints, 2026-09-29 | Retain primal/dual CDF variables simultaneously, including containment, complement identities, monotonicity, and rank-one shell bounds | Initial exact dual checks improve rank-four caps by 3.1921 bits at support 108, 4.6134 at 112, and 4.0301 at 116. A box-residual correction covers every floating-solver error | Material counting progress; extend across 96–136 and regenerate complete support covers q=20–24 |
| Joint-count support covers, 2026-09-29 | Check every support at each requested occupancy using only exactly verified count improvements | q=20 closes at 79.7107 bits, q=21 at 58.1963 bits; q=22–24 remain open. Rank-four improvement reaches 7.3728 bits at support 115 | Prefix advances from 19 to 21; first gap moves to 22. Independently regenerate at 384 bits |
| Wider joint LP, 2026-09-29 | Increase the support limit from 144 to 192 | All six selected objectives returned no usable dual under the bounded default/tighter solver attempts. Each retained its prior bound | Numerical search failure, not evidence that the larger exact system is useless. Add a relaxed-equality proposal with exact verification |
| All-ones complement constraints, 2026-09-29 | Add the rank-one identities F(u)+F(255-u)=2^128-2 for the code and its dual | An exact 0.6629-bit gain at support 120, where the earlier joint LP gave no gain. No benefit established at 119,124,128. Not yet used in a full support cover | Retain as a valid additional counting constraint; do not attribute its gain to the q=20,21 certificates |
| Independent joint-count regeneration, 2026-09-29 | Fresh counts, local census, and 384-bit all-support covers, with an additional bounded proposal strategy | q=20 and q=21 close at 81.8435 and 62.8736 bits. q=22 remains open with a 26.6139-bit floating proposal | Confirms new coverage and improves the unresolved bound; no whole-code certificate |
| Wider joint LP with repaired proposals, 2026-09-29 | Propose multipliers using relaxed numerical equality rows, then check them against the original exact equations | Exact gains of 8.3135 bits at support 112 and 8.9482 at 120. No usable gain at 119 or 128 in this bounded probe | The previous wider-model failure was numerical, not a mathematical plateau. Feed the checked bounds into a new all-support cover |
| Cheap hyperplane counting, 2026-09-29 | Every rank-h space has a hyperplane with smaller support; count such hyperplanes and their extensions in shortened codes | Exact small-code tests pass. No material gain on the BCH bottleneck: at most 0.00115 bits on the selected supports, with no gain at 116,119,120,128 | Retain as a checked negative result. Its best witnesses reduce to the existing full-hyperplane count, so do not pursue this variant further now |
| Wider-count support covers, 2026-09-29 | Combine exactly verified 144- and 192-support CDF caps, then cover all supports at selected occupancies | q=22 closes at 54.1298 bits; q=23 remains open at a 21.8377-bit floating proposal. The restricted aggregate through q=22 exceeds 44.20 bits | Prefix advances to 22 |
| Independent wider-count regeneration, 2026-09-29 | Fresh counts and 384-bit verification; increase q=23 split budget to 128 | q=22 closes at 48.2926 bits. q=23 stays open at a 3.4807-bit proposal because a useful support-114 counting witness was not found in the bounded LP search | Confirms closure, not the identical bound. The q=23 comparison is confounded by different counts, so it is not a clean split-budget plateau test. Retain exact dual witnesses for stable replay |
| Component-mass budget, 2026-09-29 | Add c<=2^520 to the greedy positive mixture, retaining the empty-atom constraint and exact shell checks | Shell domination and both budgets pass. The largest component exponent falls from 549.81 to 519.71. Both outward diagnostics remain +615.6485 and +66953.3499 | No downstream gain or coverage. Total component inflation was not the active bottleneck at these probes |
| Full-support joint LP, 2026-09-29 | Extend the joint constraints through support 256, including the endpoint total-count equalities | All five selected objectives return no usable witness: standard/relaxed solvers report numerical infeasibility and the small-coefficient solver reaches its five-second limit | Numerical failure, not evidence that the exact system is infeasible. The later 20-second coarsened control is recorded below |
| Residual-dual refinement, 2026-09-29 | Re-optimize the exact dual residual at unit scale, then recheck the sum of multipliers against the original objective | Gain beyond the existing witness is only 0.00220 bits at support 104; none at 112,120,128. Every retained result is rationally checked | No material gain. This is a control within the joint-LP approach, not a distinct mathematical success |
| Fourier-monotone comparison, 2026-09-29 | Replace nonzero labels by a Fourier-dominating lazy law, then use stochastic domination from the support CDF | Exhaustive tests over all binary maps from two two-bit packets to two output bits pass. Fresh outward diagnostics at comparison means .0005 and .032 are +2510.3930 and +125899.1624, for q>=23 | Valid lemma, but no covered point or demonstrated benefit on BCH. The mean coordinates are specific to the comparison measure, not matched message classes across models |
| Full-support numerical control, 2026-09-29 | Increase per-proposal limits to 20 seconds; explicitly budget dropped negative coefficients in a coarsened numerical proposal | The full 256-support model now yields an exact 8.2203-bit gain at support 112, weaker than the existing 8.3135-bit gain. No witness at 120. Original constraints, not the relaxed proposal, verify the result | The larger model is not wholly unusable, but this bounded control gives no new bound. Do not infer infeasibility of the exact model |
| Retained-witness search control, 2026-09-29 | Recover and save the support-114 dual; replay it and the full-support dual against fresh constraints before a 512-split q=23 search | q=23 remains open. The full-cover proposal reaches 23.20064 bits, a 1.36297-bit gain over 48 splits. From 250 to 512 subdivisions the gain is only 0.08282 bits. Both saved duals replay, the three-packet census checks 16,740,000 inputs, and all support boxes are accounted for | The small search gain does not remove the roughly 17-bit shortfall at the first unresolved occupancy. Stop at a plateau; no new coverage or dense progress |

Since the q=22 advance, **three approaches have given no material gain**:
component-mass budgeting, extending the joint counting system, and the
Fourier-monotone comparison. Residual refinement and numerical coarsening
are controls within the LP approach, not additional independent failures.
The final retained-witness q=23 search has finished. It used finer tilts,
384-bit operators, and 512 subdivisions and retained the support-114 dual.
It also found extra exact count improvements at supports 118 and 128.
These local improvements are not ignored: the resulting complete-cover
proposal improved by 1.36297 bits. However, the last 262 subdivisions
gained only 0.08282 bits, no new occupancy closed, and the dense comparison
did not improve. That is insufficient to restart the search after the
three negative approaches. A **plateau of the present methods** is now
recorded, not an impossibility result or a proof of the target.

The sparse priority is the q=23 gap, rather than extending the local
return census. The completed wider-count run checked q=22,23,24,25,26,28,32
with the componentwise minimum of the verified 144- and 192-support models.
The wider model is solved at even supports 96–136; monotonicity supplies
valid intervening caps. This does not skip odd-support message classes:
the subsequent support cover still includes every support at each q.
Its output is `tmp/shared-goal/wide-counts-q22-32.json`, with the retained
log `tmp/shared-goal-wide-counts-cover.log`. Do not treat an unfinished
run or an absent result as a verified contribution. At
support 112 the four rank contributions have log2 caps 126.85, 184.45,
245.21, and 303.62: rank four dominates. Existing iterations propagate
separate CDF caps; retaining the constraints together now demonstrably
removes some of that loss and has closed three additional occupancies.
The regional majorant remains an alternative, but its checked probes do
not justify a full-domain search yet.

The completed control is `tmp/shared-goal/q23-retained-control.json`, with
log `tmp/shared-goal-q23-retained-control.log`. The process exited normally.
Its null upper means that no requested outward certificate was obtained;
it is not a zero contribution. No further sweep was started.

## Final Audit and Next Choices

The fixed target was not reached. No full-code claim follows from the
partial sum. The first missing occupancy is q=23; no complete certificate
has been assembled for any of q=23–2048. The proof-search stopping condition
was reached by the experiments and numerical controls above.

The exact scope-and-sum audit of the retained q=1–22 receipts confirms
the same four-update setup and cutoff 209715, with no missing or repeated
occupancy in that prefix. Using the independently regenerated q=22 result,
the aggregate margin is 44.12599953929 bits. This receipt audit does not
substitute for numerical replay; the separate fresh runs are identified above.

The old archive passes all 4,113 member hashes and its recorded ZIP hash.
All 15 local links in the proof index resolve. The independent-row four-update
10% certificate still passes its scope-and-sum audit: 427 disjoint exhaustive
dense intervals, matching sparse endpoint, and total upper below 2^-48.
The updated GF16 test suite passes 243 tests; the archive helper passes four.
No production encoder, paper claim, benchmark, Git history, or old proof
receipt was changed by this investigation.

The final sources, ledger, receipts, and logs are also preserved in
`research/.proof-archives/packet-proofs-shared-gf16-plateau-20260929.zip`:
4,196 files, 195,295,905 bytes, SHA256
`f5da6b5909d78dcb4ff374732d850621c76a6a676449a6b1f6e143f9560977f0`.
All members pass the manifest hash check. This supplements the original
archive without overwriting it; both remain local and ignored by Git.

The highest-value next choices are:

- For a proved 10% code now, retain the independent-row GF16 four-update
  construction. Its matched precomputed timing was 7.506 ms versus 6.532 ms
  for this unclosed shared-shuffle version; neither was retimed here.
- For the same shared-shuffle distribution, seek a stronger bound on the
  number of rank-four BCH subspaces with union support near 100–120.
  More local return census, wider instances of the same counting LP, and
  more cover subdivisions have demonstrated low returns. Rescaling or
  reformulating the LP may merit a separate experiment, but is not a
  missing certificate implied by these results.
- For a construction change, investigate two independent pairwise coordinate
  shuffles within each four-row group. This reduces the outer counting
  problem to two rank-at-most-two pairs before their supports overlap.
  It is a proposed route to an easier proof, not a proved distance or
  performance claim. It requires its own implementation comparison and
  certificate; the present proof does not transfer automatically.

The practical recommendation is to keep the independently shuffled code
as the certified reference and explore the pairwise-shared route next,
if another construction experiment is wanted. Do not continue increasing
the current search budget without a new reason to expect a useful bound.

## Evidence

The baseline commands, exact setup, receipt names, and hashes are in
[SHARED_SHUFFLE.md](SHARED_SHUFFLE.md). The archive named in the proof index
preserves the earlier source and witness checkpoint. New experiment data
belongs under ignored `tmp/shared-goal/`; this ledger and proof sources
remain separate from numerical data.

Additional receipts in `tmp/shared-goal/` are `return-four-q19-20.json`,
`longer-shortening.json`, `zero-budget-mixture.json`, `tilted-mixture.json`,
`q20-refined-search.json`, `regional-mixture.json`, and
`extension-moments.json`. Matching `tmp/shared-goal-*.log` files retain
floating search diagnostics, including cases whose JSON upper is null.
Null means this run did not establish its requested margin; it is not a
zero contribution. The failed experiments are retained, not discarded.

The joint-count receipts are `joint-constraints.json`,
`joint-counts-q20-24-bounded.json`, `joint-constraints-192.json`, and
`joint-constraints-symmetry.json`, `joint-counts-q20-22-p384.json`,
`joint-constraints-192-relaxed.json`, and `hyperplane-counts.json` in the
same directory. Their numerical
LP statuses are diagnostic. Every adopted bound has an exact dual check.

This checkpoint also retains `wide-counts-q22-32.json`,
`wide-counts-q22-23-p384.json`, `mass-budget-mixture.json`,
`joint-constraints-256.json`, `joint-constraints-residual.json`, and
`monotone-mixture.json`. The completed numerical controls wrote
`joint-constraints-retained-witnesses.json` and
`joint-constraints-coarse-256.json`. New count receipts retain their exact
dual multipliers. `constraint_counts.py --replay` applies those multipliers
to freshly reconstructed constraints and recomputes the cap; stored caps
and solver feasibility are never accepted as premises.

The updated GF16 suite passes all **243 tests**. The archive helper passes all
four tests. The reused containment,
dual-moment, and shortening/complement self-tests pass 88, 375, and 432
exact checks, respectively. No benchmark or production change was made.

Reproduce the 384-bit q=19 check:

```text
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 19 --precision 384 --refined-counts --coupled-counts --joint-return-through 3 --lazy-density-through 6 --tilts .008 .01 .011 .012 .013 .014 .016 .02 .024 .032 --max-splits 32 --target-bits 48 --proposal-buffer-bits .125 --output tmp/shared-goal/q19-p384.json
```

Regenerate q=20,21 and the q=22 diagnostic independently:

```text
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 20 21 22 --precision 384 --refined-counts --coupled-counts --joint-counts --joint-return-through 3 --lazy-density-through 6 --tilts .008 .009 .01 .011 .012 .013 .014 .016 .02 --max-splits 64 --target-bits 48 --proposal-buffer-bits .125 --output tmp/shared-goal/joint-counts-q20-22-p384.json
```

The completed wider cover command is:

```text
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 22 23 24 25 26 28 32 --refined-counts --coupled-counts --joint-counts --wide-counts --joint-return-through 3 --lazy-density-through 6 --tilts .008 .009 .01 .011 .012 .013 .014 .016 .02 --max-splits 48 --target-bits 48 --proposal-buffer-bits .125 --output tmp/shared-goal/wide-counts-q22-32.json
```

The independent regeneration uses the same flags with occupancies `22 23`,
`--precision 384 --max-splits 128` and output
`tmp/shared-goal/wide-counts-q22-23-p384.json`.

Generate the retained numerical controls:

```text
python -B research/workstreams/permutation_locality/gf16_packets/constraint_counts.py --last 192 --supports 114 118 128 --time-limit 20 --output tmp/shared-goal/joint-constraints-retained-witnesses.json
python -B research/workstreams/permutation_locality/gf16_packets/constraint_counts.py --last 256 --supports 112 120 --time-limit 20 --output tmp/shared-goal/joint-constraints-coarse-256.json
```

The final q=23 control is:

```text
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 23 --precision 384 --refined-counts --coupled-counts --joint-counts --wide-counts --count-witnesses tmp/shared-goal/joint-constraints-retained-witnesses.json tmp/shared-goal/joint-constraints-coarse-256.json --joint-return-through 3 --lazy-density-through 6 --tilts .008 .009 .01 .011 .012 .0125 .013 .0135 .014 .016 .02 --max-splits 512 --target-bits 48 --proposal-buffer-bits .125 --output tmp/shared-goal/q23-retained-control.json
```

`--count-witnesses` rechecks rational multipliers against fresh constraints.
It never imports the recorded count or objective value as a premise.

Reproduce the two dense diagnostics (not a cover):

```text
python -B research/workstreams/permutation_locality/gf16_packets/shared_mixture.py --coupled-counts --zero-bits 64 --cost-tilt 1/4 --minimum-groups 20 --variance-bins 4 --regional-count --probe .0005 .032 --output tmp/shared-goal/regional-mixture.json
```

The first counting experiment is `tmp/shared-goal/count-refinements.json`.
Regenerate it with:

```text
python -B research/workstreams/permutation_locality/gf16_packets/count_refinements.py --iterations 8 --output tmp/shared-goal/count-refinements.json
```

Its driver authenticates the BCH inputs and regenerates every bound. The
new basis inequality passes exhaustive small-code checks. The alternating
steps reuse the exact containment and dual identities; no recorded count
is accepted as a proof input. The production implementation is unchanged.

## Counting and comparison arguments

For the basis bound, fix a binary h-dimensional subspace with union support
u. An ordered basis has mean total weight h*2^(h-1)*u/(2^h-1), when averaged
over its B=product_j(2^h-2^j) ordered bases. Convexity of x^w for 0<x<=1
lower-bounds the sum of their exponential weights. The implementation
interpolates on the lattice of attainable weights rather than rounding
the mean upward. If M(x) is the upper weight polynomial of nonzero BCH
words, M(x)^h bounds the sum over all ordered independent bases of all
subspaces. Multiplying by product_j(2^g-2^j)/B converts the resulting
subspace bound to a bound on ordered rank-h g-tuples. Rational arithmetic
verifies every selected tilt; floating optimization only proposes it.

The coupled step alternates two already verified inequalities: dual-complement
identities for shortening moments, and containment bounds relating different
subspace ranks. Reapplying containment after each dual exchange propagates
new lower-rank bounds into the higher ranks. The earlier loop exchanged
dual bounds without this intermediate propagation. Every iteration retains
the smaller of valid bounds, so stopping before convergence remains sound.

For the dense comparison, a component of mass c and Bernoulli activity p
assigns mass c*(1-p)^256 to an empty group even though it has an active
label. The initial mixture makes the sum Z of those masses exceed one.
For exactly q active groups, it then contains an all-zero-input contribution
binomial(2048,q)*Z^q. A linear inner maps that input to zero. Thus no valid
analysis of that unrestricted comparison can return a subunit total bound.
This does not imply any actual nonzero BCH message maps to zero.

The new component constraint is checked exactly together with shell
domination. It removes that obstruction without assuming a new encoder
distribution. The implemented assignment minimizes the required tilted
component mass c*(1-p+p*x)^256 separately for each shell, for a chosen
0<x<=1. This changes the choice of majorant,
not the shell inequalities that establish its validity. It is a greedy
assignment with exact verification, not a claim of a globally optimal
positive mixture.

For the extension-moment experiment, expand the binary generator over
GF16. Its codewords correspond bijectively to four binary codewords, and
their symbol support is exactly the binary union support. The extension
has 2^512 words and dual distance at least 30, since every nonzero field
dual word has a nonzero binary coordinate component in the binary dual.
Consequently, weight polynomials of degree below 30 have the same mean
as a uniform ambient GF16 word. The squared reproducing kernel has degree
28. Its mean is evaluated exactly by Krawtchouk orthogonality; its minimum
over integer weights 38 through u bounds the nonzero prefix count, after
subtracting the zero word's contribution. Polynomial recurrence,
orthogonality, and the count inequality are checked on small codes. This
valid bound is too weak at the current BCH bottleneck.

For the joint constraints, let F_h(v) count h-dimensional binary subcodes
whose union support has size at most v. The model contains these counts
for both the BCH code and its dual. For t<n, the sum of contained
h-dimensional subcodes over all t-coordinate sets is

    M_h(t) = sum_{v=0}^t F_h(v) binomial(n-v-1,t-v).

This follows by summation by parts. At t=n, M_h(n)=F_h(n). A uniform
shortening dimension bound d_t gives, for r<h,

    Gaussian(d_t,r) M_h(t) <= Gaussian(d_t,h) M_r(t).

The existing dual-complement identities give exact equalities between
primal and dual M_h values. Monotonicity and authenticated rank-one shell
bounds add further constraints. Each F_h(v) is divided by its existing
integer upper bound, so the resulting unknowns x lie in [0,1]. All row
normalizations are powers of two.

Write the resulting constraints as Ax<=b and Ex=d. For objective c.x,
any rational y>=0 and rational z give the valid upper bound

    y.b + z.d + sum_i max(0, c_i-(yA+zE)_i).

The last term uses x_i in [0,1] and repairs every residual coefficient.
The solver merely proposes y,z; the implementation rounds them to dyadic
rationals and evaluates this expression exactly. It then rounds the
resulting subcode count downward and restores the exact tuple multiplicity.
Solver failure retains the previous bound. Solver tolerances, numerical
feasibility, and approximate dual objectives are never proof premises.

The initial 60-second-per-support tighter solver sometimes cycled or
reported numerical infeasibility. Its partial log is preserved as
`tmp/shared-goal-joint-counts-cover.log`; that run was explicitly stopped
before its support cover. The revised driver tries the default and
small-coefficient solver settings, each with a five-second limit, and
retains the better exactly checked bound. This prevents a failed numerical
proposal from erasing either prior evidence or a successful alternative.

A third, bounded proposal replaces numerical equality rows by two
inequalities with small slack. Its multipliers are converted back to
equality multipliers and checked against the original exact equations.
The slack is not a proof assumption. This is a numerical fallback, not a
change to the mathematical constraints or setup distribution.

The coarsened proposal additionally discards coefficients smaller than
1e-8 in magnitude and compensates for each removed negative coefficient
using x in [0,1]. It also adds numerical slack. Signed pairs of multipliers
for the relaxed equality rows are converted back into equality multipliers.
Only verification against the original rational matrix contributes a bound.
Re-optimizing an exact residual similarly yields another set of multipliers;
adding those to the first dual is checked against the original objective.

## Fourier-Monotone Comparison

Fix all routing and inner maps, excluding the independent GF16 packet
labels. With zero initial state, the output is a binary linear function of
those labels. For 0<z<=1, write each output factor as

    z^y = (1+z)/2 + (1-z)(-1)^y/2.

The product over output coordinates is therefore a sum of characters with
nonnegative coefficients. A uniform nonzero GF16 label has character mean
1 for the trivial character and -1/15 for any nontrivial character. Replacing
each factor by its absolute value can only increase the total moment.
Those absolute values are exactly the Fourier coefficients of a label
that is zero with probability 1/8 and uniform nonzero otherwise.

In this comparison, adding an active support position replaces some factors
1 by 1/15 and leaves all others unchanged. The output-weight moment is thus
nonincreasing under support inclusion, for every fixed remaining setup.
Uniform supports of successive sizes can be coupled by taking prefixes of
one uniform permutation. The averaged moment is nonincreasing in each
group's support size as well.

Let A(u) be the actual nonzero-group support CDF and let B(u)>=A(u), with
B(0)=0 and B(256)=2^512-1. Summation by parts for the decreasing comparison
moment permits replacing A by B. The increments B(u)-B(u-1) are now masses
of a *comparison distribution*, not upper bounds on actual shell counts.
Each of its u positions survives independently with probability 7/8; hence
its mass at remaining support v is

    sum_{u>=v} (B(u)-B(u-1)) binomial(u,v) (7/8)^v (1/8)^(u-v).

`monotone_comparison.py` computes this rationally. It preserves the empty
outcome's active-group label and does not round a tiny empty mass up to one.
The positive mixture is then checked on every shell, including zero. This
argument uses linearity and zero initial state: an arbitrary affine offset
could introduce negative character coefficients. Although valid, the
implemented comparison gives no useful distance improvement on the tested
BCH instance, so it is not part of the current sparse certificate.
