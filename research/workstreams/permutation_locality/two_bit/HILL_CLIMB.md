# Hill climbing the first complete two-bit proof

The baseline was a complete 6% distance guarantee at K=2^20, with setup
failure probability below 2^-49.11. The same encoder now closes through 9.25%.
The 8%, 8.5%, 9%, and 9.25% dense certificates passed independent 256-bit replay.
Their whole-code aggregates with the prior sparse lemmas remain below
2^-49.11. See [FIRST_CLOSURE.md](FIRST_CLOSURE.md).

## Reuse proofs without inheriting an invalid claim

All three dense drivers now support `--retarget`. They reconstruct the
saved partition, recompute each old leaf at the requested cutoff, and
requeue leaves that miss the new budget. They do not trust saved scores
or upper bounds. Every requeued leaf remains part of the unresolved
partition until a new witness or subdivision covers it. Ordinary replay
still forbids changing the saved cutoff.

These runs retained all encoder parameters and closed 8.5% and 9%:

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-401-d08.json --threshold 178257 --max-cells 1500 --max-depth 48 --output tmp/two-bit-affine-dense-401-d085.json
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-401-d08.json --threshold 188743 --max-cells 2500 --max-depth 48 --output tmp/two-bit-affine-dense-401-d09.json
```

The 8.5% run completed after 854 new cells, with 820 terminal leaves and
no unresolved cells. The dense margin is above 74.62992856911 bits. Its
256-bit replay and aggregate give a full-code margin above 49.11287660853
bits. The 9% run completed after 2256 new cells, with 1453 terminal leaves.
Its dense margin is above 73.84716062499 bits. Independent 256-bit replay
and aggregation give a full-code margin above 49.11287658688 bits.
The 9% witness SHA-256 is
`57f31c6562be34f907994cbc9b470d7a7f99edd17b565b356f41e6584c53e15a`.

The unchanged-encoder 9.25% run completed after 2302 new visits, with
2235 terminal leaves and no unresolved cells. Independent 256-bit replay
gives dense margin above 71.23535577341 bits and full-code margin above
49.11287632260 bits. The witness SHA-256 is
`ff1844642df892daa79a254dc3a1c67f0fdbb213ad6e9eccd6de6a46e8324058`.
Its cutoff is 193986:

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-401-d09.json --threshold 193986 --max-cells 2500 --max-depth 51 --output tmp/two-bit-affine-dense-401-d0925.json
python -B research/workstreams/permutation_locality/two_bit/closure.py --dense tmp/two-bit-affine-dense-401-d0925.json --precision 256
```

## What limits the next climb

The first 9.5% run stopped at a depth-48 cell with log2 upper
-42.12689774975, short of its 70-bit per-cell budget. Increasing the depth
to 60 and lowering that budget to 60 left a finer cell at -61.9476258409
in the floating proposal, just short of the search's extra two-bit buffer.
Neither run covered the full domain; neither is a 9.5% certificate.

The initial output tilt was optimized for an earlier, uncoupled bound.
Optimizing it against the final coupled bound improves the first cell to
outward log2 upper -43.57222897332. The finer cell improves to
-63.66968558744. This modest improvement justifies a bounded refinement,
not a claim that witness tuning removes the remaining gap.
The optional `--polish-tilt` enables that search. `--reuse-witnesses`
tries each parent's rational witness on its children before optimizing
new witnesses; every accepted reuse receives a fresh outward bound.

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --resume tmp/two-bit-affine-dense-401-d095-deeper.json --max-cells 750 --max-depth 60 --target-bits 55 --polish-tilt --reuse-witnesses --output tmp/two-bit-affine-dense-401-d095-refined.json
```

An accepted-cell bound says nothing about the remaining cells. The
subdivision limit is a limit of this proof search, not evidence of an
actual low-distance codeword.

The bounded refinement still stopped, at a nearby cell with proposal
-52.47333641556. Inspecting this cell exposed an integer-count obstruction
in the relaxation. Its central-row fraction contracts to zero. Its double
category coordinate m2 satisfies

    5/32 <= G*m2 <= 41/256.

For a composition with no central rows, the positive contributions that
fit below the upper endpoint are 1/13 and 6/43. Thus the only attainable
sums below that endpoint are 0, 1/13, 6/43, and 2/13. All are below 5/32.
This cell contains no integer composition, despite its fractional feasibility.

The new exact projection check removes such cells. It first excludes
component types whose single contribution exceeds any cell upper bound.
For each remaining coordinate, it enumerates nonnegative integer sums
of the positive component increments, using at most G terms. When the
coordinate cap permits more than 12 positive terms, or enumeration exceeds
20,000 nodes, the check returns inconclusive and retains the cell.
Only exhaustive rational infeasibility permits discarding a cell.
Ignoring the other coordinates and the lower occupancy constraint makes
this a safe outer relaxation of integer feasibility.

Small exhaustive tests retain every actual composition and check interval
endpoints and enumeration-budget exhaustion. A 9.5% rerun used this check;
the earlier 9.25% and three-update pilots were already running and used
the weaker pruning from their start.

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --resume tmp/two-bit-affine-dense-401-d095-refined.json --max-cells 2000 --max-depth 54 --target-bits 60 --polish-tilt --reuse-witnesses --output tmp/two-bit-affine-dense-401-d095-integer.json
```

That rerun removed the original obstruction but stopped at another cell
after 19 new visits. The next check retains the coordinates jointly near
the least-active mixture vertex. Write a=min_i(f_i,1+f_i,2) over active
types. Here a=1/7, attained by the light-only pair type. Every other active
type has positive excess c_i=f_i,1+f_i,2-a. For q>=q_min,

    sum_i n_i*c_i <= G*(upper(m1)+upper(m2))-q_min*a.

When this permits at most ten nonminimal pairs, enumerate their integer
counts within that exact excess budget. For each choice, the remaining
light-type count must lie in the intersection of three integer intervals,
as well as [q_min-exceptional_count, G-exceptional_count]. Empty intersection
for every choice proves the whole cell empty. Inactive types have feature
zero. If the least-active feature is not unique, or the search budget is
exceeded, retain the cell. This check does not approximate the counting
measure or remove a low-probability event.

Both searches scale rational increments to integers before enumeration.
Their inner loops use integer arithmetic and bounded work. Exhaustive
small-composition tests check both feasible points and empty boxes. The
joint-boundary rerun is:

```sh
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --resume tmp/two-bit-affine-dense-401-d095-integer.json --max-cells 2500 --max-depth 54 --target-bits 60 --polish-tilt --reuse-witnesses --output tmp/two-bit-affine-dense-401-d095-boundary.json
```

The joint-boundary rerun stopped after 26 new visits, with 185 terminal
leaves and 702 unresolved cells. Thus exact integer pruning removes real
relaxation artifacts but has not closed 9.5%. A fresh 256-bit replay of
the complete 9% witness with both integer checks reproduces its original
aggregate margin, 49.11287658688 bits. The separate 9.25% two-update run
subsequently completed. The three-update 10% pilot exhausted its work
budget with a saved partial cover, as recorded below.

The remaining depth-limited cell has floating log2 upper +67.01823 and
really contains an integer comparison composition: 400 light-only pairs,
one heavy-heavy pair, and 3695 inactive pairs. Its mean coordinates are
((400/7+12/37)/4096, (9/37)/4096, 0). This is a comparison-measure
composition, not an exhibited BCH message or a low-distance codeword.
The next useful diagnosis is its exact-composition bound, followed by
adjusting the sparse/dense split if needed; integer pruning alone cannot
discard this obstruction.

That diagnosis now succeeds. The exact composition with a constant capped
shuffle factor still gives log2 upper +311.91237 at 9.5%. Keeping the
profile-dependent factor inside the input moment instead gives below
-181.66464. Thus the failed affine cell was not an intrinsic obstruction
even within the comparison measure.

The new [exact-composition branch](EXACT_COMPOSITIONS.md) reconstructs all
integer compositions in a small boundary cell and verifies their complete
sum. It rejects truncated lists. The first automated proposal rounded a
rare category's reference count to one and missed the useful witness;
trying zero as well repairs this search error. The next retry stopped at
a cell containing one central row. Bounding central and noncentral
exception counts separately proved that cell empty. A subsequent retry
reached 222 accepted leaves but retained 704 unresolved cells; it is not
a full certificate.

The next real composition has 398 light-only pairs, three light-light
pairs, and 3695 inactive pairs. Its first exact bound was about 2^-50.04,
short of the search's 60-bit per-cell budget. This particular miss might
have been repaired by reallocating that budget. However, nearby integer
compositions exposed a larger gap.

Let b count light-light pairs, with 401-b light-only pairs and 3695
inactive pairs. Each row of this table is a separate, complete composition
bound at cutoff 199229 (9.5%). The second column uses reference counts
estimated from the input weights. The third retunes those counts using
logarithmic derivatives of the low-output tilted envelope. All reported
bounds were recomputed outward at 192-bit precision.

| b | Initial log2 upper | Retuned log2 upper |
|---:|---:|---:|
| 4 | -14.69414 | -120.49443 |
| 6 | +65.15035 | -41.46219 |
| 10 | +225.90685 | +31.71934 |
| 24 | +274.78060 | +86.75981 |

The retuned witnesses are still integer reference counts. Numerical
derivatives only propose them; the same pointwise shuffle inequality and
outward verifier certify the result. The automatic exact-composition
branch now tries this refinement, with a bounded three-iteration search.

Thus the 9.5% gap is not just the conservative per-cell budget or a
fractional-composition artifact. Further subdivision cannot improve an
already fixed composition. This does not establish a low-distance
codeword, nor rule out a stronger bound for the same encoder. It favors
assessing a different proof split or the three-update alternative over
indefinitely subdividing these cells. The 9.25% checkpoint is complete.

## Move the boundary between the two proof methods

The split at q=401 is a proof parameter, not an encoder parameter. The
dense method could start later if the sparse method covers the additional
occupancies. At cutoff 199229, the same family with 24 light-light pairs
has these outward log2 bounds after reference-count tuning:

| Active pairs q | Light-only pairs | Log2 upper |
|---:|---:|---:|
| 401 | 377 | +86.75981 |
| 405 | 381 | -14.96833 |
| 409 | 385 | -116.99223 |
| 417 | 393 | -321.88783 |

Each composition also has 4096-q inactive pairs. These selected points
motivate a split at 409, but do not cover the dense range. The witness
for the q=409 row can be reproduced directly:

```sh
python -B research/workstreams/permutation_locality/two_bit/mixture_probe.py --theta 2/5 --threshold 199229 --composition 0,0=3687 0,3=385 3,3=24 --density-anchors 4036 58 2 --witness 45847339/500000000 239687269/1000000000 283603509/1000000000 --outward
```

The sparse driver now records its output cutoff and applies it in both
the proposal and outward Chernoff factors. Legacy records mean 209715.
Ordinary replay preserves the saved cutoff; only explicit retargeting
can change it. The local operators and the setup distribution do not
change when only the cutoff changes.

The affine driver now permits raising the minimum occupancy during an
explicit retarget. It recomputes every retained witness on the restricted
domain. Lowering that threshold requires a new cover; an old empty cell
could then contain newly admitted compositions. The aggregate verifier
still rejects a sparse/dense gap: a q>=409 dense result cannot be joined
to the old q<=400 lemmas alone.

The two complementary attempts are:

```powershell
$twoBitOccupancies = 401..408
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --threshold 199229 --updates 2 --tilts .088 .096 .104 .112 .12 --max-splits 1 --density-fold --output-degree 64 --precision 256 --certificate-dir tmp/two-bit-sparse-401-408-d095
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-401-d0925.json --minimum-groups 409 --threshold 199229 --max-cells 1000 --max-depth 54 --target-bits 60 --exact-boundary --reuse-witnesses --output tmp/two-bit-affine-dense-409-d095.json
```

The sparse extension passed all eight complete support covers at 256-bit
precision. Each used one symmetric split. Their aggregate margin exceeds
1491.07046549378 bits. The dense attempt stopped after 115 new visits,
with 1173 terminal leaves and 1119 unresolved cells.

Its depth-limited cell has log2 proposal -54.58654, but contains no
integer composition. The old ten-pair enumeration limit returned
inconclusive. Raising that limit to sixteen proves the cell empty while
retaining the 20,000-node cap and the fail-closed budget behavior.
A regression test checks this exact cell, including the old inconclusive
result.

The sparse extension has enough slack to test a later proof boundary.
The next attempt covers q=409,...,424 sparsely and retargets the saved
dense partition to q>=425. It does not change the encoder:

```powershell
$twoBitOccupancies = 409..424
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --threshold 199229 --updates 2 --tilts .096 .104 .112 .12 --max-splits 1 --density-fold --output-degree 64 --precision 256 --certificate-dir tmp/two-bit-sparse-409-424-d095
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-409-d095.json --minimum-groups 425 --threshold 199229 --max-cells 2200 --max-depth 54 --target-bits 60 --exact-boundary --reuse-witnesses --output tmp/two-bit-affine-dense-425-d095.json
```

Checkpoint inspected on 2026-09-28: all sixteen q=409,...,424 sparse
records contain complete support partitions. Their saved aggregate
margin exceeds 206.24224259823 bits. This inspection checked their scope
and partitions; it did not independently rerun their numerical bounds.
The q>=425 dense run stopped at its 2200-cell budget with 2733 accepted
leaves and 199 unresolved cells. No depth-limit failure was recorded.
The files remain in the ignored paths above, and no search is running.

The 9.5% effort is paused in favor of joint route--inner tuning.
Resuming it requires finishing the dense cover, independently replaying
the sparse extension and dense witness, and checking their aggregate.
These attempts do not yet constitute a 9.5% certificate. The earlier
complete operating points remain unchanged. The aggregate driver now
accepts explicit `--bridge` witnesses. It checks every missing occupancy
exactly once, requires matching encoder parameters, and checks cutoff
inclusion. It reconstructs local operators and replays each complete
support partition; cached upper bounds are ignored. The final sum must
meet the requested margin, which defaults to 40 bits.

## Extract the slack already present in a complete certificate

For a fixed cell witness, the cutoff d occurs only in exp(lambda*d).
The new `certificate_slack.py` replays its remaining coefficient at d=0.
It then sums these positive exponential terms and the earlier sparse
bounds, searching integer d until the aggregate reaches a requested margin.
No fresh cells, numerical derivatives, or cached upper bounds are trusted.

Before the new integer pruning, the 9% partition supported cutoff 188784
with a recomputed aggregate margin above 40.90707 bits. Independent direct
replay at that cutoff confirms the aggregate. This increases the cutoff by only 41 output bits.
The resulting fraction is about 9.00192%; the single-coordinate integer
check leaves this fixed-witness frontier unchanged. New witnesses and subdivisions
are needed for a meaningful distance improvement. This is the frontier
of those fixed witnesses, not the code's distance or the best possible bound.

```sh
python -B research/workstreams/permutation_locality/two_bit/certificate_slack.py --dense tmp/two-bit-affine-dense-401-d09.json --target-bits 40 --output tmp/two-bit-affine-dense-401-d09-frontier40.json
python -B research/workstreams/permutation_locality/two_bit/closure.py --dense tmp/two-bit-affine-dense-401-d09-frontier40.json --precision 256
```

## Where the current 10% bound is weak

Selected-composition checks put the difficult case near the lower end
of the dense range, where many pairs contain one light comparison row.
At q=401 and theta=2/5, the affine bound at that light-only endpoint gives
log2 upper +1113.28294 at the 10% target. This is a failed upper bound,
not a lower bound on the actual code's failure probability.

Changing the comparison bias alone did not repair that endpoint:

| Comparison bias theta | Outward log2 upper at 10% |
|---:|---:|
| 1/3 | +2344.41578 |
| 0.36 | +1480.12185 |
| 0.38 | +1203.44500 |
| 0.40 | +1113.28295 |

These are different positive comparison measures for the same code.
They do not describe different outer distributions. Other selected
compositions had substantial slack: the direct composition calculation
for two central rows per active pair at q=401 gave log2 upper below
-6469.27978, even using the coarser universal shuffle loss.
Selected compositions never replace a full mixture cover.

## One additional update is a concrete alternative

The iid-reference kernel now permits an explicit update count. After r
independent transvections, its lazy probability is
2^-r and its remaining mass refreshes uniformly among nonzero states.
The fixed-maps, output-before-update convention is unchanged. Small exact
state-chain tests check the envelope for r=1,2,3,4.

At the same q=401 endpoint and 10% target, outward bounds are:

| Updates per step | Log2 upper |
|---:|---:|
| 2 | +1113.28295 |
| 3 | -646.03718 |
| 4 | -2218.15262 |
| 6 | -4326.90432 |

```sh
python -B research/workstreams/permutation_locality/two_bit/endpoint_probe.py --theta 1/3 9/25 19/50 2/5 --updates 2
python -B research/workstreams/permutation_locality/two_bit/endpoint_probe.py --theta 2/5 --updates 2 3 4 6
```

Thus three updates repair this particular bottleneck without changing the
state size. This is not a whole-code three-update certificate. The existing
sparse lemmas concern two updates and cannot be reused merely by assuming
that more mixing improves distance. A three-update attempt must rederive
or recheck the sparse operators, cover all occupancies, and aggregate again.
The dense and sparse drivers now accept `--updates`. New dense records
store the count; legacy records mean two. Sparse records include it in
their ensemble identifier. Replay and resume reject a changed count.
Explicit retargeting recomputes every retained witness for the new inner.
The aggregate driver still rejects three-update dense records because
its documented sparse lemmas cover only two updates.

The sparse local operators now accept an explicit update count as well.
Their lazy and uniform-refresh terms use 2^-r and 1-2^-r. The translated
density cache was computed with lazy factor 1/4, so those cached terms
are rescaled by 2^(2-r); the underlying character and overlap counts do
not depend on r.

A fresh 256-bit one-pair run at r=3 and 10% distance gives margin
49.1911704615769 bits. It covers every nonzero pair message, every support,
and all 4096 pair locations. This result replaces, rather than assumes,
the old q=1 lemma for the three-update candidate. A separate complete q=2
calculation gives 83.12117219910 bits, using census cutoff two and tilts
.001, .0025, and .0064. Its saved witness was also retargeted back to two
updates: fresh replay gives 81.58829858595 bits, demonstrating that the
driver changes the actual operators rather than only their metadata.
These q=2 bounds use a small test grid, not the earlier optimized grid.
The subsequent q=3,...,32 batch passed every complete support cover at
256-bit precision, with no support split. Its aggregate margin exceeds
145.73784546884 bits. Thus q=1,...,32 now has three-update coverage;
q>=33 still requires coverage before a whole-code claim.

```sh
python -B research/workstreams/permutation_locality/two_bit/probe.py --one-group --updates 3 --precision 256
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies 2 --updates 3 --cutoff 2 --tilts .001 .0025 .0064 --max-splits 0 --precision 256 --certificate-dir tmp/two-bit-r3-covers
python -B research/workstreams/permutation_locality/two_bit/affine_mixture.py --retarget tmp/two-bit-affine-dense-401-d085.json --threshold 209715 --updates 3 --max-cells 1800 --max-depth 48 --output tmp/two-bit-affine-dense-401-r3-d10.json
```

The complete-domain three-update dense pilot stopped at its 1800-cell
budget with 1036 terminal leaves and 633 unresolved cells. No depth-limit
failure was reached. This is saved progress, not a complete dense proof
or evidence that the construction fails. If that cover later closes,
recheck the remaining sparse range rather than assuming monotonicity in r.
There is no new performance measurement for that variant yet.

```powershell
$twoBitOccupancies = 3..32
python -B research/workstreams/permutation_locality/two_bit/full_cover.py --occupancies $twoBitOccupancies --updates 3 --tilts .00064 .001 .0016 .0025 .004 .0064 .01 .016 .025 --max-splits 5 --precision 256 --certificate-dir tmp/two-bit-r3-covers-3-32
```

## An optional geometric refinement

[Clipping impossible mixture vertices](CLIPPED_DOMAIN.md) is implemented
and tested, but remains opt-in with `--clip-domain`. Its first diagnostic
was slightly worse than the original proposal, so the running full covers
retain the established rectangular search. No production code or original
four-bit proof was changed. Generated witnesses remain in ignored `tmp/`.

There are 89 passing tests, including safe retargeting, rational polytope
vertices, the clipped convex-moment bound, ensemble guards, witness reuse,
the fixed-witness cutoff search, exact integer pruning, and complete
composition enumeration and replay. The new density-anchor tests cover
serialized rational witnesses and categories with zero comparison mass.
Sparse-cutoff tests check the exact exponential factor, serialized scope,
and unchanged legacy defaults. Occupancy-retarget tests reject silent
scope changes and expansion of a previously covered domain.
Bridge tests reject gaps, duplicates, incomplete support partitions,
incompatible update counts, and insufficient cutoffs. A shared-driver
smoke test recomputes the q=2 bound despite a deliberately corrupted
cached upper bound.
