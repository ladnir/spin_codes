# Bounded diagnostics of the whole-state bridge

These floating-point checks evaluate the bound in
[the counting bridge](SMALL_STATE_RANK_BRIDGE.md), not the encoder's
minimum distance. They cover selected pair profiles only. No row is a
certificate for an entire source size or the whole code.

**Later checkpoint:** replacing the IID-conditioned route majorant with
exact regional coefficient extraction closes source counts g=7,...,16.
Those are separate directed checks, recorded in
[the exact-route gate](SMALL_STATE_EXACT_ROUTE_GATE.md). The failed scores
below remain useful diagnostics of the older majorant; they are not a
limitation of the construction or of the complete bridge interface.

The construction remains the historical 64-input/16-state geometry with
independent full-GF32 nonzero DHD diagonals. The target is K=2^20,
rate 1/2, and output weight at most 734003 (35% of N=2^21, rounded down).
The direct bounds for actual source groups 1 through 6 are separate;
their combined failure probability is below 2^-44.

## Corrected paired layout

The initial diagnostic incorrectly placed both members of a physical
pair in the same RS16 block. All numerical results from that version
are withdrawn. The construction and analytic bridge were unchanged.

The retained encoder stripes a pair across two different RS16 blocks.
For each half-batch, its occupancy polynomial is
`(p0+p1*(x+y)/2+p2*x*y)^16`. The two half-batches are independent under
the IID proposal. Conditioning on the exact pair histogram recovers the
actual regional distribution. The implementation sums ordered occupancies
before using the symmetric local bound.

The corrected schema names end in `v2-striped`. Tests check 15 small
oriented-pair polynomials, 36 ordered-matrix operations, 3 profile masses,
and all 289 cells of the singleton law. An explicit negative control
rejects the former same-block layout. Independent source review confirms
the corrected layout and positive matrix operations.

## Retained-marker probe

The probe borrows seven rational output markers from the older GF256
calculation. For each source size it tries three nearby fractional powers,
four zero-state penalties, and five positive IID tilts. These are parameter
choices only; the older certificate does not transfer to this operator.

There are 50 selected pair profiles and 600 profile/marker evaluations.
The source-set factor is the binomial overcount `binom(8192,g)`.
Both rows and columns of each regional matrix use the actual paired
distribution. Positive matrix powers retain their order. The final vector
is bounded by `M*ones`; no actual final refresh is charged.

The following numbers use the cutoff-source version of the bound.
For each profile, take its best tested parameter choice. The columns
give the minimum and maximum of these log2 upper-bound scores over
the selected profiles. A positive score is an inconclusive upper bound,
not a low-distance word.

| Actual source groups | Best selected profile | Worst selected profile |
|---:|---:|---:|
| 7 | 3352.50 | 10118.64 |
| 16 | 4833.64 | 18893.69 |
| 64 | 1293.19 | 4323.49 |
| 256 | 3214.49 | 8711.40 |
| 1024 | 900.20 | 19989.52 |
| 4096 | -62268.27 | -9540.32 |
| 8192 | -101712.54 | -101712.54 |

The run completed in 5.15 seconds with about 50 MiB peak resident memory.
This is checker runtime, not encoding performance. It wrote no files.
The dense cases justify further proof work; the sparse and intermediate
cases do not close under these choices.

## Joint-weighting comparison

The joint Jensen refinement in Section 8 of the bridge removes the
missing-source cutoff. It is proved and implemented as an optional mode;
the default cutoff calculation is preserved. Its identical-grid probe
also completed all 600 evaluations and all local checks. Peak resident
memory was about 50 MiB; checker runtime was 5.16 seconds.

| Actual source groups | Best selected profile | Worst selected profile |
|---:|---:|---:|
| 7 | 3374.86 | 10141.00 |
| 16 | 4882.04 | 18942.09 |
| 64 | 1336.88 | 4367.19 |
| 256 | 3214.37 | 8711.28 |
| 1024 | 899.92 | 19989.25 |
| 4096 | -62269.28 | -9541.12 |
| 8192 | -101715.98 | -101715.98 |

There is no qualitative change in this probe. Every selected best result
uses omega=1, so penalizing auxiliary zero states did not help within
the tested grid. Joint weighting slightly improves the dense scores but
increases the source prefactor at small sizes. Retain both analytic bounds;
do not expand this marker sweep merely to tune the small correction.

## Reproduction and scope

Run serially from the repository root, with NumPy available:

```text
python -B research/workstreams/one_route_coupled_scan/theory/test_small_state_bridge.py
python -B research/workstreams/one_route_coupled_scan/theory/probe_small_state_anchors.py --seconds 50 --memory-mib 1024
python -B research/workstreams/one_route_coupled_scan/theory/probe_small_state_anchors.py --source-mode joint --seconds 50 --memory-mib 1024
```

Use an external 60-second timeout and a 1-GiB address-space cap for the
probe. The helper also checks the free-space floor and its own resource
limits. It emits compact JSON to stdout, without creating caches or data.

The exact GF4/GF16 toy test checks 288 fixed zero/path existence bounds,
288 zero-source fibers, 35 retained-support cases, one empty-input path,
135 root-before-marker inequalities, 180 projective-count inequalities,
7172 exact source-stratum comparisons, and 35 joint-weighting inequalities.
It passes, but does not replace a directed finite K20 calculation.

Before a whole-code claim, cover every remaining source size and actual
pair profile with directed arithmetic. A failed IID-conditioned envelope
does not by itself identify whether route averaging, source-density
relaxation, or the auxiliary-state relaxation is responsible.

## Refinement proposed before exact-route closure

Retain the RS source constraints before replacing the outer distribution
by its density bound. Start with seven actual source groups, immediately
beyond the separately proved direct gates.

Fix an actual source set S and a whole-state path. A quiet-to-quiet batch
requires `C X_t=0`. For each RS column j, let R_j sum the ranks of these
constraints over its four regions. Their input coordinates are disjoint;
the invertible column mixer preserves their joint rank. Any 32 RS columns
are independently surjective. Thus the source constraints have rank at
least the sum of the largest 32 values among R_0,...,R_63. If this sum
reaches `128*|S|`, the path has no nonzero source and must be removed
before density domination. The direct sparse gates already retain this
cutoff; the present two-state diagnostic does not.

This is a targeted proof refinement, not a proposed change to the code.
It remains to compose the cutoff with the new whole-state argument and
evaluate its complete seven-group bound. No closure is asserted here.
The subsequent exact-route checkpoint closes seven groups without this
cutoff. Preserve the cutoff idea for later source classes rather than
repeating the failed seven-group rank-potential search.

The at-mean histogram divisor alone costs fewer than 1107 bits at g=7
and 3235 bits uniformly over all g, across 256 regions. These follow from
Stirling bounds for the trinomial probability at its mean. They are much
smaller than several failed scores. This observation does not bound the
entire benefit of exact coefficient extraction: tilting and off-histogram
terms incur additional losses. Exact small-g extraction remains a useful
control, but removing the divisor alone is not the main missing step.

## Complete dense-anchor control

The later [bounded anchor helper](theory/probe_small_state_dense_anchor.py)
covers every feasible matching count, rather than selected profiles.
It partitions the counts into at most 128 intervals. Each interval uses
a fixed positive rational IID proposal and the smaller endpoint histogram
probability. Log-concavity makes this a valid interval bound. Up to five
proposals can be combined by entrywise minima of their regional upper
matrices. The [bridge, Sections 9 and 12](SMALL_STATE_RANK_BRIDGE.md),
justifies this interface and raw-operator interpolation across source sizes.
These computations are FLOAT diagnostics, not new certificates.

At full support g8192, z=1/500 and alpha=1 give a log2 bound near -526601
with proposal pseudocount 1/65536. The default pseudocount 1/2 instead
retains enough off-histogram empty batches to make the bound fail. A tilt
of 256 also removes that artifact. This does not certify a dense interval.

At the inherited fractional marker, the complete g4096 anchor fails.
At z=370391/5732639,
alpha=2939512/14963605, pseudocount 1/65536, and tilts 1 and 5/3,
its score is about +14748.49. The worst interval is f=2040..2048 in
both matchings. The other two retained alpha choices and a g2048
interpolated marker also fail. No further broad sweep was run.

A discriminating check explains why the earlier selected-profile scores
were encouraging. At the same g4096 marker, the old five-tilt helper gives
about -9543.45 for (f_even,f_odd)=(2048,1024), but +29273.96 for the
realizable profile (2048,2047). A single occupied run of 4096 groups can
have the latter counts. It was absent from the earlier selected grid.
Changing its pruned matrix to the raw matrix costs only about 0.012 bits;
empty-input pruning is not the main loss at this anchor. The complete
interval helper's bound is better because it uses different proposal
parameters, but still fails.

### Separating coefficient slack from marker choice

At the same retained z and alpha, optimizing each regional entry's IID
tilt gives a one-run log2 bound near +4146.05 after charging exactly 8192
source sets. The new
[targeted coefficient control](theory/probe_small_state_cluster_route.py)
extracts the two regional coefficients (0,2048) and (2,2047) with positive
logarithmic convolution. Both profile orders are included. It improves
the one-run score to +1611.75, so exact conditioning alone does not close
this retained witness. Neither positive score is a low-distance word.

The coefficient control agrees with the shared engine at g=8 and g=16
within 9.4e-13 log2 units, and with a constant-matrix fixture within 1.8e-15.
The large control takes 0.51 seconds and 51 MB. It remains FLOAT arithmetic.
A conservative two-, four-, or eight-batch Jensen envelope gives no gain
at the retained witness; the original bound wins each entry.

The fractional parameter was not yet bracketed. Holding z fixed and
lowering alpha from 2939512/14963605 to 3/20 changes the optimized one-run
score to approximately -16934.91. More importantly, a complete 121-cell
anchor with the single tilt 3/2 gives -3247.12, including all source sets
at g=4096. Thus this failed witness did not require an encoder change.

The subsequent directed replay closes the full source interval
4096..4108 below 2^-100; see the
[dense gate](SMALL_STATE_EXACT_ROUTE_GATE.md#a-directed-dense-interval-without-changing-the-encoder).
It complements, rather than fills the gap after, the g<=128 checkpoint.
At g=2048, z=1/5, alpha=1/10 and 1/20 still fail in bounded controls. Retune
markers and inspect the responsible profiles there before claiming
interpolation or adding construction work.

## Intermediate anchors and two-batch grouping

The g2048 diagnosis separates three sources of proof slack: the output
marker, conditioning on the exact regional histogram, and splitting
state paths before taking a fractional power. None changes the encoder.

At z=23/100 and alpha=1/10, exact regional coefficient extraction gives
a one-run score near -118.90, counting all 8192 source sets in that class.
The corresponding complete cell census still fails. Adding exact
coefficients for the nine highest full-pair counts and IID tilts
1, 4/3, and 3/2 gives +3386.50. Its dominant run counts are near 120;
the passing one-run example is not a complete-size certificate.
The retained controls are
[exact tails](theory/probe_small_state_exact_tail.py) and
[complete cell census](theory/probe_small_state_census_anchor.py).

Grouping two ordered batches before the fractional power retains the
sum over their intermediate state. On the full-pairs-only one-run
control, this changes the score to -3478.61, even after conservatively
deleting the two singleton inputs in the other regional matching.
The actual source-size factor remains g2048 after deletion.

The [grouped control](theory/probe_small_state_grouped_route.py) implements
the ratio-bin bound in bridge Section 14 for all 4845 local occupancies.
At the same z and alpha, a complete 157-cell census with tilt 4/3 and
32 singleton tail cells gives +13.19 with bin width 1, but -106.34 with
width 1/4. Thus finer proof bins close this FLOAT diagnostic without
changing the local operations or route.

A short marker bracket improves the complete bound further. With
z=23/100, alpha=19/200, tilt=33/25, pseudocount=1/65536, and width 1/4,
the g2048 score is -423.37. Averaged raw monotonicity extends this
FLOAT bound to 2048..2051 with combined score -125.59. With width 1,
the shorter interval 2048..2050 has score -92.12. These are diagnostic
values, not directed certificates; any accepted replay is recorded
separately in the exact-route gate.

The complete grouped calculation takes about 1.5 seconds and 51 MB.
These are proof-evaluation costs, not encoding times. A direct control
enumerates all 23,474,025 ordered local-profile pairs for one proposal.
Every binned entry exceeds the corresponding direct grouped value,
by 0.00155..0.00349 in log2 units. At alpha=1, the grouped operator
matches the previous one-batch square within 1.43e-13. All 1809 exact
ratio-bin boundary checks pass. Earlier small-g controls also agree
with exact-tail extraction within 9.4e-13, and the g8 complete census
agrees with direct profile enumeration. These FLOAT checks support
implementation correctness but do not replace directed arithmetic.

## Actual-size reuse and longer windows

Direct conditioning at each actual size retains a route improvement that
averaged-subset transport discards. The
[direct interval control](theory/probe_small_state_grouped_interval.py)
fixes one IID proposal family for each parity and caches every ordered
unnormalized regional product. It then recomputes complete matching-count
cells, endpoint histogram divisors, and the source census at every size.
The retained 2048 witness covers 2048..2073 in this calculation. Its
directed replay is recorded separately in the exact-route gate.
Keeping those same proposals through 2074 fails, because the extended
lowest cell receives an increasingly poor histogram bound. This is
proposal-reuse slack, not a demonstrated bad codeword.

At g1024 the inherited marker gives a complete two-batch score near
+13570.97. A bounded bracket improves it to +2505.27 at z=1/2,
alpha=3/40, tilt=6/5, pseudocount=1/65536, and ratio-bin width1/4.
Multiplying the singleton proposal weight by additional factors does
not help the preceding alpha=2/25 witness. The next change therefore
retains longer state paths, rather than continuing that proposal sweep.

The [window control](theory/probe_small_state_window_route.py) represents
a positive row by its scale and coordinate ratio. It rounds that ratio
upward on a rational geometric grid. Each occupancy sends each grid
state to one next state; the calculation does not split its quiet/live
paths before the fractional power. All 4845 occupancies remain present.
Bridge Section 16 proves the positive majorant. FLOAT ceiling assignments
and arithmetic make these diagnostic values, not certificates.

At the preceding g1024 markers, using grid ratio9/8, all 153 matching-count
cells and all 1024 source-run classes are covered:

| Consecutive batches retained together | Complete log2 bound |
|---|---:|
| 2, exact ratio-bin grouping | +2505.27 |
| 4, projective grid | +812.67 |
| 8, projective grid | -35.57 |
| 16, projective grid | -450.52 |
| 64, projective grid | -758.88 |

The 16-step calculation uses 762 grid states, takes 3.48 seconds, and uses
183 MB. These are proof-computation costs, not encoder timings. It leaves
enough diagnostic slack to justify directed replay. The source hash is
`0a7f941e8c0e94c937bc7f5982622eab1ac85a93da9dec00baaf313da0e9ac8f`.

Reproduce the retained diagnostic, under the usual external limits:

```text
python -B theory/probe_small_state_window_route.py --groups 1024 --z 1/2 --alpha 3/40 --tilt 6/5 --window 16 --rho 9/8 --seconds 20
```

Bounded fixtures compare the helper against direct ordered products.
A rank-one matrix with an exactly represented invariant ratio agrees
through 32 steps within 3.82e-12 relative error; repeated FLOAT sums account
for that discrepancy. The initial 1e-12 fixture tolerance was too strict.
For a noncommuting two-matrix mixture, exhaustive four- and eight-step
enumeration lies below the projective bound and above that bound divided
by rho^(alpha*(window-1)). This also holds at alpha=1. These checks test
normalization and chronology, not rigorous numerical enclosure.

A [common-functional control](theory/probe_small_state_functional_route.py)
tries to replace the ratio distribution by one weighted linear functional.
Its concave scalar moment has a global tangent upper bound, but the complete
g1024 scores stay near +12582 for log2 weights 64..96. The worst scalar ratio
can recur at every step of this bound, even when the actual ratio evolution
would not keep returning there. Preserve the failed simplification; it is
not evidence that extra arithmetic is needed in the construction.

Next: replay the 16-step positive calculation with directed arithmetic.
Do not promote its negative FLOAT score to coverage or extrapolate it to
other source sizes. No construction change or new native timing is claimed.
