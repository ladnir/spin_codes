# Reset timings and a routing-aware weaker reset

Update: [45% certificates](GROUPED_45_CERTIFICATE.md) now close group32 and
the simpler running-sum reset, with both uniform and two-diagonal MDS
feedback. The new complete inner is untimed; the reset-only measurements
below remain unchanged. The remaining group32 proof gap refers to49.1%.

Follow-up: the [joint count and paired-reset alternative](GROUPED_RESET_JOINT.md)
retains the cost of inactive interior coordinates. Group32 is still open;
group64 passes bounded selected-class screens but has no full certificate
or new timing. The measurements below are unchanged.

## What the timing pass settled

At K=2^20, the implemented MDS reset alone consumes the full 3--5ms
end-to-end target. This justifies testing a weaker reset before further
micro-optimization. It is a result about these kernels, not a lower bound
on every implementation of the certified construction.

Measurements process N=2^21 consecutive 128-bit output elements. They
include reading those elements, the reset transforms, and the persistent
state update. They exclude feedback, the outer, routing, setup, and
allocation. A 128-bit element is represented in four coordinates over
GF(2^32); multiplying each coordinate by the same subfield scalar is the
required scalar multiplication in GF(2^128). This does not replace uniform
GF(2^128) feedback with small-field feedback; feedback is simply not timed.

Host: Peach, Ryzen7950X, GCC15.2, `-O3 -std=c++20 -march=znver4`, CPU2.
Runs were serial, with three warmups. The first comparison used15 samples
per variant; control and small-group replays used31. Values below are
medians, or ranges of medians, not minima. Machine noise remains possible.

| Reset | RS[256,128], state512 | RS[512,256], state256 |
|---|---:|---:|
| Direct shifted transpose FFT, first pass | 11.467ms | 11.507ms |
| Gray, first pass | 6.515ms | 6.032ms |
| Fully normalized Gray, first pass | 6.698ms | 6.439ms |
| Partially normalized Gray, first pass | 6.934ms | 6.241ms |
| Gray,31-sample replays | 6.45--6.59ms | 5.90--6.24ms |
| Independent groups of32, two replays | 2.795--2.824ms | 2.797--2.820ms |
| Independent groups of8, two replays | 1.993--2.006ms | 2.012--2.064ms |

The first four variants retain the [two certificates](IDEAL_RESET_LADDER.md).
The two small-group variants do not. All timings are reset-only. In
particular, the2ms entry is not a2ms encoder or a certified operating point.
No outer-plus-route comparison was needed to reject the present MDS reset
implementation from the3--5ms total budget: reset alone already exceeds it.

The direct shifted FFT also reads about106MB of prepared coefficients;
the Gray variants reuse small tables. Its timing therefore includes that
real precomputed-table traffic, not just a comparison of arithmetic counts.
Normalization does not win here, despite reducing product sites. The
extra scaling and wider coefficients are not free.

`native/reset_native.cpp` uses the existing GFNI field arithmetic, SIMD
butterflies with compile-time sizes, fixed aligned scratch, and no hot-loop
allocations or runtime callbacks. Direct checks compare160 coefficient
fixtures with the scalar field oracle, compare selected transform columns
with polynomial evaluation at sizes8/32/256/512, and couple the coordinate
variants through complete reset streams. Grouped streams are also checked
against separate group streams. Optimized and UBSan builds pass. The
final source and optimized binary hashes are:

- source: `95f77f80d4d38050069acc30ed7eda5fd7e08624c8705c656a0c363702df8dc6`
- binary: `c5f588e29e8bcf31d6dc908e7ad3e37cca4cf242057074e7dcbcf7bc675d9f82`

## Cheaper construction: preserve memory, reduce reset coupling

Keep the outer, route, batch size b, and state size b of either profile.
Split each output batch and state into G=b/d contiguous groups of d
symbols. The candidates use d=8 or32. Each group maintains its own
distinct-point MDS syndrome, using the same time schedule for all groups.
Use the d-point reset recurrence from
[the coordinate construction](CHEAP_RESET_COORDINATES.md) independently
on each group. Retain one full uniform b-by-b feedback matrix per batch.

Thus the feedback still couples all groups; the encoder is not a direct
sum of independent codes. Setup is sampled once and fixed for all messages.
There is still one global routing stage and no new global permutation.
Each reset injection is invertible, so the full state-reset probability
from a nonzero state remains q^-b, where q=2^128.

The deterministic reset property is now weaker. During a closed excursion,
let a be the number of groups with any nonzero output. Each active group
has at least d+1 nonzero symbols, and its output support imposes rank d.
The whole closed support has rank da, not necessarily b. A single active
group permits a weight-(d+1) closed output segment. This deliberately
abandons the full-state MDS guarantee. It does not itself establish that
the inverse input lies in the routed outer code.

The previous small-memory obstruction is not bypassed by notation: the
state really retains b independent coordinates, not d. Uniform feedback
can activate all of them. The saving comes from local reset arithmetic,
not from claiming a small reachable state has become large.

## How the outer and routing can pay for weaker resets

Fix c excursions occupying v batches. Leave the final excursion open,
and put r=c-1. Fix the nonempty set of active groups in each closed
excursion, and let A be the sum of their cardinalities. The symbol A in
this note is a count, not the inverse proof's coordinate set.

Two consequences must be counted together:

1. Closed supports impose dA independent reset equations.
2. At the first batch of an excursion the state is zero, so X=Y. Every
   inactive output group is therefore an input hole at that activation.

Let S contain every noninitial excursion batch, the selected groups of
each closed excursion's first batch, and the entire final first batch.
Any compatible outer input is supported in S, and

    |S| = b(v-c+1)+dA.

The full excursion support had size bv. The activation holes remove
b(c-1)-dA coordinates from that support. This is especially valuable when
few groups are active, exactly where the reset rank is weak.

Batches divide regions in both profiles. A first-batch hole stays in
the starting endpoint's region. All partially selected regions of S
therefore belong to the at most2c endpoint regions. The existing regional
shortening moment applies to S despite the within-batch holes. It does
not require S to be a union of complete intervals.

For a fixed route, let D(S) be the dimension of the outer code shortened
to S. The forward argument counts this smaller source subspace and charges
q^[-b(c-1)] for the required live-state resets. Projective counting removes
the zero message and scalar duplicates, as in the retained proof. Its
regional moment must be evaluated at |S|, not bv.

The inverse argument instead charges the actual rank dA. Inactive groups
are zero throughout their closed excursions, each of which spans at
least two batches. The number of available output coordinates is at most

    E = b(v-2(c-1))+2dA.

This gives an output-volume bound on E coordinates with dA dimensions
removed. The inverse randomized coordinate set still has size b(v-c).
Use the smaller complete forward or inverse bound for each class; do not
multiply them as if their routing events were independent.

Group patterns cost something. For fixed A, their number is the coefficient
of z^A in ((1+z)^G-1)^r, bounded by binom(Gr,A). The current diagnostic
charges this factor in full, in addition to excursion placements. It also
enforces A>=r and (d+1)A+1<=W, where W=floor(.491N).

## What was checked, and what remains open

The coarse alternative of merely replacing reset rank b by d does not
close the complete numerical bound, even in the tested d64/d128 cases.
Those failures are failed upper bounds, not counterexamples. We therefore
added the activation-support restriction above rather than tuning that
coarse bound indefinitely.

For d32 and the larger-outer profile, the troublesome selected class
v4096,c64,A252 changes to a FLOAT log2 upper bound below-249.29. Sparse
activity uses the smaller outer support; dense activity uses larger reset
rank. That is a concrete benefit of using the construction's components
together, rather than requiring all output patterns to be MDS.

The broader screen does NOT close. Across105 selected classes for the
larger-outer/d32 profile, a class at v4224,c128,A254 remains near+851.80.
Across100 selected classes for original-outer/d32, a class at
v2176,c128,A1016 remains near+1443.94. The d8 critical examples are weaker
still. None of these FLOAT values is a directed certificate or a proved
lower bound on failure. No new whole-code distance certificate is claimed.

The remaining gap is intermediate activity. There are enough active-group
patterns to incur a large counting cost, but not enough reset rank to
recover the full-MDS bound. We also separately maximize output capacity
and route moments over episode lengths and placements. That can charge
the same class for incompatible worst cases.

Next test: retain the actual closed lengths and activation holes jointly
for a failing intermediate class. First compare its possible shortened
dimension with the required reset/output-zero constraints. Then determine
whether joint counting closes it. If it does not, consider inexpensive
mixing between groups rather than repeatedly tuning the same coarse bound.
Do not shrink the total state or change feedback distribution in this test.

Exact finite checks pass4096 grouped-support ranks,24 active-pattern
counts, and755 activation-support/capacity geometries. They include a
local weight-three reset at d2 to distinguish the weaker construction
from the full MDS control. These are local interface tests, not a K20 proof.

## Reproduction and limits

Compile `native/reset_native.cpp` with the native directory,
`spin/experiments/char2_field_spin`, and its `field128_lift` directory on
the include path. Use the GCC flags above. `--check` runs exact checks;
`--bench 256 1 31` selects the larger-outer profile's Gray reset.
Replace256 with512 for the original-outer profile. Modes0/1/2/3 select
direct/Gray/full-normalization/partial-normalization; modes4/5 select
uncertified groups8/32. Keep runs serial and pinned to CPU2.

Builds used a2GiB memory cap and60-second timeout. Runs used a1GiB cap
and60--90-second timeout. Only source and compact results were retained.
No raw timing dumps, local binaries, caches, production/paper edits,
commits, or promotions were made.

Run the exact Python tests with `python -B .../coded_state/test_striped_reset.py`.
The bounded FLOAT diagnostics are in `striped_reset_probe.py`; use
`--outer 512 --group 32 --screen` or `--outer 256 --group 32 --screen`.
The `--cover` option tests only the coarse rank-floor envelope and never
claims a directed certificate. Frozen helpers and prior proofs are unchanged.
