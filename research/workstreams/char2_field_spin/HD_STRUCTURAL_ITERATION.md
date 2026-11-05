# Structural optimization of the grouped GF256 code

## Objective and retained reference

Iterate on the construction, its proof, and its implementation together.
The current target is K=2^20, rate 1/2, at least 35% symbol distance, and
at least 40 bits of whole-code setup-failure margin. Time complete forward
encoding of 128-bit records with setup and allocation outside the clock.
Seek a repeatable material reduction from about 4.0 ms, rather than another
tile-size adjustment. Do not replace a certified checkpoint with an
uncertified timing candidate.

The retained [grouped-outer checkpoint](HD_GROUPED_OUTER.md) uses four
RS[64,32] rows, independent invertible Toeplitz column mixers over GF256,
256 independently permuted regions, and the h4 common-refresh inner.
Its complete certificate gives 35% distance and 72.3975395430 bits of
margin. The frozen native binary has SHA256
`6c0f42d94ef06b987dec159f092ea30ee286b0544f15e757fdf502fb43154264`.
Its output checksum at K20 is 3966879751720036268.

No production profile, paper result, or existing certificate is changed.
The experimental selectors default to the previous construction. Run
benchmarks serially, under both benchmark locks on Peach CPU 2. Another
user's CPU 0 workload shares L3; do not interfere with it. Compile and run
correctness checks before, not during, timed comparisons.

The latest [closed paired-route checkpoint](HD_SHIFTED_CLOSED_CERTIFICATE.md)
now gives 35% distance with 55.0927 bits of margin. A matched three-way
comparison measures 3.798 ms versus 4.049 ms for the reference and 3.871 ms
for its layout-only refinement. The previous checkpoint remains preserved.

### Latest continuation: exact-constraint controls remain far from closure

This structural-search pass ends at the certified paired-route checkpoint.
Its matched complete-encoding time is 3.798160 ms versus 4.048848 ms,
a 6.19% reduction. The unchanged target is rate1/2, 35% symbol distance,
and at least40 bits of ideal-setup margin; the winner gives55.0927 bits.
The original72.3975-bit checkpoint remains preserved. This does not
close or promote the faster G16 candidate.

The preceding status-only turn made no research progress. This continuation
executed both pending controls and changed the next action: neither justifies
another witness grid or a larger version of the same relaxed state model.
The complete 3.798-ms/55.0927-bit checkpoint and experimental
3.397943-ms binary are unchanged. The latter remains uncertified.

The [fixed-source partition](HD_TWO_SUM_GF16_SOURCE_STRATA.md) retains
the maximum full-pair count F of each indexed source set. Its 320 controls
reduce the complete source unions at r=18,24,32 only to +2697.88,
+3785.56, and +6035.83. One permitted full pair already restores the
problematic small-occupancy intercept. These are failed upper bounds,
not low-weight codewords.

Two new six-control evaluations use only F=1 at those same three weights,
with two frozen marker templates each. They replay the old baseline before
adding a constraint, rebuild no regional providers, and fit no parameters.
The [positive coverage recursion](HD_TWO_SUM_GF16_STATE_COVER.md) restores
the exact off-diagonal zeros at empty occupancy. It improves the old-zero
bounds by 0.2641, 0.0000085, and zero bits. Its additional Minkowski
relaxation is controlled by taking a pathwise minimum with the old bound.

The [positive activity recursion](HD_TWO_SUM_GF16_GROUP_SEEN.md) instead
excludes wholly empty groups or imposes at least 9r active group-column
incidences. These are independent necessary conditions, not their joint
intersection or each group's individual nine-column floor. The aggregate
constraint improves r=18 by 307.95 bits to +2449.84. Neither condition
improves r=24 or32. The activity patterns remain fixed through the type
Hölder step; a counterexample in the note explains why a coefficient tail
cannot simply be moved inside that norm.

The controls completed serially in 4.974 and 0.651 seconds. Their negative
results concern these envelopes and witnesses, not an impossibility theorem
for the construction. They do not combine into a whole-code certificate.
Further progress on this candidate needs a different way to retain
row/type correlations or a changed outer, not additional bins in these
same comparisons. All earlier partial proofs remain indexed and intact.

A separate [exact mixer probe](HD_LOWRANK_OUTER_PROBE.md) tests
H_d=I+V_d V_d^T for d=2,3,4. Extra moments preserve full invertibility
but introduce 256 singular principal subsets and weight-14 pair outputs.
The current rank2 matrix has neither defect. This is not a promising
drop-in proof repair, and no native work followed from it.

The fresh lightweight regression passes 257 two-sum tests and four
low-rank tests. Native source and the complete paired-route manifest
still match their frozen hashes. No timing run, production edit,
paper change, commit, push, or promotion occurred in this continuation.

The final [coverage checker](theory/hd_shifted_verified_cover.py) passes
at384-bit aggregation precision. Its new receipt,
`theory/hd_structural_final_coverage_check384.json`, is byte-identical to
the pinned384-bit receipt: SHA256
`446c2ce09488fcaf86c2070d86b8890d510a767a21cd134a2ef9c5c20672ba8e`.
It authenticates previously observed256/384-bit numerical replays and
freshly rebuilds exact profile coverage and directed aggregation. It
does not rerun every coefficient calculation. All11,465,830,401 feasible
profiles and every active-group weight1..8192 remain covered.

Read-only Peach checks confirm the retained reference binary SHA256
`6c0f42d94ef06b987dec159f092ea30ee286b0544f15e757fdf502fb43154264`,
paired-route winner `103c38e52809ab4a1030150280e55080dc1fdba844b4ee17a7ca95e791ed7312`,
and experimental G16 binary `239278430a28e0a1f6b5e563c64cb0eb7e8f9a03e3de34dfb06f5912f18d28d5`.
Fresh serial winner checks pass the30720 input/state basis directions,
pair masks and cyclic boundary cases, small dense-oracle matrix, and
K20 oracle with checksum7701295210980077712. No new timing is inferred
from those correctness runs. The published matched timing triples are
retained; the G16 remote timing log was also inspected and matches its pin.

The [larger-Toeplitz review](HD_LARGER_TOEPLITZ_COST.md) supplies another
negative structural result. Its proof interface transfers, but both the
recursive circuits and straightforward bilinear evaluation circuits cost
more than the current certified outer. Fixed transforms cannot disappear
merely by moving them across the row RS encoder. A new counted circuit
must pass the thresholds in that note before native work is justified.

The stopping evidence is therefore a complete measured winner, preserved
checkpoints, and several independently reviewed structural alternatives
that either remain uncertified or add work. This is not a claim that the
encoder is optimal or that G16 has low distance. Reopen G16 for a genuinely
different correlation-preserving argument, or the stronger-mixer route
for an explicit cheaper circuit. Another grid over the current relaxed
envelopes is not the recommended next step.

### Previous continuation: high/low arity and the quiet-state obstruction

The experimental binary is unchanged. A new
[high-total-arity comparison](HD_TWO_SUM_GF16_HIGH_ARITY_CERTIFICATE.md)
now has fresh directed receipts at three anchor weights. The classes
(r,A0)=(18,247),(24,317),(32,463), with total active-row count A>=A0,
have a combined 62.24877612169-bit bound. All 555 ordered-profile
ceilings agree at 256/384 bits. Independent 512- and 1024-bit rounded
unions agree. These are selected classes, not complete weight coverage.

The complementary [low-arity control](HD_TWO_SUM_GF16_LOW_ARITY.md)
fails at all three anchors. Its dominant term uses the quiet-to-quiet
state path and intermediate arities 8 or 9, already below the high-side
cutoffs. A new positive-polynomial calculation retains total arity
exactly, but therefore does not remove those dominant terms.

A separate provider composes two good packet matrices before taking
the fractional power. It is entrywise no worse than the previous
provider and agrees at alpha=1. Rebuilding the three banks and repeating
all 45 identical numerical witnesses leaves the selected bounds
unchanged: the same regional envelopes dominate.

The [joint input/output counter](HD_TWO_SUM_INPUT_WEIGHT.md) now retains
the stronger RS support floor: a active rows have at least 9a nonzero
pre-mixer symbols across the sixteen columns. Exact integer Mobius
inversion gives positive joint-weight coefficients. Fifteen bounded
saved-bank controls still fail with the existing regional envelopes.
This is a useful counting interface, not a new certificate.

The next structural loss is explicit. At r=18, the fractional regional
quiet-to-quiet entry has log2 value 0 at empty occupancy, approximately
-28.89 at occupancy 2, and -42.05 at occupancy 18. One exponential
envelope forced through the empty value loses roughly 28 bits per
nonempty region. The new [empty-separated control](HD_TWO_SUM_GF16_EMPTY_CUT.md)
keeps that branch separate while preserving the two inner states and
fixed types across columns. Twenty-four fixed controls reduce the best
bounds to +2700.251233, +3789.765729, and +6039.627030 at the three
anchors, but none closes. At r=24/32, the dominant term now alternates
activation and return eight times rather than staying in state zero.
The next structural proof target is dependence in those paths, not
another fit of the same regional slopes. This differs from the earlier
scalar common-potential cut control.

Next priority: retain the fixed source's matching stratum instead of
maximizing every region over all full-pair counts, and preserve the exact
zero off-diagonal entries at empty occupancy. These are untested
structural proof changes. Repeating the current witness grids is deferred.

This iteration is classified as progress: three new directed high-arity
classes close, exact joint support counts are available, and a measured
proof-envelope loss is removed. Whole-code coverage has not advanced
to a new complete certificate, and no new speedup is claimed this turn.

The full lightweight regression currently passes 220 tests. The native
outer, encoder, and closed paired-route manifest hashes match the
retained pins. No native timing run, production change, or promotion is
part of this continuation.

### Previous continuation: expand proof coverage for the 3.398-ms candidate

The preceding goal turn made progress: it changed the experimental outer
and its dataflow, measured a repeatable reduction, and derived the exact
closed support kernel. This turn keeps that binary unchanged and tests
whether its proof scales beyond a selected source arrangement.

Fresh directed replay closes the fourteen one-dense-group arities at
r=17/source `(8,0)`: the rounded-class union is
-585.542121955663 bits. All arity ceilings agree at 256/384 bits; separate
512- and 1024-bit sums agree. This is a scoped certificate, not the
construction's complete margin.

More importantly, a new source-independent envelope covers all fixed
source supports for r=2..17 with zero or one dense group, plus both sparse
r=1 classes. Its 322 cells agree at 256/384 bits, with a rounded-cell
union of -563.328194368258. The exact combined dense moment, source
count, row choices, and support floors are rebuilt, not transferred from
the older candidate. The [small-cover record](HD_TWO_SUM_GF16_SMALL_COVER.md)
and [directed receipt](HD_TWO_SUM_GF16_SMALL_CERTIFICATE.md) specify scope.

The two-dense-group weight DP also passes all 105 separate arity pairs at
r=17/source `(8,0)`, with FLOAT union -476.62769692618843. A deliberately
merged-arity control instead gives +258.67213841832563. Empty columns
prefer small arity while occupied columns prefer large arity in that
relaxation; the actual fixed row set cannot switch between them. Both
results are preserved in the [two-group record](HD_TWO_SUM_GF16_DENSE_TWO.md).

The remaining one-group dense class initially fails a common-potential
bound. Keeping the complete two-state matrix through all sixteen columns
and using three fixed rational marks closes it. Fresh 256/384-bit replays
agree on all fourteen arities. Independent 512- and 1024-bit rounded sums
give -53.223053203418. Combining this disjoint class with the 322 cells
above gives the same rounded upper bound for every source support with
r=1..17 and at most one dense group. The
[single-group certificate](HD_TWO_SUM_GF16_SINGLE_DENSE_CERTIFICATE.md)
records the positive matrix calculation. No change to the mixer, route,
or inner was needed for this improvement.

The first complete all-cohort envelope for r=18..32 fails. It replaces
each regional contraction by A B^m; the contact at empty occupancy loses
too much for positive occupancies. Nine fixed mark/potential pairs confirm
the weakness of this envelope, not a low-weight codeword. The
[coverage plan](HD_TWO_SUM_GF16_COVERAGE_PLAN.md) records these negatives.
A nonempty-region comparison improves the complete-band FLOAT bound by
336 bits but still fails by a wide margin. A subsequent control retains
the inner's two states while summing all fixed group types. All twelve
fixed mark/floor schedules still fail: the best complete-band FLOAT
log2 bound is +16494.53668044997. The 274-class path calculation took
7.453 seconds and does not change the encoder. Further fitting of this
first-moment envelope is deferred.

The new [fractional comparison](HD_TWO_SUM_GF16_FRACTIONAL_MASKS.md)
averages outer and inner setup conditional on the route, then applies a
fractional moment before averaging permutations. Its packet powers precede
placement averaging; its fixed-row-set moments and sparse branch weights
are rebuilt. Twenty-five fixed witness/floor controls fail. Diagnosing
the live-state envelope exposes up to 460 bits of avoidable slack per
region. A bounded saved-bank repair removes much of this slack but still
does not close the band. Both controls and their limitations are retained;
neither is evidence of an actual low-weight codeword.

A stronger algebraic simplification now gives the next broad partition.
For a fixed target, exposing D2 and applying H inverse determines every
active constituent row before any D1 key is exposed. Thus the union of
all row sets with total active-row count at least A0 has pointwise outer
measure at most 255^(-8 A0), with no row-set or cohort multiplicity.
This admits a uniform-input comparison for the high-row-count part.
The complementary mask comparison can retain total row count with one
marker, instead of enumerating every cohort. Two independent reviews
confirm both steps; their combined numerical coverage remains to be run.

Classify this turn as progress: a formerly failing small class has a
directed certificate, the source-independent small cover is complete in
its stated scope, and the new row-count partition removes an unnecessary
combinatorial penalty. The next iteration should test that partition,
starting from the earlier working sparse-route witnesses near z=0.98 and
alpha=1/5, rather than continue the rough fractional grid. There has been
no new native benchmark, production change, paper edit, commit, or push.
The final focused regression passes all 166 `hd_two_sum*_tests.py` tests
in 6.679 seconds. The native header, native driver, and closed paired-route
manifest retain their previously frozen hashes. No numerical job remains
active at this checkpoint.

### Previous iteration: remove cancellation rather than only tighten its bound

The next implementation/proof cycle now has a concrete structural win.
Sixteen shorter RS[16,8] constituents reduce outer arithmetic, but initially
lose to Group8 because the compiler outlines their final store helper.
That interface adds 32,768 calls and 64 MiB of intermediate stack traffic
at K20. An explicit inlined sixteen-store path removes the materialization
without changing the matrix, diagonals, RS schedule, or route. A fresh
serial three-round comparison gives 3.398 ms versus 3.616 ms for the
eight-row affine candidate and 4.067 ms for the original sixteen-row path.
All oracle checks pass and existing control binaries remain byte-identical.
The [native record](HD_TWO_SUM_GF16_NATIVE.md) retains both the negative
and the repair, rather than reporting only the winner.

The new [sixteen-row sparse screen](HD_TWO_SUM_GF16_PROPOSAL.md) covers
all 81 source classes and 18 sparse cohorts at r=17, with a FLOAT log2
union of -511.43278211723936. It uses the changed RS support floor and
normalizer, not an old certificate. Dense cohorts and other active-group
counts remain unproved. The faster implementation is not a certified winner.

For the eight-row affine candidate, fresh 256/384-bit replay now confirms
the entire fixed-r17 sparse family at 246.921306449012 bits. Independent
512- and 1024-bit rounded-cell sums agree. Its first mixed dense gate also
improves substantially with an existing dense witness; that work is recorded
separately and does not fill the whole-code coverage gap.

The larger proof improvement is algebraic. For either inverse-disjoint
affine label set, all principal mixer submatrices are invertible. Uniform
input on a fixed coordinate space makes its inside outputs independent
uniform values. Every outside coordinate is a full-support linear form
in those values, and two outside zeros force both moments to zero.
This gives the [exact closed support kernel](HD_TWO_SUM_OUTER.md#closed-kernel-for-the-affine-variants).
Summing fixed row sets of one arity yields a weight depending only on
output cardinality. An ordered coefficient DP replaces exponential mask
enumeration without commuting matrices or resampling row sets.

The first [Group16 arity gate](HD_TWO_SUM_GF16_ARITY_DENSE.md) then covers
all 65,399 dense row sets of arities 3..16, with sixteen two-row groups,
at r=17/source `(8,0)`. One existing witness gives a FLOAT union of
-585.5421219556642 in 2.78 seconds, including reproduction of the older
Group8 labelled calculation. The narrower full-arity two-dense-group
control gives -476.80214895720337. These are not directed certificates
or complete source/cohort coverage.

For Group8, the one-dense-group directed calculation now gives
-184.772727559433; separate 512/1024-bit rounded-class sums agree. A
bounded two-dense-group witness fit followed by one frozen all-21-type
screen gives -2001.873508517472 in FLOAT arithmetic. The winning literal
lies on both fit-box boundaries; neither optimality nor a construction
limit is inferred. The [record](HD_TWO_SUM_AFFINE_DENSE.md) preserves the
initial failed witness and all once-only counting factors.

The same record rules out width64/state8 at 35% distance: a fixed whole-group
message subspace contains a nonzero word whose feedback stays zero and
whose output weight is at most N/4. Stronger within-group mixing does not
evade that dimension argument. Width32/state8 avoids this particular
obstruction but requires a new proof and deliberate SIMD packing.

Next: certify the new Group16 cardinality calculation, extend it to two
dense groups and broader source/group counts, then reassess structural
inner changes. Preserve the 3.398-ms experimental binary and all fully
certified checkpoints; do not promote the sparse/dense gates to whole-code
distance claims.

Final focused regression for this cycle: all 78 `hd_two_sum*_tests.py`
tests pass. Native benchmark jobs are stopped. Research changes remain
uncommitted; no production profile or paper result is changed.

The fixed-r17 sparse proof has made broader progress. One bounded joint
witness refinement gives a complete sparse-family proposal near 50.89 bits,
preserving the stronger homogeneous and sibling schedules. Its remaining
sixteen sparse cohorts have a FLOAT union below 2^-156.78. The proposal and
shared directed-replay plan are in
[the sparse-cohort record](HD_XOR_SPARSE_COHORT_REFINE.md). The sibling replay
has now completed both precisions: all 324 cell ceilings agree, and separate
512- and 1024-bit unions give -70.354148046674 for that cohort.

The twelve-type Hölder compression now has a directed 56.177151934929-bit
receipt for its one-source, one-dense-group class. Fresh 256/384-bit
calculations agree on all twelve cells; independent 512- and 1024-bit sums
of the rounded cells agree. This is not a whole-code certificate.

For two dense groups, replacing labelled masks by occupancy totals loses
too much: even the one-dense-group control loses 467 bits. The subsequent
[labelled occupancy algorithm](HD_XOR_DENSE_TWO_LABELLED.md) retains all
6,561 two-group occupancies through shared suffixes. It reproduces the
one-group control and improves the two-group bound by 539.85 bits. The
single inherited-witness two-group screen is nevertheless +316.24, not
closed. That is a failed witness/envelope control, not a codeword example.

A new [two-sum outer](HD_TWO_SUM_OUTER.md) attacks the construction-level
cause. Instead of only S=sum(x_i), also compute T=sum(u_i*x_i) and emit
x_i+S+u_i*T. For u=0,...,7, four fixed GF operations replace the problematic
equal-pair behavior: every one- or two-input column has at least seven
outputs. An affine set u=8,...,15 adds two more fixed operations but makes
every singleton full and every input pair's support distribution identical.
Neither variant adds random keys. The latter admits an exact positive
hole-histogram calculation without the old pair-multiplicity partition.

The first bounded native comparison passed local basis and full K20 oracle
checks. Serial three-run medians were 3.663 ms for the certified Toeplitz
control, 3.465 ms for XOR2, and 3.537 ms for the u=0,...,7 two-sum candidate.
These matched measurements put the new candidate 2.06% above XOR2 and
3.45% below the certified control. A second matched comparison gives 3.574 ms
for the affine candidate against 3.668 ms certified and 3.471 ms XOR2 controls.
It adds 2.97% to XOR2's time while simplifying the sparse support distribution.
All basis and full-size oracle checks pass; the old builds remain byte-identical.
The [native record](HD_TWO_SUM_NATIVE.md) preserves timings and binary pins.
No local mixing theorem is being promoted to a distance claim.

The affine sparse proof interface subsequently passes every r=17 source
geometry and sparse cohort in one FLOAT screen: 1,458 cells, log2 union
-246.921306449020, 20.03 seconds, one fixed inherited witness. It averages
the exact joint hole distribution before the potential quotient. Regional
subset envelopes remain, but the old column-occupancy resource maximum and
pair-multiplicity strata are absent. Directed replay, other group counts,
and dense cohorts remain open. This is a structural simplification, not a
claim that the full 35% certificate now follows.

## Iteration 1: retain the SIMD layout through routing

The outer writes four adjacent 128-bit records per cache line, one from
each of four source groups. The current route loads these records from
independently selected positions. The inner reconstructs SIMD vectors and
then transposes its output back to the previous coordinate order.

Two changes must be tested separately:

1. **Native output order.** Within each 64-record inner step, store logical
   coordinate `16*lane+j` at physical coordinate `4*j+lane`. This is a
   fixed output permutation. It leaves the sampled map unchanged up to
   that permutation, so the existing distance certificate transfers
   exactly. The intermediate output transpose is no longer required.
2. **Four-record packet routing.** Independently in each region, permute
   the 2048 cache-line packets rather than the 8192 individual records.
   Retain the lane inside each packet. Sixteen source packets fill one
   inner step, with logical coordinate `16*lane+j` taken from lane `lane`
   of source packet `pi[16*step+j]`. The outer and inner algebra stay fixed,
   but the route distribution changes. This candidate needs a new proof.

The second change replaces four indexed record loads and SIMD assembly
with one indexed vector load. It also needs one route index per four
records. The proof must account for active groups sharing source packets;
the old independent scalar-route certificate is not evidence for it.

The test matrix compares the previous schedule, native output alone,
packet routing alone, and their combination. For each sampled setup,
the independent dense oracle uses the expanded logical route and the
explicit output permutation. Packet basis tests must compare both emitted
output and next state, including arbitrary incoming states.

## Iteration 2 candidate: fewer operations in each outer column

In parallel, examine eight RS[32,16] rows with independent nonzero scalar
randomization followed by the XOR map `y_i=x_i+sum_j x_j`. With eight rows
over a characteristic-two field, this map is its own inverse. A singleton
input spreads to seven outputs, but two equal active coordinates have zero
sum and receive no additional spreading. Invertibility is not a distance
proof. The current pointwise outer-measure bound cannot simply be reused.

The potential saving is structural: shorter RS constituents and one
variable field multiplication per output, instead of nine per four
outputs in the current column mixer. A useful candidate must bound the
correlated cancellation patterns over all messages, not only a fixed
random input. Its implementation and proof status will be recorded
separately from the packet-routing experiment.

## Progress log

| Stage | Performance evidence | Proof / validation status | Decision |
|---|---|---|---|
| Retained reference | About 4.0 ms full; outer about 1.85 ms, routed inner about 2.0 ms | Complete 35% / 72.3975-bit certificate; independent native and sanitizer checks | Preserve |
| Eleven outer tile sizes | Best typical difference about 1--2%; reversed runs overlap | Exact same code | Stop tile hill-climbing |
| Native output permutation | 3.823--3.826 ms full | Exact distance-preserving permutation; basis and dense-oracle checks pass | Retain as certified-layout option |
| Four-record packet routing, native output | 3.473--3.500 ms full | Correct implementation; current proof fails on medium fully occupied packet classes | Preserve timing candidate; do not certify |
| Two-record packet routing, native output | 3.499--3.509 ms full | h4 bounds improve, but medium packet classes remain open | Strengthen proof and break fixed pairing |
| Two-record routing with global h16 | 4.217--4.258 ms full | All native tests pass; exact four-step FLOAT closes the tested full-pair 1024-group class | Reject as performance winner; preserve proof lesson |
| Shifted two-record routing, random parity | 3.514--3.519 ms full | Native checks pass; new support-count bound, no complete distance proof | Continue proof screen |
| Shifted two-record routing, alternating parity | 3.557 ms initial repeated runs; 3.607--3.634 ms in later comparison | Native and sanitizer checks pass; three formerly obstructing subclasses close in FLOAT screens, not a complete distance proof | Broaden support coverage, validate arithmetic |
| Eight-row XOR column mixer, H D1 | 3.148--3.182 ms full with packet4; 3.596 ms scalar route | Correct implementation; exact membership formula and cancellation obstruction identified; no certificate | Preserve as speed ceiling, study richer randomization |
| Eight-row XOR column mixer, D2 H D1 | 3.252 ms packet4; 3.331 ms packet2; 3.686 ms scalar route, first samples | Conditional value distribution simplifies proof; full correlated-mask bound remains open | Screen before promotion |
| D2 H D1 outer plus alternating pair route | 3.296--3.332 ms full | Native checks pass; both changed interfaces still require joint proof | Preserve combined timing candidate |
| Toeplitz outer, alternating pairs with independent lane swaps | 3.798 ms matched median; earlier 3.660--3.778 ms | Complete 35% / 55.0927-bit certificate; native, independent oracle, and sanitizer checks | Preserve certified research checkpoint |

After each meaningful result, update the table, identify the remaining
dominant work, and choose another construction-level change. Keep failed
bounds distinct from actual low-distance examples, and keep FLOAT screens
distinct from directed complete certificates.

## Measurement and validation receipts (2026-10-08)

These are serial Peach CPU 2 measurements, with setup, allocation, and
input generation excluded. Each full run has 501 timed calls after three
warmups. The unrelated CPU 0 job still shares L3. Matched frozen-reference
medians were 3.907, 3.863, 3.949, 3.954, and 3.892 ms; a freshly rebuilt
reference gave 3.889 ms. Report repeated ranges, not an isolated minimum.

Four-record routing without the output-order change measured 3.560 ms.
Thus most of the combined improvement is the route, not removal of the
output transpose. Independent phase runs gave 1.329 ms for packet4/native
inner versus 1.831 ms for scalar/native inner. The H D1 outer took 1.563 ms
versus 1.871 ms for the certified Toeplitz outer. Phase timings are not
additive estimates of full timing: cache state differs between experiments.

The original control and the new builds pass small full-encoder checks
against the independent dense oracle, including guard regions, input/output
aliases, unaligned input, and output/scratch offsets 0 through 3. Each
optimized route build also checks 30,720 input/state basis directions
against the previous inner. K20 full dense-oracle comparisons pass for all
tested builds. Native-output-only, alternating-shifted-route, and combined
XOR2/alternating-route ASAN/UBSAN checks also pass. These tests establish
implementation consistency, not distance.

The final combined theory regression run passes 72 tests, including the
prior grouped/Toeplitz/35% checks and the new packet, shifted-pair, and XOR
screens. A new Fourier extraction test initially exposed poor relative
accuracy for an extremely small coefficient. That fixture is now explicitly
rejected by the FLOAT screen's normalized cancellation-floor check; it was
not hidden by a relaxed global tolerance. Well-scaled coefficients still
pass a meaningful relative check. This diagnostic is not a certified
rounding-error enclosure.

## Proof-guided next iteration

The [packet-route analysis](HD_GROUPED_PACKET_PROOF.md) retains packet-mask
correlations. At width four, a fully occupied packet activates all four
local inner blocks together. Permuting its lanes does not change that
fact. The FLOAT bound is positive for 64 and 512 active groups even after
two-step composition. These are failures of the available bound, not
demonstrated low-distance codewords.

Width two is materially more promising. At 512 active groups, the two-step
h4 bound is about +83 bits, compared with thousands of positive bits for
width four. Replacing the four RS[16,4] inner maps by one global RS[64,16]
map makes that sampled width-two class negative, about -159 bits. Other
classes, including 1024 active groups, still need better bounds. The
global map keeps width 64, state 16, the existing refresh, and CA=0; its
larger transform is an implementation cost to measure, not assume away.

That cost is now measured: global h16 takes 4.217--4.258 ms full and about
2.324 ms for the inner alone. It is slower than the retained reference,
so it is not the next performance candidate. Exact four-step composition
reduces the previously positive n2=512 h16 FLOAT bound to about -985 bits;
this identifies real state-path fractionalization slack, not a reason to
accept slower encoding. A corresponding h4 screen must average additional
within-step half profiles; it closes n2=256 at about -201 bits but still
fails at n2=512. The partial proofs and implementation remain preserved.

The next structural change is [shifted adjacent pairs](HD_SHIFTED_PAIRS.md).
Each region uses either even pairs or odd pairs on the cyclic group array.
The same unaligned two-record loads suffice, with one duplicated boundary
record per region. Random parity measures 3.514--3.519 ms; deterministically
alternating parity measures about 3.557 ms. Matched baseline runs were
3.869--3.877 ms. Both pass 30,720 basis checks, full guard/alignment/alias
checks at small sizes, and the K20 dense oracle. The h4 inner is unchanged.
The proof benefit is a lower count of supports that remain fully paired
under both matchings: such supports have few runs. Three exact counting
tests pass, but regional mixed-mask operators and complete coverage remain
open. The deterministic schedule removes unnecessary parity randomness
and forces both pairings to appear exactly 128 times.

The refined FLOAT screen now covers three difficult subclasses of the
1024-active-group / 512-full-even-pair class. At 128, 256, and 448 contracted
runs, a positive root-of-unity alias extraction of the odd-region operator
gives approximately -756.60, -1647.41, and -3793.41 bits with a 64-cubed grid.
A 128-cubed crosscheck improves the first to -818.25, consistent with
removal of nonnegative aliases. The even-region operator uses exact regional
extraction with two-step composition. This is a meaningful proof-screen
improvement at unchanged encoding cost. These numbers are not whole-code
margins: other support triples and general mixed orientations remain open.
The numerical cushion in the FLOAT Fourier calculation is diagnostic,
not a rigorous enclosure. Directed complex arithmetic/error bounds are
required before certification.

The [XOR-column analysis](HD_XOR_COLUMN.md) records why invertibility alone
is insufficient. Two constituent rows can make their scaled outputs equal
on 31 columns, cancelling the shared XOR sum. The outer can therefore have
weight at most 70 in 256 coordinates for every setup. This is not a
whole-code counterexample, but a proof must charge that correlated family.
The D2 H D1 option supplies independent uniform nonzero amplitudes
conditional on its support mask, which is a better proof interface than
H D1. Its masks across regions must still be treated jointly.

Retaining the exact RS support counts and the fact that a single present
input is nonzero materially strengthens its sparse proof screen. The
one-active-group classes with one and two active constituent rows give
FLOAT union bounds around -1044 and -301 bits, including the choices of
group and rows. This does not cover more active groups or all row mixtures,
and does not transfer automatically to the changed packet route.

## General-support and directed-proof continuation

The next proof iteration makes no encoder change and runs no new native
timing experiment. Its aim is to recover the speed candidates' distance
bounds before paying for more mixing.

The [general cyclic count](HD_SHIFTED_GENERAL_COUNTS.md) classifies every
support by `(E,O,a,b)`: active even positions, active odd positions, and
full-pair counts in the two adjacent matchings. It supplies an exact
integer count, a constant-time feasibility test, and a positive transfer
matrix bound for boxes of profiles. Exhaustive small-cycle tests validate
the counts and their symmetries. Count symmetries are not used to exchange
the ordered regional operators.

The generalized screen retains separate `10`, `01`, and `11` packet types.
The two histograms are

    H0 = (4096-E-O+a, E-a, O-a, a),
    H1 = (4096-E-O+b, O-b, E-b, b).

The regional product is always `(R0 R1)^128`, starting from zero state.
The following are representative class contributions at the unchanged
35% threshold, not setup margins for the whole construction:

| Profile `(E,O,a,b)` | Method | log2 upper bound | Status |
|---|---|---:|---|
| `(1,0,0,0)` | Exact univariate regional extraction, paired steps | -115.08 | FLOAT |
| `(8,8,0,0)` | Three-type positive aliases, grid 64 | -187.11 | FLOAT |
| `(64,64,1,1)` | Three-type positive aliases, grid 64 | -1604.80 | FLOAT |
| `(512,512,64,64)` | Positive-marker regional bound | -4360.74 | FLOAT |
| `(512,512,512,384)` | Paired even region, odd aliases at grid 64 | -756.5967856 | Directed, 256/384 Arb bits |
| `(512,512,400,400)` | Both regional aliases at grid 128, exact support count | -208.2461070 | Directed, 256 Arb bits |
| `(512,480,450,420)` | One-step regional aliases at grid 128, exact support count | +227.04 | FLOAT; old bound fails |
| `(512,480,450,420)` | Selective paired regional aliases at grid 64, exact support count | -394.2067275 | Directed, 256 Arb bits |
| `(896,128,100,100)` | Three-type positive aliases, grid 64 | -2977.21 | FLOAT |
| `(2048,2048,1024,1024)` | Positive-marker regional bound | -26694.09 | FLOAT |

The [directed alias replay](HD_SHIFTED_ALIAS_RIGOR.md) starts from exact
rational packet caps. Arb and outward binary64 complex intervals enclose
the polynomial coefficients, transforms, matrix powers, and selected
inverse-transform coefficient. Positivity makes extra aliases a sound
upper bound. The proof does not rely on a guessed FFT error or the FLOAT
screen's diagnostic cushion. An independent audit and exact integer
noncommuting-polynomial tests check the arithmetic and schedule order.

The mixed profile `(512,512,400,400)` has singleton packets under both
matchings. Its directed grid-64 bound is -39.82744, so that replay correctly
fails a -40 threshold. Grid 128 gives -208.24610 without changing the
encoder or probability distribution. This is a numerical extraction
improvement, not additional code distance.

The asymmetric profile was the next discriminating test. Its exact support
count improves the cheap count by 30.21 bits. Refining the positive markers
and increasing the grid reduces the bound, but still leaves it positive.
This is a limitation of the current bound, not a low-distance codeword.
The new reduction composes two ordinary packet steps before applying the
fractional power, while retaining a separate looser bound for steps with
multiple singleton packets. At the same witness, it reduces the grid-64
bound from +425.18 to -394.21. It changes proof work rather than encoding
work.

More precisely, let G collect step patterns with at most one singleton
pair and D collect the remaining patterns. Both use the prior entrywise
fractional-power cap. For good-good pattern pairs, compose the original
two-state matrices first, apply the fractional power, and then sum their
weighted type monomials. Call the resulting polynomial C. The two-step
majorant is

    C + G D + D G + D D.

Every coefficient is nonnegative. Entrywise subadditivity of a fractional
power shows C <= G G, coefficientwise. No subtraction of numerical bounds
is needed. The good geometries reduce to 153 patterns without a singleton
and 272 patterns with one specified singleton orientation, so exact
good-good composition needs only 425^2 ordered pairs. The two singleton
orientations remain separate through the polynomial variable x10+x01.
The regional power is now 64 ordered two-step blocks.

The discriminating FLOAT witnesses are saved in
`theory/hd_shifted_general_sweep.py`. From the theory directory, replay
them serially with:

    python -B hd_shifted_general_sweep.py --replay singleton mixed asymmetric --grid 128

These commands reproduce the current bounds without searching again.
Use `hd_shifted_general_certificate.py --grid 128` for the directed mixed
profile receipt; the FLOAT replay is not its numerical proof.

### Directed neighborhood coverage

`hd_shifted_box_screen.py` bounds neighborhoods instead of listing their
individual profiles. It combines the exact rational transfer-matrix count
from `box_count_upper`, a regional coefficient envelope, and the directed
scalar helper `hd_shifted_box_normalizer.py`. The separate
`hd_shifted_box_certificate.py` now replays the complete calculation with
directed arithmetic, including the alias envelope.

For a profile h, write each regional denominator as

    D_j(h) = multinomial(H_j(h)) * product_i marker_(j,i)^H_(j,i)(h).

The regional alias matrices have not yet been divided by these quantities.
The profile-dependent scalar is the outer fractional-moment factor times
`(D_0(h) D_1(h))^-128`. Its logarithm is an affine term plus positive
multiples of `log2(H_(j,i)(h)!)`. Log factorial is convex on integer counts;
a chord between the minimum and maximum possible counts is an affine
upper bound. Summing these chords and maximizing the affine function over
the 16 box corners gives a directed scalar envelope. Infeasible corners
can weaken that bound but cannot invalidate it.

For the operator envelope, retain all coefficients up to the largest
histogram target in the box, compute the positive aliases, and take an
entrywise maximum over every residue that the box can request. Keep the
two regional operators in their original order. Multiplying their
envelopes and the scalar envelope by the box's total support count bounds
the union over that box, not merely its center.

At the asymmetric witness above, the FLOAT grid-64 box bounds are:

| Inclusive radius in each of E,O,a,b | Number of profiles | log2 box union upper |
|---|---:|---:|
| 0 | 1 | -375.52 |
| 1 | 81 | -278.21 |
| 2 | 625 | -120.92 |

All displayed profiles in these small boxes are feasible. The radius-zero
result uses a transfer-count upper bound, which loses 18.69 bits relative
to the exact count. The prototype checks containments and the point-bound
reduction. Its scalar helper is separately checked against a 512-bit
direct formula on exhaustive small boxes.

The radius-two box now has directed replays at 256 and 384 Arb bits. Both
give the exact rational log2 ceiling

    -12092295652059/100000000000.

Thus the entire 625-profile neighborhood contributes less than
`2^-120.92295652059`, at the unchanged 35% weight threshold. The replay
uses the full interval inverse transform, divides by the exact binary
normalization factor, and takes upper endpoints for every requested
residue. It does not certify selected samples and interpolate between
them. An independent audit checked the count, scalar, and operator
envelopes together and found no gap. In particular, the support count
already sums over the box; there is no additional factor of 625.

Replay from the theory directory:

    python -B hd_shifted_good_pair_certificate.py
    python -B hd_shifted_box_certificate.py --radius 2
    python -B hd_shifted_box_certificate.py --radius 2 --precision 384

These are certified contributions, not a disjoint cover of all supports
or a whole-code distance certificate. A broader [frontier screen](HD_SHIFTED_FRONTIER.md)
now has 13 saved FLOAT witnesses, all below -40 at grid 64. Its initially
positive imbalanced profile `(600,424,400,380)` closes after witness
refinement at -392.51, without a new operator. The weakest of these
13 samples is about -291.72. Unsampled profiles remain unproved.

Independently, the [XOR2 multigroup analysis](HD_XOR_MULTIGROUP.md) derives
packet operators for independent nonzero amplitudes, including feedback
cancellation at occupancies five and above. A common positive potential
controls each regional operator conditional on its support count. This
retains arbitrary correlations among regional masks and pays the endpoint
potential factor only once. It gives a directed -97.40067-bit union for
all group counts 1 through 8192 when every active group has exactly one
active RS constituent, including all group and constituent choices.
The full 512-bit replay checks and sums 84 disjoint interval witnesses;
the earlier 29-interval sparse replay is preserved separately.
Groups with two or more active constituents remain open. In particular,
the coarse potential loses an already established one-group two-row
bound, so its failure does not yet justify more expensive randomization.

The combined local regression run passes 137 tests across the prior 35%
and grouped-outer helpers and the new counts, general operators, interval
alias replay, selective paired composition, box normalizers, and XOR
multigroup bound. This run includes the original
packet and shifted-pair tests. A separate read-only audit of the XOR
argument found no gap in nonzero conditioning, regional ordering,
correlated-mask handling, or projective normalization. No benchmark was
run concurrently with this proof work.

## Structural proof reductions and reusable coverage

The [coverage engine](HD_SHIFTED_COVERAGE.md) now reuses each untruncated
directed regional transform across support boxes. Exact support counts,
normalizers, and an explicit bisection tree remain profile-dependent.
It certifies all 6,561 profiles in the radius-four neighborhood of
`(512,480,450,420)`, using 76 disjoint leaves. Their sum is at most
`2^-79.675797981654`; both 256- and 384-bit replays give the same rational
ceiling. The original radius-two receipt is preserved. This is complete
coverage of a declared local box, not coverage of every support.

A second structural reduction [removes lane balance](HD_SHIFTED_UNORIENTED.md)
from the operator bound. It replaces separate singleton-orientation
variables by one singleton variable, after dividing out the orientation
choice count and taking an entrywise maximum. Good-good two-step
coefficients are unchanged. At `(r,a,b)=(992,450,420)`, the directed
grid-128 replay gives a log2 upper bound below -508.65, including all
even/odd splits and their support counts. The regional transform is now
two-dimensional, and the support search has three parameters. No encoder
work or randomness was added for either proof improvement.

The reduced global bound also has a convex profile exponent. Directed
eight-corner evaluations cover two disjoint boxes containing 2,146,689
and 16,974,593 triples, with log2 union bounds -337.92079 and -30008.47234.
Each triple includes every even/odd split. Both rational ceilings agree
at 256 and 384 bits. An optimized route-independent ball bound covers
all active-group counts `r>=7649` at -293.68492, without a separate sum
over supports. The unresolved domain below that cutoff is still large;
these bulk bounds do not establish complete coverage.

Independently, the [XOR2 column potential](HD_XOR_COLUMN_POTENTIAL.md)
averages a column's cancellation mixture before applying a potential.
At two active groups with two active RS constituents each, its directed
log2 failure bound is at most -108.4887608648336599. Independent replay
and source-case coverage checks pass. This closes an interacting-group
case that the coarser regional potential could not handle. It does not
close every constituent-count class or transfer automatically to packet
routing.

Permuting group labels leaves regional occupancy counts unchanged. Using
that symmetry reduces the universal column checks to multisets, while
preserving the order of the eight labelled regions. The three- and
four-group cases now give directed log2 bounds -76.44222 and -44.76798.
The bounded five-group screen reaches only about -13.62. This is a limit
of the current factorized column potential, not a low-distance witness.
Its marks are selected by particular five-active-group column patterns
and then applied to all patterns. Retaining their joint support counts
is the next proof-side refinement to test before adding encoder work.

The final combined local regression run passes 171 tests. Independent
audits checked the orientation-binomial cancellation, the two-dimensional
directed transform against exact integer convolution, the convex bulk
bounds, and the XOR multiset/factorization argument. This round changed
proof and validation tools only. The native performance measurements,
certified reference, production library, and paper are unchanged.

## Domain coverage and a polynomial-time XOR envelope

The new [domain ledger](HD_SHIFTED_DOMAIN_COVERAGE.md) caches directed
regional operators and clips support boxes to the feasible histogram
domain. Exact counts and convex scalar bounds cover boundary boxes
without enumerating their profiles. Three frozen proposal plans pass
replay at both 256 and 384 bits. Their combined completed upper range
is every r>=6016, including all lane balances, with failure probability
at most 2^-104.530220525833. Lower-weight classes remain explicit pending
work; this does not certify the whole faster construction.

The remaining failures are not all losses in the regional operator.
The initial bulk search kept alpha near 0.114. For a difficult
low-empty-pair profile, optimizing alpha as well as the weight marker
changes a positive bound to a strongly negative one. A nearby witness
with alpha about 0.536 closes the formerly worst box
`r=5120..5247, a=1024..1278, b=1024..1231` at more than 244.60 bits.
This is proof slack, not evidence for adding another mixing layer.
The box overlaps other explored domains and is not added to their
reported coverage count without reconciliation.

The [joint XOR column count](HD_XOR_JOINT_COLUMNS.md) retains the exact
category contraction before relaxing the RS support constraints. A
positive polynomial of degree 64r counts the enlarged set, retaining
the total support floor across all active constituents. It certifies
the five-group, two-constituent-per-group class at 129.92 bits. The old
factorized count's five-group failure was therefore avoidable slack.

A second reduction bounds the eight regional occupancies by their
coordinate interval and total count, conditioned on the number of
cancellations. Nonnegative matrix propagation and componentwise maxima
give a resource dynamic program. Averaging the exact binomial
cancellation law then bounds every source configuration. Independent
review checked the chronology, relaxed constraints, directed maxima,
and multiplicities. At six groups it certifies 127.85 bits in under a
second; the literal source replay required about three minutes.
These are proof-check runtimes, not new encoding timings.

The same resource envelope now certifies each group count from six
through sixteen. The saved literal witnesses pass both 256-bit and
384-bit replay; their summed failure bound has more than 63.04 bits of
margin. This bound still concerns exactly two active RS constituents
in every active group. Larger counts and mixed constituent classes
remain unresolved. No exponential enumeration was used for the extension.
Adding the earlier one-through-four-group union and five-group receipt
gives more than 44.76 bits for the homogeneous two-constituent family
through sixteen groups. The combined scope remains a strict subclass
of messages, not a certificate for the XOR construction.

The combined regression run passes 233 tests; a subsequent 28-test run
also checks the saved alpha-correction proposals and ledger edge cases.
The final 13-test joint-column run passes, including directed replay
of the entire six-through-sixteen-group envelope family.

The runtime improvement candidates and the certified reference are
unchanged. The XOR proof here uses scalar routing and cannot be
combined automatically with the alternating-pair proof.

## Coverage extension and a structural obstruction

Adaptive alpha proposals complete the band r=5120..6143. Directed replay
at 256 and 384 bits accepts all 542 leaves in that plan. Joining its
assigned band with the earlier disjoint bands and route-independent tail
covers every r>=5120. The combined upper bound is 2^-104.530220525834.
The ledger records 664 accepted boxes, covering 2,420,639,745 feasible
unsigned profiles; 9,045,190,656 profiles below that band remain unresolved.
These counts describe support geometry, not individual messages.

Two proof improvements help smaller boxes without changing the encoder.
Cached coefficient aliases recover thousands of bits lost by evaluating
regional polynomials only at positive markers. A two-state transfer matrix
counts indexed cyclic supports by their weight and number of runs.
Combining this count with the scalar tilt before maximizing avoids a
second profile-count factor. These improvements close selected local
boxes. Larger run-aligned boxes still have positive bounds; their negative
results are retained in the orientation-independent proof note.

The fixed-lane route has a more fundamental difficulty. Activate every
other outer group. Each routed pair then contains one active lane, and
fixed lane order restricts every inner input to two of its four RS cosets.
The always-zero-state branch gives a positive lower bound on the current
marked comparator for every weight marker. Further marker optimization
cannot close this case with that comparator. This is not a proof that the
actual code has distance below 35%.

This diagnosis suggests a [construction change](HD_SHIFTED_PAIR_SWAPS.md) with little additional data
movement: independently swap the two records of every routed pair during
setup. Encoding can retain two adjacent-record loads, then swap their SIMD
lanes. The proof must average over these independent swap coins, not treat
the two lanes as independently routed. The fixed-lane measurements and
receipts remain separate research checkpoints.

The endpoint now has a directed log2 bound below -95969.87, corroborated
at 256 and 384 bits. This bound covers exactly the two alternating
supports. For general supports, the swap-averaged unsigned coefficient is
a convex combination of the old normalized oriented coefficients. It is
therefore bounded by their former maximum. Good-good two-step caps are
unchanged, and positive composition preserves domination. The previous
unsigned-envelope coverage, including r>=5120, consequently transfers.
Fixed-oriented receipts do not transfer through this argument.

The implementation adds `HDGROUP_PAIR_SWAPS`, default zero. Its enabled
path packs one swap bit into each 16-bit route entry and uses one masked
shuffle per four records. All 30,720 packet basis checks, four swap masks,
odd wrapping pairs, full K20 oracle checks, and ASAN/UBSAN checks pass.
With the flag disabled, GCC reproduces the frozen fixed-lane binary's
SHA256 exactly. The existing theory regression passes 1,200 tests; the
subsequently added endpoint and inheritance suite passes another 11.

Three serial 501-call full-encoding runs on Peach CPU 2 give these medians,
in microseconds. Setup, allocation, and input generation are excluded.
Run orders were fixed/swapped/reference, reference/swapped/fixed, and
fixed/reference/swapped. An unrelated CPU 0 benchmark still shared L3.

| Path | First run | Second run | Third run |
|---|---:|---:|---:|
| Certified reference | 3945.337 | 3998.175 | 3908.687 |
| Fixed-lane paired route | 3576.507 | 3555.498 | 3880.665 |
| Independently swapped pairs | 3750.142 | 3659.662 | 3777.653 |

The third fixed-lane run demonstrates material cross-run noise. The swap
candidate retains a speed advantage over the reference in these runs,
but the exact overhead of swaps is not resolved to a stable percentage.
This round's gain is a tractable endpoint at roughly 3.7 ms, not a new
fully certified speed record.

The XOR-outer branch has a separate outcome. Its current resource envelope
does not close unrestricted r=17, and an actual repeated-pair geometry
attains the troublesome category bound. The subclass in which all 17
groups use one shared constituent pair closes at 67.5436393309 bits,
agreeing at 256 and 384 bits. That smaller subclass has 28 pair choices,
not 28^17. It does not extend the unrestricted certificate beyond r=16.

## Swap-route proof iteration

The [exact sparse calculation](HD_SHIFTED_SWAP_SPARSE.md) now covers every
support with 1 through 64 active outer groups. Its directed union bound is
`2^-107.393558356229`. All 64 per-weight rational ceilings agree at 256
and 384 bits. This closes 23,936 unsigned profile triples, including every
even/odd lane split. It does not cover the remaining intermediate supports.

The calculation retains low-degree regional coefficients exactly rather
than using a broad coefficient bound. Nonnegative degrees permit truncation
after each matrix multiplication without losing any requested coefficient.
A carry-free embedding lets Arb compute these bivariate products with
directed rounding. The successful fractional exponent is one for the
smallest supports and decreases along the saved rational witness schedule.

A second improvement normalizes each regional coefficient before taking
the worst case over a box. This retains cancellations lost by maximizing
the coefficient and its divisor separately. At both arithmetic precisions,
the box `r=500..524, a,b=188..212` has log2 bound at most
`-385.442619507605`; the box `r=1012..1036, a,b=388..412` has bound at most
`-96.292234390170`. Each box contains 15,625 unsigned profiles.
`theory/hd_shifted_projected_bound.py` freezes their exact rational witnesses.

The first general swap proposal has 118 accepted boxes in `r=3072..5119`.
Directed replay at 256 and 384 bits accepts every box, covering
1,912,995,840 profile triples with union logarithm at most
`-554.006034989175`. This is a partial cover. The combined 384-bit replay
with the earlier unsigned bounds, dense tail, two projected boxes, and
alternating endpoint records 4,333,666,836 covered triples and
7,132,163,565 unresolved triples. Its accepted-event union is at most
`2^-96.287463764292`; it is not a whole-code certificate. These combined
counts exclude the separately replayed sparse band.

A subsequent complete plan covers every support profile with
`4096<=r<=4351`. All 3,937 leaves pass at both 256 and 384 bits, with
identical rational ceilings and aggregate log2 bound
`-96.594485017714`. This adds a complete intermediate band, not complete
intermediate coverage. Its 1,009,047,040 triples overlap the preceding
partial plan; a combined ledger must replace that band's earlier
assignments rather than add the profile counts. The plan hash and replay
command are in the [swap proof note](HD_SHIFTED_PAIR_SWAPS.md).

These partial support classes refer to GF(256) messages. A complete
scalar-code distance statement extends to larger alphabets coordinatewise;
projecting a block message can shrink its active-group count, so partial
scalar classes alone do not classify block-message supports.

Construction identity is now checked by the coverage ledger. Importing an
old bound requires the explicit unsigned-domination or route-independent
tail lemma. Saved FLOAT plans supply witnesses only; replay recomputes
every accepted bound. No count of successful probes substitutes for an
exhaustive partition of the support domain.

## Structural candidate: remove the column-mixer layer

Replace four RS[64,32] constituents and their Toeplitz column mixers by
one RS[256,128] over GF(256) per existing 128-symbol input group. Multiply
each of its 256 outputs by an independent nonzero GF(256) scalar sampled
during setup. Keep the 256-region route, independent pair swaps, and h4
inner unchanged. This is a separate outer ensemble, not an implementation
of the previous map.

The new ensemble satisfies the same outer comparison. For a nonempty
support `S` of size `w`, the shortened RS code is zero if `w<129` and has
dimension `w-128` otherwise. An information set injects its exact-support
words into `(GF(256)*)^(w-128)`. There are therefore at most
`255^(w-128)` such words. The independent nonzero diagonal maps any one
of them to a prescribed vector on `S` with probability `255^-w`.
Consequently the expected number mapping to that vector is at most
`255^-128`, exactly the existing nonzero-vector envelope. Independent
groups retain `beta=256^256/255^128`. Injectivity and the projective factor
are unchanged. Existing route-and-inner bounds transfer through this
envelope; gaps in those bounds remain gaps.

At K20, nominal circuit counts predict 655,360 fewer variable GFNI products
and at most 524,288 additional fixed-coefficient butterflies. Compact
coefficient storage decreases from 4.5 MiB to 2 MiB. Both schedules use
an 8-KiB local vector array and the same payload routing. Additional FFT
work could erase the savings; a matched native comparison is required.
The earlier long-RS measurements used GF65536 and older inner schedules,
so they did not settle this byte-field comparison.

The matched native experiment is now complete and loses. The following
medians are microseconds from serial, CPU-2, 501-call runs. Full runs use
Toeplitz/RS256/RS256/Toeplitz order. Setup and allocation are excluded;
the unrelated CPU-0 job still shares L3.

| Measurement | Four RS64 constituents + Toeplitz | Diagonal-randomized RS256 |
|---|---:|---:|
| Full encode, first comparison | 3645.887 | 4229.948 |
| Full encode, reversed comparison | 3702.362 | 4227.072 |
| Outer only | 1881.433 | 2393.518 |
| Inner only | 1769.224 | 1780.675 |

The loss is in the outer. Counting the actual instruction classes gives
only about 6.5% fewer GFNI operations (1096 versus 1025 per four-group
batch). Fixed-affine operations increase from 392 to 769; FFT dependency
depth increases from five to seven. The long circuit has sixteen distinct
leaf bodies instead of four reused across constituents. Explicit scratch
traffic decreases, approximately 96 KiB to 32 KiB per batch, so additional
scratch traffic is not the explanation.

Whole-command `perf stat` measurements, including setup and 501 outer
calls, show fewer retired instructions but lower IPC (2.10 versus 1.63)
and more cache misses (17.78M versus 77.00M). These counters suggest that
operation count misses important locality and dependency costs. They do
not isolate the cause inside the encoding kernel. Keep the smaller-RS
Toeplitz outer as the performance candidate; park long-RS as a correct,
slower control rather than promoting it.

The default-off `HDGROUP_LONG_RS` flag preserves the old optimized binary
bit-for-bit. The new path passes all 128 polynomial-basis inputs across
four SIMD groups, the existing 30,720 inner basis cases, small/natural-size
and K20 dense-oracle checks, and ASAN/UBSAN checks. Its sampling uses
independent nonzero byte diagonals and no column-mixer keys. The native
setup stream remains a reproducible PRNG, distinct from ideal randomness
in the theorem.

Build with the common flags below, replace `HDGROUP_TOEPLITZ=1` by
`HDGROUP_TOEPLITZ=0 -DHDGROUP_LONG_RS=1`, retain pair swaps, and use
`--bench 501 full 20 nt 1 4`. The `1` denotes one RS256 constituent, not
one outer group in the full code. Binary SHA256 is
`f0a8be74a41a6d7988bf493353c5df54f2c2012fa7c8aef1219570d6a5610737`;
the K20 output checksum is `15312363359593041432`.

## Support coverage and proof reduction, 9 October

This iteration changes the proof computation, not the native encoder.
The independent-swap route remains an approximately 3.7-ms candidate;
the certified approximately 4.0-ms checkpoint is unchanged.

The complete sparse range is now `r=1..256`, with directed union ceiling
`-55.360239472934` at both 256- and 384-bit precision. Its 1,431,040
feasible matching profiles include all lane splits. Exact regional
coefficient calculations are reused across nearby weights, without
interpolating probability bounds. The last batch initially missed the
40-bit target only at its left endpoint; moving that weight into the
preceding batch's existing witness closes it. No encoder parameter changed.
The literal schedules and failed first partition are documented in
[the sparse proof](HD_SHIFTED_SWAP_SPARSE.md).

The next dense band, `r=4352..5119`, also closes at both precisions:
435 leaves cover all 2,305,327,616 feasible profiles, with exact union
ceiling `-12080248013127/100000000000`. Its frozen plan4 SHA-256 is
`f023bfd9c50f0ddb9a1bf14d86630d667cacc5ec751233a22e9f511b53c3c972`.
It extends the previously completed `4096..4351` band toward the inherited
`5120..8192` coverage. Plan5 for `3072..4095` remains FLOAT-only and partial:
72.73% of its profile domain passes, with 903,322,494 profiles unresolved.
Its SHA-256 is
`d86046f447f5839c482ce8276145ea9f0d81c26e5eaf6c3bff002fd4d1362bf5`.
Do not relabel this proposal as a certificate or overwrite the frozen plans.

A more scalable [run-convolution bound](HD_SHIFTED_PAIR_SWAPS.md#common-state-potential-and-cyclic-runs)
now replaces ordered matrix powers for every matching-count pair by one
common state potential and one positive polynomial convolution. Cyclic
run counts retain the dependence between the two alternating matchings.
Directed checks close every profile at `r=128,512,1024` with margins
360.76, 831.90, and 81.57 bits. The two intermediate rows are complete
weight classes, not isolated matching configurations. They still do not
provide coverage between those weights.

The plain matching-marginal Cauchy bound was insufficient: it included
impossible or overcounted simultaneous high-collision configurations.
Keeping cyclic run counts fixes that particular loss. A fixed bank of
regional aliases can then be reused across different weights because the
target-specific histogram normalization is performed only after selecting
the coefficient. The stronger frozen frontier witness now closes the
entire `1008..1040` band in directed arithmetic, with exact union ceiling
`-105547875859763/1000000000000`. All 33 per-weight ceilings and the union
agree at 256/384 bits. Each replay builds only 22 directed alias tables,
then evaluates every integer weight; it does not interpolate from samples.
This covers 8,677,104 profiles, including the older 1024 checkpoint.

A bounded six-candidate witness check at `r=1280,1536` finds FLOAT
single-weight bounds of `-81.5180588075` and `-41.3155331617`. The latter
has insufficient slack to infer a wide band or a full-code union. Both
remain proposals, not directed certificates. Exact banks and all six
outcomes are frozen in `theory/hd_shifted_run_tuning_1280_1536.json`, SHA-256
`198c80bf5203a3658cf83114e5f7bbc1b057339162e3af8835db35799f581383`.
These checks tune proof witnesses only, not the encoder or its setup law.

The remaining proof work is still substantial: complete the unclaimed
weights between 257 and 4095, retaining the complete sparse, intermediate,
and dense receipts above. No whole-code certificate or native performance
improvement is claimed for this proof-only iteration. Both native source
hashes from the preceding measured experiment remain unchanged.

### Locating the remaining proof loss

The apparent bottleneck at 1536 active groups was largely a poor proof
witness, not evidence that the encoder needed more mixing. Keeping the
original regional bounds, exact ordered products improve the FLOAT bound
from -41.3155 to -44.7687 bits. Exact correlated support counts instead
improve it to -45.0809 bits. These are separate controls, not gains that
can simply be added. A grid increase from 128 to 256 gives -78.4926 bits
at the center, but does not close the surrounding 33-weight band.

Three bounded changes to the proof's weight marker, with all encoder
parameters unchanged, give a much larger improvement. Multiplying the
original marker by `2^-0.04` and fitting the coefficient markers only at
the center gives a FLOAT bound of -177.6480 bits. The resulting frozen
bank closes every weight in `1520..1552`. Independent directed replays at
256 and 384 bits agree on every per-weight ceiling and on the union
ceiling `-31247068409433/200000000000` (-156.235342047165 bits).
The proposal is `hd_shifted_potential_band_1536_zminus004.json`; its
SHA-256 is `2a7151b798aa08adf244c8fcbaf9be43c23c23928e62dc13af9266baca0a1fe7`.

The existing 512 bank also closes the full `458..566` interval without
new fits. Both directed precisions give the union ceiling
`-49106949211489/500000000000` (-98.213898422978 bits), covering 7,212,399
profiles. This supersedes the narrower `496..528` receipt for coverage
accounting. One new central fit at 992 closes `975..1007` with union
ceiling `-143355167233523/1000000000000`; both precisions agree.
These bands remain partial certificates, not a certificate for all codes
in the changed ensemble.

The directed band helper now accepts an explicit `--max-weights` budget.
Its default remains 33; the natural domain permits at most 8192 weights.
Increasing the budget only allows more individually replayed weights.
It does not interpolate, trust saved FLOAT values, or allow a proposal to
raise the caller's budget. For example, replay the 109-weight proposal
with `--max-weights 109` at each precision.

The next proof step is therefore broad coverage using fixed banks, not
another expensive change to the encoder. Native timings and source hashes
are unchanged throughout this proof iteration.

Subsequent broad checks also close `280..349` at 86.242628866867 bits and
`700..809` at 82.297346435715 bits. Both directed precisions agree on every
weight and the band union. The 1536 bank extends to `1500..1571` at
85.153916334037 bits, again with matching 256/384-bit replays. This wider
interval replaces, rather than adds to, the earlier `1520..1552` interval.

An independent audit identified an additional structural proof improvement:
the outer is already Toeplitz, but the numerical prefix still uses the
older diagonal-Cauchy bound. Its sharper information-set bound is
`beta_T=256^160/255^32`, saving `0.54207006155*alpha*r` bits without changing
the code. See [the derivation](HD_SHIFTED_PAIR_SWAPS.md#use-the-actual-toeplitz-outer-bound).
Keep this as an explicit proof mode; do not silently change old receipts.
This is preferable to changing the column mixer solely for a stronger
atom bound before using the strength of the current mixer.

The first fresh directed replay with this sharper prefix closes every
weight in `1256..1325`, using the frozen 1280 bank. Both precisions give
the exact union ceiling `-710866278539/7812500000`
(-90.990883652992 bits), covering 29,219,225 profiles. The proposal SHA-256
is `af186c43e45728767c3e5616d0450303a16c9914211449164aac1094aa6b9a47`.
The separate 2048 bank closes `1882..2153` under the legacy bound, with
union ceiling `-23841229745461/250000000000` (-95.364918981844 bits) over
277,611,632 profiles. All 272 per-weight ceilings agree at both precisions.

The fixed 3072 bank also closes every weight in `2828..3285`, using the
explicit sharper Toeplitz bound. Both precisions agree on all 458
per-weight ceilings and the union ceiling
`-113060545953299/1000000000000` (-113.060545953299 bits).
This covers 1,072,732,409 feasible profiles. The
[receipt record](HD_SHIFTED_PAIR_SWAPS.md#the-3072-bank-directed-receipt-and-float-edge-loss)
contains both saved receipts and their SHA-256 hashes. The disjoint
coverage table now contains 5,580 of 8,192 nonzero weight classes;
2,612 weights remain in nine explicit gaps. This is not a whole-code
certificate.

The saved FLOAT screen substantially understates the lower-edge margin.
At `r=2828`, its log2 bound is -125.74103779126781, whereas directed
replay gives -593.414269637913: a 467.6732318466-bit difference with the
same witnesses. Center and upper-edge differences are below 0.000016 bits.
Before fitting more banks, reuse logarithms of directed alias upper
tables to check whether numerical screening artifacts create apparent gaps.
No neighboring weights are certified by this comparison alone.

Two exact changes reduce replay work. A singleton residue query now reuses
the same directed alias numerator instead of reconstructing a general
range query. Cyclic-run counts use the integer recurrence
`N[1]=n`, `N[c+1]=N[c]*(r-c)*(n-r-c)/(c*(c+1))`; every division is checked
for exactness. Neither change modifies an operator or count. The 2048
cross-precision replay also compares the old path with these optimizations
and retains every rational ceiling. The final focused suite currently
contains 160 passing tests, including the directed-source screening adapter.
Native sources remain unchanged.

The prefix-only proposal converter preserves all support sums and exact
witnesses; it does not certify saved FLOAT values. Its output still needs
fresh directed replay. The saved-witness CLI now preserves exact literals,
fixing a fail-closed path that previously produced unreplayable proposals.

The new screening adapter builds the directed alias tables once, then uses
their logarithms for fast proposal evaluation. It restores the stored
scaling factors in log space and preserves exact zero entries. This avoids
the coarse FFT floor without changing a proof witness. Its final support
sums remain FLOAT; a selected interval still requires fresh directed replay.
An independent review checked scaling, exact witness preservation, and zero
handling.

The first production control uses the frozen 1280 bank on all 301 weights
in `1150..1450`. At a 90-bit per-weight threshold, the passing interval is
still `1256..1325`. Across these 70 weights, the FLOAT values differ from
the existing directed ceilings by less than `1.5e-9` bits. The screen is
`theory/hd_shifted_potential_band_1280_directed_source_1150_1450.json`,
SHA-256 `3ff09b10ed774e1542ea9dcb013b8369afd3dffbff959e910a624b4216a75706`.
Thus the neighboring gaps are not removed by correcting this numerical
floor. Do not rebuild or refit this bank merely to repeat that control.

A prefix-only inventory of the saved screens adds only 68 currently
uncovered weights at a per-weight 85-bit threshold, or 96 at 60 bits.
Neither selection closes a complete gap. These are proposals, not new
certificates. After the 3072 lower-edge diagnostic, prioritize new
witnesses near 1728, between the 2048 and 3072 banks, and above 3285.
Do not infer a similar numerical gain at the 3072 bank's upper edge:
its FLOAT loss at 3285 is only about `2e-6` bits.

The complete 3072 control then evaluates `2600..3285` using the same 59
markers and directed source tables at 256-bit precision. Its passing range
at the 85-bit per-weight threshold is `2820..3285`; the passing-only FLOAT
union is `-113.06054595188293` bits. The other 220 requested weights fail
that threshold, so the union over the entire requested screen does not pass.
The new values agree with the existing 458 directed ceilings within `8e-9`
bits. Correcting the floor adds only the eight proposed weights
`2820..2827`, not a broad bridge. No fresh directed replay of these eight
weights was performed in that control; its coverage checkpoint remained
5,580 classes. The later intermediate-bank round below supersedes this count.

This control is saved as
`theory/hd_shifted_potential_band_3072_2600_3285_directed_source256_screen.json`,
SHA-256 `01768be5777a2e7af374b3c4ea895a29d0653dafefc10f079da2e85b33573eba`.
Fold these eight weights into a later coverage replay. Prioritize a new
intermediate witness over another old-bank scan or an isolated eight-weight
replay. The two controls distinguish numerical slack from the remaining
parameter-coverage work without changing the construction.

### Smaller outer groups: a bounded next structural question

A distinct candidate uses four RS[32,16] constituents per group: 64 input
symbols and 128 output symbols, rather than 128 and 256. At K20 this
doubles the group count to 16384 and halves the number of regions to 128.
The 64-input, 16-state inner and total number of inner steps stay fixed.
The same information-set argument gives
`beta_128=q^64*(1-1/q)^(-16)` for independent Toeplitz4 column mixers.

This removes one RS stage but does not remove mixer or payload traffic.
A source-level count estimates about 10.95% fewer outer GFNI products;
it is not an end-to-end speedup estimate. Prior g8/RS[32,16] tests kept
256-symbol groups and do not test this design.

The proof change is not a mechanical transfer. There are 8192 routed pairs,
256 inner packets per region, and 128 regions. At matched active-coordinate
density, doubling the number of active groups increases the support-count
entropy by about 4452 bits near the old 1024-group class. The outer prefix
is unchanged at that matched density, but the regional matrices also
change, so subtracting this entropy from an old margin is not a proof.

The [bounded B128 screen](HD_SHIFTED_B128.md) now completes that check:
the FLOAT bounds for one and two active groups are -47.5761 and -96.9903
bits, but the complete new r=2048 screen with the old 1024 witness gives
+4845.9257 bits. No witnesses were fitted and no native work was done.
This is a failed fixed-witness bound, not a construction counterexample.
Finish the current route's certificate first; if B128 is revisited,
improve this middle-weight gate before a native implementation or broad
proof search. B128 is a separate construction, not an upgrade to the
current route's receipts.

### Intermediate-bank closure round (historical checkpoint)

The unchanged paired encoder now has 31 disjoint directed proof bands,
covering 7,651 of 8,192 nonzero active-group weights. All new bands have
matching 256- and 384-bit ceilings, including each individual weight.
The [coverage table](HD_SHIFTED_PAIR_SWAPS.md#disjoint-coverage-checkpoint)
records their domains and union margins. Three gaps remain:
`1326..1499`, `3286..3325`, and `3769..4095`.

Each new bank uses one central coefficient fit followed by at most one
optimization of a common positive state potential. Its exact witnesses
then stay fixed throughout the assigned interval. Denser coefficient
markers resolve the failed coarse screens near 1064--1232; they do not
change the construction. The weakest new directed band is `1200..1255`,
with a 58.8263-bit union margin. This is sufficient for a 40-bit target
if the other classes also close with adequate margin.

The 3584 bank illustrates the numerical distinction: its FLOAT screen
gives about 143 bits over `3326..3768`, whereas directed replay gives
1168.6895 bits for the same witnesses. Accurate-source screening outside
that interval still fails at `3286..3311`. Similarly, the 3936 bank's
accurate-source screen reaches 4050 but fails at `4051..4095`.
These remaining failures require different proof witnesses or a different
bound; improving numerical precision alone is not enough for these banks.

The new helpers preserve exact source exponents and potentials across
central fitting. A potential retune explicitly records its changed value.
The focused suite contains 174 passing tests, including the B128 control,
fixed-bank selection, directed atlas conversion, and exact cyclic counts.
No native source changed in that proof round, and no new timings were claimed.

### Complete closure and matched comparison

Seven further disjoint rows close the remaining 541 weights. The final
cover has 38 rows, covering all 8192 active-group weights and all
11465830401 feasible matching profiles. Each coefficient calculation was
replayed at 256 and 384 bits. The final accounting checks exact witnesses,
dual per-weight or per-leaf ceilings, all assigned domains, and a fresh
directed probability sum. Its log2 union ceiling is
`-55092738508223/1000000000000`, giving a 55.0927-bit whole-code margin.
An independent 512-bit raw sum agrees. The
[closed certificate](HD_SHIFTED_CLOSED_CERTIFICATE.md) records the manifest
and the precise ideal-setup event.

The final 4068--4095 gap reused the existing dense coefficient bank and
subdivided its boundary boxes. It needed no new coefficient fit or encoder
change. The 4051--4067 prefix uses the separate 4080 whole-weight bank.
This distinguishes a proof-partition limitation from a need for more mixing.

The broader focused theory suite passes 181 tests. The final exporter,
manifest-checker, and boundary-helper suite separately passes 27 tests.
Native basis checks and small/full K20 independent-oracle comparisons were
rerun successfully before the serial timings.

Three matched cycles, each with 501 complete encodings per build, give:

| Build | Run medians, ms | Median, ms |
|---|---|---:|
| Retained reference | 3.996491, 4.118699, 4.048848 | 4.048848 |
| Native output only | 3.870926, 3.881285, 3.870415 | 3.870926 |
| Independently swapped paired route | 3.798160, 3.800274, 3.781779 | 3.798160 |

The candidate reduces time by 6.19% against the retained reference, or
1.88% against native output alone. The shared-host noise does not justify
replacing these matched values with the best older isolated run. This is
a certified improvement, not the much larger structural gain still sought.

## Next iteration

Preserve the now-certified paired-route checkpoint. Stop extending its
coefficient banks: its whole-code target is closed. The next structural
candidate is the eight-row RS[32,16] outer with normalized D2(I+11^T)D1
column mixing, using the same independent swaps and h4 inner.
It removes about 19.7% of outer GFNI products and 983040 mixer XORs from
the current static schedule; hot keys fall from 4.5 to 3.75 MiB. These
counts do not predict an equal reduction in complete encoding time.
The larger-RS substitution is now a measured negative result; do not
repeat it solely because it uses fewer variable multiplications.
Preserve old receipts; transfer only bounds whose operator domination has
been established for the new distribution. Do not spend another large
coverage search on the old alternating endpoint's failing comparator.

The [paired XOR2 gate](HD_XOR_PAIRED_GATE.md) now closes the complete
17-active-group, two-active-constituent class at 50.89 bits.
It retains both source matching counts and constituent-pair multiplicity
across all32 columns. Source-dependent state potentials resolve the failed
concentrated classes without changing the encoder. Fresh 256/384-bit
directed replays agree on all 405 cell ceilings and the aggregate
-50.891481584656. Independent 512-bit summation of rounded cells gives
-50.891481584655. No whole-code claim follows from this class gate.

A matched serial comparison now gives3.456162 ms for XOR2 and3.762143 ms
for the certified paired-route control, an8.13% reduction. The route,
swaps, h4 inner, and native output layout match. Separate outer times are
1.628049 versus1.864741 ms; routed-inner times are essentially unchanged.
Basis and independent K20 oracle checks pass. Older unswapped XOR timings
are not used in this comparison. Keep the certified control preserved.

Next extend the gate to neighboring mixed constituent counts. A denser
group can change its active pair between columns, so its multiplicity
cannot be counted as a fixed two-row choice. The mixed proof must use
fixed row-set containment or a proved dense-group comparison.

The [exact one-group mask calculation](HD_XOR_FIXED_MASK_FLOOR.md)
removes an identified source of slack. Keeping the exact output-support
mixture, the fixed constituent set, and the 17-nonempty-column floor
gives a directed 164.796-bit bound for all 219 row sets of sizes 3 through 8.
Both precision replays agree. The product-density comparison failed even
this one-group gate; do not carry that avoidable loss into mixed profiles.

The [rank-one refresh study](HD_RANKONE_STATE_MASKS.md) starts the next
structural inner iteration. It could delete20 of the current54 GFNI calls
per packet, but the initial three-component cap provably cannot close
medium weights. Retaining all 15 nonempty state-block masks improves the
tested bound, but it remains positive. Grouping masks before taking
fractional powers also fails the frozen witness. The next bounded test
allows the new cap its own proof parameters; it does not revisit the
three-component cap's proved obstruction. A bounded joint fit still leaves
the tested 17-state medium bound above +92000 bits. At that witness, two
independent rank-one updates leave a bound above +65000 bits. Neither is
a construction counterexample, but neither justifies native implementation
or broad coverage work. Further progress needs a stronger state argument
or a different refresh distribution, not another nearby marker sweep.

The next distribution-level option is
[diagonal-XOR-diagonal state refresh](HD_XOR_STATE_REFRESH.md). Its fixed
map I+11^T would replace the current interpolation map while retaining
both independent random diagonals. Five exact local tests confirm the
support law and conditional independent nonzero amplitudes. This differs
from a uniform rank-one update, whose rare branch preserves a scalar line.
The [native prototype](HD_XOR_STATE_NATIVE.md) now deletes those 20 fixed
GFNI products, with nine XORs and two lane shuffles. It passes basis,
independent K20 oracle, and sanitizer checks. Isolated inner time falls
about 20%, but matched complete time falls only 5.14% with the Toeplitz
outer and 1.23% with XOR2; the latter is small relative to run variation.
Neither changed-inner build has a distance certificate. The isolated
phase times must not be added to predict complete time.

An [exact local kernel](HD_XOR_STATE_LOCAL.md) handles up to two potential
input positions in one four-state-coordinate RS block. Seven tests check
actual subset ranks, integer counts, and independent enumeration. It does
not assume selected RS columns are MDS. Twelve local and refresh tests
pass in total. The input comparison is uniform-field, not forced-nonzero.
Keep these tools, but do not prioritize a full 625-state sweep for this
modest measured end-to-end improvement.

The cheaper outer's [mixed-group analysis](HD_XOR_MIXED_MASK_FLOOR.md)
now closes a bounded directed gate at 109.86 bits. It covers source matching
counts (8,0), sixteen two-row groups with multiplicity16, and one dense
group with any of 219 fixed row sets. Three fixed witness choices suffice;
an independent audit reproduces all assignments and the union. Fresh
256/384-bit replays agree on all 219 ceilings and the raw union
-109.863950510801. Independent 512-bit summation of rounded cells gives
-109.863950510800. This certifies only the stated mixed source class,
not the whole construction.

The resource constraints used in that mixed gate do not depend on the
common-pair assumption: the pairwise restriction was discarded. A separate
512-bit count transfer therefore replaces 28 common-pair choices by 28^16
arbitrary assignments. Applied to the verified rounded-cell bound, it gives
an upper log2 bound of -37.753626679935 for the same source class and all
constituent geometries. Two subsequent bounded fits target its largest
remaining fixed row-set terms, without adding encoder work. Fresh full
256/384-bit replays now agree on all 219 cell ceilings; a separate 512-bit
sum gives **-61.000916136179** after the exact `28^15` count transfer.
Five fixed witness classes suffice. This is a verified bound for source
matching counts (8,0), sixteen two-row groups, and one dense group with any
fixed row set. It is not a whole-code certificate.

The adjacent [sparse mixed cohort](HD_XOR_SPARSE_MIXED.md), with one
single-row group and sixteen two-row groups, also needed a different
witness from the homogeneous two-row case. A bounded 32-evaluation fit
followed by a complete 324-cell fixed-bank screen gives -70.35414804665
in floating point. Its fresh directed replay is in progress. This family
covers all 81 ordered source classes, unlike the one-dense-group gate.

The [coverage audit](HD_XOR_COVERAGE_PLAN.md) yields a separate structural
proof simplification: invert the fixed inner, route, output diagonal and
column mixer, and leave the input diagonals random. For each possible
low-weight output, its uniquely determined preimage has a fixed number of
active RS rows. The per-row membership cap gives a route-independent
terminal bound for every cohort mixture with at least 61,191 active rows.
Fresh 256/384-bit scalar arithmetic gives -165.775266208552. This removes
a dense tail from the remaining enumeration; it does not settle the
middle-density cases or allow the old Toeplitz cover to transfer wholesale.

The [whole-word Holder envelope](HD_XOR_DENSE_HOLDER.md) gives a candidate
way to avoid enumerating 219 row sets independently for every dense group.
It is a proved domination of fixed-row-set counting measures, not an
assumption that a row set is resampled per column. Its first one-type
control fails. Five witness classes also fail, at -14.99 bits: the
class containing heterogeneous row counts loses about 130 bits internally.
Splitting the five classes by active-row count leaves only twelve nonempty
types and gives **-56.17715193493034 FLOAT**, within 4.83 bits of the exact
219-mask directed result. Eleven tests cover exact toy-field counts and
the Holder inequality, subset moments and the type partition. The new
interface still needs a directed replay and a joint multi-group bound.
For two dense groups there are 78 type histograms, but each requires one
global Chernoff mark and potential; separate per-group witnesses cannot
simply be multiplied.

A second reuse opportunity concerns all sparse cohorts at a fixed source
size. The expensive column-cap bank depends on source geometry and the
witness, not on which groups contain one or two active constituents. Two
fixed banks cover all eighteen such cohort splits for source (8,0) at r=17
with a **-49.60373901560979 FLOAT** union. This takes about 5.8 seconds and
requires no new fitting. The full 81-source extension checks 2673 cells in
396.5 seconds. Its two-bank union fails, so the concentrated source alone
was not sufficient. Failures are confined to n1=0,1,2, where n1 counts
single-row groups. The n1=0 class already has its separate 50.89-bit directed
schedule; n1=1 has the 70.35-bit proposal undergoing directed replay.
All fifteen cohorts with n1>=3 jointly give **-41.4914794112 FLOAT** over
every source class. The n1=2 cohort gives only -35.4960930, dominated by
one cell: source (0,7), multiplicity15, at -35.51892616. These are proof-bound
failures, not low-weight codewords. The raw failed controls remain preserved.
Next target that cell, then increase the fifteen-cohort union's reserve
before attempting a whole-code union. Do not describe the full sparse
family as certified on the strength of these floating-point results.

The [two-coset scheduling experiment](HD_COSET_PIPELINE.md) keeps exactly
the old map and transform count, but consumes the first half of the
outer planes before producing the second. Native oracle, overlap, and
sanitizer checks pass. The matched control remains byte-identical to its
previous binary. NT producer stores regress from 3.784825 to 4.025123 ms;
cached stores take 5.866731 ms. Preserve the default-off research option
and negative result; do not promote it or expand a nearby tile sweep.
The next proof task is broader coverage for the cheaper XOR2 outer.

### Changing the RS message basis does not yet remove work

An invertible change of message basis preserves the outer code image for
every fixed setup. It therefore preserves the spectrum used here, but
changes the map from a particular message to its codeword. Making the
first RS coset systematic is a legitimate research option, not an
identical-message encoder replacement.

The existing zero-coefficient-eliminated transforms use 49 GF products
and 113 XORs per dimension-16 constituent, or 129 products and 289 XORs
per dimension-32 constituent. Each operation acts on a ZMM containing
four independent groups. These counts exclude column randomization.
The naive systematic inverse-first-coset/forward-second-coset composition
has the same counts, not half the work.

The translation map P between cosets satisfies P^2=I, but P-I has rank
L/2 at dimension L, not a constant rank. Its natural recursive realization
also recovers the current counts. For nonzero coset shift a, write

    f(x) = g(x^2+a*x) + (x/a)*h(x^2+a*x),
    f(x+a) = f(x) + h(x^2+a*x),

where g and h have degree below L/2. Two half-size rate-half transforms
and the final combination give M(L)=2M(L/2)+L-1 GF products and
X(L)=2X(L/2)+2L-1 XORs, with M(1)=X(1)=0. At L=16 and 32 these are
exactly the counts above. The center butterfly fusion likewise saves
nothing beyond the existing zero-coefficient specialization. GFNI affine
maps act within each byte; they do not perform the needed cross-record sums.

No new native experiment follows from rank alone. Revisit this direction
only with an explicit factorization that beats these counts after routing
and layout costs. The assessment did not modify kernels or setup.

### Absorbing the outer diagonal does not yet remove work

For the XOR2 candidate, let D be the outer output diagonal after routing
one 64-record packet. The exact coordinate change

    Y' = v + D^-1 A Q,
    Q_next = T(Q + C D Y')

replaces input X=Dv and output Y=X+AQ by Y'=D^-1Y. Here T is the unchanged
sampled state refresh. Coordinate scaling preserves output weight. The
coupled maps A'=D^-1A and C'=CD retain C'A'=CA=0; however, C' is not the
transpose of A'. Reusing that transpose would break the cancellation rule.

The current XOR2 kernel already fuses D2 into its final producer store.
There is no separate payload pass to delete. The certified Toeplitz path
does not apply this outer D2 at all. Source-level operation counts for the
XOR2 path, per 64-record packet and excluding unchanged refresh/routing,
are:

| Realization | Variable GFNI | Fixed GFNI | XORs | Sampled-key expansions |
|---|---:|---:|---:|---:|
| Current outer D2 plus inner A/C | 16 | 26 | 86 | 16 |
| Explicit coordinate conjugation | 32 | 26 | 86 | 32 |
| Direct preweighted A/C | 118 | 0 | 114 | 118 |

Straightforward conjugation therefore adds 524,288 GFNI operations at K20.
Prepacking routed D and its inverse also increases these hot key tables
from 2 to 4 MiB. The random input loads, swaps, and output stores remain.
These are source-level counts, not new native timing measurements.

State-only absorption has an algebraic obstruction. For one RS[16,4]
emission block, DA=AB implies det(B)=product(i in S)d_i for every four-row
subset S, since all such RS minors are nonzero. Comparing subsets forces
all sixteen d_i to agree. Only four blockwise scalar gauges, rather than
64 independent output scales, can move wholly into the state. This does
not remove a vector pass in the present producer layout.

The read-only audit therefore does not justify another implementation.
Revisit only with an explicit weighted A/C circuit beating the current
42 GFNI operations per packet, including displaced outer work. Deleting
D2 outright or independently sampling new weighted maps is a different
distribution and does not retain the existing proof interface.

## Native reproduction identifiers

Build the existing `spin_char2_hd_grouped_forward` target with GCC 15.2,
Release, and these common compile definitions:

    -DHDGROUP_SUBFIELD=1 -DHDGROUP_FUSE_LAST=1
    -DHDGROUP_BYTE_KEYS=1 -DHD_SUMRANK_BYTE_KEYS=1
    -DHDGROUP_INNER=5 -DHDGROUP_FUSE_REFRESH=1
    -DHDGROUP_NATIVE_OUTPUT=1

For the alternating-pair candidate add

    -DHDGROUP_TOEPLITZ=1 -DHDGROUP_PACKET_ROUTE=1
    -DHDGROUP_ROUTE_WIDTH=2 -DHDGROUP_SHIFTED_PAIRS=2

For the certified independently swapped candidate, also add

    -DHDGROUP_PAIR_SWAPS=1

Omitting that definition selects the earlier unswapped distribution.
The next XOR2 proof candidate retains this definition as well.

For the combined XOR2 candidate replace `HDGROUP_TOEPLITZ=1` by
`HDGROUP_XOR_MIXER=2`; retain the default gauge normalization. Run `--check`
and `--large-check` before timing. Full timing arguments are
`--bench 501 full 20 nt 4 4` for Toeplitz and
`--bench 501 full 20 nt 8 4` for XOR2. Both are h4, not h16.
Run only one benchmark at a time. The sanitizer checks use
`-fsanitize=address,undefined -fno-omit-frame-pointer` and are not timed.

Frozen native hashes and K20 checksums for this round:

| Candidate | SHA256 | Output checksum |
|---|---|---:|
| Retained reference | `6c0f42d94ef06b987dec159f092ea30ee286b0544f15e757fdf502fb43154264` | 3966879751720036268 |
| Toeplitz, alternating pairs | `6f0b2210a510ae1e57976573da75b4912253af08f1ccd2eac3031baf3581e6c4` | 18267286937415057400 |
| XOR2, alternating pairs | `2e9d6a8bcd0909d8b5f9c425b2256246297d9e0a6afd3b462b49bc6e4c8a7527` | 8789024993925240057 |
| Toeplitz, alternating pairs with independent swaps | `103c38e52809ab4a1030150280e55080dc1fdba844b4ee17a7ca95e791ed7312` | 7701295210980077712 |

For the swapped route, add `-DHDGROUP_PAIR_SWAPS=1` to the Toeplitz
alternating-pair definitions above. The prototype reserves route bit 15
and accepts at most 32768 rows; the K20 target uses 8192. The default-off
path retains its previous range and code generation. Native setup uses
a separate reproducible stream for swaps; the proof concerns ideal
independent fair setup bits. The swaps are sampled once, not at each call.
