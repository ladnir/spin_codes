# High-distance SPIN: proof and latency goal

Started 2026-10-06. The retained RS baseline is 16.622 ms for one forward
K20 codeword on one Peach core. It has rate 1/2 and deterministic symbol
distance slightly above 50%. See the [RS measurements](../rs_distance_encode/K20_RESULTS.md).

Current direction: the user approved [independent randomized outer outputs](HD_RANDOMIZED_OUTER.md)
after reviewing the binary proof. The new checked K20 prototype takes7.02ms,
with directed sparsequiet and densewhole-class bounds but no complete certificate.
It continues the [blockwise SPIN restart](HD_BLOCKWISE_RESTART.md)
of 2026-10-06. The [accumulator-based candidate](ACCUMULATOR_DETOUR_CHECKPOINT.md)
is preserved, not promoted. Keep the acceptance targets below; use an outer,
one routing stage, and a batched stateful inner without a pre-accumulator.
The entries below remain a chronological research record.

## Acceptance targets

Keep requested K=2^20, rate 1/2, 128-bit application records, precomputed
setup, and one global routing stage. A regional permutation is allowed;
an additional global routing pass is not assumed.

1. Close a finite-length certificate for at least 40% symbol distance,
   aiming for a 40-bit whole-code setup-failure margin. Measure the same
   construction's forward encoder below 5 ms on the RS measurement host.
2. Improve toward 3 ms and 45--50% distance. Record the achieved frontier
   if these targets conflict; neither time nor distance may be silently relaxed.

Setup and allocation are excluded from precomputed encoding time. Required
data conversion, input movement, routing, arithmetic, and output writes
remain included. Compare one codeword with one codeword. Run all benchmarks
serially and retain independent correctness checks outside timing.

## Starting point

The [4K internal-expansion candidate](HIGH_DISTANCE_LOWER_EXPANSION.md)
uses GRS[256,64], regional routing, and a four-input/four-state Cauchy inner
over GF(2^16), retaining two outputs per step. Its sampled 40% bounds are
positive, but they do not constitute a complete certificate. The existing
15.6 ms measurement is transposed, not directly comparable to forward RS.

The [high-rate winner](GF8_HIGH_RATE_WINNER.md), older proof controls,
and production profiles remain unchanged.

## Initial work split

- Implement a forward control for the unchanged 4K candidate. Fuse routing
  with the input diagonal and inner; avoid a separate regional temporary.
  Use an explicit SIMD outer and verify it against the independent scalar map.
- Search lower-expansion constructions and cheaper inners. Apply the
  [window obstruction](HIGH_DISTANCE_WINDOW_OBSTRUCTION.md) and quiet-gap
  dimension bound before interpreting a failed numerical upper bound.
- For a plausible cost/proof candidate, cover all source classes with
  directed arithmetic and verify joins. Point screens remain labeled as such.

Progress entries below must identify a concrete result, a rejected design,
or a remaining obstruction. A speed estimate is not a timing receipt.

## Progress

- Goal created. Forward implementation and lower-expansion proof search
  started independently.
- The exact-map forward control now passes independent scalar, extension-field
  adjoint, alias, offset, guard, Release, and ASAN/UBSAN checks. Its selected
  K20 medians are 16.863 and 16.840 ms. This does not beat the RS baseline
  materially, and its 40% certificate is still incomplete.
- Simpler 2K/3K geometries and smaller state failed tested upper bounds.
  The [proof screen](HD_FAST_PROOF.md) separates those failures from actual
  window obstructions. An expanding-inner alternative is being checked
  before implementation; it cannot use an arbitrarily high-rate outer.

## First direction-matched implementation control

The new `HdFastForward.h` and `HdFastOuter.h` implement the same [256,64]
outer, regional route, GF(2^16) diagonals, and keep-two Cauchy inner as the
old transposed control. The source contains K 128-bit records and output
contains 2K records. The outer materializes 4K internal records. No distance
claim changes with the execution schedule.

The first fused route/inner implementation took 27.2 ms. Explicit 128-bit
record loads and folded setup indices reduced that to 23.5 ms. Splitting
routing from the inner reduced it to 21.1 ms. Sequentially staging the
256-KiB source region before its random reads reduced the routing phase
from about 10.2 to 5.5 ms, including the copy. The inner itself took about
3.9 ms. This identifies a source-access problem in the fused implementation;
the experiment alone does not distinguish cache misses from translation
or access-scheduling effects.

The selected outer processes 64 rows per tile, with coset-major passes.
Its diagnostic time is about 6.6 ms. A 256-row tile was inconsistent across
processes and was not selected. Explicit hardware gathers regressed the
routing phase to about 17.7 ms. Region prefetching alone did not reproduce
the staging improvement. These screening phase medians need not sum to
the final full-call medians.

The final two launches each used three warmups and 101 timed calls on
Peach CPU2, with a 4.5-GHz userspace request and boost disabled. GCC15.2
compiled the explicit AVX512/GFNI kernels with `-O3 -mtune=znver4 -fno-lto`.
Input is regenerated with non-temporal stores and a fence outside timing.
Timing includes the outer, regional copy, routing, both diagonals, inner,
and cached output stores. Setup, allocation, validation, and hashing are
excluded. Payload allocations request huge pages above16MiB; physical
backing was not verified. All benchmark processes ran serially under the
shared host locks. The observed 16--17 ms range is not an improvement over
16.6-ms RS and does not justify promoting this construction.

The output checksum is `6203579142741173778` for every tested schedule.
Tests compare against independent scalar Cauchy matrices and direct
additive-polynomial evaluation, including shortened sources and dirty
destinations. A K16 extension-field adjoint test checks the existing
transposed map. The registered Release test and ASAN/UBSAN run both pass.

The selected binary SHA256 is
`bd8acc55ee6cca59145ff1d18dfddd639bbb829549d43aaa49fcfed44fc3b3e0`.
The ignored source snapshot `tmp/hd-fast-selected-source.tgz` has SHA256
`37d8fa8cff16548ac322450eda97365bc3f735df89dcd819b86ad3570d2c5ac7`.
Raw receipts remain ignored under `receipts/hd-fast/receipts/`.
The retained remote workspace is `/tmp/spin-hd-fast-20261006-PTWML1`.

To build and validate the selected research control:

```text
cmake -S spin/experiments/char2_field_spin -B build/hd-fast -DCMAKE_BUILD_TYPE=Release
cmake --build build/hd-fast --target spin_char2_hd_fast_forward -j 3
ctest --test-dir build/hd-fast -R char2_hd_fast_forward --output-on-failure
build/hd-fast/spin_char2_hd_fast_forward --bench 101 full produced 20
```

Acquire the shared benchmark locks and pin CPU2 for comparative timings.
The compile-time `HDF_*` switches retain exact-map research controls;
they do not select different code distributions.

Next prioritize a cheaper construction or a stronger lower-expansion proof.
This first implementation round supplies a checked forward control, not
completion of either acceptance target.

## Second round: sequential precode and one global permutation

The [striped two-scan control](HD_STRIPED_NATIVE.md) reduces raw expansion
to2K and forward time to roughly12ms, but it is now withdrawn. A deterministic
grouped-check argument gives distance at most832321 atK20 for its[128,64]
outer, below the838861 required for40%. This rules out every setup,
not just the tested numerical bounds. Smaller outers also have explicit
band-to-suffix failure events; see [the proof work](HD_BAA_PROOF.md).

The complete map passes scalar, alias, offset, guard, and sanitizer checks.
Prefetch, split-scan, and exact bucket-routing controls do not approach5ms.
The global gather alone is too costly in the measured schedule. These
results change the next action: use a stronger recurrence and investigate
cache-local routing, rather than pursue the two-scalar-scan certificate.

The proposed follow-up was GRS[64,32], a true order-two sequential
precode, one regional routing stage, and a final weighted scan. Although
the degree-six ceiling did not reject40%, the stronger cuts below do.

## Third round: reject local-scan shortcuts before implementation

The [common-row-order proof](HD_REGIONAL_THREE_PROOF.md) gives distance
at most699169, or33.339%, for every setup of that three-scan construction.
It even allows arbitrary global routing. Preserving common SIMD row blocks
while changing their orders does not repair it: regrouping their checks
recovers a degree-four obstruction. Nine exact regression tests pass.

A rational scalar post-inner has a dense inverse and cheap batched products,
but its inverse has rank-one tails. The [run-cut proof](HD_RATIONAL_POST_OBSTRUCTION.md)
shows that its proposed regional construction fails40% distance with
probability at least0.9638277172631955. The bound is an exact rational
calculation with a conservative variance estimate, not a failed numerical
upper bound. Seven independent tests pass. Neither rejected encoder was
implemented or benchmarked.

The [wide-packet screen](HD_WIDE_PACKET_SCREEN.md) increases state enough
to escape the earlier fixed-window test. Its sampled middle-class bounds
still fail; fifteen frozen witnesses and45 operator regressions are retained.
No proof or implementation win is claimed by increasing width alone.

The remaining new direction uses [two additive states and a varying observation](HD_ADDITIVE_TWO_STATE.md).
Its exact four-state law is checked over small fields. Its inverse genuinely
has rank-two tails, so the scalar run cut does not carry over unchanged.
Both state updates permit batched products and XOR prefixes. A local rank
envelope also survives an independent adaptive-quotient audit, but the first
global-route numerical bounds remain too loose and are not valid substitutes
for a regional-route proof.

Next analyze regional occupancy together with the retained source/state
rank. Require positive proof evidence before a native tuning campaign.
The40%/40-bit certificate and sub5ms forward target remain open. No native
benchmarks ran in this proof round; prior timings and frozen winners are
unchanged.

## Fourth round: an actual two-state obstruction

The [mixed-interval cut](HD_TWO_STATE_CUTS.md) rejects the two-state
candidate with setup-failure probability at least0.9985378145964314.
Unlike the earlier safe-run cut, it permits unsafe holes within a selected
interval. Each hole costs one rank dimension, but merging runs can save
two endpoint dimensions. Exact finite-block moments account for dependence
inside each regional permutation. The resulting bound is independent of
the field coefficients. No native implementation was started.

The [regional contraction](HD_TWO_STATE_REGIONAL.md) is retained as reusable
proof infrastructure. It normalizes every region separately and carries
the inner state between regions. A sharper live-scalar pivot factor is
also checked: max(1/(q-1),1/tau), rather than the previous sum. It does not
rescue the rejected candidate, and it must not be applied to a state known
only to be nonzero as a vector.

We also tested smaller outers in the original4K control, leaving width4,
state4, keep2, and GF65536 unchanged. These are floating point-class bounds,
not full certificates or impossibility results:

| Outer | A=64 margin | Input weight near .03N margin |
|---|---:|---:|
| [64,16] | -380.516 bits | -7783.749 bits |
| [128,32] | +38.798 bits | -1121.590 bits |

The selected fixed witnesses replay with `hd_fast_screens.py --check --name small`.
They were selected with distance2/5, K20, and35 optimizer iterations.
The exact density-class input weight is125829. All four fits reported
convergence. Shrinking this outer does not yet provide positive evidence
for the desired complete bound.

The next two changes address different causes. Distinct fixed row orders
could remove the shared source-row geometry, at the cost of less-local
outer stores. A direct outer with a larger additive state removes the two
pre-scans, but the small[64,32] outer already fails the fixed-window test.
Neither direction is a proved or timed replacement. Continue proof gates
before another native tuning campaign. The acceptance targets and all
frozen winners remain unchanged; no benchmarks ran in this round.

## Next checkpoint: fixed affine placement, with fewer scans

The selected fixed64-stride placement has a directed96-bit spectral
receipt: its largest nonconstant positive Cayley eigenvalue is below.323.
With one pre-prefix, every set of8192 source rows contains at most258144
safe check columns, below its262144 checks. Thus the quarter-row cut that
rejected the common-order two-state code cannot apply at that row-set size.
The result covers every set of that size, not just a sampled cut.
It does not rule out other sizes or more general check subspaces.

The primary next candidate is therefore [64,32]GRS with fixed affine
placement, one input-weighted prefix, one sampled regional route, and
the two-additive-state inner. The fixed placement can be fused into the
outer stores. It is an additional address map, not a second full-record
permutation pass. Four adjacent source records may nevertheless land on
four cache lines, so the store cost needs measurement before extensive
implementation tuning. [The proof note](HD_TWO_STATE_REGIONAL.md) includes
the stride list, arithmetic receipt, coefficient-cancellation treatment
for the two-prefix control, and13 exact regression tests.

The direct-outer alternative uses a larger vector state and independent
uniform field vectors for its update and output observation. Its
[local laws and rank envelope](HD_VECTOR_OBSERVATION_PROOF.md) are useful,
but the bounded r8 screens fail: [128,64] passes the one-active-row class
at about187bits, while tested A16 and A64 bounds remain negative. The
minimum-weight mark does not repair those two gates. These failures do
not prove that the code has low distance. They give no reason to prioritize
its native implementation over the changed-placement candidate.

Next test remaining cuts for fixed affine placement, then measure the
outer-store cost and the full one-prefix/two-state pipeline if those gates
survive. Require a complete finite certificate before claiming success.
The16.6msRS baseline,16.8ms high-distance control, and frozen3.2ms high-rate
winner are unchanged.

Validation for this checkpoint: all72 `hd_*tests.py` tests pass, including
the exact finite cut, directed stride spectrum, small-field rank checks,
and regional normalization. The four smaller-outer witnesses and all nine
vector-observation witnesses also replay. Existing implementation files
were not changed, and no commit or push was made in this proof round.

## Affine placement: measured cost and exact source interface

The [new native checkpoint](HD_AFFINE_NATIVE.md) measures the complete
GRS[64,32], one-prefix, regional-route, two-state encoder. Individual-row
affine placement costs about12.17ms after batching. Preserving four-row
groups changes the construction but permits full-cache-line stores.
Streaming stores and a16-row outer tile reduce full forward time to
8.327--8.360ms over two101-call launches. This is approximately twice as
fast as the retained RS baseline, but it has no full distance certificate.
The current split is2.980ms outer,1.711ms prefix, and3.658ms route/post.

The [cut audit](HD_AFFINE_CUTS.md) improves the individual-row safe-count
exclusion to every source-row set of size1--8567. Its bounded menus find
no retained counterexample among92 individual-row,69 two-row-group, and
129 four-row-group sets. These stationary numerical screens do not prove
a finite failure bound for grouped placement.

The [exact source law](HD_AFFINE_PROOF.md) averages out prefix weights,
then expresses a first-moment upper bound as a positive rank-defect moment
of the remaining map. It applies to either fixed placement. Contracting
inactive-row prefix equalities gives a concrete repeated-component model
for the next all-support bound. Ten exact tests and an independent audit
check this interface, not a completed K20 certificate.

The [vector-source refinement](HD_VECTOR_SOURCE_REFINEMENT.md) also rules
out missing row-minimum constraints as the cause of the older failed
vector bounds. The useful missing variable is state-image rank; local
quiet and singleton-return refinements could improve the cheaper Cauchy
control. No new full certificate follows from those local lemmas.

Validation: all96 high-distance helper tests pass. The grouped native
encoder passes scalar-reference, alias, alignment, guard, and ASAN/UBSAN
checks at four sizes through K4096. Benchmarks ran serially. Production,
certified winners, and the RS baseline are unchanged; no commit or push.

Next pursue the contracted component-rank bound, with a finite test of
which rare profiles dominate. Keep the cheaper Cauchy rank refinement as
the alternate route toward5ms. The40%/40-bit and sub5ms goal remains active;
the8.3ms result is an uncertified performance checkpoint, not completion.

## Component-class proof progress

The [component contraction](HD_COMPONENT_RANK.md) exposes a concrete issue
in the previous proof target: no single common Chernoff mark can bound its
sparse and dense terms simultaneously. The positive expansion can instead
use a separate mark for each permitted-row count a. These auxiliary row
sets are not the exact active outer rows of a codeword.

The complete a=1 term is now bounded above by 2^-84.93786249 uniformly over
all fixed regionwise placements. The [new pair analysis](HD_COMPONENT_PAIRS.md)
also bounds the complete a=2 term above by 2^-41.11393101. Both use directed
arithmetic and include all row choices. Independent audits confirmed the
normalization, global carries, zero subtraction, and positive upper bounds.
These terms together consume less than 47% of the 40-bit total budget.
Every a>=3 term remains to be bounded; there is no complete certificate.

The same positive post operator covers a [three-product inner](HD_THREE_PRODUCT.md):
set the feedthrough coefficient equal to the first update coefficient,
and observe the updated first state plus a random multiple of the old
second state. Its local map remains invertible. Five exact tests and an
independent audit validate the new support law. This is a changed ensemble,
not an exact-map optimization. No native performance claim is made for it.

The alternate [Cauchy cap refinement](HD_CAUCHY_RANK.md) is composable and
tested, but adds only80 bits at A16 and7 bits in the low-middle witness.
Its large negative margins remain. A bounded [heterogeneous grouped-cut
search](HD_GROUPED_CUT_FOLLOWUP.md) also found no new obstruction; this is
not a finite failure upper bound. Those results prioritize the component
proof over more cut sweeps or Cauchy singleton tuning.

Validation: all 126 high-distance helper tests pass. The independent pair
audit includes 636 exact local comparisons and a 160-bit directed replay.

Next extend the component bound to a>=3, retaining the number and lengths
of live components. Benchmark the three-product distribution only after
stronger proof evidence. The8.327--8.360ms checkpoint and16.622msRS baseline
are unchanged. No native benchmark, production change, commit, or push
was performed in this proof round. The full goal remains active.

## Rare-route multiplicity: change the proof objective, not the goal

The next round identifies a genuine obstruction to the uncapped codeword
first moment for both measured source placements. The
[frozen group-four witness](HD_COMPONENT_DENSE.md) has a route-event
probability near 2^-1587966. On that event, every admissible coefficient
setup contains a low-weight actual-code subspace of dimension at least
126841 over GF(65536). Consequently the expected number of low-weight
projective codewords exceeds 2^441457. The actual MDS equations are already
imposed, so refining their enumerator cannot eliminate this contribution.

This is not a rejection of the 40-bit setup-failure target. It rules out
using the uncapped first moment to prove that target. Even a separate
Chernoff mark for each permitted row set cannot repair its middle term.

There are two positive results. An information-set argument bounds the
entire dense tail a>=29943 below 2^-322.58291595. More importantly, the
[new conditioning gate](HD_ROUTING_CONDITION.md) counts each bad route
once. Over every row set with a>=18300 and every output support of weight
at most 838860 with at most 64 runs, the contained-component dimension
obstruction has probability below 2^-155.23606481. Regional occupancies
need not be balanced. Independent review confirmed both the probability
argument and a 256-bit directed replay.

The remaining proof must control low-weight subspaces on the complement
of suitable routing events. This includes smaller row sets, output supports
with more runs, and cancellations outside those supports. The gate above
does not characterize all low-weight vectors and is not a full certificate.

Bounded sparse refinements remain indexed. Minimum-occupancy and joint-cut
bounds fail at a=3; an independent cyclic-group-offset proposal improves
the tested upper but does not pass. The fixed group-zero triple also remains
unclosed under the inexpensive total-degree MDS filter. No new distribution
was implemented or benchmarked on the strength of these failed bounds.

Next pursue a routing-conditioned rank-tail bound. The prior 8.327--8.360ms
native checkpoint, 16.622ms RS comparison, and certified high-rate winner
are unchanged. No production change, remote run, commit, or push was made.
The complete 40%/40-bit and sub-5ms goal remains active.

Validation for this checkpoint: all 158 high-distance helper tests pass,
including exact small-field actual-code witnesses and directed arithmetic.
No native benchmarks ran in this proof round.

## User-directed restart: blockwise SPIN

The user requested preservation of the accumulated work and a restart from
SPIN's blockwise structure. The [checkpoint index](ACCUMULATOR_DETOUR_CHECKPOINT.md)
records the measured/proved frontier and a source-only archive whose
578 files were verified after extraction. Existing sources, proofs, and
certified winners were preserved; there was no commit or push.

The active [restart](HD_BLOCKWISE_RESTART.md) uses a fast outer, one route,
and a batched stateful inner. It omits the extra pre-accumulator. The
targets remain K20, rate 1/2, 40% distance/40-bit margin and sub-5ms
forward time, followed by the higher-distance/lower-latency frontier.

The first finite geometry screen distinguishes outer sizes quantitatively.
At 40%, balanced-window lower bounds require at least 46 state symbols
for GRS[64,32] merely to escape rejection, but only three for GRS[128,64].
Thus the initial shortlist is GRS[128,64], t=32/64, and s=16/32. These
are unexcluded candidates, not certified operating points.

A cheap field-valued rank-one state mixer now has a checked exact local
law and a literal packet recurrence, inverse, and transpose. At t=64,s=16,
one refresh requires 0.5 variable field products per output record,
excluding fixed maps and routing. It cannot be treated as a full-rank
uniform-state-map refresh in a quotient proof. The fixed-map joint ranks
and complete distance proof remain open.

All 176 high-distance helper tests pass, including the 18 new restart
tests. No restarted native timing is claimed. Next choose and analyze the
rectangular fixed maps, then measure a validated full forward control.

## Direct-SPIN maps, capped route event, and native control

The restart now has exact rectangular interpolation maps with C=A^T and
CA=0. A positive three-component local bound retains zero-state mass,
arbitrary live-state mass, and uniform-live density. It handles correlated
emission and returns without treating rank-one refresh as a fresh dense
matrix. See [the local analysis](HD_BLOCKWISE_MAPS.md).

The [global check](HD_BLOCKWISE_GLOBAL.md) proves that an uncapped actual
codeword first moment also fails for this direct construction. It then
bounds the complete-row, connected-window trap family below2^-1511.82,
counting each bad route once. Partial-row shortening and general low-weight
output supports remain open; this is not a complete distance certificate.

The [native forward control](HD_BLOCKWISE_NATIVE.md) uses the same outer,
maps, route, and rank-one mixer. Its final K20 medians are8.958--8.977ms,
with about3.056ms outer and5.902ms routed inner. Exact pre-routed inner
inputs cost3.047ms in a diagnostic that excludes routing materialization.
Prefetching is retained; hardware gather was slower in its measured control.
The construction has no pre-accumulator and no extra routing stage.

All190 high-distance helper tests pass. Release and ASAN/UBSAN native checks
pass against independent dense scalar maps, including nonzero state, both
index widths, exact aliasing, offsets, guards, and K=12288. The old archive
and accumulator source hashes remain unchanged. No production changes,
commit, or push were made. Both full certification and sub-5ms remain open.

Next combine the local operator with source-aware capped routing bounds,
and assess cheaper fixed maps/state sizes and routed-load overlap against
this correctness-checked cost control.

## Quiet-subcode rejection and smaller packets

The [quiet-subcode theorem](HD_BLOCKWISE_QUIET_SUBCODE.md) now rejects
width64/state4 and width64/state8 for every setup. With identity feedthrough,
their distance is at most12.5% and25%, independently of outer length or
state refresh. The proof restricts messages to a small set of outer rows
and annihilates every preterminal feedback syndrome. A larger outer alone
cannot repair this failure.

The next controls shrink width together with state. GRS[128,64], t32/s8
and GRS[256,128], t16/s4 both have a50% quiet-subcode ceiling and escape
the tested paired-window obstruction. The complete one-active-row bounds
are143 and223 bits at40% distance. These restricted-class bounds do not
certify either whole code; the many-row classes remain open.

The [native comparison](HD_BLOCKWISE_NATIVE.md) selects the second control
at6.348--6.352ms forward. Its measured phases are about2.301ms outer and
4.043ms staged routing plus inner. The full time includes the sequential
region copy; setup and allocation alone are excluded. Cached stores take
7.36ms, and separating gather from the inner takes7.41ms. The original
9-ms direct baseline and archived8.3-ms accumulator detour remain unchanged.

All219 high-distance helper tests and the selected control's Release and
ASAN/UBSAN checks pass. There is no complete high-distance certificate,
production change, commit, or push. Next prioritize a source-aware capped
many-row proof for the fixed6.35-ms geometry, while preserving the sub-5-ms
target rather than relabeling this checkpoint as completion.

The subsequent [interval-span bound](HD_BLOCKWISE_INTERVAL_SPANS.md) controls
the dimension-forced complete-row, one-interval-per-region trap family for
that exact geometry below 2^-528.39. Exact order statistics replace the
loose union over window endpoints, and a minimal-row witness removes the
unnecessary sum over larger row sets. The actual-code first moment still
exceeds 2^740 on a much rarer contained-row event; counting that route once
avoids this source of slack. Partial-row shortening and arbitrary output
histories remain open. This is concrete proof progress, not a full certificate.

## User-directed structural co-design

The user explicitly redirected work away from timing/layout tuning. The
6.35-ms implementation is now a fixed reference, not the next optimization
target. New operations must address a distance obstruction at low cost.

The first [shared-transform proposal](HD_BLOCKWISE_SHARED_MAPS.md) uses
eight-point MDS feedthrough, a four-symbol state, and 64 inputs per packet.
One inverse FFT supplies both input coefficients and feedback; state
coefficients are added before the forward FFT emits the output. This
spreads quiet-state inputs instead of passing them through unchanged.
The exact map identities, all 12,869 tile minors, and scalar basis cases
pass. No whole-code bound or native timing is claimed.

The [cost model](HD_BLOCKWISE_SHARED_COST.md) gives 80 GFNI operations per
64 outputs, versus 88 in the reference, and one refresh boundary instead
of four. That modest count reduction alone is not a claim of reaching
5 ms. Its value is obtaining substantial input mixing without a separate
matrix application, enabling a different distance/cost tradeoff.

A [global-prefix cut](HD_FAST_PROOF.md#global-prefix-strengthening) rules
out the proposed 2/3- and 3/4-rate outers with a proportionally expanding
small-state inner: their distance ceilings approach the outer redundancy,
regardless of the route or local MDS feedthrough. The active shared circuit
therefore retains the half-rate outer. A cheaper polynomial state update
is parked as a separate conditional idea, with no selected polynomial or
claim across nonzero inputs.

Next derive the shared circuit's joint tile support and feedback-rank bound,
retaining actual outer constraints and capped route events. Benchmark only
after the new map survives that structural analysis.

All 224 high-distance helper tests pass after the shared-map checks. The
archived accumulator source and archive hashes remain unchanged. This
structural round made no production changes, native implementation of the
new circuit, commit, or push.

## Structural progress: shared sums and joint ranks

The next round made progress without native tuning. It preserves the6.35-ms
reference and develops an exact cheaper circuit for a stronger local map.

- [Joint support analysis](HD_BLOCKWISE_SHARED_SUPPORT.md) identifies each
  tile as a GRS[16,12,5] map when its four state coefficients are included.
  After eliminating tile inputs, output zeros impose one shared four-column
  state system. Its ranks must be stacked, not added. Shift8 has a concrete
  paired zero-feedback swap; shift64 removes that duplication. Some remaining
  rank-three examples show why independence still cannot be assumed.
- [The custom cube circuit](HD_BLOCKWISE_SHARED_FACTOR.md) computes the same
  shift64 map using59 packed GFNI arithmetic calls per64 outputs, versus82
  for the two-FFT circuit. It reuses19 local cube sums per tile and22 global
  moments. Exact basis and extension-field checks pass. The extra shuffles
  and live moments are unmeasured, so this is not an end-to-end speed claim.
- [The tiled quiet bound](HD_BLOCKWISE_SHARED_BUDGET.md) holds for every
  route. Eight-point mixing with width64/state1 cannot exceed22.439%
  distance. Width64/state2 has only a40.342% ceiling; width128/state4 is
  similarly marginal. The active width64/state4 geometry survives this gate.
- [A source-aware quiet-class bound](HD_BLOCKWISE_SHARED_CUTS.md) covers
  every nonzero quiet word on at most256 outer rows, regardless of output
  weight. Its directed route probability is below2^-2236.70950574. The
  local quiet packet has combined input/output distance13. Neither result
  covers larger quiet supports or nonquiet histories.

Quiet distance also descends exactly fromGF65536 toGF256, since the quiet
constraints contain no extension-field mixers. This halves the field-entropy
cost of counting that class. It does not remove the rare-route problem:
an explicit640-row all-quiet trap still gives an uncapped base-field first
moment above2^4110.59. The event itself has probability about2^-231841.
Use base-field counting inside a capped or route-conditioned proof, not as
a reason to repeat the failed unconditional count.

The concrete next step is to extend the source-aware bound beyond256 rows,
retaining feedback-image ranks and outer constraints. In parallel, an
explicit register/shuffle schedule can determine whether the59-call circuit
deserves a native experiment. The objective remains the same complete finite
certificate and measured sub-5-ms forward call. No production change,
native benchmark, commit, or push accompanies this structural checkpoint.

Validation: all254 high-distance helper tests pass. Independent reviews
checked the joint-rank interface, sparse route bound, state-budget theorem,
and59-call arithmetic accounting. The generic state-budget helper's
rejection flag was corrected to include the ordinary Singleton bound.
The6.35-ms source hashes, accumulator source, and frozen archive are unchanged.
The documentation-link check found207 valid local links and no broken links.

## Native cost correction and remaining proof structure

The [native shared-map controls](HD_BLOCKWISE_SHARED_NATIVE.md) pass independent
scalar and sanitizer checks. At K20, the same beta64 map takes6.541ms with
FFT,7.651ms with the59-product cube, and7.004ms with a63-product whole-vector
rewrite. The frozen identity-feedthrough reference reruns at6.346ms. Lower
GFNI product counts have not lowered complete-call latency; this tuning path
is stopped in favor of construction changes that remove shared-map work.

Two new proof interfaces are saved and independently reviewed:

- The [quiet composition bound](HD_BLOCKWISE_SHARED_QUIET_COMPOSE.md) closes
  all base-field quiet witnesses whose top128 regional type-information sum
  is at least7,693,904bits, with failure below2^-80.8712. It uses actual outer
  value consistency. The low-information complement remains open.
- The [joint-image transfer](HD_BLOCKWISE_SHARED_TRANSFER.md) gives the exact
  rank-up probability, tight rank-down bounds, and a capped whole-subspace
  return event under rank-one refresh. It exposes deterministic return when
  the joint image contains a two-dimensional diagonal subspace. It does not
  make ranks alone a complete whole-code proof state.

All265 high-distance helper tests pass at this checkpoint. Native checks
cover all packet bases, genuine extension-field states, route widths, aliasing,
offsets and guards, and complete outer integration at three sizes. Both cube
schedules pass ASAN/UBSAN. The archived sources and hashes are unchanged.
The construction still lacks a complete40%/40-bit certificate and a sub-5ms
implementation. No production change, commit, or push was made.

The next structural questions are whether one rank-four local factor can
serve both input mixing and feedback, and whether a stronger low-cost state
randomizer can eliminate guaranteed subspace cancellations. These are new
construction questions, not optimizations with inherited distance proofs.
