# One-permutation field-SPIN: progress record

The new objective is a useful finite proof and fast precomputed transpose
with one outer-to-inner routing stage. No second global shuffle is allowed.
The two-permutation construction remains a proved reference, not a candidate
for this objective. General large-field multiplications remain excluded.
The regional permutation is one structured permutation of all coordinates;
it is not uniform over all permutations. Both distributions must be labeled.

## Starting checkpoint

The GF256 GRS[16,8], t16, two-symbol SL2 inner measures approximately
0.133/0.575/6.78 ms at K16/K18/K20 on 128-bit records. Its K16 one- and
two-active-row bounds jointly contribute below 2^-43. A separate dense
bound contributes below 2^-64. Active-row classes 3 through 2878 remain
uncovered. These numbers are not a whole-code certificate.

Ordinary bad-word counting cannot complete that proof. Rare regional routes
admit large subspaces of bad messages, making the expected count too large.
The rank-based continuation bounds in [SUBFIELD_PROOF.md](SUBFIELD_PROOF.md)
instead count independent constraints. Their local inequalities have been
audited; complete composition across regional routing is the next proof task.

## Iteration 1: reject an unstructured small-outer route

A uniform global permutation has a stronger problem: late activation gives
an appreciable probability that the code actually has small distance.
This is not the bad-word overcounting issue above.

Fix a repeated MDS[b,k] outer with L rows, so K=kL and N=bL. Let the causal
linear inner start in zero state. Choose a suffix of W output coordinates,
aligned to its step boundary, with W below the desired distance. If at least
b-k+1 coordinates of one outer row route into this suffix, shortening the
MDS row to those coordinates leaves a nonzero message. Its entire routed
input lies in the suffix. Causality confines its output to the suffix too.
Nonzero diagonal coefficients and arbitrary causal state updates cannot
remove this event.

Let p be the exact hypergeometric probability of this event for one row,
and p2 its joint probability for two disjoint rows. For any a<=L rows,
Bonferroni gives an actual code-failure lower bound

    Pr[bad code] >= a*p - binom(a,2)*p2.

The independent test enumerates all small suffix sets to check both
hypergeometric probabilities and the union inequality. At K=65536,
N=131072, t=16, and W=13104<13108, exact rational calculations give:

| Half-rate outer | Code-failure lower bound, rounded down |
|---|---:|
| [8,4] | 0.5009 |
| [16,8] | 0.04717 |
| [32,16] | 0.000005192 |
| [64,32] | 0.0000000000001488 |

The last row leaves at most about 42.61 bits before considering any other
failure. None of these lower bounds is a positive distance certificate.
The regional layout avoids this particular event: a short terminal window
cannot contain enough distinct regions from an MDS outer row.

Reproduce with `theory/one_route_tail_obstruction.py`. The rational lower
bounds are exact; decimal probabilities and displayed bit margins are
diagnostics. The three tests are in `test_one_route_tail_obstruction.py`.

## Iteration 2: exact-map implementation changes

`spin/experiments/char2_field_spin/OneRoute.h` implements eleven execution
schedules for the same GF256 code. The new successful schedule processes
regions in reverse order, carries the inner state across region boundaries,
and reorders one region immediately. Thus the random reads access a
region-sized temporary rather than a full-code-sized temporary. Non-temporal
stores avoid fetching destination cache lines before overwriting them.
The outer still executes the same fixed butterfly circuit.

Serial measurements on the Ryzen 7950X, GCC 15.2, CPU 2, use precomputed
setup and 128-bit records. Two independent setup seeds, 51 full-encode calls
each, give K20 medians of 4.081 and 4.096 ms for `region-nt-all`, compared
with about 6.76 ms for the original scatter. Setup is excluded; routing,
arithmetic, all output stores, and the store fence are included.
Huge-page-preferred buffers are used when the payload exceeds the existing
allocation threshold. The new schedule needs one additional region-sized
temporary (2 MiB at K20).

The same change loses at smaller sizes. The original scatter remains best:
K16 is 0.132 ms and K18 is 0.585 ms in the holdout. Streaming stores at these
sizes evict useful locality. Full-code gather, a separate full-code reorder,
source-sorted buckets, explicit prefetch, and hardware gathers also lost.
These negative experiments remain available in the research frontend.
No production dispatch or installed API has changed.

## Iteration 3: rank analysis and its limitation

A sharper return-event bound separates constraints on the random state map
from constraints on the fresh input diagonal. This recovered roughly 130
bits in one previously failing GF256 screen, but did not close the bound.
Exact local small-field tests check the separation.

An independent audit also identified a composition limitation: a maximum
over continuation ranks cannot in general be moved outside a sum over
future outer-support types. Those types determine the initial shortened
source space. Results using that interchange are heuristic diagnostics, not
certificates. The linear sum-envelope remains valid, but looser. No complete
certificate relies on the problematic maximum.

## Iteration 4: smaller coefficient field, wider state

The current proof candidate uses GRS[16,8] over GF16, four state symbols,
and groups of 16 symbols. It retains the one regional routing stage and
fresh nonzero input scalars. Its state update is a uniform GL4(GF16) map.
The expansion rows are explicitly A_j=(1,j,j^2,j^3), for j in GF16, with
feedback C=A^T. This is a construction change, not merely a kernel change.

GF16 with four state symbols and GF256 with two both have 65,536 states.
However, the GF16 outer has far fewer words per active row, reducing the
cost of the whole-code union bound. This also removes the earlier rare-window
first-moment obstruction at 10%: its leading threshold q^(-1/2) is now
25%, rather than 6.25%. This observation alone is not a distance proof.

For exactly one active input in an inner step, the exact joint moment of
nonzero incoming state, return to zero, and output weight is

    J_1(z) = (13*z^15 + z^16 + 2*z^14) / (16*(16^4-1)).

The first generic-return screen, summed over active-row classes 4 through
80, reached only 38.66 bits at K16 and 10% distance. The exact singleton
return, exact source-support counts for classes 1 through 4, and a stronger
multi-input return bound close this gap. The latter uses CA=0, so I+AC is
invertible, to bound how often a return can simultaneously erase outputs.

All active-row classes are now rigorously covered:

| Classes | Method | Combined failure upper bound |
|---|---|---:|
| 1--127 | Directed local/region transfer and source enumeration | 1.763459244866e-13 |
| 128--5997 | Positive-coefficient Chernoff bound, eight endpoint witnesses | 3.740908268730e-22 |
| 5998--8192 | Exact-integer geometric tail, z=1/52 | 2^-64 |

The decimals in this table round the certified bounds upward. Their sum
is below 2^-42. Thus the ideal GF16 ensemble at K=65536, N=131072 has
minimum symbol distance at least 13108 except with probability below
2^-42. This is a complete finite certificate, independently audited; it is
not a claim that each individual benchmark seed has been distance-tested.
It meets the requested 10%/40-bit target. K18 and K20 have not yet received
this full replay.

The first native implementation passes all six Release and ASAN/UBSAN test
targets, including scalar/SIMD adjoints over all four field widths, exact-map
schedule comparisons, and natural-size/alignment guards. Serial holdouts
measure about 0.183 ms at K16, 0.760--0.773 ms at K18, and 5.36--5.38 ms
at K20 on 128-bit records. These use precomputed setup and include all
encoding work; setup itself is about 2.3/9.4/39 ms respectively.

GF16 arithmetic will be embedded into the AES GF256 representation for
GFNI. Raw byte labels 0 through 15 are not the required subfield. Scalar
extension then provides encoders over GF(2^m) for m=8,32,64,128 without
general multiplication in those larger fields. The resulting guaranteed
symbol distance is inherited from GF16; it does not improve merely by
increasing m.

## Iteration 5: exact inner FFT, rejected on performance

Factoring each inner map through the basis (1,X,X^2+X,X^3+X^2) halves its
GFNI multiplication count. The basis conversions are XOR-only, but the
implementation also requires two four-by-four record transposes. Isolated
basis-vector checks and full scalar adjoints confirm the map is unchanged.
The full-encoder A/B loses at every size: approximately 209 versus 185 us
at K16, 882 versus 803 us at K18, and 5.88 versus 5.17 ms at K20 in that
matched sweep. The direct map stays the default; the optional FFT path is
retained only as a measured research alternative.

A subsequent full-mode holdout of the direct kernel, with two setup seeds
and 101 measured calls per seed, gives 0.1857--0.1861 ms at K16,
0.7643--0.7645 ms at K18, and 5.0506--5.0638 ms at K20. Earlier campaigns
are retained above rather than attributed to a construction change.

## Iteration 6: packed coefficients, same sampled code

The direct four-state kernel now has an optional `packed` schedule. It loads
the sixteen diagonal coefficients once per step and the sixteen GL4 entries
once per update, then expands their rows with fixed byte shuffles. The
previous direct kernel repeats those loads and broadcasts. Setup storage,
the sampled code, and its certificate are unchanged.

Two reversed-order, full-mode holdouts use 101 calls per seed:

| Size | Packed, seeds 113 / 7919 | Direct, same seeds |
|---|---:|---:|
| K16 | 177.191 / 176.299 us | 184.274 / 183.984 us |
| K18 | 739.249 / 738.999 us | 767.813 / 767.943 us |
| K20 | 5.284411 / 5.281827 ms | 5.401560 / 5.375081 ms |

These matched comparisons show about 4% improvement at the smaller sizes
and 2% at K20. They do not establish improvement over the separate earlier
5.05 ms campaign. Small sizes use scatter; K20 uses `region-nt-all`.
Both versions remain selectable, and `direct` remains the explicit baseline.

## Iteration 7: close K18 by retaining the spectrum and widening the state

The four-state proof does not carry over automatically to K18. Its sparse
middle becomes the limiting region, even though one-row classes retain
over forty bits. This is failure of the upper-bound method, not evidence
that the code itself has poor distance.

The same outer and route, with six state symbols and uniform GL6 refresh,
now has a complete [K18 certificate](GF16_SIZES.md). The proof also retains
more of the exact outer spectrum through a positive exponential envelope.
At K=262144 and N=524288 it proves minimum symbol distance at least 52429
with setup-failure probability at most 8.272540769181098e-14 < 2^-43.

Coverage is complete: classes 1--223 use directed regional coefficients,
224--23991 use eight endpoint witnesses, and 23992--32768 use an exact
integer dense-tail comparison below 2^-80. An independent audit and full
replay passed. The receipt verifier also checks the displayed endpoint
bounds and their consistency with interval contributions. The original
K16 certificate and its witnesses remain unchanged.

The separate six-state native kernel passes all seven Release and ASAN/UBSAN
tests. Its first K18 phase measurement is 1.040 ms: 0.877 ms for inner plus
routing and 0.163 ms for the outer. The matched four-state packed run is
0.725 ms, but lacks a K18 certificate. Six-state `region-nt-all` loses at
1.493 ms. Setup is 14.3 ms and is excluded from encoding time.

An exact six-state additive-FFT variant is implemented and passes the same
tests, including all input/state basis vectors and minimum-alignment cases.
It saves twenty-four GFNI512 products across the two fixed maps. Its timing
was temporarily blocked by the benchmark host refusing SSH connections;
no performance advantage is claimed from operation counts alone.

## Iteration 8: distinguish useful proof refinements from dead ends

The [fresh-feedback alternative](RESET_FEEDBACK.md) has a valid exact
state-weight transfer, but its first-moment screen still fails at K18/K20.
It is parked without a native implementation or performance claim.

Shortening the inner from sixteen to eight symbols preserves CA=0 on an
eight-point additive subspace. However, the four-state K20 screen still
fails badly in the sparse middle. `theory/gf16_width_screen.py` records
the experiment; three tests check the unchanged t16 limit and local counts.

For the existing t16 construction, `theory/gf16_return_spectrum.py` enumerates
the exact return-output spectrum for input weights through four. If w<=s,
the selected MDS feedback columns are independent, so CX is nonzero. A
return fixes the refreshed state to CX, and emits (I+AC)X. Enumerating
projective input values gives 1,800 / 126,000 / 6,142,500 cases for w=2/3/4.
The source file computes integer histograms rather than storing raw dumps.

These exact transfers improve the bound but do not close K18 with four
states: the a=64 screen remains about -27.56 bits. This is another negative
upper-bound screen, not a distance counterexample. Three tests verify
singleton agreement, projective normalization, and complete output moments.
The existing certificates do not depend on either new floating screen.

## Iteration 9: measure the six-state FFT and change the outer

The benchmark host is available again. Two-seed full-mode holdouts, with
101 calls per seed, compare the exact six-state FFT against its direct map:

| Size | FFT, seeds 113 / 7919 | Direct, same seeds |
|---|---:|---:|
| K16 | 247.891 / 248.052 us | 255.886 / 253.592 us |
| K18 | 1013.356 / 1014.797 us | 1069.560 / 1073.166 us |
| K20 | 6.499350 / 6.479625 ms | 6.528735 / 6.499592 ms |

Small sizes use scatter; K20 uses `region-nt-all`. Setup is excluded.
K18 improves by about 3--5.5% across the matched phase and full campaigns;
K20's improvement is below 1%. The six-state K20 timing is not a certificate.
Four states remain cheaper and certified at K16.

Three construction screens followed. A [16,7] outer with four states loses
rate yet still misses forty bits in the sparse union; [16,6] has more slack
but writes one-third more symbols. [OUTER_RATE.md](OUTER_RATE.md) preserves
those partial screens. Scaling the fixed feedback permits states above
eight without singular return maps. Fourteen states look promising in
sampled classes, but have no complete certificate or native timing; see
[LARGE_STATE.md](LARGE_STATE.md).

The successful change is a longer, independently sampled systematic
Toeplitz[64,32] outer over GF16. It retains rate 1/2 and needs only three
state symbols. The [complete K20 proof](K20_RANK_CONTINUATION.md) gives
minimum distance at least 209716 with failure below 2^-80. Fifteen sparse
classes and eleven endpoint intervals cover all 32768 active-row counts.
A separate 120-digit directed replay confirms the bound. Two states fail
the present bound in the middle, despite excellent sparse screens.

The proof averages the complete outer-value counting measure before
analyzing the inner. That stronger identity also removes the independent
input-scaling layer. This is an ensemble first-moment statement, not a claim
that each fixed message has random systematic coordinates. One 64-region
routing stage and independent uniform GL3 updates remain. The native
implementation passes all eight Release and ASAN/UBSAN test targets and
full-size dense-reference comparisons. Initial full timings are about
0.278/1.13/6.81 ms; the longer outer is costly, so this is not yet a speed
improvement over the smaller-outer candidates. Separate directed replays
also close the same construction at K16/K18 below 2^-84/2^-82. Both received
independent 120-digit audits. Every reported Toeplitz size now has a proof.

## Iteration 10: stop reproducing the binary baseline

The previous discussion corrected the objective but did not itself execute
a new experiment. This continuation revalidated the files and moved to
higher-rate and higher-distance gates. The half-rate/10% native work remains
parked; no encoder benchmark or production change accompanies this iteration.

The higher-rate rectangular Toeplitz and GRS screens remain unclosed after
replacing separate return/stay estimates with an exact joint local transfer.
The transfer now has a positive integer input/output enumerator, eliminating
signed floating cancellation. Exact tests cover GF16 with three state
symbols and independently enumerate small-field cases. This is a reusable
proof calculation, not a new whole-code certificate.

More importantly, the full-uniform-permutation high-distance proposal is
actually wrong for a small fixed outer. A causal inner preserves every
suffix-shortened subspace. With GRS[16,8], the 30% target fails with probability
at least 1-2^-282 already at K16, independently of width or coefficient
choices. The exact finite proof and tests are in
[HIGH_DISTANCE_GLOBAL.md](HIGH_DISTANCE_GLOBAL.md). The proposal is withdrawn.
At K20, even MDS[256,128] cannot obtain 40 bits against this suffix event;
the necessary-only half-rate length screen first permits length 378.

Regional routing remains essential for the small-outer design. Removing
all per-step state refresh initially looked attractive, and GF256 direction
cosets repair the one-active-row issue. However, two rows can cancel through
repeated injection labels and matching diagonal ratios. The first two
Bonferroni terms give an actual failure lower bound above 0.0017 at K20 for
the GRS[16,12], period-16 proposal. It is rejected rather than promoted.

The next candidate used GRS[32,24], one 32-region route, sixteen-symbol
inner steps, four state symbols, independent nonzero GF256 input scalars,
and a uniform GL4 state-basis change only at region boundaries. It removes
the per-step random matrix work. A one-active-row return-count bound gives
about 45.56 / 43.56 / 41.56 bits at K16/K18/K20. Those are contributions from
one class, not full certificates. The fixed outer transpose also has a
verified [34-product circuit](GRS32_24_COST.md), with 2297 exact comparisons.

The multi-row gate then rejected this period-16 version. On seventeen
matching regional labels, two shortened outer rows have eighteen free
dimensions and only seventeen cancellation equations. A nonzero confined
pair exists for every diagonal and every regional matrix. Exact Bonferroni
bounds give setup failure at least 6.8736e-9 / 1.3202e-7 / 1.9696e-6 across
the three sizes. Period thirty-two still fails forty bits at K18/K20.
These lower bounds received an independent dimension and arithmetic audit.

The next repair candidate expanded the label pattern to sixty-four
positions using four local basis choices. It keeps the same one-row
argument and random region-boundary mixing. The complete tested family of
matching-label paired-return windows, summed over all outer weights, has
upper bound 7.20223e-16 at K20. This is one event family, not the full
two-row contribution or a distance certificate. The exact local/global
and adjoint reference also passes: phase changes need four fixed GF256
products and four XORs rather than five. No native kernel was built.

The broader multi-row gate rejected this repair as well. At K20, choose
37 outer rows shortened to nine coordinates. Their source dimension is 37.
If the row labels fall in nine aligned windows of 15520 positions each,
zero state at all nine window ends imposes at most 36 equations. A nonzero
word remains inside 139680 positions, below the forbidden threshold 139878.
Exact box-event Bonferroni bounds give setup failure above 2^-38, with
independent arithmetic and dimension audits. This witness is independent
of diagonals, regional matrices, and injection period. Increasing the
period cannot fix it.

This turn therefore establishes actual construction limits, not another
certificate for a guarantee already available from binary SPIN. The
verified circuits and narrow bounds are preserved as building blocks,
but no new one-route profile meets the useful complete-certificate gate.

## Iteration 11: separate dense membership from sparse geometry

The preceding turn made progress by rejecting a construction. This turn
compared concrete parity/state repairs and proved a complete dense class,
instead of attempting another half-rate/10% certificate.

The new membership argument first fixes the inner and a bad output, then
averages the regional route and independent coordinate diagonals. Tilted
Bernoulli conditioning and two arithmetic--geometric-mean bounds remove
the regional weight profile. Rational interval checks cover all inverse
input densities above the chosen cutoff, not a sampled grid of weights.
For [32,22], density at least .24 is covered with conservative class margins
578/3381/14745 bits. For [34,24], density at least .25 is covered with
527/3244/14318 bits. The proof does not depend on inner-state mixing, but
the inner must remain independent of the route and external diagonals.
See [REGIONAL_GLOBAL_ROUTE.md](REGIONAL_GLOBAL_ROUTE.md).

The structural gate now covers unequal window lengths and every starting
position. For [32,23], equal-length windows gave 65.78 bits at K20, but
unequal lengths reduce the bound to 27.17 bits. Neither number is an
actual failure lower bound. One more parity dimension changes this
comparison materially: [32,22] gives 90.85 bits, with a 70.64-bit one-row
bound. [34,24] gives 74.44 and 46.68 bits. All these classes remain separate;
different row supports and more general output patterns are not covered.
See [REGIONAL_CAPACITY_COMPARISON.md](REGIONAL_CAPACITY_COMPARISON.md).

Fixed-circuit analysis makes this a meaningful construction/implementation
comparison. [32,23] adds only one product to the prior outer, while [32,22]
adds seven. [34,24] adds fourteen but spreads that work over 24 source
symbols. Including the unchanged inner, [32,22] and [34,24] have almost
equal arithmetic volume; stride and memory effects remain unmeasured.

A separate Cauchy packet-map reference now has an exact positive local
enumerator for arbitrary tested state sizes, without enumerating every
state vector. Exhaustive GF8/GF16 checks, branch-weight bounds, and a
normalization test pass. It is not promoted: sampled ordinary source-word
bounds still fail for sparse/intermediate classes, and the elementary
native budget is higher. Its per-packet refresh is not interchangeable
with the regional candidate's once-per-region refresh.

No native benchmark, production edit, commit, or push was performed.
The complete useful certificate remains open. The next work is the
sparse covering theorem, not a new performance campaign.

## Iteration 12: repair the emission geometry and expose the missing rank

The preceding conversational reply corrected the objective but added no
evidence. This continuation revalidated the checkout and pursued useful
higher-rate constructions. It did not resume the redundant half-rate/10%
benchmark work or change production code.

A fixed reversal J of the four state coordinates makes the cubic
translation maps preserve the same pairing in both directions. Injection
and emission now use 64 distinct evaluation labels; both maps are MDS.
An unchanged nonzero state emits at least 61 symbols per full period,
instead of the archived map's 58. The revised local transpose places J
on the moment output and is not the forward step. The separate scalar
reference and negative adjoint regression preserve this distinction.
No additional field products are required, but native elapsed time is
unmeasured. See [DUAL_MDS_PHASE.md](DUAL_MDS_PHASE.md).

The sparse analysis now includes a complete support cover and an alive
rank potential that does not count a branch after a required input form
vanishes. The previous independent-packet relaxation still fails.
Even preserving the exact periodic emission of unchanged states does not
close it, so that numerical relaxation is parked rather than retuned.

The new packet argument first imposes all zero-output equations. The
remaining state update uses only diagonals at allowed output positions;
those diagonals have not yet been revealed. Deterministic root equations
also consume rank, which the previous envelope omitted. Exact zero/one-
output transitions and a designated-fresh-variable property of shortened
MDS rows are now available for a joint finite operator. This is verified
local progress, not a composed whole-code bound. See
[REGIONAL_SPARSE_COVER.md](REGIONAL_SPARSE_COVER.md).

Fresh random Toeplitz packet maps were also tested as a proof-friendly
alternative. Their exact source counter removes a large middle-class
overcount, but the tested rate-2/3 profiles remain unclosed. A rate-5/8
backup has positive isolated classes; its current arithmetic estimate is
about 2.17 times the regional [32,22] estimate, even after two further
Hankel-circuit reductions. A scalable source envelope
was proved valid and tested, but is much looser than the exact counter.
Neither finding justifies a new kernel or a lengthy backup certification.
See [RANDOM_PACKET_GATE.md](RANDOM_PACKET_GATE.md) and
[NONEMPTY_SOURCE_ENVELOPE.md](NONEMPTY_SOURCE_ENVELOPE.md).

The next action is the joint rank/root/freshness operator for [32,22] with
the dual-MDS repair. It must cover every sparse class and then be combined
with the existing dense theorem. The full goal remains open; no new
performance or complete distance claim was made.

Root replay passed 33 unit tests across the dual-MDS, source-envelope,
fresh-packet, and existing capacity/window suites. It also passed 7,290
reported sparse algebra cases, the archived 72 global/local/adjoint cases,
and the fresh-packet circuit/inverse/adjoint reference. An independent
read-only reviewer found no defect in the new local proof statements and
confirmed that the unproved global dependency issue remains explicit.
The production include/source/test directories have no changes.

## Iteration 13: return equations, sparse-envelope limits, and cheap refresh

The previous conversational reply aligned the objective but made no
technical change. This continuation made proof and construction progress
at rate 11/16, not on the redundant half-rate/10% control.

The preserved eta-scaled regional construction gained a directed complete
input-weight class. Adding independent zero-state return equations closes
[.19N,.24N), with margins 646/3616/15647 bits at K16/K18/K20. The earlier
dense theorem covers the rest above .24N. This leaves weights below .19N
unresolved for that construction; see INTERMEDIATE_RETURN_CERTIFICATE.md.

The five-state packet operator and exact row-freshness automaton were
implemented and independently audited. A nonlinear continuation action
avoids duplicating probability mass among unknown successor ranks. It
improves the evaluated bounds but does not close sparse inputs. An active-row
comparison with exact modal-conditioning costs and Bernstein bounds was
also executed. It remains too loose: independent regional maxima can use
an empty local profile in each region, losing the outer row's mandatory
span. These failed bounds are not evidence that the construction fails.

Full GL4 refreshes every16 or64 symbols were screened as changed ensembles.
Reoptimizing the tilts, rather than retaining the old witnesses, matters:
refresh every16 gives negative exponents at densities .10 and .18 but not
.01. Its cost motivates a more structured state randomizer.

The preferred new ensemble uses Cauchy-form emission and four independent
nonzero state scalars every packet. For an arbitrary incoming rank-r image,
the first moment of the surviving state-subspace size gives a concave
rank-potential bound. This has almost the root-equation contraction of a
uniform GL4 matrix, without sampling one or assuming a uniform image.
Exact checks cover 50,421 small-field subspace/diagonal pairs.

The new fixed maps decompose into four GF4 XOR-circulants, each with five
scalar products. Packing four record lanes gives five GFNI512 products
per map. No dense basis conversions or phase-state translations are needed.
The complete arithmetic count is 169/88 per source record, versus181/88
for the previous regional candidate. This is not a runtime result.

Directed replay for the new ensemble covers every W in [.11N,.24N), with
class margins 202/1833/8535. Combining the independent dense class covers all
W>=.11N. The unresolved range is now (0,.11N), not the old .19 cutoff.
The new construction and conditional proof are in DIAGONAL_REFRESH_PROOF.md;
the circuits and actual field nodes are in CAUCHY_DIAGONAL_PACKET_COST.md.
The first overly sparse set of numerical witnesses failed a replay cell;
the final frozen one-percent witness grid passes every one of the1,300
rational density cells. No optimization result is trusted as a certificate.

A new fixed-word analysis closes exactly one active outer row with
107/105/103 bits at K16/K18/K20. Jointly exposing the input and state
diagonals gives return probability at most 255^-4; a one-input return has
all sixteen outputs nonzero. Exact GRS support counting and positive
regional transfer are necessary: the looser coupled first moments failed
even at one row. The multirow versions still fail, so the remaining task
is sparse multirow coverage, not further singleton polishing.

All earlier references and binary production profiles remain unchanged.
No native benchmark, production edit, commit, or push was performed.
The full goal remains active: a complete sparse certificate and matched
performance measurements are still required.

## Iteration 14: sparse classes and a rank-aware zero-set transfer

The immediately preceding conversational turn restated the objective and
made no technical change. This continuation made proof progress on the
higher-rate candidate; it did not resume the redundant half-rate control.

Exactly two active outer rows now have a directed certificate with margins
65/72/77 bits at K16/K18/K20. Exact per-row MDS support counts, distinct
regional route slots, and packet collisions are retained. Eleven tests
cover the source quotient, route normalization, positive recurrence, and
directed arithmetic. See CAUCHY_Q2_TRANSFER.md. An exact rare-window
calculation also gives expected bad-word count above 2^111 across at most
sixteen rows. This is not bad-code probability: it explains why an ordinary
word-counting first moment cannot finish that class.

The joint return analysis enumerates all 480 pair and 2,240 triple packet
placements. Two-/three-input returns have real 255^-3 exceptions, all with
sixteen outputs. A separate rank argument imposes selected output zeros
and return constraints jointly, consuming input and state diagonals once.
It sums incoming support branches and uses projective strata only to check
feasibility, not as a word multiplicity. An independent audit and eight
exact test groups passed; see CAUCHY_JOINT_RETURNS.md.

The new global transfer uses the positive identity obtained by summing
over subsets of a witness's zero outputs. Its multiplier is tau^(D-1),
where D is the shortened outer dimension. It therefore retains rank
clipping instead of reverting to q^D message counting. Exact regional
coefficients retain both the maximum of A active slots per region and
the factor binom(A,h)/binom(L,h). The homogeneous-slot relaxation and
the exact nonempty-row truncation were tested and remain too loose.

At K20, the exact-regional rank screen with joint returns gives margins
102.21, 71.98, -19.29, -143.75, -363.71, and -752.85 bits for A equal to
3, 4, 8, 16, 32, and 64. These are floating upper-bound screens, not
directed certificates or construction-failure lower bounds. The first
version without joint returns gave only 36.66 bits at A=4 and -163.22
at A=16. An independent implementation reproduced those earlier values.

Larger outers at the same rate were also screened. GRS[64,44] gives 80.56
bits at A=8 but -67.02 at A=16; GRS[128,88] gives 64.38 at A=16 but
-275.46 at A=32. Halving the inner width to eight, with the original
conservative multi-input return terms, gives -154.20 at A=16. None fixes
the complete gap. Structured circuit estimates, not timings, are recorded
in OUTER_LENGTH_TRADEOFF.md. No such alternative is promoted.

The next promising refinement tracks actual state-coordinate support.
Empty-input packets preserve it, and one-input births have full support.
The current two-state transfer forgets the cost of reaching the sparse
states responsible for exceptional returns. A support-aware operator
must preserve that information without duplicating transition mass.

The bounded sixteen-mask kernel is now implemented for zero and one input.
It retains exact masks on empty packets, sends a one-input birth to full
support, and reproduces the y^16/tau^4 return from full support. Exact
finite-field rank/feasibility tests and ninety-six small-field dependent-form
pivot cases pass. Packet weights two and above are explicitly unsupported;
no global mask-transfer certificate is claimed. The next turn should extend
this kernel and compose it with the existing regional transfer.

The improved zero-set transfer was also tested in the homogeneous-density
rank bound. Its normalized log-upper values at densities .0001, .001,
.01, .04, and .08 are .0001075, .0008324, .0061428, .0138064, and .0037106.
They are substantially smaller than the older support cover, but all remain
positive even before finite conditioning costs. That relaxation is not a
shortcut around the unresolved multirow class.

No native encoder, benchmark, production edit, commit, or push was made.
All original references remain intact. The full goal is still open.

## Iteration 15: a verified obstruction, not just remaining proof slack

The preceding turn completed the local two-input mask reference while
answering the user's objective correction. This continuation evaluated
that reference globally and changed the next construction decision.
Half-rate/10% remains a parked control, not the target.

The 17-state regional transfer was implemented with exact support-mask
kernels for packet weights zero through two and a pooled fallback.
Independent audits checked the common-history pivot argument and the
hypergeometric convolution. Seven composition tests and 144 exhaustive
dependent-form witness checks passed. Optimized K20 floating margins at
4/8/16 active rows are 70.870/-49.025/-249.904 bits. Branch splitting
outweighs the extra state information. See MASK_REGIONAL_SCREEN.md.

A separate exact total-weight truncation retains 11A<=W<=32A in the
old two-state regional operator. At eight active rows it changes the
margin from -19.293908 to -19.290279 bits; at sixteen, from -143.748251
to -143.748240. The sixteen-row weighted mean W is 268.466, above its
minimum 176. Thus impossible low total weights do not explain this gap.
The new cauchy_weight_cut_screen.py has four passing polynomial tests.
These are floating diagnostics, not certificates.

The decisive check allows different outer supports per source row and
unequal packet-aligned windows in the 32 coordinate regions. For hit
counts c_i, the shortened source dimension is sum_i max(c_i-10,0).
Forcing zero state at every window exit costs at most 128 equations.
Any dimension above 128 therefore forces a nonzero output confined to
the windows, for every fixed assignment of the inner coefficients.

At K20, bounded searches on route seeds 0/1/2 found dimensions 123/136/133
within 152064 output positions, below the 152576 required distance.
Excluding the final 18 complete source rows removes all row-major padding;
the two witnessed dimensions remain 136/132. Their kernel dimensions are
at least 8/4. Nine tests pass. Literal fixtures replay without optimization,
and an independent inverse-permutation count confirms both witnesses.
See WINDOW_DIMENSION_SEARCH.md. No ideal-ensemble failure-probability
lower bound is claimed from these three sampled routes.

This parks the unchanged GRS[32,22]/four-state profile. Its already-proved
dense, one-row, and two-row classes remain valid. Stronger four-state
refresh, including the separately recorded full-matrix local lemma, cannot
repair the witnessed geometry. The old common-support window union did
not cover the differing supports used here.

The next rate-preserving candidate is GRS[64,44] with the existing inner.
Its recorded product count is 5.62% larger, but Boolean work and register
pressure are material unresolved costs. One bounded seed-one window
screen with 24 starts and 128 length exchanges stayed at D=0 against 256
exits. It was trapped at a flat threshold objective, so this gives no
safety evidence. A rate-2/3 GRS[33,22] alternative adds outer parity with
a smaller geometry change. State five or more needs a new circuit audit.
None of these candidates has a full certificate or a runtime result.

No native benchmark, production edit, commit, or push was made. All
reference proofs and binary production profiles remain intact. The goal
is active: no complete useful higher-rate one-route result is claimed.

## Iteration 16: complete higher-rate K20 proof

The objective correction remains binding: half-rate/10% is already proved
by binary scalar extension. This continuation closes a different point,
GRS[256,176] at rate 11/16, with one regional routing stage, width sixteen,
state eight, and independent GF256 nonzero coordinate diagonals.
The [complete certificate](FIELD256_K20_CERTIFICATE.md) gives 10% symbol
distance and at least 119 whole-code setup-failure bits at K20. This is
concrete proof progress, not an improvement to the old redundant target.

The first state-eight repair of GRS[64,44] was insufficient in the sampled
64/128-active-row classes. Tracking state-rank caps improved those bounds
by only about ten bits. A connected-window capacity moment gave a
157-bit directed bound for its specified obstruction family, but did not
constitute a distance proof. These attempts and their tests are retained
in STATE_CAPACITY_REPAIR.md, LARGER_OUTER_STATE8_SCREEN.md and
WINDOW_CAPACITY_MOMENT.md.

The larger GRS[256,176] outer preserved rate and closed the sparse screen.
A frozen 80-interval cover includes every active-row count 1 through 600.
Its Arb replay uses rational marks and outward sums, returning an upper
endpoint near 1.156409247462385e-36, below 2^-119. Independent replay,
five tests and 168 additional exact source-max checks passed. A typed
cache fix prevents a floating rank-moment result from contaminating the
rational replay.

The old pure-return envelope still failed intermediate input densities.
The successful refinement jointly ranks the selected output-zero and
state-return equations. It applies one fresh-diagonal pivot bound to the
combined matrix, rather than multiplying two bounds that reuse randomness.
Superregularity, CA=0 and the weighted incoming state flats supply the
rank bound. Independent finite-field checks passed.

The middle replay then covered all 323421 weights from 46203 through
369623 in 216 gap-free cells, with 344 bits. It uses twelve frozen
rational anchors and exact rational source-odds retargeting. Eight tests
and an independent full replay passed. See STATE8_DENSITY_CERTIFICATE.md.

The dense theorem covers weights 369624 through N with 12072 bits.
Its old GRS[32,22] rational witnesses lift by an eighth-power inequality;
all 152 intervals are accepted by exact comparisons. The common geometry
is L=6016, Kpad=1058816, N=1540096, and forbidden output weight e=154009.
For more than 600 active rows, the MDS minimum gives W>=601*81=48681,
so the three classes cover every nonzero word. Their exact upper-endpoint
sum remains below 2^-119. The aggregate checks geometry and coverage
explicitly; it does not reuse saved success flags. Independent review
confirmed the same construction and conditioning in all three bounds.

The complete generator is GF256-linear. Scalar extension preserves its
minimum symbol distance for each setup, so the bound holds simultaneously
for compatible GF(2^m), m=8,32,64,128. No partial density claim is lifted
in isolation and no union over extension-field words is needed.

Width thirty-two offers a lower arithmetic count but its present density
bound still fails, including near 14%. That is a failed upper bound, not
a realized bad-code witness. It is not covered by the width-sixteen
certificate. The large outer also has substantial Boolean and scratch
costs; the complete width-sixteen product estimate is 1959/704 per source
record. No arithmetic count is reported as a runtime win.

No native benchmark, production edit, commit or push was made. Existing
proofs and production profiles remain intact. The goal stays active:
the higher-rate proof is now complete at K20, but useful native performance
and the smaller-size points remain to be established.

## Iteration 17: native validation and exact-map optimization

Continuation audit: iteration 16 made progress by closing the full K20
certificate. The intervening conversational response only restated the
objective (no progress). This pass
revalidated the saved certificate, implemented and checked its native map,
and measured competing exact-map factorizations. It makes concrete progress
toward the higher-rate target, but does not pass the performance gate.

The isolated Cauchy256 implementation now executes the certified
GRS[256,176], width-sixteen/state-eight construction. An independent dense
Lagrange oracle checks its field adjoint, generic scalar widths, boundary
sizes, padding, alignment, scratch guards, aliases and repeated calls.
K16/K18/K20 native checks and ASAN/UBSAN pass for all five variants. The
registered CTest passes, as do the five aggregate-certificate regression
tests. The proof scripts and production kernels are unchanged.

The first field implementation took about 0.61/2.50/11.85 ms in phase
diagnostics at K16/K18/K20. Minimizing field-product counts produced too
much Boolean work and scratch traffic. Three exact-map changes improve it:
reused packed broadcasts in the dense inner; an additive-FFT outer; and
parity-only outer evaluation with systematic terms fused into final stores.
Independent scalar generation checks pass for the FFT (672) and parity
wrapper (128), followed by native tests against the separate scalar oracle.
No distribution, map, or certificate is transferred to a different code.

Final same-buffer, 101-call comparisons at seeds 113 and 7919 give
0.279/1.227--1.234/7.322--7.383 ms for the optimized field encoder.
Matched binary controls give 0.0984--0.0987/0.4608--0.4617/3.154--3.174 ms.
Setup and restores are excluded, batches alternate order, and all
benchmarks run serially. See FIELD256_NATIVE.md for source, variants,
page policy, phase diagnostics, executable hashes and reproduction.

The native map is now validated and substantially faster than its first
implementation, but still 2.3--2.8 times slower than binary SPIN. Its
higher rate and 119-bit K20 certificate remain useful proof evidence, not
a performance success. K16/K18 certificates remain open. Nothing was
promoted, committed, pushed or removed from production or prior references.

The goal remains active. Both inner and outer costs must fall; improving
the already sufficient failure margin is not the next task.

## Iteration 18: smaller outer closed and measured

Continuation audit: iteration 17 improved the exact native implementation.
The intervening scope clarification was not technical progress. This pass
revalidated pending work, rejected an exact-map optimization by measurement,
closed the smaller outer's full proof, and implemented its matching encoder.

The exact residue-coordinate FFT inner passes scalar/native, alias, and
sanitizer checks. It reduces products from 70 to 34 per packet but regresses
full encoding by about 8% at K16/K18 and 5% at K20. Two seeds and reversed
variant order agree. Keep the dense inner; the certificate is unchanged.
Cheaper state-coordinate bases have explicit deficient-rank witnesses,
so their old proof cannot be reused. State-six [128,86] still has a middle
gap under the audited joint-return cover. Higher-distance screens also
fail some multirow classes; these are failed bounds, not bad-code witnesses.
See CHEAPER_DESIGN_STATUS.md for scripts, tests, and exact mixer obstructions.

GRS[128,87], the unchanged width16/state8 inner, and one regional route
now have a full K20 certificate. L=12096, Kpad=1052352, N=1548288, and
forbidden output weight e=154828. Sparse A=1..1105 is covered by 129
directed intervals; the middle covers every W=46449..371589 in 210 cells;
the dense bound covers W=371590..N with 152 rational cells. Their bounds
sum to less than 2^-67. The MDS minimum 42 joins all A>=1106 to the weight
cover. Independent audits/replays agree; the aggregate checks geometry,
coverage and the exact sum. No partial class is lifted in isolation.
See FIELD128_K20_CERTIFICATE.md and its replay.

The new outer circuit and explicit four-row-interleaved source order
avoid the final record transpose. Arbitrary zero padding restricts the
full code and retains the distance guarantee. Independent dense interpolation,
scalar widths, adjoints, boundary/tail, route, alignment, and alias tests
pass. Large native schedule comparisons, ASAN/UBSAN, and both registered
Cauchy CTests also pass. The old native implementation and proof are preserved.

The final local regression run passes 52 tests across the new sparse,
middle, dense, and aggregate replays, the variable-distance and reduced-state
screens, joint-return and cheap-outer checks, and the older aggregate tests.

Matched 101-call, same-buffer, rotating-batch comparisons at two seeds
measure 0.288/1.13--1.17/6.01--6.10 ms for K16/K18/K20. The previous
higher-rate reference measures 0.278/1.24--1.29/7.25--7.28 ms; binary SPIN
measures 0.099/0.46/3.17 ms. Thus K20 improves by 16--17%, K18 improves
with visible timing noise, and K16 mildly regresses. The new effective
K20 rate is 128/189; output is 26.171875% shorter than half-rate SPIN.
The phase-only 5.63 ms number is not substituted for the matched result.
See FIELD128_NATIVE.md for hashes, conventions, reproduction and caveats.

This is genuine higher-rate proof/performance progress, not another
half-rate/10% proof. It still fails the competitive-encoding gate and
does not close the high-distance track. K16/K18 are timed but not certified
for their own padded geometries. Production and previous proofs remain
unchanged. No commit, push, deletion, or promotion was made.

## Iteration 19: smaller sizes closed; exact-map variants evaluated

Continuation audit: iteration 18 closed and measured the smaller K20 outer.
This round completes the naturally padded K16/K18 certificates and tests
two load/routing changes plus an alternative SIMD layout. It makes proof
progress beyond binary scalar extension, but does not pass the performance gate.

The new size replay keeps GRS[128,87], the symmetric width16/state8 maps,
fresh nonzero diagonals, and one regional route. For each column's prescribed
count h, the comparison samples Bernoulli entries at h/L before conditioning.
Its exact modal probability is at least C(L,L/2)/2^L. Replacing the loose
(L+1)^128 conditioning factor by the reciprocal central-mode factor removes
finite-size slack without changing the ensemble or discarding weight profiles.
The dense bound separately rechecks all 152 rational witnesses at exponent
1.30. The sparse transfer and K20 scripts are unchanged.

K16: L768, Kpad66816, N98304, minimum distance9831, requested rate2/3.
Sparse A1..122, middle W5161..23592, and dense W23593..98304 have component
margins67/63/1669 bits; their exact sum is below2^-63. K18: L3072,
Kpad267264, N393216, minimum distance39322, requested rate2/3. Sparse
A1..280, middle W11797..94371, and dense W94372..393216 have component
margins99/79/10260 bits; their exact sum is below2^-79. The MDS minimum42
joins the classes. All10 new tests and independent full replays pass.
See FIELD128_SIZE_CERTIFICATES.md for the lemma, scope and receipts.

The native plan now exposes six exact-map controls. Compact regional
indices fit uint16 when L<=65536, with a uint32 fallback. A separate option
batches sixteen input-diagonal coefficients. The first four-way build shows
roughly4--5% K20 improvement for their combination; K16/K18 changes are tiny.
Both tables remain in setup for controls, increasing total setup storage
rather than claiming a memory reduction.

The lane-local FFT inner transposes 32-bit chunks inside each ZMM and
keeps the two state vectors in that representation. Its diagonals transform
consistently, so the code and certificate are unchanged. It halves dense
field-product work, but adds local shuffles and boundary permutations.
Despite exact packet/layout tests, it remains about5% slower at small sizes
and loses at K20. It is retained as an optional control, not selected.

All six variants pass independent scalar/full-map, basis-state, both-direction,
alignment, alias, boundary, guard, and large-size tests. Release and
ASAN/UBSAN pass; both Cauchy CTests pass. No benchmark ran concurrently
with another benchmark. Full timings include routing, diagonals, layout
changes, outer and stores; setup and untimed restoration are excluded.

Final public comparisons measure approximately0.288/1.128 ms at K16/K18.
K20 has unexplained process-level variation: seed113 gives5.73477--5.75997 ms,
while seed7919 gives7.14506--7.46337 ms in the same final build. Matched
binary controls remain about3.17 ms. The faster isolated diagnostics and
earlier build are preserved, not substituted for these slow launches.
See FIELD128_NATIVE.md for every retained comparison and executable hash.
There is no dependable5.7 ms promotion or speed win over binary SPIN.

The user requested a status discussion after this round. The round is now
complete; no further experiments are queued. The goal remains unresolved:
we have a complete higher-rate proof at all three sizes, but not competitive
encoding performance or the higher-distance half-rate construction.
Production, old proofs, and references remain unchanged. No commit, push,
deletion, or promotion was made.

## Iteration 20: medium-alphabet proof and arithmetic comparison

The user requested GF256 versus GF32 exploration for the 10% distance
track, retaining 128-bit application records. This authorizes a new
bounded exploration after the preceding status discussion. The useful
target remains higher rate and competitive cost, not another half-rate
distance proof.

Added a genuinely field-parameterized packet/rank screen with seven
passing exact/regression tests. The fixed GF256 maps retain their ranks
under an embedding into GF32; fresh diagonals now range over all nonzero
GF32 elements. Paired screens seed the GF256 witness into the GF32 search
and also compare the same witness directly. Singleton margins improve,
but GRS64/48 states4/8 and GRS128/87 states4/6 still have unresolved middle
classes. At density.10 the latter state6 comparison improves from about
-26089 to -26087 bits. This is a failed upper bound, not code rejection.
See MEDIUM_FIELD_PROOF.md for the construction, rank argument and results.

Implemented isolated arithmetic controls: bytewise GF256, full GF32 via
VPCLMUL, full GF32 via a GF256 tower, and restricted polynomial-byte
coefficients in GF32. All128-bit records use sixteen or four field lanes,
sharing one coefficient per record. The tower preserves the cheap bytewise
action of existing GF256 fixed maps. The native polynomial representation
does not claim that embedding for free. The restricted-coefficient control
does not receive the uniform-GF32 probability bound.

All native checks pass in Release and ASAN/UBSAN, including independent
scalar/quartic oracles, basis products, tails and aliases. The arithmetic
CTest passes. Serial CPU2 measurements on Peach show K18 coefficient
passes of84.5--85.2us forGF256 versus337.7us forVPCLMUL and340.1--341.3us
for the tower. The restricted coefficient control takes145.2--145.6us.
MEDIUM_ALPHABETS.md records all sizes, methods and the executable hash.
These measurements do not include a full encoder and do not establish
a lower bound on optimized GF32 arithmetic.

Keep GF256 as the implementation baseline. The next useful task is a
middle-class analysis or changed source/packet geometry that can justify
less total work. Merely enlarging the coefficient alphabet does not pass
that gate. No two benchmarks ran concurrently; no production, frozen
certificate, commit, push, deletion, or promotion was made.

## Iteration 21: GF256 rank history and coefficient-basis outer

The user requested a goal to explore the eight-bit case and report a
meaningful result. The goal retains one regional route, 128-bit application
records, higher rate at 10% symbol distance, and competitive full cost.

Generalized the rank-cap/joint-return transfer to state counts four, six
and eight, with independent continuation-cone and rank audits. Exact
return enumeration gives six-state minima 15/13 at input weights one/two.
Feasible MDS source truncation is justified only when beta>=2; the screen
falls back outside that domain. Sharper modal conditioning contributes
827 bits at K20. Thirteen transfer/cost tests and three return-geometry tests pass.
GRS[128,86] state six still fails the sampled .10 class by 17151 bits;
the .14 class reaches +124 bits. This is not a complete certificate or a
construction counterexample. See GF8_RANK_HISTORY_SCREEN.md and its audit.

An independent opportunity removes systematic outer interpolation. Taking
normalized additive-polynomial coefficients as source coordinates gives
E_c=E_s B with B invertible. The entire padded codeword image is identical
for every setup, so existing 10% certificates with 63/79/67 bits transfer
without loss. Requested-prefix subcodes may change but inherit the same
full-space guarantees. This is a changed source interface, not byte-level
equivalence with old outputs. Production remains unchanged.

The generated coefficient transform uses 321 GFNI and 769 XOR operations
per four rows, versus 516 and 1275 for the old systematic transpose. An
additional fused leaf loader removes 16 KiB of local scratch traffic per
tile. All 4147 generated-code/basis checks pass, as do native Release,
ASAN/UBSAN, CTest, large-size, alias, guard and irregular-stride checks.

Serial matched full-encoder runs on Peach CPU2 give about .244/.956 ms at
K16/K18, versus .286/1.13 ms systematic: about 15% less time. Both new
load schedules lose at K20: gather6.06--6.08 ms and fused6.77--6.87 ms,
versus systematic5.76--5.83 ms. Phase measurements isolate the regression
to the outer; the precise cache/latency mechanism remains unproved.
The systematic K20 path remains the preferred reference. No hidden
size-dependent source-basis change is introduced. See GF8_EXPLORATION.md
for methods, both rounds, hashes, and reproducible commands.

A bounded lower-rate check finds positive sampled middle bounds for
GRS[128,80] state six at densities .06/.08/.10/.12/.14/.18. This is not a
whole-code certificate. The retained SIMD packing still uses two state
vectors. Relative to the systematic K20 reference, the audited combined
GFNI count falls only 4.84%, while output/route traffic grows 8.47%.
That is not a compelling reason to start another native kernel or full
certificate campaign. The next worthwhile inner change should reduce
actual packed work, not merely scalar state dimension.

This exploration milestone is complete. It yields a checked higher-rate
implementation improvement at two sizes, not a speed win over binary
half-rate SPIN or a new rate--distance certificate. No concurrent native
benchmarks, commit, push, deletion or production promotion occurred.

## Iteration 22: outer-dimension sweep with the inner fixed

The user authorized changing only the outer dimension to seek a higher
rate at 10% distance and forty bits. The regional route, GF256 diagonals,
width16/state8 maps, and natural row multiple64 remain unchanged.
No new implementation or benchmark was started during the proof sweep.

The GRS[128,88] K18 profile closes. It has L3008, paddedK264704, N385024,
minimumdistance38503, designrate11/16, and requestedrate32/47. SparseA1..375,
middleW15401..92405, and denseW92406..385024 give component margins
72/129/15031. The exact sum is below1.253e-22 and retains72bits.
The minimumouterdistance41 joins the sparse range to the middle range.
An independent complete replay verified geometry, arithmetic and coverage.

New parameterized helpers preserve the old frozen certificate files.
Sparse uses outward Arb matrix-polynomial arithmetic; middle retains the
old two-state transfer with new rational witnesses and modal conditioning;
dense uses exact endpoint brackets and outward log bounds. Tests cover
the component replays and aggregator corruption/coverage cases. See
GF8_OUTER_RATE_SWEEP.md and its three component notes for commands and
receipts.

K20 dimension 88 is not closed: optimized two-state sparse bounds dip to
about 4.59 bits at A=256. An audited nine-state rank-history composition
improves that point to 30 directed bits at the optimized marks, still not forty. Its middle and
dense components are complete but cannot be presented as a whole code.
Dimension89 has refined middle point bounds near-7295/-7391bits at
densities.13/.14. These are proof gaps, not construction counterexamples.

NaturalpaddinglimitsK16: dimensions87..93 allretainN98304. The first
smalleroutput, dimension94, fails sampled middle bounds by thousands
ofbits. AtK18 the closed dimension88 reducesNfrom393216to385024, only2.083%.
Do not attach dimension87 timing numbers to this new certificate. A
materially stronger inner or proof is the next useful direction; a major
implementation campaign for this small rate change is not justified yet.
No commit, push, native benchmark, or production change occurred.

## Iteration 23: repaired cubic width-eight inner

The user authorized an inner-design round after the fixed-inner rate
sweep. Three parallel investigations covered packet width/state geometry,
explicit compatible maps, and cheap local mixing. No native benchmark or
production change was made.

The useful result is the repaired embedded-GF4 cubic packet: width 8,
state 4, eta 2, A=(u+v)^3+1, C=(u+v)^3+3. The old injection constant had
singular coordinate minors. The new maps pass all 7904 square-minor checks,
satisfy CA=0, and have exact return minima [8,6,5,4,3,2,2,1]. A five-product
paired circuit yields 9 GFNI per eight outputs including diagonals and eta,
with 15 forward or 16 adjoint Boolean operations and at most 16 record shuffles. Arithmetic and
adjoint checks pass; these are not native timings.

GRS[128,83] at K20 closes 10% distance with 65 bits of whole-code margin.
L=12672, Kpad=1051776, N=1622016, requested rate 64/99, dmin=162202.
Sparse A=1..1057 has 65 bits; middle W=48661..389283 has 3303;
dense W=389284..N has 212483. The MDS minimum 46 joins the ranges because
46*1058=48668. The exact union retains 65 bits. Complete replay includes
the local map audit, geometry checks, and both joins; all 46 new tests
pass. See [the checkpoint](GF8_PACKET_CERTIFICATE.md).

Higher dimensions remain bounded research screens: repaired-map dimension 84
gives only 22.70 floating sparse bits at A=128; dimension 85 gives -359.82. The
normalized symmetric fallback is weaker. Wider packets without more
state fail the middle screen; extra output-pair mixing costs more and
does not repair return-rank deficiencies. These failures are bound
failures, not bad-code counterexamples.

The new certified point emits 4.762% more records than [128,87], but its
proposed total inner GFNI count is about 73% lower. The next action is one
faithful native implementation and serial matched whole-call benchmarks.
Do not transfer the old encoder timings to this changed construction,
and do not promote it before that gate. Old proofs and implementations
are unchanged; no commit, push, or raw-data addition occurred.

## Iteration 24: native gate and genuinely higher-rate screens

The user authorized the native cubic8 test and explicitly asked how far
the rate can rise beyond two thirds. The certified [128,83] point remains
a control, not the target endpoint. A fresh full replay preserves its
65-bit whole-code margin at K20.

The native implementation passes scalar maps, direct outer evaluation,
shortened adjoints, aliases, boundaries, all schedules, and Release plus
ASAN/UBSAN large-size checks. The first serial matched timing round does
not win: about .277 / 1.09 / 6.86 ms at K16/K18/K20. Phase measurements
locate the regression in inner plus route. The [native receipt](GF8_CUBIC_NATIVE.md)
retains all controls, including noisy systematic K20 launches. A subsequent
exact-map paired schedule removes a physical dependency between input-only
moments and incoming state. Release and sanitizer checks pass; matched
times fall14--16% to .237 / .94 / 5.71--5.86ms. K20 inner plus route drops
from4.19 to3.13ms. The first negative result remains recorded. No promotion.

Three parallel proof investigations target nominal rate3/4 at10% distance:

- [Larger outers](GF8_HIGH_RATE_OUTER.md): length128 to256 gains only
  3--4k bits against middle-bound deficits of60--90k. Removing finite
  row-selection entropy optimistically does not close the gap. Actual
  short-window obstructions are demonstrated for [32,24], not for the
  [128,96] or [256,192] candidates.
- [Rank history and returns](GF8_HIGH_RATE_INNER.md): five-state caps
  improve the A64 sparse bound from about-2131 to-2058 bits. Full GL4
  refresh barely improves marginal factors. An explicitly invalid
  no-return diagnostic reaches +8863 bits after retuning, identifying
  short return paths as a major proof cost.
- [Pair premixing](GF8_PAIR_PREMIX.md): a safe mixed-zero-flat cover
  preserves raw source rank and fresh-diagonal chronology, but its
  conservative sparse and middle screens are worse. Singleton birth
  spreading alone does not repair returns. The proof interface and
  independently audited small-field tests are retained.

These negative upper-bound screens do not prove the large-outer codes
have low distance. The GF256 random-code reference at10% is about84.1%
rate, so3/4 remains a useful target. A final bounded gate tests width8/state6
without requiring CA=0. Two superregular Cauchy pairs have checked
invertible extended maps and exact return minima. Sparse A64 improves to
about-1311 bits; the best tested middle points improve to about-18242 at
density.10 and-7927 at.14, but no candidate closes. Both profiles' algebra
and 6004 minors pass tests. This opens a valid larger-state proof interface,
not a higher-rate certificate. Next, improve activation and return behavior
jointly; the mu4 feedthrough map is a concrete unclosed candidate.
No existing proof or production code is replaced, and no commit or push
was requested.

## Preserved decisions for the prior status discussion

The user corrected the objective after iteration 9. The half-rate/10%
guarantee already follows from existing binary SPIN, and the new implementations
are slower. Extra failure margin does not make that a useful field-code result.
The previous kernel and smaller-half-rate-outer tasks are stopped.

1. Decide whether a roughly25% shorter output justifies further work on a
   currently much slower encoder. GRS[128,87] is now certified at all three
   sizes; more confidence bits alone do not improve the tradeoff.
2. If continuing this kernel, first explain the K20 process-level variation
   and obtain a stable comparison. The tested local changes are too small;
   both fixed maps and data movement remain important. A construction change
   must pass the multirow/window proof gate before implementation.
3. Alternatively, prioritize substantially higher distance at rate1/2,
   retaining the one-route and inexpensive-coefficient constraints. That
   track remains open. Use the existing binary proof and implementation as
   the baseline throughout; see [ACTIVE_TARGETS.md](ACTIVE_TARGETS.md).

All performance experiments run serially. Preserve the existing source and
measure exact-map alternatives before changing a construction. The binary
production profiles and prior proof receipts remain unchanged.
