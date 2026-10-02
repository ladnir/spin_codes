# Designing for K = 2^16

**Latest matched result:** the larger-outer K20 design reaches **3.259 ms**
with a complete **10% distance / 68.104-bit setup-failure certificate**.
See [the matched K20 result](#matched-k20-result-10-distance-in-326-ms) for
the exact construction, reproduction commands, and retained hashes.
The earlier sections below preserve the K16 starting point and the search history.

The target is a new rate-1/2 binary code at K = 65,536, with at least 10%
relative distance except with probability at most 2^-40 over setup. The
performance workload is precomputed transposed encoding of 128-bit elements.
Setup is sampled once and reused. Outer code, routing, local mixing, and inner
parameters are all design variables; the K = 2^20 configuration is a control,
not a constraint.

This is a research branch of the design. No public default or existing proof
has been replaced. The preceding common-API promotion is committed and pushed
as `5f5ef5d0` on `codex/packet-inner-tuning`. On the Ryzen 7950X, the new
RS16/GF16 prototype takes about 0.126 ms, versus 0.205 ms for that baseline.
The implementation and proof calculations are isolated here; the public
library and its defaults are unchanged.

The complete 256-bit outward replay closes the binary distance target with
49.7217337682 bits of setup-failure margin. This covers every nonzero message,
not only the sparse screening cases.

## Sparse occupancy bounds

Let q be the number of outer groups with nonzero message input. The table
bounds the expected number of low-weight outputs contributed by q = 1 or 2,
including all choices of active groups and all their nonzero messages.
The reported bits are minus the base-two logarithm of this upper bound.
They are not whole-code margins or measured distances.

| Outer and local mixing | Groups | Regions | Physical steps per region | q = 1 bits | q = 2 bits |
|---|---:|---:|---:|---:|---:|
| Current four-row BCH256 with GL32 mixing | 128 | 256 | 8 | 79.1387 | Not run here |
| Independent random-systematic rows with GL4 packet labels | 128 | 256 | 8 | 62.7931 | Not run here |
| Four parallel RS[8,4] rows over GF256 with GL32 mixing | 512 | 64 | 32 | 47.8166 | 118.3216 |
| Four parallel RS[16,8] rows over GF16 with GL16 mixing | 512 | 64 | 32 | 49.7217 | 118.3251 |

All use the selected t = 64, s = 16 inner maps, independent uniform GL16 state
updates, independent uniform column and regional shuffles, zero initial
state, and no final flush. N = 131,072 and the integer bad-weight cutoff is
13,107. The initial calculations use 192-bit outward arithmetic; the complete
RS16 replay uses 256 bits. No K20 numerical endpoint is reused. The
random-systematic and RS counts are exact expectations
over their local setup. The BCH calculation uses authenticated count bounds.

The current BCH packet design does not exhibit a q = 1 obstruction at K16.
The RS q = 2 calculation covers every support pair, retains state across all
regions, and samples regional packet positions without replacement. Its
RS8 result agrees at 192- and 256-bit precision. The complete RS union combines
these components with the higher-occupancy bound described below.

## Byte symbol RS construction

The RS candidate changes both outer arithmetic and routing geometry:

1. Partition the message into 512 groups of 128 bits. Within each group,
   arrange the bits as four rows of four GF256 symbols.
2. Encode each row with systematic RS[8,4]. Use polynomial-basis GF256 with
   modulus 0x11b; interpolate at 0,1,2,3 and evaluate at 0 through 7.
3. At each of the eight symbol positions, collect the four row bytes and
   apply an independently sampled uniform invertible 32-by-32 binary map.
4. Regard the result as four rows of 64 bits. Independently shuffle the
   64 columns of each group. Each column is a four-bit packet sent to its
   corresponding region; independently shuffle the 512 packets in each region.
5. Apply the existing t64/s16 recursive inner as the initial control. Each
   region occupies 32 physical inner steps; state persists across regions.

All random maps and shuffles above are mutually independent and sampled once
at setup. Field operations act on binary code coordinates, not on the
128-bit elements carried through the transposed encoder.

Four parallel RS rows have the symbol-support enumerator of an MDS[8,4] code
over an alphabet of size Q = 2^32. This follows by extending the GF256
generator matrix to a degree-four extension field. If A_h is its exact number
of words with h nonzero symbols, then the expected nonzero packet-support
enumerator of one group is

    sum_{h=5}^8 A_h * (((1 + 15 z)^8 - 1) / (2^32 - 1))^h.

The coefficient of z^v counts, in expectation, nonzero group messages yielding
v active packets. The sum of these coefficients is 2^128 - 1. A uniform GL32
map sends each fixed nonzero input to a uniform nonzero 32-bit string. This
gives the displayed factor and independent uniform nonzero packet labels
conditional on support. The independent column shuffle makes that support
uniform conditional on its size. These are first-moment statements; no
independence between different messages or concentration of a sampled
enumerator is assumed.

The parity map has an exact five-byte-map factorization, versus sixteen maps
for direct matrix evaluation. `rs_maps.py` checks the generator, every
four-symbol erasure pattern, the factorization, and its full binary adjoint.
The transposed encoder needs adjoint 8-by-8 binary maps, not ordinary GF256
multiplication: symmetry of the field matrix is not binary symmetry in the
polynomial basis. These adjoints can be expressed as constant GFNI affine
maps. The prototype uses this factorization together with the retained packed
inner kernel. The measurements below include packing, GL32 mixing, routing,
and the inner.

There is also a concrete geometric difference. Conditional on two active
packets reaching the same region, their positions after its uniform shuffle
are distinct and uniform. A physical inner step holds 16 packets, so the
chance they share a step is 15/(L-1), where L is the number of groups. This
is about 11.81% for the current K16 packet layout and 2.94% for the RS layout.
The q = 2 calculation above includes these collisions explicitly. The
collision probability alone is not a whole-code distance bound.

## Four bit symbol RS construction

The faster variant keeps the same 128-to-256-bit group and routing geometry.
Each row instead contains eight GF16 message symbols and uses RS[16,8]. The
field modulus is 0x13, with systematic evaluation points 0 through 7 and
parity evaluation points 8 through 15. Each aligned four-row symbol has
16 bits and receives its own independent uniform GL16 map. Its four
four-bit packets then enter the same 64-column group shuffle.

The effective MDS alphabet for counting is now Q = 2^16. Replace the symbol
length and dimension by 16 and 8, and the packets per symbol by four, in the
exact count formula. The minimum possible packet support rises from five to
nine. The exact expected counts, rather than this minimum alone, enter the
sparse bounds.

`rs16_maps.py` derives the parity map through Lagrange interpolation. Its
additive-coset structure gives a factored circuit with 19 nontrivial GF16
products and 24 basis-change XORs. Binary adjoints again replace ordinary
field multiplication in the transpose. Smaller symbol maps reduce the
mixing cost enough to offset the additional parity arithmetic in the measured
implementation.

## Covering every occupancy

The useful simplification is a pointwise outer bound, not an assumption that
the inner receives independent random bits. Let mu(x) be the expected number
of nonzero group messages whose randomized outer output equals a particular
256-bit vector x. For the RS16 construction,

    mu(x) <= (2^16 - 1)^(-8) = beta * Pr[U = x],
    beta = 2^256 / (2^16 - 1)^8,

where U is uniform on all 256-bit vectors. This includes x = 0, for which
the actual count is zero. The uniform comparison deliberately adds zero-vector
mass; it does not turn an active message group into an inactive one.

To see the bound, view the four parallel RS rows as an MDS[16,8] code over
Q = 2^16. Fix a symbol support of size h >= 9. Shortening to that support
gives dimension h - 8. Projection onto information coordinates is injective,
so at most (Q - 1)^(h - 8) codewords have exactly that support. Independent
uniform GL16 maps send each such word to any particular labeling of the
support with probability (Q - 1)^(-h). Their product is (Q - 1)^(-8).
Smaller supports have no nonzero codewords. This argument concerns expected
counts; it requires no independence between different messages.

Independent group setups let us multiply this bound across q active groups.
In the comparison distribution, each four-bit packet is active independently
with probability 15/16, with a uniform nonzero label when active. Regional
shuffling then places j active packets into a uniform j-subset of the 512
slots, without replacement. If R_j is the corresponding regional upper-bound
operator, the averaged operator is

    R(q) = sum_{j=0}^q binom(q,j) (15/16)^j (1/16)^(q-j) R_j.

For each positive Chernoff tilt lambda, the contribution of occupancy q is
at most

    binom(512,q) beta^q exp(lambda * 13107) e_zero R(q)^64 terminal_mass.

The operator keeps state between physical steps and regions; terminal_mass
sums the retained state coordinates without a final flush. The calculation
chooses a valid tilt separately for each q and then sums every q from 3
through 512. It combines that tail with the sharper exact-support q = 1
and q = 2 calculations. The independently reviewed implementation checks
full, disjoint occupancy coverage, identical inner maps and geometry, source
hashes, and outward numerical endpoints before forming the final union.

The complete fresh 256-bit replay gives 49.7217337682 bits for RS16; its
q >= 3 tail alone gives 66.1541237479 bits. In particular, for the stated
independent-uniform setup distribution,

    Pr[d_min <= 13107] <= 2^(-49.7217) < 2^(-40).

Thus the binary [131072,65536] code has minimum distance at least 13108
(just over 10%) except with that probability. This is an ensemble guarantee,
not a measurement of a sampled code's distance. The final receipt records
`fresh_replay=true`, `all_occupancies_covered=true`, and
`whole_code_certificate=true`. The replay took 196 seconds locally and agrees
with the 192-bit union. Receipt combination without fresh replay deliberately
does not set the certificate flag. A second run through `reproduce_rs16.py`,
starting without saved numerical receipts, reproduced the same endpoint in
196 seconds.

## Full transposed encoding performance

Measured on Peach's Ryzen 7950X, CPU 15, GCC 15.2, Release, with
`SPIN_TUNE=znver4`. All buffers use normal 64-byte-aligned storage. Setup,
allocation, and input generation are outside timing; every measured call
includes the full inner, routing, outer, and output stores. Calls reuse the
in-place buffer, matching the retained packet benchmark.

| K = 65,536, 128-bit elements | Median ms | Range of process medians ms |
|---|---:|---:|
| Current BCH256 packet code through `spin::Code` | 0.205158 | 0.204722–0.206054 |
| RS[8,4] over GF256 with GL32 | 0.147466 | 0.147065–0.147976 |
| RS[16,8] over GF16 with GL16 | 0.126061 | 0.124803–0.131625 |

Each entry is the median of eight process medians: four seeds and two run
orders, with five warmups and 501 measured calls per process. The serial
order for each seed was BCH, RS8, RS16, RS16, RS8, BCH. No benchmarks ran
concurrently. The RS16 prototype uses about 38.6% less time than the current
packet baseline, a 1.63x throughput improvement.

Isolated RS16 phases at seed 1, with 1001 calls, took 0.063118 ms for the
reverse inner and routing and 0.062537 ms for the outer. These isolated phase
medians describe where work remains; the full-pipeline measurements determine
the speed comparison.

Both variants pass independent scalar-reference, forward/transpose dot-product,
in-place, suffix-preservation, alignment, route-bijection, and matrix-rank
checks at K = 4096 and 65,536 with two seed pairs. ASan/UBSan runs also pass.
The generated packing and outer circuits pass exhaustive binary-coordinate
checks. The C++ API here is a research prototype; its scalar forward encoder
is a validation oracle, not an optimized public forward implementation.

## K20 performance and proof status

**Baseline correction:** the following sweep compared against the initial
promoted library, which had a K20 performance regression. It does not establish an
RS16 win over our retained best BCH implementation. The historical mode-19
winner measured 4.7706 ms; a new same-host serial replay at seeds 1 and 17
measured 4.717013, 4.771895, 4.778337, and 4.759052 ms. The interleaved library
control measured 7.051826 and 7.048058 ms. RS16's current 5.398090 ms has not
beaten that retained winner.

The source-level discrepancy was explicit: the historical K20 route uses
`_mm512_stream_si512` followed by one `_mm_sfence`, while the initial promoted
BCH kernel and this RS prototype used cached `_mm512_store_si512` at every
size. The small-size tuning report explicitly required preserving K20's
separate routing policy. Promotion timings covered K16 and K18, not K20.
Restoring only streaming stores recovered approximately 5.42 ms with prepared
huge-page buffers but did not recover the winner. The public library now restores
the complete retained mode-19 inner, routing, and tile-mode-1 BCH outer at K20.
A four-seed, reversed-order replay measures 4.92885 ms with ordinary owned buffers
and 4.68464 ms with prepared huge-page buffers. K16/K18 retain their measured
implementation. See the [library restoration results](../../../spin/README.md).
The original RS sweep below used the cached-route path. The subsequent retained-inner
experiment compares both RS paths against the restored public implementation.
See [the retained K20 campaign](../permutation_locality/encoder_cost_REPORT.md)
and [the size-selection recommendation](../permutation_locality/inner_size_REPORT.md).

The same construction runs at K = 1,048,576: 8,192 groups, 64 regions, and
512 physical t64 steps per region. Both RS variants passed the complete
scalar-reference, adjoint, route, matrix-rank, alignment, and in-place checks
at this size with route seed 1 and inner seed 18.

The serial K20 sweep used the same Ryzen 7950X, CPU 15, compiler, and flags
as K16. Each cell below is the median of eight process medians: four seeds,
two run orders, five warmups, and 501 timed calls. For each seed and memory
policy, the order was BCH, RS8, RS16, RS16, RS8, BCH. Setup, allocation,
initialization, and checksums were excluded; every call included the complete
in-place transposed encoder on 128-bit elements.

| Encoder | Default allocation ms | Matched PreferHugePages ms |
|---|---:|---:|
| Promoted BCH256 packet code with cached routing | 7.151353 | 7.035247 |
| RS[8,4] over GF256 with GL32 | 6.729851 | 6.245041 |
| RS[16,8] over GF16 with GL16 | 5.743289 | 5.398090 |

Default allocation is not identical page treatment: the BCH library advises
its owned scratch buffer even with `MemoryPolicy::Normal`; the research RS
buffer does not. The matched column explicitly selects `PreferHugePages` for
both input and scratch in every encoder. This aligns and rounds allocations
to 2 MiB and requests huge pages during setup; the request remains best effort.
The benchmark does not change global page settings.

Against this regressed library baseline only, RS16 uses 23.27% less time
with matched policies, giving 1.30x throughput. This is not an improvement
over the retained 4.77-ms result.
Its process medians range from 5.380566 to 5.409591 ms; BCH ranges from
7.021381 to 7.066455 ms. Under the existing defaults, RS16 uses 19.69% less
time. Every seed's output checksum agrees across run orders and memory policies.

Separate phase runs at seed 1, using 501 calls and `PreferHugePages`, measured
3.637322 ms for RS16's reverse inner plus routing and 2.050178 ms for its
mixed outer. RS8 measured 3.636370 and 2.848569 ms, respectively. These phase
measurements have different cache states from consecutive phases in a full
call, so they are diagnostic and must not be added to predict full latency.

K20 has a 32 MiB input and 32.5 MiB RS scratch,
before route, update, and outer-map storage; this exceeds one CCD's 32 MiB L3.
That footprint makes the missing size-specific store policy important to
test; footprint alone does not account for the lost historical result.
The same-session K16 control
still measured BCH at 0.205454–0.205834 ms and RS16 at 0.125905–0.132407 ms.

These are implementation results, not a K20 whole-code distance certificate.
The first K20 q = 1 calculation gave 58.9695337534 bits. The subsequent
all-occupancy investigation below improves that sparse bound but finds an
unresolved intermediate-occupancy gap; the K20 10%/40-bit target remains open.

## Retained streaming inner in RS16

The RS prototype now has an explicit `rs16-stream` candidate. It uses the complete
retained mode-19 reverse inner, writes whole packet cache lines with non-temporal
stores, and fences before the RS outer reads scratch. The RS outer itself is unchanged.
This candidate includes `PacketLargeInner.h` directly; it does not reconstruct or
generalize the retained inner's unrolled schedule.

The sampled code does not change. Setup samples the same update matrices and only
reverses their coefficient-byte layout for the retained GFNI kernel. The route,
64-region geometry, 260-record group stride, symbol maps, and continuous inner state
remain unchanged. Tests compare these sampled objects with the cached implementation,
then check each kernel against the scalar map and forward/transpose identity.
The candidate remains research-only and explicitly selected; no public family or
automatic size crossover has changed.

The following full-encoder replay uses the same Ryzen host, CPU 15, compiler, and
flags as above. K16 uses explicit Normal buffers; K20 uses explicit PreferHugePages
for both input and scratch in every encoder. Each result is the median of eight
process medians: four seeds, both execution orders, five warmups, and 501 timed calls.
The order was frozen RS16, rebuilt cached RS16, streaming RS16, corrected BCH,
followed by the reverse order. All benchmarks ran serially.

| Encoder | K16 ms | K20 ms |
|---|---:|---:|
| Corrected public BCH256 packet code | 0.205553 | 4.653647 |
| Frozen RS16 cached implementation | 0.125074 | 5.403656 |
| Rebuilt RS16 cached implementation | 0.125905 | 5.449126 |
| RS16 with retained streaming inner | 0.175297 | 3.309507 |

At K20, streaming RS16 uses 28.88% less time than corrected BCH, giving 1.406x
throughput. Its eight process medians range from 3.304027 to 3.319015 ms; corrected
BCH ranges from 4.646794 to 4.671659 ms. This comparison does not use the regressed
7 ms BCH implementation. Every RS seed's checksum agrees across both kernels and
the frozen binary. The frozen binary's SHA256 is
`a8dab215faabf5c4ae5c761d707a8f2ccbbd20390df09bb56c888ea202508708`.
The measured streaming-candidate binary is also preserved, with SHA256
`9b11dc91c654a9f69237c0023440a14c8a921841a6729f3686187e469701f02b`.

At K16, streaming is slower. Keep the cached implementation there: its rebuilt
median differs from the frozen control by 0.66%, and the existing 49.72-bit whole-code
certificate still applies to that unchanged construction. K20 still needs its
own all-occupancy certificate; the faster implementation does not supply one.

At K20, isolated reverse-inner/routing time falls from 3.635255 to 1.503537 ms.
The unchanged outer takes 2.083709 and 2.083990 ms, respectively. These phase runs
use seed 1, 501 calls, and explicit PreferHugePages. Their cache states differ from
the full encoder; do not add the isolated medians to predict end-to-end time.

Both RS variants and both kernel selections pass Linux release and ASan/UBSan
checks at K4096, K16, and K20. Checks include sampled-map identity, table packing,
scalar agreement, binary adjoints, in-place suffixes, guards, and legal input/output
alignments. The Windows build succeeds; its runtime checks are skipped on the
available machine because the required ISA is unavailable. The research driver
checks AVX512DQ as well as the packet kernel's ISA requirements before execution.

To select the large-size candidate explicitly:

```sh
taskset -c 15 tmp/k16-rs-build/spin_k16_rs check 1048576 1
taskset -c 15 tmp/k16-rs-build/spin_k16_rs rs16-stream 1048576 1 501 huge
taskset -c 15 tmp/k16-rs-build/spin_k16_rs bch 1048576 1 501 huge
taskset -c 15 tmp/k16-rs-build/spin_k16_rs rs16-stream-phases 1048576 1 501 huge
```

## K20 certificate investigation: target not closed

The target is the unchanged RS16 construction at K = 1,048,576 and
N = 2,097,152: 8,192 outer groups, 64 regions, and 256 two-step macros per
region. The state starts at zero, persists across all 32,768 physical steps,
and is not flushed. The bad event is a nonzero message producing weight at
most 209,715. The desired setup-failure probability is at most 2^-40.
All probabilities concern the independent uniform setup distributions stated
above. No implementation parameter was changed during this investigation.

The existing K16 49.7217337682-bit certificate remains intact. Both retained
K16 whole-replay receipts still match every pinned mathematical source.
The K20 timing of 3.309507 ms must not be advertised as a 10%/40-bit certified
operating point.

The following are negative base-two logarithms of valid first-moment upper
bounds for individual occupancies, not measured distances. The q = 1 and
q = 2 results use 256-bit Arb arithmetic; q >= 3 uses 192-bit Arb arithmetic
and exact without-replacement regional placement. Each bound includes the
choice of active groups and their nonzero messages.

| Active groups q | Bound in bits | Counting method |
|---|---:|---|
| 1 | 59.6310 | Exact expected RS16 support shells |
| 2 | 117.9008 | Exact shells and every ordered support pair |
| 19 | 81.474 | Uniform outer majorant, exact regional placement |
| 20 | 69.290 | Same |
| 21 | 58.813 | Same |
| 22 | 46.018 | Same |
| 23 | 36.581 | Same |
| 25 | 14.752 | Same |
| 32 | -56.8689 | Same |

The partial union q = 1 through 22 has 45.5383 bits. It does not cover all
nonzero messages. A negative entry means that this upper bound exceeds one;
it is not a lower bound on failure probability and does not exhibit a bad code.
The q = 32 value uses output tilt .0072. Refining the coarse factor-two tilt
grid improved it by about 23 bits, but did not reach the target. A separate
fresh 256-bit evaluation reproduces -56.86893414998522028 bits at that tilt.

Three checks isolate that gap:

- Replacing one uniform-majorized group by its exact RS16 support shells
  improves q = 32 by only 0.000176 bits. The dominant support sizes are
  46 through 51, not the excluded tiny supports.
- Selecting birth-density rows using sparse activity instead of 1/2 gives
  no material gain at the checked points. Disabling the capped-density
  refinement altogether changes the q = 32 result by less than 3e-9 bits.
- A floating diagnostic of the regional upper-envelope matrix identifies
  alternating zero/nonzero-state paths across region boundaries. Paths that
  are zero at every region boundary contribute a negligible fraction of the
  computed moment. This is a diagnostic of the bound, not a probability
  measurement or a proof that the sampled code has low distance.

Empty physical steps cannot turn a nonzero state into zero: their update is
invertible. Increasing K nevertheless increases the number of active-group
choices. At fixed q, the term binom(L,q) alone grows by about 4q bits when
L grows from 512 to 8,192. Thus the K16 certificate need not transfer to K20
even though larger regions can reduce some collision probabilities.

### Dense comparison and replay safeguards

The dense evaluator supplies a rigorous alternative to building every regional
coefficient. In the uniform outer comparison, a region contains q uniform
four-bit packets at distinct uniform slots. Instead mark each of its L slots
independently with probability p and put a uniform packet in each marked slot.
Conditioning on exactly q markers recovers the required regional distribution.
For a nonnegative path functional, removing this conditioning costs at most

    D(q,p) = 1 / [binom(L,q) p^q (1-p)^(L-q)].

The 64 independent regional conditioning events cost D(q,p)^64, without
resetting the inner state. In the unconditioned experiment, packet activity
is 15p/16. A valid iid physical-step operator H therefore gives the bound

    binom(L,q) beta^q exp(lambda * 209715)
      * D(q,p)^64 * e_zero H^(32768) terminal_mass.

Any 0 < p < 1 is valid. The endpoint p = 1 is used only for q = L.
Choosing p = q/L maximizes the conditioning probability but can give a very
poor whole bound: its unconditioned experiment admits excessively sparse
inputs. The shared-p evaluator instead searches p as a coefficient-bound
parameter. One outward matrix power per (lambda,p) serves every occupancy.
The coarse complete dense sweep has gaps and is not a passing certificate
recipe; its sampled successes do not fill the sparse q = 23 onward gap.
Independent two-parameter optimization followed by fresh 256-bit outward
checks does close individual large-occupancy points at the 10% cutoff:

| q | Output tilt lambda | Marker probability p | Bound in bits |
|---|---:|---:|---:|
| 1,024 | .247480623888 | .394890289669 | 1,116.1321 |
| 2,048 | .488911459571 | .612692069402 | 18,745.8500 |

These use the existing conditional-occupancy local operators, averaged under
the iid comparison. They cover only the two stated occupancies. In particular,
the large holes in the first coarse dense sweep are not evidence that those
occupancies intrinsically fail, and the pointwise successes are not an
interval certificate. A direct iid envelope was also checked; it did not
remove the intermediate-occupancy obstruction.

`packet_rs_k20_whole.py` checks exact coverage of q = 1 through 8,192,
the correct finite geometry, matching regenerated maps, outer-count hashes,
the required mathematical source pins, and positive outward endpoints.
A saved-receipt sum never claims a fresh certificate. A fresh replay recomputes
all components at at least 256-bit precision and sets its certificate flag
only if the complete union meets 40 bits. The current default witness grid
is explicitly exploratory and does not meet the target.

Representative reproduction commands (use fresh output paths):

```text
python -B research/workstreams/k16_design/packet_rs_k20_sparse.py q1 --precision 256 --tilts .00016 .00024 .00032 --output tmp/rs16-k20-new/q1.json
python -B research/workstreams/k16_design/packet_rs_k20_sparse.py q2 --precision 256 --tilts .00032 .00048 .00064 --output tmp/rs16-k20-new/q2.json
python -B research/workstreams/k16_design/packet_rs_k20_sparse.py tail --q-min 19 --q-max 32 --precision 192 --tilts .0036 .004 .0044 .0048 .0052 .0056 .006 .0072 .008 --output tmp/rs16-k20-new/prefix.json
python -B research/workstreams/k16_design/packet_rs_k20_whole.py --dry-run --output tmp/rs16-k20-new/whole.json
```

The local screening receipts are under `tmp/rs16-k20-close-*`; they remain
experimental data, outside the source artifact. The next useful step is to
analyze or suppress the intermediate-occupancy return cycles, not to repeat
the same tilt sweep. A larger inner state is a candidate construction change,
but needs a fresh certificate and performance measurements. The current work
does not establish that such a change is necessary.

## K20 trial: twenty-coordinate inner

This trial keeps the RS16 outer, the four-bit packets, all routing distributions,
the 64-bit physical step, and the 10% distance / 40-bit margin target. It changes
only the persistent inner state from sixteen to twenty coordinates. In particular,
the outer symbol mixers remain independent uniform GL(16,2); the inner updates
are now independent uniform GL(20,2). The initial state is zero, state persists
between all steps and regions, and there is no final flush.

`packet_inner_s20.py` preserves the sixteen selected expansion vectors and appends
the first independent quadratic truth tables in lexicographic variable-pair
order: x0*x1, x0*x2, x0*x3, and x0*x4 on the 64 six-bit coordinate labels.
The feedback map is again the transpose of this expansion. Fresh enumeration
of all 1,048,576 states verifies expansion rank 20, minimum nonzero expansion
weight 16, feedback-kernel minimum weight 6, and rank four for each packet.
The product of feedback and expansion is zero. The derived-map identity is
`b703d01842724952a2f818b58083c9f4af017590b9f7ac3ebcd11adea18f1218`.
No sixteen-bit spectrum or numerical endpoint is reused.

For a fixed nonzero entering state, the refreshed state is uniform over
1,048,575 nonzero vectors, rather than 65,535. This improves the eligible
return-to-zero atom by approximately a factor of sixteen. The output-weight
and feedback profiles also change, so that factor is not by itself a distance
proof. They are regenerated for the new maps before any transfer calculation.

The following component bounds include the active-group choices and message
counts. They use the same bad-output cutoff 209,715 as the sixteen-coordinate
investigation. Negative entries mean that the computed upper bound exceeds
one; they do not establish a bad code or an actual distance limitation.

| Active groups q | Bound in bits | Calculation |
|---|---:|---|
| 1 | 59.9683 | Exact expected RS16 shells, 256-bit Arb |
| 32 | 86.0990 | Exact regional placement; independent fresh 256-bit replay |
| 48 | -52.2882 | Exact regional placement, refined tilt grid, 256-bit Arb |
| 64 | -179.8540 | Exact regional placement, coarse tilt grid, 192-bit Arb |
| 96 | -424.8810 | Same |
| 128 | -578.5568 | Same |

At q = 32 and output tilt .0072, this improves the sixteen-coordinate bound
by about 143 bits. It does **not** close the K20 whole-code certificate: the
intermediate-occupancy bounds still leave a gap, and most occupancies have
not been checked for the new construction. The existing K16 certificate
and its mathematical sources remain unchanged.
The q = 48 refinement checks tilts .0104, .0108, .0112, and .0116; .0108 gives
the best bound in that grid. It remains about 92 bits short of the component
target. This is a gap in the present bound, not a demonstrated failure event.

The alternative fugacity calculation was also retuned using a floating search
over marker probability and output tilt, followed by fresh 256-bit outward
evaluation. It closes individual q = 1,024 and 2,048 components with about
2,437 and 22,005 bits, respectively, but leaves the intermediate range open.
Those are isolated points, not a covered interval. Fugacity removes 64 separate
regional occupancy conditions; its failure must not substitute for an exact
regional check.

No twenty-coordinate encoding kernel has been implemented or benchmarked.
The retained 3.309507 ms K20 time belongs to the sixteen-coordinate version.
A twenty-two-coordinate quadratic extension is a possible next screen: adding
x0*x5 and x1*x4 reaches full rank 22 while preserving rank-four packets and
zero feedback-times-expansion. It fits the same three-byte GFNI update layout
as the twenty-coordinate candidate, but has more ordinary state words and
four times as many states to enumerate. Neither a full certificate nor equal
end-to-end encoding cost follows from that layout observation. Check the
remaining difficult occupancies before changing the encoder.

New isolated proof tools are `packet_inner_s20.py`, `packet_rs_state_sparse.py`
(q1/q2 with an explicit state dimension), `packet_rs_s20.py` (exact regional or
fugacity components), and `packet_rs_s20_proposal.py` (floating witnesses only).
Component receipts explicitly distinguish the GL16 outer and GL20 inner and
never assert a whole-code certificate. Local numerical receipts remain under
`tmp/rs16-s20-*`, outside the source artifact.

Representative fresh replay:

```text
python -B research/workstreams/k16_design/packet_rs_s20.py --points 32 --tilts .0072 --precision 256 --output tmp/rs16-s20-new-q32.json
python -m unittest discover -s research/workstreams/k16_design -p "test_*.py"
```

The complete suite passes 123 tests, including 36 new state/map, algebra,
conditioning, and proposal checks. The existing import-graph test now runs its
graph probe in a fresh interpreter, so optional helpers loaded by earlier
toy-map tests do not create an order-dependent failure. No old mathematical
source or certificate endpoint was changed to accommodate these tests.

## Length hill climb with the unchanged sixteen-coordinate inner

This investigation changes only the message length K. The outer remains four
parallel RS[16,8] words over GF16, with independent uniform GL16 symbol maps,
64 four-bit packets per group, independent group packet shuffles, and 64
independently shuffled routing regions. The inner remains the selected t64/s16
map with independent uniform GL16 updates, zero initial state, continuous state
between regions, and no final flush. Rate is 1/2 and the target distance is 10%.
No encoder kernel, mixing depth, packet pooling, or state dimension was changed.

`packet_rs_lengths.py` specializes the existing proof routines to an explicit K.
For now it requires positive multiples of 4,096: every region then has a whole
number of the existing two-step proof macros. This is a proof-tool condition,
not a new restriction on the encoding API. `packet_rs_length_whole.py` checks
all occupancies from 1 through K/128, minimizes overlapping component bounds,
and sums their upward dyadic endpoints exactly. Saved-component assembly is
marked separately from a fresh complete replay.

At K = 98,304, a fresh 256-bit Arb replay gives **50.908994 bits** of setup-failure
margin. It covers q = 1 and 2 using the exact expected RS support counts,
q = 3 through 129 using exact regional placement, and q = 130 through 768
using the conditioned-iid fugacity bound. The bad-output cutoff is 19,660,
so the resulting minimum-distance guarantee is 19,661 out of 196,608 bits.
The one-active-group contribution gives 51.003348 bits and the combined
q >= 3 contribution gives 54.890451 bits. This is a pointwise-in-length
certificate; no claim about every intermediate K follows by interpolation.

The next completed point is **K = 114,688**, with **42.916564 bits** from a
fresh complete 256-bit replay. Its cutoff is 22,937, giving minimum distance
at least 22,938 out of 229,376 bits. Exact regional placement covers q = 3
through 164; the dense method covers every q from 165 through 896. Adding
the three output tilts .098, .110, and .122 improved the preliminary tail
union from 39.771194 to 43.574141 bits. This was a numerical-witness refinement,
not a construction change. The final q1 bound is 44.366420 bits: this replay
uses only .0032 and .0048 low-occupancy witnesses. The earlier
screen also used .0024, which improves small-support contributions and gives
55.604455 bits for q1. The narrower witness set still closes the full target;
the discrepancy is not a precision effect or a change to the code.

At the immediately next proof-macro size, K = 118,784, a fresh 256-bit check
of q = 46 through 64 gives only 36.201634 bits for that partial union, using
tilts .092 through .128 in increments of .004. Several individual bounds near
q = 52--55 are below 40 bits. This does not establish failure of the ensemble,
or optimality of the present bound. It identifies the first sampled size at
which a modest proof or construction improvement is worth testing. No full
certificate or interpolation claim is made at this length.

At K = 131,072, the refined full-coverage screen gives only 12.826922 bits,
not the requested 40. The difficult classes are near 57 active groups;
the smallest refined individual bound is about 15.55 bits. Exact regional
placement and finer output tilts improve substantially over the cheap initial
screen but do not close this target. This is a limitation of the current bound,
not an upper bound on the code's actual distance. It motivates smaller length
steps before changing the construction.

To reproduce the K = 98,304 certificate from the repository root, run this
Python code in the research environment; it reads no saved numerical receipt:

```python
import sys
sys.path.insert(0, "research/workstreams/k16_design")
from packet_rs_length_whole import replay

replay("tmp/rs16-length-k98304-new.json", K=98304,
       exact_occupancies=list(range(3, 130)), dense_min=130,
       q1_tilts=[".0032", ".0048", ".0064"],
       q2_tilts=[".0032", ".0048", ".0064"],
       exact_tilts=[".0064", ".01", ".016", ".0256", ".04", ".064",
                    ".08", ".1024", ".128", ".1536", ".18", ".2048",
                    ".256", ".3072"],
       dense_tilts=[".4096", ".6144", ".8192", "1.2288", "1.6384", "2.1972246"],
       marker_probabilities=[".55", ".7", ".75", ".8", ".85", ".9", ".95", ".975", "1"],
       precision=256)
```

For the K = 114,688 replay, use the same entry point with the following recipe.
The dense witness grid contains every retained winning witness; its extra
valid candidates may only improve the bound.

```python
from math import exp

dense_tilts = []
x = .0512
while x < 2.1972246:
    dense_tilts.append(format(x, ".12g"))
    x *= 1.25
dense_tilts.append("2.1972246")
markers = [format(1 / (1 + exp(-(-5 + .2*i))), ".12g") for i in range(52)] + ["1"]
replay("tmp/rs16-length-k114688-new.json", K=114688,
       exact_occupancies=list(range(3, 165)), dense_min=165,
       q1_tilts=[".0032", ".0048"], q2_tilts=[".0032", ".0048"],
       exact_tilts=[".00512", ".008", ".0128", ".02048", ".032", ".048",
                    ".064", ".08", ".092", ".098", ".104", ".110", ".116",
                    ".122", ".128", ".144", ".16", ".192", ".24", ".28",
                    ".32", ".4096"],
       dense_tilts=dense_tilts, marker_probabilities=markers, precision=256)
```

The displayed floating calculations only propose rational witness values,
which the evaluator parses and recomputes with outward arithmetic. They are
not numerical upper endpoints. For a stronger low-occupancy proposal, retain
.0024 as well; do not discard it based on the cumulative trial margin, since
different support classes choose different tilts.

The replay prepares the maps once, recomputes every component, and checks that
the complete source snapshot stays fixed. Preparation may import an additional
proof helper, but may not change or remove any earlier source pin. Only a
complete fresh replay below 2^-40 sets `whole_code_certificate` to true.
Both retained K16 certificate source snapshots still match their original
91/92 pins. The then-current research suite passed 163 tests. Numerical
receipts remain local under `tmp/rs16-length-*`, outside the source artifact.

## One additional state coordinate for the length hill climb

The unchanged s16 bound did not close the next sampled lengths after
K = 114,688. At K = 118,784, a finer output-tilt grid gives 36.616364 bits
for the partial union q = 50--58. Changing the valid birth-row selection
activity improves the q = 54 bound at tilt .11 by only about 0.1 bit.
Neither experiment closes the 40-bit target. These are refinements of the
bound for the same construction, not evidence of poor realized codes.

`packet_inner_small_extension.py` makes a separate, explicit construction
change. It preserves the selected sixteen expansion rows and appends the
quadratic truth table x0*x1 to obtain s17. The optional s18 control also
appends x0*x2. In each case A has 64 rows, C=A^T, and each physical state
update is independent uniform GL(s,2). The RS outer, four-bit packets,
64 routing regions, zero initial state, state continuity, and absent final
flush are unchanged. Fresh preparation enumerates all 2^s states; the
expansion distance remains 16 and the feedback-kernel distance remains 6.

This targets the repeated-return contribution directly. For a fixed
nonzero entering state and a fixed input, a uniform GL(s,2) update returns
to zero with probability 1/(2^s-1) when the feedback is nonzero, and with
probability zero otherwise. Increasing s from 16 to 17 roughly halves
each such return probability. The new coordinate also changes the fixed
maps, so the numerical comparison regenerates their full profiles instead
of merely replacing this denominator in an old bound.

The s17 construction closes **K = 131,072**, at rate 1/2 and 10% distance,
with **51.498372 bits** of setup-failure margin from a fresh 256-bit Arb
replay. It covers every occupancy q = 1 through 1,024. The bad-output cutoff
is 26,214, so the minimum-distance guarantee is 26,215 out of 262,144 bits.
This is a certificate at this specific length and for this new inner, not
an interpolation claim or a result for the unchanged s16 construction.

| Component | Occupancies | Margin in bits |
| --- | --- | ---: |
| Exact expected one-group support counts | 1 | 51.535851 |
| Exact expected two-group support counts | 2 | 124.970509 |
| Exact regional placement | 3--188 | 56.783598 |
| Conditioned-iid fugacity bound | 189--1,024 | 169.338811 |
| Complete first-moment union | 1--1,024 | **51.498372** |

The one-group contribution now dominates again; the intermediate-occupancy
obstruction has been removed at this length. All 100 replay source pins
match, as do the retained 91/92-pin K16 and 96-pin K98,304/K114,688 records.
A separate exact-rational calculation checked all 1,024 selected component
endpoints and their sum against the final outward endpoint and 2^-40.
The then-current research suite passed 163 tests, including the new map identities,
dimension-specific arithmetic, source stability, and complete-coverage checks.

To reproduce this point from the repository root without saved numerical data:

```python
import sys
sys.path.insert(0, "research/workstreams/k16_design")
from packet_rs_small_state import replay

replay("tmp/rs16-s17-k131072-new.json", K=131072, bits=17,
       exact_occupancies=range(3, 189), dense_min=189,
       q1_tilts=[".0024", ".0032", ".0048", ".0064"],
       q2_tilts=[".0024", ".0032", ".0048", ".0064"],
       exact_tilts=[".005", ".01", ".02", ".04", ".06", ".07", ".08",
                    ".09", ".10", ".104", ".11", ".12", ".14", ".16",
                    ".20", ".24", ".32"],
       dense_tilts=[".32", ".48", ".64", "1", "1.5", "2.2"],
       marker_probabilities=[".3", ".4", ".5", ".6", ".7", ".75", ".8",
                             ".85", ".9", ".95", ".975", "1"],
       precision=256)
```

The retained local receipt is `tmp/rs16-s17-k131072-whole-fresh-p256.json`,
with four adjacent component receipts. Numerical receipts remain outside
the source artifact. The s18 constructor was initially a control; the next
section records its subsequent full certificate at a larger length.

The proof-side change is not an implemented or benchmarked encoder change.
A straightforward padded implementation would increase the dense state
update from eight to eighteen wide GFNI operations. A more economical
candidate keeps the existing sixteen packed coordinates plus one ordinary
128-bit word. Splitting the sampled GL17 matrix into a 16-by-16 block and
its border gives an eight-GFNI core plus two wide GFNI operations for the
border row; the border column uses masked XORs. Payload-bit ordering needs
an additional narrow correction. This remains a fully general GL17 matrix,
not a restricted distribution. The added expansion coefficient costs one
XOR, and its feedback reuses an existing quadratic moment. These are
operation counts and an implementation plan, not measured latency claims.

For the transposed encoder, split the matrix only after forming
T=P17*M^T*P17^-1, with P17=diag(P16,1) and P16 the existing packed-state
basis. Splitting M itself swaps the border roles incorrectly. Keep all
seventeen bits of each row until extracting the sixteen-bit core; preserve
the core kernel's existing GFNI row-byte ordering.

## Larger lengths with the same outer and routing

The s18 construction closes **K = 163,840** at rate 1/2 and 10% distance,
with **44.195767 bits** of setup-failure margin. A fresh 256-bit Arb replay
covers every occupancy from 1 through 1,280. Its bad-output cutoff is 32,768,
so the minimum-distance guarantee is 32,769 out of 327,680 bits. The outer,
four-bit packets, routing distribution, and physical step length t = 64 are
unchanged. Relative to s17, the fixed maps gain the x0*x2 state coordinate,
and each physical update is independently uniform in GL(18,2).

| Component | Occupancies | Margin in bits |
| --- | --- | ---: |
| Exact expected one-group support counts | 1 | 57.531828 |
| Exact expected two-group support counts | 2 | 125.632745 |
| Exact regional placement | 3--256 | 44.195907 |
| Conditioned-iid fugacity bound | 257--1,280 | 1,156.056756 |
| Complete first-moment union | 1--1,280 | **44.195767** |

Here the intermediate occupancies dominate, rather than the one-group
contribution. The weakest individual bounds in the replay are near q = 66
and q = 69, at about 47.1 bits before summing all occupancies. This is a
certificate at one length, not a certified interval or a claim about the
unchanged s16 implementation.

The new `packet_rs_state_ladder.py` selects an explicit map constructor and
supports both the retained linear placement recurrence and a new binary
polynomial-powering evaluator. For a region of E ordered macros, each with
32 packet slots, let T_j be the local transition envelope for j occupied
slots. Define B(x) = sum_j binom(32,j) T_j x^j. The regional envelope for q
occupied slots is the coefficient of x^q in B(x)^E, divided by binom(32E,q).
Matrix multiplication retains chronological order; the T_j need not commute.
Truncating degrees above the requested q range cannot affect the retained
coefficients. `packet_regional_power.py` rounds each completed coefficient
sum and final division upward. Binary powering changes the proof calculation,
not the code, its distribution, or the quantity being bounded.

The replay records its backend and pins its sources before computation.
Independent checks verified all 1,280 component endpoints, their exact
rational sum, the final 256-bit upward rounding, and all 103 source pins.
The five retained earlier whole-code receipts still match their source pins.
The research suite passes 198 tests, including exact noncommuting-matrix
comparisons and outward enclosure checks for the new evaluator.

To reproduce this point from the repository root:

```python
import sys
sys.path.insert(0, "research/workstreams/k16_design")
from packet_rs_state_ladder import replay

replay("tmp/rs16-s18-k163840-new.json", K=163840, bits=18,
       placement_backend="binary", precision=256,
       exact_occupancies=range(3, 257), dense_min=257,
       q1_tilts=[".0016", ".0024", ".0032", ".0048", ".0064"],
       q2_tilts=[".0016", ".0024", ".0032", ".0048", ".0064"],
       exact_tilts=[".004", ".008", ".016", ".03", ".05", ".065", ".075",
                    ".085", ".09", ".095", ".10", ".105", ".11", ".12",
                    ".14", ".16", ".20", ".24", ".32", ".4"],
       dense_tilts=[".32", ".48", ".64", "1", "1.5", "2.2"],
       marker_probabilities=[".3", ".4", ".5", ".6", ".7", ".75", ".8",
                             ".85", ".9", ".95", ".975", "1"])
```

The retained local receipt is `tmp/rs16-s18-k163840-whole-fresh-p256.json`,
with four adjacent component receipts. No numerical receipts enter the
source artifact. Encoder implementation and benchmarking remain deferred.

For further length exploration, `packet_inner_quadratic_extension.py`
preserves the selected sixteen rows and extends them by independent
quadratic truth tables to s = 19, 20, 21, or 22. The s20 map matches the
retained twenty-state map exactly; s22 spans all of RM(2,6). Each candidate
uses its own full state census and independent uniform GL(s,2) updates.
`packet_rs_length_proposal.py` supplies floating tilt proposals, never
certificate endpoints, and rejects arithmetic underflow instead of silently
discarding positive mass. Larger-length screens are not substitutes for
fresh, complete occupancy replays.

### Power-of-two checkpoint: K = 2^18

The s22 construction closes **K = 262,144**, again at rate 1/2 and 10%
distance, with **42.090840 bits** of setup-failure margin. This fresh
256-bit replay covers every occupancy from 1 through 2,048. Its bad-output
cutoff is 52,428, giving minimum distance at least 52,429 out of 524,288.
The fixed expansion spans RM(2,6), with minimum weight 16; its transpose
feedback kernel has minimum weight 8. The physical step length remains 64,
and each state update is independently uniform in GL(22,2).

| Component | Occupancies | Margin in bits |
| --- | --- | ---: |
| Exact expected one-group support counts | 1 | 60.751047 |
| Exact expected two-group support counts | 2 | 126.546358 |
| Exact regional placement | 3--384 | 42.090844 |
| Conditioned-iid fugacity bound | 385--2,048 | 1,198.721804 |
| Complete first-moment union | 1--2,048 | **42.090840** |

The intermediate range still dominates. The weakest selected individual
bound is at q = 89, with 44.841583 bits before the occupancy union. This is
not a claim that s22 is necessary, or that the bound is tight. In particular,
the witness grid could be refined; the present grid already closes the
stated 40-bit target.

An independent audit checked all 2,048 endpoints, component hashes, the
exact rational union, and the final 256-bit upward rounding. It also checked
rank 22, preservation of the original sixteen rows, C=A^T, CA=0, every
four-bit packet's rank, and equality to the full RM(2,6) row span. All 103
source pins match. Across the seven retained and new whole-code receipts,
all 681 source-pin checks still pass. These remain ideal-ensemble results,
not guarantees for every setup seed or a certified interval of lengths.

To reproduce the new point from the repository root:

```python
import sys
sys.path.insert(0, "research/workstreams/k16_design")
from packet_rs_state_ladder import replay

replay("tmp/rs16-s22-k262144-new.json", K=262144, bits=22,
       placement_backend="binary", precision=256,
       exact_occupancies=range(3, 385), dense_min=385,
       q1_tilts=[".0008", ".0012", ".0016", ".0024", ".0032"],
       q2_tilts=[".0008", ".0012", ".0016", ".0024", ".0032"],
       exact_tilts=[".003", ".006", ".012", ".024", ".04", ".055", ".065",
                    ".075", ".08", ".085", ".0875", ".09", ".0925", ".095",
                    ".0975", ".10", ".105", ".11", ".12", ".14", ".16",
                    ".20", ".24", ".32", ".4"],
       dense_tilts=[".32", ".48", ".64", "1", "1.5", "2.2"],
       marker_probabilities=[".3", ".4", ".5", ".6", ".7", ".75", ".8",
                             ".85", ".9", ".95", ".975", "1"])
```

The local receipt is `tmp/rs16-s22-k262144-whole-fresh-p256.json`, with its
four adjacent components. The source and reproduction recipe are retained;
raw numerical receipts remain outside the source artifact. No encoder
implementation or performance measurement was changed for either new point.

The unchanged 128-to-256 outer at K = 2^20 remains open. Coarse intermediate-occupancy
screens did not close K = 196,608 with s19 or K = 262,144 with s20/s21.
At K = 524,288 and s22, the checked floating q = 160 bound had an estimated
margin near -318 bits. These are screening observations, not proofs of bad
codes. The s22 map exhausts the quadratic family at t = 64, so more state
within that family is not available. The next investigation separated
occupied zero-to-zero transitions from
nonzero-to-zero transitions. Its outcome and the larger-outer alternative
are recorded below; none of these changes replaces the retained K18 proof.

### Larger RS outers close K19 and K20

Keeping the **same t64/s22 inner**, larger RS constituents close both
K = 524,288 and K = 1,048,576 at rate 1/2 and 10% distance. These are
complete **derived certificates**: fresh q1/q2 calculations are combined
with a proved transfer of the authenticated K18 tail, not an all-occupancy
numerical replay. All 2,048 possible active-group counts are covered.

| K | Four parallel RS rows | Group input/output bits | Routing regions | Full margin |
| ---: | --- | ---: | ---: | ---: |
| 262,144 | [16,8] over GF(16), retained | 128 / 256 | 64 | 42.090840 bits |
| 524,288 | [16,8] over GF(256) | 256 / 512 | 128 | **120.359711 bits** |
| 1,048,576 | [32,16] over GF(256) | 512 / 1,024 | 256 | **255.335060 bits** |

The new outers use independent uniform GL32 maps on their aligned four-byte
symbols. Each group then has an independent packet shuffle; each routing
region has an independent permutation of its 2,048 slots. The retained
inner has fixed A and C, independent uniform GL22 updates at every physical
step, zero initial state, continuous state, and no final flush.

The new proof fact is a zero-state comparison under the **full uniform-input
majorant**. For fixed future routing and linear maps, a segment's output
from an entering state is uniform on a coset of a binary linear subspace.
The function z^wt has nonnegative Fourier coefficients for 0 < z <= 1,
so the unshifted subspace maximizes its average. Consequently, the actual
moment of consecutive independent segments is at most the product of their
zero-start moments, even though the encoder never resets its state.
This argument does not apply after conditioning marked packets to be
nonzero. The [transfer proof](OUTER_BLOCK_TRANSFER.md) states the precise
comparison and numerical conversion.

| Component | K = 524,288 | K = 1,048,576 |
| --- | ---: | ---: |
| Fresh uniform-majorant q1 | 120.359711 | 255.335060 |
| Fresh uniform-majorant q2 | 242.975241 | 511.581457 |
| Derived q3--2,048 tail | 270.458250 | 571.325095 |
| Complete first-moment union | **120.359711** | **255.335060** |

Margins are bits of setup-failure probability, over the stated ideal
ensemble. The bad-weight cutoffs are 104,857 and 209,715, respectively;
the guaranteed minimum distances are therefore at least 104,858 and
209,716. These are specific constructions at specific lengths, not claims
about every seed or every intervening size.

Local whole receipts are
`tmp/rs-s22-k524288-larger-outer-whole-p256.json` and
`tmp/rs-s22-k1048576-larger-outer-whole-p256.json`.
They retain the original K18 receipt hash, fresh sparse receipt hashes,
map identity, mathematical source pins, and the transfer proof/source hashes.
The union sums all disjoint dyadic endpoints exactly before outward rounding.
Raw numerical receipts stay outside the source artifact.

An independent audit recomputed both unions with integer arithmetic,
verified that rounding adds less than one endpoint ulp, and checked every
source and receipt hash. Each new sparse component has 95 source pins;
the retained base has 103. The final whole-receipt SHA256 values are:

- K19: `2001ef7ef62edcc6ee8c11d95020ae27a5d5ad0bfcf95503a826cbe41a28d980`.
- K20: `08949a8601cc58a928daa464c762f8938ad627ea4ec5f894f9f1a6c67e18fd7c`.

To reproduce the new sparse computations and derive the complete bounds
from the retained, audited K18 receipt:

```python
import json, sys
from pathlib import Path
sys.path.insert(0, "research/workstreams/k16_design")
import packet_inner_quadratic_extension as maps
import packet_larger_outer_sparse as sparse
from packet_outer_block_whole import assemble

data, record = maps.prepare(22, birth_density="capped")
components = sparse.run(data, record, precision=256,
    tilts=[".0008", ".0012", ".0016", ".0024", ".0032"])
source = "tmp/rs16-s22-k262144-whole-fresh-p256.json"
for multiplier, component in components.items():
    partial = Path(f"tmp/rs-outer-x{multiplier}-new-sparse.json")
    with partial.open("x") as handle:
        json.dump(component, handle, indent=2)
    assemble(source, partial, multiplier=multiplier, precision=256,
             output=f"tmp/rs-outer-x{multiplier}-new-whole.json")
```

The transfer deliberately accepts only the audited K18 receipt hash. Its
fresh-replay recipe remains above; a newly generated base receipt must be
audited before adding a new trust anchor. No existing certificate source
was modified. The larger outers have **not** been implemented or benchmarked;
timings for the old outer do not carry over. This is a proof anchor, not an
optimized final parameter choice.

### What the cancellation and mixer investigations showed

At K19 with the unchanged outer, q = 160, and tilt 75/1024, the baseline
bound has margin -318.290353 bits. Deleting occupied zero-state self-loops
changes this diagnostic to -185.449338 bits; deleting only two-packet
self-loops gives -186.002572 bits. Deleting live-state returns to zero
instead gives +6,406.976430 bits. The zero-only path contributes at the
+6,525.944421-bit scale. These artificial deletions are not probability
bounds for any claimed construction. They identify repeated returns to
zero as the dominant issue in this envelope, not merely the initial
two-packet cancellation.

Two separately analyzed mixer candidates did not close that witness.
Moving the existing GL update after feedback gives -364.718453 bits at
q = 160 on the tested tilt grid; adding one post-feedback transvection
to the original update gives -317.211338 bits. Neither result demonstrates
that the code itself is bad. The moved-update candidate has a new exact
two-state invariant; its return/emission correlation differs from the old
recurrence and is bounded separately. Small-occupancy return histograms
are exact through three active packets.

An isolated t128/s20-or-s22 adapter is also retained for a possible larger
inner step. Its rank, transpose, high-half handling, and one-step semantics
are tested. It has no new full-size numerical certificate. In both maps,
all two-packet restrictions still have rank seven, so this change alone
does not remove the two-packet feedback relation. The larger-outer proof
currently provides the clearer route forward.

## Alternative: random-systematic rows

`random_systematic.py` gives exact expected packet-support counts for four
independent binary [256,128] systematic rows, each with its own uniformly
sampled 128-by-128 parity matrix. Independent uniform GL4 maps on the
four-bit columns supply uniform nonzero packet labels. This removes GL32
symbol mixing but introduces random outer maps. Its setup and encoding costs
have not been measured. Structured parity implementations need separate
distribution and implementation checks before using these counts.

## Code and checks

- `rs_maps.py`: scalar RS map, five-map parity factor, and binary adjoint.
- `rs_outer.py`: exact MDS counts and GL32-to-packet support transport.
- `rs16_maps.py` and `rs_variants.py`: GF16 maps, factorization, and variant counts.
- `random_systematic.py`: exact expected support counts by active-row class.
- `packet_q1.py`: finite-geometry occupancy-one screening using retained
  inner-kernel and count-bound routines without changing them.
- `packet_q2.py`: every pair of active groups, including regional collisions.
- `rs_uniform_envelope.py`: the exact pointwise outer-domination lemma.
- `packet_uniform_tail.py`: every occupancy in a requested q >= 3 interval.
- `packet_rs_whole.py`: source-checked union and fresh complete replay.
- `reproduce_rs16.py`: from-scratch 256-bit replay with fixed tilt proposals;
  no saved experimental receipt is required.
- `packet_rs_k20_sparse.py`: explicit K20 q = 1, q = 2, and exact regional-prefix
  bounds, without modifying the K16 proof sources.
- `packet_rs_k20_dense.py`: conditioned-iid and shared-p dense bounds, with
  positive-kernel tests and optional exploratory distance thresholds.
- `packet_rs_k20_whole.py`: fail-closed whole-union validation and fresh replay;
  the K20 target is not yet closed.
- `packet_tail.py`: an alternative density-mixture screen, not used in the
  reported complete RS bound.
- `implementation/`: separate research prototype, scalar oracles, generated
  binary-adjoint kernels, correctness checks, and serial benchmark driver.
- `packet_rs_lengths.py`: explicit-length q1/q2 and exact regional components
  for the unchanged RS16/t64/s16 ensemble.
- `packet_rs_length_whole.py`: complete exact-dyadic occupancy union and fresh
  replay at those lengths; saved assembly cannot claim a fresh certificate.
- `packet_inner_small_extension.py`: isolated s17/s18 quadratic-prefix maps,
  with a fresh full state census and authenticated actual map identity.
- `packet_rs_small_state.py`: exact-placement and fugacity components plus
  fresh complete replays for those maps; the RS outer remains unchanged.
- `packet_inner_quadratic_extension.py`: isolated s19--s22 extensions of the
  same selected base, with dimension-specific full censuses and source pins.
- `packet_rs_state_ladder.py`: explicit s17--s22 full replay, recording the
  actual state dimension and selected exact-placement backend.
- `packet_regional_power.py`: outward truncated matrix-polynomial powering,
  checked against the retained linear placement recurrence.
- `packet_rs_length_proposal.py`: scaled floating length/tilt proposals;
  no numerical endpoint or certificate is produced.
- `packet_return_ablation.py`: noncertificate fixed-witness path diagnostics.
- `packet_post_feedback.py`: isolated moved-GL and post-transvection models.
- `packet_regional_log.py`: log-domain floating proposals without absolute
  entry underflow; selected witnesses still require interval replay.
- `packet_outer_geometry_proposal.py`: explicit larger-outer floating screens.
- `packet_inner_t128_extension.py`: a separate one-step t128 map adapter;
  not a wrapped pair of t64 steps and not a certified replacement.
- `packet_larger_outer_sparse.py`: fresh uniform-majorant q1/q2 for the two
  larger GF256 RS outers, using the retained t64/s22 inner.
- `packet_outer_block_transfer.py` and `packet_outer_block_whole.py`:
  authenticated tail transfer and complete derived-certificate assembly.
- `test_*.py`: 253 passing regression tests, including exhaustive small-field and small-map
  cases, exact total counts, binary adjoints, continuity across regions,
  outward arithmetic, pointwise domination, and complete occupancy coverage.

From the repository root:

```text
python -B -m unittest discover -s research/workstreams/k16_design -p "test_*.py" -v
python -B research/workstreams/k16_design/implementation/rs_outer_codegen.py --check
python -B research/workstreams/k16_design/implementation/rs16_outer_codegen.py --check
python -B research/workstreams/k16_design/reproduce_rs16.py --output tmp/rs16-reproduction/whole-p256.json
python -B research/workstreams/k16_design/packet_q1.py --outer rs --exponent 16 --precision 192 --tilts .00512 .01024 .0256 .0512
python -B research/workstreams/k16_design/packet_q1.py --outer bch --exponent 16 --precision 192 --tilts .00256 .00512 .01024 .0256 .0512 .1024
python -B research/workstreams/k16_design/packet_q1.py --outer random-systematic --exponent 16 --precision 192 --tilts .00256 .00512 .01024 .0256 .0512 .1024
```

The screening driver requires the existing research Python environment
(including python-flint and the retained kernel's dependencies), selected-map
inputs, and BCH evidence for the BCH run. It is not a standalone artifact.
Optional `--output` writes a local diagnostic receipt; experimental receipts
remain outside the committed source. Component receipts distinguish covered
occupancies from the remaining range and set `whole_code_certificate` to false.
Only the complete wrapper's successful fresh replay sets that flag to true.

On the measured Linux host, build and check the research encoder with:

```sh
cmake -S research/workstreams/k16_design/implementation -B tmp/k16-rs-build \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15 -DSPIN_TUNE=znver4
cmake --build tmp/k16-rs-build -j4
taskset -c 15 tmp/k16-rs-build/spin_k16_rs check
for seed in 1 2 3 4; do
  for mode in bch rs rs16 rs16 rs bch; do
    taskset -c 15 tmp/k16-rs-build/spin_k16_rs "$mode" 65536 "$seed" 501
  done
done
```

The optimized path needs AVX512, DQ, VBMI, and GFNI support. These commands run
one benchmark process at a time; do not run another benchmark concurrently.

To reproduce the K20 sweep after building:

```sh
taskset -c 15 tmp/k16-rs-build/spin_k16_rs check 1048576 1
for seed in 1 2 3 4; do
  for memory in normal huge; do
    for mode in bch rs rs16 rs16 rs bch; do
      taskset -c 15 tmp/k16-rs-build/spin_k16_rs "$mode" 1048576 "$seed" 501 "$memory"
    done
  done
done
taskset -c 15 tmp/k16-rs-build/spin_k16_rs rs16-phases 1048576 1 501 huge
```

## Next decisions

1. Preserve the proved and measured RS16 configuration as the K16 candidate.
   Its unchanged proof now also closes K = 98,304 and 114,688 at 10%/40 bits.
   These are individual certified lengths, not a certified interval. Keep the
   byte-symbol variant as an independently checked control.
2. The larger-outer proof now reaches K = 2^20. Preserve these complete
   certificates as anchors. Next test whether the 256-to-512 outer used at
   K19 also suffices at K20, before retaining the larger 512-to-1,024 outer.
   Another option is reducing s while retaining the larger outer. Both need
   new complete occupancy checks; the present slack is not a certificate
   for smaller parameters. The joint goal below reopens encoder work.
3. Preserve the fast s16 K16 implementation and the bordered-s17 implementation
   plan without changing production code. The retained s16 implementation has
   a measured K20 win over corrected BCH, but its full 10%/40-bit claim there
   remains open. The new larger-outer certificates do not certify that
   implementation. Revisit costs after selecting a proved larger-length design.
4. Decide separately whether to promote this family into the common public
   interface. Natural sizes, generic payload fallback, and certificates at
   other lengths require their own checks; the K16 result does not settle them.

Prior BCH/IMT and packet proofs remain in their existing workstreams. Their
certificates do not automatically transfer to this new ensemble.

## Joint Proof/Performance Search: K20 near 3.3 ms

The target is K = 1,048,576, rate 1/2, at least 10% distance and 40 bits of
setup-failure margin, with precomputed transposed encoding of 128-bit
elements near 3.3 ms on Peach's Ryzen 7950X. A timing and a certificate must
describe the same construction. Existing frozen proofs and measured kernels
remain controls, not automatic replacements for one another.

A fresh serial baseline check on CPU 15 reproduced the accepted s16
RS16-stream binary at **3.303528 ms** (seed 1, 201 calls, huge-page policy).
Its SHA256 remains
`9b11dc91c654a9f69237c0023440a14c8a921841a6729f3686187e469701f02b`.
This is a baseline check, not a new multi-seed campaign or a K20 certificate.
The benchmark acquired all existing shared benchmark locks; no other
benchmark ran concurrently.

The first implementation candidate retains that exact streaming inner and
uses eight GF16 RS[16,8] rows with independent GL32 symbol maps. Its group
dimension is 256, output length is 512, and it has 128 four-bit regions.
The averaged outer law equals the four-GF256-row RS[16,8]/GL32 law already
used in the larger-outer proof, but the nibble RS circuit can be reused.
This equality preserves the outer proof interface; changing s or K still
requires a matching full certificate. The initial implementation uses s16
to establish its cost before introducing a larger state.

Work proceeds on three fronts: screen lower-state/larger-outer bounds,
implement candidates in isolated research targets, and compare full
encoders serially against the accepted binary. Exact transitive field
multipliers are another candidate for cheaper symbol randomization; unlike
a heuristic mixer, their uniform-nonzero image law can match the proof
exactly. No candidate is promoted on a phase timing or an incomplete bound.

### Initial exact-family implementation results

`implementation/rs16x8` keeps the retained physical t64/s16 inner and its
non-temporal reverse-route stores. The outer has 256-to-512 groups. Three
independently sampled symbol-map families implement the same expected outer
measure: dense GL32, nonzero GF(2^32) multipliers, and their binary adjoints.
The adjoint family makes the transposed encoder use ordinary field
multiplication. Its nine byte coefficients occupy 0.5625 MiB at K20,
compared with 4.5 MiB for the nine-map affine tower circuit. The proof and
the distinction between field multiplication and its binary adjoint are
in `TRANSITIVE_SYMBOL_MAPS.md`.

Initial serial Ryzen 7950X measurements, CPU 15, GCC 15.2.0, `znver4`,
huge-page policy, five warmups and setup excluded:

| Candidate | Calls | Median ms | Meaning |
|---|---:|---:|---|
| Retained smaller outer, s16 | 501, seeds 1/2 | 3.293 / 3.295 | Unchanged performance reference |
| Larger outer, dense GL32, s16 | 301, seed 1 | 3.779 | First implementation control |
| Larger outer, affine field32, s16 | 501, seeds 1/2 | 3.668 / 3.661 | Exact transitive family |
| Larger outer, byte-MUL field32 adjoints, s16 | 501, seeds 1/2 | 3.450 / 3.431 | Best initial larger-outer variant |

These are development measurements, not the final matched proof/performance
result. The 256-to-512 outer with s16 does not close in the present K20
bound. A 19-bit extension and a 512-to-1024 outer with s16 are being checked.
The SIMD implementation passed the independent scalar-transpose and
forward/transpose adjoint tests for all three outer families, at K = 4096,
12288, 65536 with two seeds and K20 with another seed. The tests also check
unaligned buffers, in-place operation, guards, routing bijections, and ranks.

An optional `InnerRandomizer::TowerByte16` experiment chooses a conjugated
adjoint-field family for the state update. It preserves the fixed-state
uniform-nonzero transition law, but measured 3.512 / 3.556 ms versus a
3.392 ms dense-inner control in the same binary. The reduced coefficient
storage and instruction count did not improve the dependent inner loop.
It remains an explicitly selected experiment, not the default. Its seven
scalar/adjoint checks passed, including K20. No production sources changed.

Versioned local archives `tmp/spin-k20-goal-{dense,byte,field}-source.tar.gz`
retain the corresponding source snapshots; remote binaries live under
`/tmp/spin-k20-goal-3gmCfg/build-{dense,byte,field}`. Raw measurements and
archives remain outside the committed artifact. The public library and
all existing certificates remain unchanged.

## Matched K20 Result: 10% Distance in 3.26 ms

The joint proof/performance target closes for K = 1,048,576 and N = 2,097,152.
Under the independent setup distribution specified below, the probability
that any nonzero message has output weight at most 209,715 is bounded by
3.153017e-21, or less than 2^-68.1037. Outside that event, the minimum
distance is at least 209,716, slightly above 10% of N.
This is a whole-code setup-failure bound, not a per-message probability or
an estimate of the realized minimum distance.

The same construction encodes transposed 128-bit payloads in **3.258973 ms**
on the retained Ryzen 7950X host. The target was about 3.3 ms with at least
40 bits of margin. The result remains an isolated research implementation;
it does not change public defaults or certify other lengths.

### Construction and What Changed

Each of 4,096 groups maps 256 input bits to 512 output bits:

1. Encode eight parallel systematic RS[16,8] rows over GF(16), using the
   retained nibble-field parity circuit and common evaluation points.
2. At each of the 16 aligned 32-bit symbols, independently sample a nonzero
   scalar c in GF(2^32). The forward symbol map is the binary adjoint of
   multiplication by c. The transposed encoder therefore uses ordinary
   field multiplication, implemented with nine byte-GFNI multiplications
   and fifteen XORs per payload half.
3. Split each group into 128 four-bit packets. Independently shuffle each
   group's packets into 128 regions, then independently shuffle the 4,096
   group slots within each region.
4. Process the resulting stream with the physical t = 64, s = 20 inner.
   The state starts at zero, continues across regions, and has no final flush.
   At each step, an independently sampled uniform GL(20,2) map updates the
   old state. The fixed expansion and feedback maps are the recorded
   selected-16 maps extended by the four quadratic rows (0,1), (0,2),
   (0,3), and (0,4), with feedback C = A^T.

More explicitly, a physical step consumes x in F2^64 and state a in F2^20.
It emits x + A a and sets the next state to M a + C x, where M is that
step's sampled GL(20,2) matrix. A and C remain fixed throughout the code.

Setup is sampled once and the resulting linear code is reused. The proof
uses this ideal independent distribution. Benchmarks use reproducible
seeded setup draws of the implementation.

Relative to the 3.30 ms smaller-outer reference, the group dimension doubles
and the state grows from 16 to 20 bits. There is no extra mixing stage or
global permutation. The local symbol mixer grows from 16 to 32 bits;
the field-adjoint family replaces a dense random matrix without changing
the fixed-input law required by the first-moment proof.
[The symbol-map argument](TRANSITIVE_SYMBOL_MAPS.md) establishes that transfer.

The measured kernel uses the row-fused inner refresh, fused two-plane RS
parity/output processing, and non-temporal output stores. The complete
memory fences are inside the timed call. The implementation retains an
ordinary-store fallback when output alignment does not permit streaming stores.

### Complete Certificate and Independent Checks

The fresh replay uses 256-bit outward Arb arithmetic throughout:

| Component | Active groups q | Weakest individual contribution, bits |
|---|---:|---:|
| Exact regional placement | 1--442 | 100.25949 |
| Conditioned independent-slot bound | 443--4,096 | 68.10376 |
| Exact sum over all occupancies | 1--4,096 | **68.103757217935** |

The assembler sums all positive dyadic endpoints exactly, then rounds
the total upward once. No floating-point screening endpoint enters the
certificate. Rational tilt choices are search hints only; all bounds are
recomputed from the fixed maps and construction.

Independent audits checked the 4,096-entry coverage, the exact union,
all 104 mathematical source hashes, both component hashes, and the proof
note. They also checked the implementation's 20 expansion rows and all
64 feedback columns against the certificate. The map rank, C A = 0,
state spectrum, dual spectrum, and field-adjoint distribution agree.
The focused regression suite passes all 39 tests.

The retained first full replay is
`tmp/rs-k20-outer256-s20-whole-fresh-p256.json`, with SHA256
`97396aee84ff16c05924df80ce6b4c7722bff9ebd16643e4c15bcaf69d508e84`.
Its exact-prefix and dense-suffix component hashes are, respectively,
`d1c76dff27e89324ed5bd51a3188547e984eea3aa8d7a0e5325a29395196ac95` and
`6ca3df4814c8d17cceed264470a5e31d4a0a9be2f55ec47049c026937a6124c1`.
These raw receipts remain local; the source and fixed rational recipe
reproduce the computation without them.

### Serial Performance Campaign

Measurements use one pinned core of the Ryzen 7950X, GCC 15.2.0,
`SPIN_TUNE=znver4`, 128-bit payloads, and the large-page buffer policy.
Setup, allocation, and input preparation are outside the timing.
Every timed call includes reverse inner/routing, the outer, and store fences.
The input has N payload elements; the transposed result has K elements.

For each of seeds 1 through 4, fresh processes run in A/B/B/A order.
A is the frozen smaller-outer s16 reference; B is the new s20 code.
Each process performs five warmups and 501 measured calls. No other
benchmark runs concurrently. The table reports medians across the eight
process medians, not a selected fastest trial.

| Code | Median ms | Range of process medians ms |
|---|---:|---:|
| Retained smaller outer, s16; not certified at K20 | 3.301277 | 3.289965--3.309713 |
| New 256-to-512 outer, s20; 10% / 68.10-bit bound | **3.258973** | **3.222891--3.327525** |

Correctness checks cover an independent scalar transpose, the forward/transpose
adjoint identity, route bijections, map ranks, in-place operation, unaligned
buffers, and guard regions. The checks include K20 as well as smaller natural
sizes. Downstream application performance and setup latency are separate
measurements, not included in this result.

### Retained Implementation and Reproduction

Reproduce the complete bound from sources and fixed rational choices:

```sh
python -B research/workstreams/k16_design/reproduce_rs16_k20.py \
  --output tmp/rs16-k20-reproduction/whole-p256.json --dry-run
python -B research/workstreams/k16_design/reproduce_rs16_k20.py \
  --output tmp/rs16-k20-reproduction/whole-p256.json
```

Python, python-flint, NumPy, and SciPy are required for the real replay.
The recorded run uses Python 3.14.4, python-flint 0.9.0, NumPy 2.4.4,
and SciPy 1.18.0.
The dry run needs only the standard library and writes nothing. The real
run refuses existing output files and freshly computes every bound. Its
16 prefix tilt intervals and 15 dense tilts contain no saved endpoints.
Four reproduction tests check the recipe, dependency-free dry run, and
output protection. Receipt hashes can differ with invocation paths and
timing metadata; compare the exact numerical union and construction records.

A second fresh run through this standalone entry point reproduced all
4,096 dyadic endpoints, the final union, and the map record exactly in
516 seconds. Its receipt, `tmp/rs16-k20-reproduction/whole-p256.json`, has
SHA256 `08915a0e2340843999fceff2a94e2ba6927fa1d3375f2408f333c88f4d571642`.
Its 105 source pins add the reproducer itself; all 104 original pins are unchanged.

Build from the repository root on the supported x86 host:

```sh
cmake -S research/workstreams/k16_design/implementation/rs16x8_border \
  -B build-rs-k20 -DCMAKE_BUILD_TYPE=Release -DSPIN_TUNE=znver4
cmake --build build-rs-k20 -j4 --target spin_rs16x8_border
build-rs-k20/spin_rs16x8_border check
build-rs-k20/spin_rs16x8_border check 1048576 11 20
taskset -c 15 build-rs-k20/spin_rs16x8_border fused-nt 1048576 20 1 501 huge
```

Run benchmarks serially; use the shared benchmark locks on the shared host.
The measured executable is
`/tmp/spin-k20-goal-3gmCfg/build-border/spin_rs16x8_border`, SHA256
`acbee1ee7b332607c5b86f4ad084eadc7eba63e0060016deee1e96691aafb548`.
The complete source archive is `tmp/spin-k20-goal-campaign-source.tar.gz`,
SHA256 `3e03aa3308264a5267e0469ca41ed35104301313aa207916f07d7b7da57431d5`.
All 36 archived experiment C++ source files still match the working sources.
The local raw campaign is `tmp/spin-k20-campaign.json`; it is not a library dependency.

The next step is controlled promotion of this exact kernel. Preserve it as
the reference, add the common API around it, and benchmark each generalization.
Natural sizes, generic payload fallback, and certificates at other lengths
remain separate tasks. The old BCH, packet, and smaller RS proofs remain intact.
