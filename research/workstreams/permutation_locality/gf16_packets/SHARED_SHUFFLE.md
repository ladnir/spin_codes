# Shared Row Shuffles with GF16 Packet Mixing

This research variant recovers most of the routing advantage of the earlier
shared-shuffle encoder. Its **two-update** version now has a complete
certificate above **7% distance**, with **44.61 bits** of setup-failure
margin, and a retained precomputed transpose time of **5.863 ms** at K=2^20.
The four-update version takes **6.532 ms**, versus **7.506 ms** for the
certified independent-row version in a matched run; its 10% proof remains open.
Production defaults, the independent-row certificate, and the old lane-only
proof track are unchanged.

## Lower-Distance Search (2026-09-29)

The current follow-on targets the measured **5.863 ms two-update** version.
The requested whole-code setup-failure margin is now **20 bits**; distance
is a search variable. The four-update 10% investigation below remains a
separate, preserved partial result. Changing the update count does not
transfer its certificate.

The first complete result excludes every nonzero message whose output has
weight at most 104857. Thus the minimum distance is at least 104858, strictly
above 5% of N=2097152. The bad-setup probability is below 2^-73.38 for the
ideal setup distribution specified below. Setup is sampled once for the
whole code; this is not a per-message probability claim.

The fresh 384-bit assembler rechecked all 68 dense intervals for q=33..2048,
regenerated the all-support covers for every q=1..32, and summed their
outward endpoints. Its aggregate margin is 73.38284642743238 bits. The
separate scope-and-sum audit also passes with a requested 73-bit margin.
The complete receipt is `tmp/shared-relaxed/shared-r2-d05-complete-p384.json`,
SHA256 `463c237f1b9da08cdc8196d34cafcdf8424ed6a8f6c89532099009aaf9ef7ba7`.
Its dense input is `dense-d05-parallel.json`, SHA256
`19eb890bbcc5b53c76cde4b773a6e35492bb8cd8cf0e2e305e7d3a70356f6d65`.

The next complete result excludes weights through 125829. It proves minimum
distance at least 125830, strictly above 6%, with aggregate margin
40.68968248832404 bits. Its fresh 384-bit replay covers all 114 dense intervals
and regenerates every sparse occupancy q=1..32. The separate audit passes
with a requested 40-bit margin. Thus this operating point does not require
the proposed relaxation to 20 bits.

| Relative distance | Whole-code margin | Dense intervals | Status |
|---|---:|---:|---|
| >5% | >73.38 bits | 68 | Complete, fresh 384-bit replay |
| >6% | >40.68 bits | 114 | Complete, fresh 384-bit replay |
| >7% | >44.61 bits | 303 | Complete, fresh 384-bit replay, 2026-09-30 |

All rows describe the same measured 5.863 ms encoder, not different parameter
sets. The bounds concern the ideal setup distribution and cover all nonzero
messages. They do not estimate the realized code's actual minimum distance.
The 7% result also implies its 44.61-bit bound at every lower threshold.
The 6% row preserves the earlier, now dominated proof; these rows are not
an optimized margin curve.

The 7% certificate excludes weights through 146800, proving minimum distance
at least 146801 out of N=2097152. Its whole-code aggregate margin is
44.6143782060373124 bits. Fresh 384-bit verification rebuilt all 303 dense
intervals for q=33..2048, then regenerated every sparse occupancy q=1..32.
Their separate margins are 48.72857190017 and 44.70018670023 bits; the
certificate sums their probability bounds, rather than taking their minimum.
The independent scope/hash/sum audit passes with a requested 44-bit margin.

The complete receipt is `tmp/shared-relaxed/shared-r2-d07-complete-p384.json`,
SHA256 `60c86cf9ded774c81889a4375e317ef643460e869ecbe4e89e3be189f214c80e`.
Its dense input is `dense-monotone-d07-tight.json`, SHA256
`5006749b5d022113eb7a6c0ebe8b6a074a3b9c1aa0f3a9fe6acef74b41a615b4`.
The fresh sparse receipt is `shared-r2-d07-complete-p384.sparse.json`, SHA256
`92f33dfdc676fc3b67fe042b99c6128261887c016a57df671b5e44f3475542cb`.
All six counting-witness hashes were separately checked after verification.
The final replay used eight processes, with search caches disabled.
Every leaf was recomputed, including the search's retained hints.
The proof-tool suite passes all 405 tests.

The 6% receipt is `tmp/shared-relaxed/shared-r2-d06-complete-p384.json`, SHA256
`53278aa5f1f948bbbdbdcddbc01cab9ca2ff46c9319c43e7034b12b14b54b2d5`.
Its dense input is `dense-pruned-d06-parallel.json`, SHA256
`1b3245c73688a6f098f15e7a79d10296263a2f68953f1cffd4fd09fb34fcabde`.

The fresh 256-bit sparse run covers every q=1..32 at cutoff 104857. Summing
its exact dyadic upper endpoints gives a restricted margin above 125.57
bits. Its receipt is `tmp/shared-relaxed/sparse-r2-d05-q1-32.json`, SHA256
`7bfec72dcb6b17f0ba151849ad24f7a5d999e19b2a912da35546379a15d93b62`.
The completed assembler regenerates that sparse result rather than trusting
the saved endpoints. The first search retained the census builder's 192-bit
context instead of its requested 256 bits. That metadata issue was corrected;
all bounds in the final assembly were replayed at 384 bits. Aggregation uses
exact dyadic sums followed by explicit upward rounding for compact receipts.

The sparse follow-on also closes q=1..32 at 7% with 44.70 aggregate bits,
and at 8% with 34.80 aggregate bits. A separate 6% run closes q=33..64
with over 1236 bits. These are bounds on restricted message classes;
none supplies the dense complement. At 9%, q=16 closes, but q=32 remains
unresolved in the bounded support search.

The original dense comparison fails the requested 20-bit test at 6% and
comparison mean .032: its directed log2 bound is -8.1596. Reducing its
positive component masses, with exact verification against every shell
cap, improves that point to -669.6609. This changes the proof bound only.
The verified pruning retains 11 of the original 21 components. Its LP
success flag is not a proof premise; the assembler checks the rational
coefficients against freshly regenerated caps. With the initial counting
inputs at 7%, pruning fixes the mean-.008 diagnostic but not the mean-.032
diagnostic. At 6%, the exhaustive
dense search closed 114 intervals with no unresolved intervals, followed by
the successful fresh whole-code assembly above.

Two independent improvements target higher distance. Extending the sparse
range through q=256 would strengthen the dense occupancy constraint. Selected
all-support covers at q=128,192,256 close at 6%, but the entire prefix is not
yet covered. Larger Chernoff tilts were necessary at q=192 and 256; the
initial capped grid's failures were not proof obstructions. Separately, new
exact joint-count duals improve the BCH support bounds near weights 103--128.
Tighter counts alone do not guarantee a better
scalar comparison: rebuilding its mixture can move mass to poorly matched
activities. The new search can instead reduce the existing component masses,
preserving its geometry. Those probes alone did not close 7%; the exhaustive
replay above supplies the whole-code certificate.

With the six retained count-witness files, all eleven sampled dense means
pass at 7%. The weakest sampled bounds are below 2^-102 at mean .008 and
below 2^-111 at mean .032. The receipt is
`tmp/shared-relaxed/alternative/canonical-monotone-d07-gaps.json`, SHA256
`5ad769c6d67222cbd2305696f9e56e0555a7638c7b1b73e1ea62cec193dc8da4`.
The completed interval cover now also covers the gaps between those points.
`shared_relaxed_warm_start.py` transfers the old 6% partition and rational
witnesses to the stronger comparison. It marks every interval unresolved;
no saved bound or acceptance decision transfers to the new distance.
Five additional gap means (.012, .024, .040, .048, .096) also pass under
exactly that comparison. Their weakest bound is below 2^-86 at .096.
The separate receipt `alternative/canonical-saved-d07-gap-probes.json`
under `tmp/shared-relaxed/` has SHA256
`9953de1aa4ad6abe46bf91b94936665e8616a15c24750e93388c7ecedcd0f00f`.

The search can rebase these point witnesses onto finite intervals, but
must recompute each interval's bound. A passing midpoint can trigger
bisection; it cannot accept the interval. Near mean .032, two bisections
give a finite interval bound below 2^-50 with the reused point witness.
The old witness fails on that same interval. An optional worker cache
reuses freshly computed regional polynomials for identical inner data,
precision, output tilt, and local refinements. It still recomputes the
cell-dependent bounds. Full verification rejects this search-cache context
and rebuilds the operators independently.

An independent 384-bit regeneration of every sparse occupancy q=1..32 at
7% gives 44.7001867002305 aggregate bits. Its receipt is
`tmp/shared-relaxed/bridge/sparse-r2-d07-q1-32-p384-independent.json`, SHA256
`411430c52476cc7eed73d4db140596069af2f59f047fa89c79c5eb2533cb6e38`.
The weakest class is q=7. The assembler's search grid now includes .0064
and .012, matching the additional choices used in this independent run.
These are search parameters, not changes to the encoder or setup distribution.

The independent 384-bit sparse replay at 8% covers every q=1..32 with
34.8057764260 aggregate bits. Its receipt is
`tmp/shared-relaxed/bridge/sparse-r2-d08-q1-32-p384-independent.json`, SHA256
`8163fe2a7247ecea1319ded18fbc8526ddcb99b7c861540c6915e514a1cbe066`.
Occupancy q=1 dominates this restricted sum. The dense complement remains
open: at 7.5%, the sampled mean .008 passes, but .032 does not. A successful
point at a larger distance can have a stronger reported bound when its
search invokes additional refinements. This reflects different witnesses,
not an improvement caused by increasing the requested distance.

The bounded frontier screen uses the same six-file comparison and q>=33,
without repruning. Fresh 256-bit evaluations give the following log2 upper
bounds. Positive entries do not establish low-weight codewords; they leave
the requested probability bound unresolved.

| Requested distance | Mean .008 | Mean .032 | Mean .096 |
|---|---:|---:|---:|
| 7.5% | -206.9142 | +1130.4155 | +3512.6862 |
| 8% | +132.4370 | +2509.0717 | +7603.9758 |

The receipts are `alternative/canonical-saved-d075-frontier.json`, SHA256
`66783ecebfe5be725b49ffd10d8f350c4b78b6f30ba4730f9ec65be1e52d3ef7`,
and `alternative/canonical-saved-d08-frontier.json`, SHA256
`b8b8e85656c7a857501a1696eb1bbf902951567bdff2342715adb367dd585925`,
under `tmp/shared-relaxed/`. These screens motivated finishing 7% before
extending the distance again. Neither higher-distance screen is a whole-code
certificate; the 7.5% and 8% targets remain open.

At 7.5% and mean .032, one variance bin contributes almost the entire
failed bound. Its log2 terms are +6108.23 from the outer comparison,
-25648.56 from the inner, and +20670.75 from the output cutoff. Splitting
the two dominant variance bins four ways and regenerating their duals
improves the upper bound by only 13.91 bits, to +1116.50. Thus finer mean
cells cannot resolve this singleton gap, and the tested variance refinement
is also insufficient. The diagnostic receipt is
`tmp/shared-relaxed/strategy/d075-m032-obstruction.json`, SHA256
`6341e747e4702c9e6bf9ce2d35dff06f698e314e96758f1df6a5d6b28ccde7c5`.
The outer dual is dominated by about 168 active groups from one comparison
component. Raising the sparse/dense boundary only to 65 would not remove
that obstruction; even a larger bridge would leave the separate .096 gap.

Further bounded count experiments reached a practical plateau. Adding
complement symmetry, including at the full support limit, did not improve
the retained bounds. Exact dual recombination gained only 0.000338 bits
at support 115 and 0.0000194 bits at 116. Those diagnostic witnesses are
preserved separately; the 7% comparison's counting inputs are unchanged.

Independent reviews checked shared-support conditioning, positive
majorization, the true-zero distinction, chronological regional operators,
and exhaustive sparse/dense coverage. They found no blocking issue. The
numerical certificates still rely on the authenticated BCH premises and
outward arithmetic; these reviews are not proof-assistant verification.

`shared_relaxed_dense.py` regenerates the shared support bounds, verifies
the positive mixture exactly, and screens or covers the dense comparison.
`shared_relaxed_parallel.py` continues a full partition with separate
arithmetic contexts and, by default, rechecks every retained leaf. Its optional
`--completion-driven` scheduler refills individual worker slots and retains
in-flight intervals in every checkpoint. The scheduler passed an independent
review, unit tests, and a real one-cell Windows spawn test. An additional
`--retain-search-leaves` option carries prior leaf labels only as unverified
search hints. It discards saved endpoints, records the retained paths and
source digest, and defers their numerical checks to the final verifier.
It cannot be combined with distance retargeting. The separate
`shared_relaxed_strategy.py` assembler regenerates the model and sparse
prefix, replays the complete dense partition, and sums exact dyadic upper
endpoints. Only that exhaustive assembly can establish the whole-code claim.
`shared_relaxed_audit.py` independently checks the final scope, source hash,
and exact endpoint aggregation; it does not replace numerical replay.
Its receipts live under ignored `tmp/shared-relaxed/`. New tests and proof
sources are separate from experiment data. The encoder is unchanged, and
the old independent-row certificates remain valid for their own ensemble.

With the archived dense witnesses available, reproduce the complete 6%
result from the repository root with the commands below. Set
OPENBLAS_NUM_THREADS, OMP_NUM_THREADS, and MKL_NUM_THREADS to 1, and choose a
new output path: verification deliberately refuses to overwrite old receipts.

```text
python -B research/workstreams/permutation_locality/gf16_packets/shared_relaxed_strategy.py tmp/shared-relaxed/dense-pruned-d06-parallel.json --precision 384 --workers 4 --bits 40 --sparse-target-bits 32 --count-witnesses tmp/shared-goal/joint-constraints-retained-witnesses.json tmp/shared-goal/joint-constraints-coarse-256.json --output tmp/shared-relaxed/shared-r2-d06-fresh-reproduction.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_relaxed_audit.py tmp/shared-relaxed/shared-r2-d06-fresh-reproduction.json tmp/shared-relaxed/dense-pruned-d06-parallel.json --distance .06 --bits 40
```

Reproduce the stronger 7% certificate with the six retained count-witness
files. The worker count affects runtime, not the claimed distribution or
bound. Again, use a new output path.

```text
python -B research/workstreams/permutation_locality/gf16_packets/shared_relaxed_strategy.py tmp/shared-relaxed/dense-monotone-d07-tight.json --precision 384 --workers 8 --bits 40 --sparse-target-bits 40 --count-witnesses tmp/shared-goal/joint-constraints-retained-witnesses.json tmp/shared-goal/joint-constraints-coarse-256.json tmp/shared-relaxed/alternative/joint-115-116.json tmp/shared-relaxed/alternative/joint-117-128.json tmp/shared-relaxed/alternative/joint-108-113-gaps.json tmp/shared-relaxed/alternative/joint-103-107.json --output tmp/shared-relaxed/shared-r2-d07-fresh-reproduction.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_relaxed_audit.py tmp/shared-relaxed/shared-r2-d07-fresh-reproduction.json tmp/shared-relaxed/dense-monotone-d07-tight.json --distance .07 --bits 44
```

The [proof index](../PROOF_INDEX.md) records all retained constructions and
the local source/witness snapshot. The new contiguous sparse result below
covers q=1–22 at the 10% target with a combined restricted-class margin
above 44.20 bits. It is not a whole-code certificate.
The [progress ledger](PROGRESS.md) records the coupled-count improvement,
joint-CDF constraints, and dense-comparison experiments. The bounded search
stopped at a plateau on 2026-09-29: q=23 and the dense range remain open.
The stages below retain their original scope and reproduction commands.

## Construction and Implementation

Keep BCH[256,128], groups of four adjacent outer rows, and IMT(128,19).
The 8192 outer rows form 2048 groups. For each group, sample one permutation
of its 256 columns and apply it to all four rows. Each column becomes a
four-bit packet in one of 256 regions. Independently permute the 2048 packets
within each region. Multiply each packet by an independent uniform nonzero
element of GF16, in the polynomial basis modulo x^4+x+1. The inner uses
independent transvections, either two or four per step, and no final flush.
These distributions describe the ideal setup analyzed below.

The transpose applies the transposed field multiplier while emitting four
128-bit lanes. One 64-byte streaming store places them directly in canonical
four-row BCH order. The independent-row route instead needs a local copy
and coordinate rearrangement before BCH. Shared coordinates remove both.
The fixed-width BCH circuit and unrolled inner remain unchanged.

Research modes in `../joint.cpp` are 7/8 for existing independent-row GF16,
9/10 for shared-row GF16, and 11 for shared-row lane permutations without GF
mixing. Modes 7/9 use streaming stores; 8/10 use cached stores. Mode 11 uses
streaming stores. All choices are compile-time specializations.

## Conditional Distribution

Fix four binary outer words c_0,...,c_3. Before shuffling, column j is the
field element v_j=sum_i c_i[j] alpha^i. Let u count the nonzero v_j: it is
the size of the union of the four binary supports.

The shared uniform column permutation places those u packets in a uniform
u-subset of the 256 regions. For fixed v_j != 0, multiplying by uniform
a_j in GF16* makes a_j v_j uniform on GF16*. Independent multipliers make
the values independent conditional on the support. Thus the group has the
same distribution as a uniform u-subset filled with independent uniform
nonzero field elements. Its distribution depends only on u, not the rank
or the original nonzero packet patterns.

Different groups have independent column permutations. The independent
regional permutations give the same conditional packet placement used by
the existing GF16 inner operators. This permits reuse of the *conditional
inner bounds*, not the independent-row outer counts or complete certificate.
It is a statement for each fixed message; the first-moment argument does
not require independent events across messages. GF multiplication preserves
zero packets, so it does not remove correlations in the union support.

## Counting Shared Supports

Let C be the implemented binary BCH code. For h=1,...,4, let B_h(u) bound
the number of ordered four-tuples in C^4 of binary rank h whose union support
has size at most u. The shared-group CDF is bounded by sum_h B_h(u), with
final value 2^512-1. The old shared-route support bounds therefore apply
directly, without a penalty for the nonzero packet-value patterns.

`shared_support.py` regenerates authenticated BCH caps, shortening dimension
bounds, and basis-count bounds. `shared_sparse.py --refined-counts` also
regenerates positive-polynomial shortening, complementary dual shortening,
containment moments, and three iterations of the dual-moment bounds. These
are existing exact refinements, not assumptions about an unknown spectrum.

For one active group, `single_group.support_moments` averages the ordered
inner operators over all supports of size u. We minimize the resulting
Chernoff bounds over a tilt grid and take a nonincreasing majorant before
folding against the upper CDF. Differences of CDF caps are not treated as
actual shell counts. Multiplication by 2048 covers the possible active groups.
Reported probability uppers use Arb outward arithmetic at 256-bit precision.

At rank one, each tuple is c times a nonzero four-bit coefficient vector,
giving exactly 15 A_u tuples of union support u. The dominant one-row term
in the independent-shuffle calculation has four choices of row. Their ratio
15/4 explains about **1.907 bits** of the observed one-group margin loss.
Higher-rank contributions to that bound are tiny. These are contributions
to an upper bound, not observed minimum distances.

## First Proof Screens

The cutoff is 209715, at N=2^21. Excluding every nonzero message with output
weight at most that cutoff would imply distance at least 209716, strictly
above 10%. Each entry below bounds only its specified occupancy class.
An active group contains at least one nonzero outer row.

| Active groups | Two updates: margin | Four updates: margin |
|---:|---:|---:|
| 1 | 40.216 bits | 46.748 bits |
| 2 | 58.557 bits | 73.803 bits |
| 4 | unresolved in initial screen | 142.012 bits |
| 8 | unresolved in initial screen | 52.461 bits |
| 16 | unresolved in initial screen | 40.160 bits (refined counts) |

Each numerical entry covers **every union support** at that occupancy.
Occupancies 3, 5--7, 9--15, and 17 onward do not follow from this table. Differently
loose bounds are used for each entry; the margins need not vary smoothly.

The basic four-update 16-group cover failed. Floating singleton probes
localized its difficulty to medium union supports: the support-160 probe
had log2 upper about +97.88. Replaying the older stronger support counts
improved the full-domain proposal from log2 +268.46 to -25.82, using the
same 64-split budget. That still did not yield a 40-bit outward certificate.
A finer tilt grid and cover are a separate check. These screens identify
slack in our bounds, not bad codewords.

With 512 splits and the finer grid, the full-domain floating proposal reaches
41.6519 bits. The usual search buffer demands 42 bits before attempting a
40-bit outward replay. The optional `--proposal-buffer-bits` changes only
that stopping rule; it does not change the final Arb comparison with
2^(-target_bits). Existing searches keep the two-bit default.

With a 0.125-bit search buffer, fresh 256-bit outward replay succeeds after
nine splits and 112 leaves. It bounds the contribution of **all messages
with exactly 16 active groups**, at cutoff 209715, by
8.139751051264960e-13, giving **40.1600805622134 bits**. This proves the
table's 16-group entry, not the floating 41.6519-bit proposal and not a
full-code distance claim. The run receipt is
`tmp/shared-gf16-q16-outward-40.json`; regenerate it with the command below.

The independent-row dense mixture cannot be imported: its four independently
shuffled row supports have a different distribution. A complete certificate
needs shared-support bounds for every occupancy, including a suitable dense
comparison measure.

A sufficient dense interface is a positive Bernoulli mixture. Write S_u
for the unknown number of group messages with union support exactly u.
Choose c_j>=0 and 0<p_j<=1 such that, for each u>=1,

    S_u <= binomial(256,u) sum_j c_j p_j^u (1-p_j)^(256-u).

Then the shared-shuffle group measure is pointwise dominated by that mixture
of Bernoulli supports, with independent uniform nonzero GF16 labels. For a
particular labeled support, both sides have the same factor 15^(-u);
uniform support placement supplies the binomial denominator on the left.
Include the zero message as a separate point mass. Products preserve the
domination because all coefficients are nonnegative. Thus a checked mixture
could feed the existing positive-mixture inner analysis. Using the group
CDF cap at u as an upper for S_u is sufficient but potentially loose.
The checked mixture below implements this interface; no complete dense
cover is claimed.

## Contiguous Sparse Coverage and the Next Gap

A fresh refined-count run covered every occupancy q=2,...,15 at cutoff
209715. Combining those outward dyadic endpoints with the independently
replayed q=1 result and the prior q=16 result gives

    U_{1:16} < 8.307638e-13 < 2^(-40.13).

Thus every nonzero message supported on at most sixteen groups is covered.
This does not bound the remaining messages. The aggregate margin is
40.1306269568 bits; q=16 dominates, with q=1 and q=15 each contributing
about 8.4e-15. The q=2,...,15 sum uses the freshly regenerated bounds even
where an older individual bound is stronger.

The new receipt is `tmp/shared-gf16-r4-d10-q2-16-fill.json`, SHA256
`af95441802dfb9e90b93e4ebf12ec686e27825833b10d63c8128503c4d43a147`.
Its q=16 entry is unresolved at the requested 44-bit search target; the
aggregate above instead uses the separate 40.16008-bit q=16 receipt.
No failed entry is counted as a proved bound.

A bounded follow-up tested q in {17,18,20,24,32}, using refined outer
counts, a finer tilt grid, and 32 support splits per occupancy. None closed.
Their full-domain floating log2 upper proposals were respectively about
7.92, 70.88, 197.23, 453.17, and 997.35. The largest leaves mix union
supports around 110–140. These are unsuccessful upper bounds, not evidence
of actual low-weight codewords. The receipt and log use the prefix
`tmp/shared-gf16-r4-d10-q17-32-probe`.

### Conditional Returns Improve the Prefix

The next run retained the weighted joint-return bounds through local
packet occupancy three and lazy-state density bounds through occupancy
six. Both are conditional GF16 inner bounds, so they apply to the new
shared-support counts without importing an independent-row certificate.
These are proof refinements only: the encoder and setup distribution do
not change. The existing defaults leave both options off.

Fresh 256-bit outward verification gives these all-support bounds:

| Active groups | Margin, rounded down |
|---:|---:|
| 16 | 113.706350 bits |
| 17 | 44.811070 bits |
| 18 | 59.021393 bits |

Replacing the old q=16 endpoint and adding q=17,18 yields

    U_{1:18} < 4.918876e-14 < 2^(-44.20).

The dyadic sum has an outward-verified margin of 44.2086647664 bits.
It includes the q=1 result and every entry q=2,...,15 from the earlier
fill run. The new receipt is `tmp/shared-gf16-r4-d10-return-q16-32.json`,
SHA256 `1ae290a3061b0809344bdece4d42d3635e8df02c398a5e0fba23d28bd2cbfc00`.
The same run leaves q=20,24,32 unresolved, with floating log2 upper
proposals about +42.53, +268.01, and +748.76. Occupancies above eighteen
remain outside the contiguous proof. Better return bounds recover real
slack, but do not yet close the intermediate-support gap.

A separate 384-bit run reproduced q=17 and improved q=18 to 61.895078
bits using an additional tilt .0105. It also tested q=19 with seven tilts
and 64 splits. The q=19 floating full-domain proposal reached only
12.7567 bits, so no outward 40-bit certificate was returned. This identifies
19 as the first checked gap, not merely an untested occupancy. The receipt
is `tmp/shared-gf16-r4-d10-q17-19-p384.json`. The aggregate above retains
the weaker 256-bit q=18 bound and remains valid without this improvement.

## Checked Shared-Support Mixture

`shared_mixture.py` now constructs the positive comparison measure above.
For each support size u, it assigns the cap B(u) to a Bernoulli component
whose binomial probability at u is largest among the chosen centers. The
coefficient of a component is the maximum ratio
B(u)/(binomial(256,u) p^u (1-p)^(256-u)) over its assigned supports.
It then checks every shell inequality over the rationals. The zero message
has its own point mass; all other components keep their active-group label
even if a comparison draw happens to be zero.

With support-grid step eight, the construction gives 29 positive components
and verifies all 256 nonzero shells. Exhaustive toy-code tests also check
the labeled packet measure directly. The proof is pointwise domination
followed by multiplication of nonnegative measures, not an independence
assumption about the BCH rows. It does not treat differences of upper
CDF caps as shell counts.

This first majorant is too loose for the desired distance proof. Selected
mean-activity diagnostics under the existing basic scalar inner bound are
strongly positive over much of the domain. These points do not constitute
a cover, and a negative endpoint is not a dense certificate. The diagnostic
receipt is `tmp/shared-gf16-r4-mixture-probe.json`. This establishes a checked
interface for future stronger counts, not closure of the dense range.

## Matched Performance and Validation

Peach Ryzen 7950X, GCC 15.2, core 15, K=2^20, 128-bit elements. Each sample
is the median of 101 in-place precomputed transpose calls. Setup, allocation,
and reference checks are excluded. Two seeds and balanced ABCDDCBA order
give four samples per configuration. The table reports their median and range.
Runs are serial and acquire all three existing benchmark locks.

| Route and inner | Median | Sample range |
|---|---:|---:|
| Shared, lane-only, two updates | 5.440 ms | 5.426--5.446 ms |
| Shared, GF16, two updates | 5.863 ms | 5.858--5.884 ms |
| Shared, GF16, four updates | 6.532 ms | 6.527--6.544 ms |
| Independent, GF16, four updates | 7.506 ms | 7.476--7.538 ms |

Shared routing saves 0.975 ms, or 13.0%, for four updates. GF16 costs about
0.423 ms in the two-update comparison. Four updates add about 0.669 ms to
the shared GF16 implementation. The older roughly 5.4 ms result is reproduced,
but is not a timing for the GF16-randomized construction.

Separate phase-instrumented samples locate the four-update saving in BCH
preparation: independent versus shared take 3.668 versus 2.700 ms for
rearrangement plus BCH. Inner plus routing take 3.815 versus 3.828 ms.
These phase means are diagnostic runs, separate from the uninstrumented
headline medians.

Normal checks pass for shared GF16 streaming/cached routes and the lane-only
control, updates 2--4 and seeds 1/17 at K=2^14. They also pass for shared GF16
at K=2^20, the independent-row control, and the two-bit control. Each case
checks three input patterns against a separately materialized encoder, the
dense inner and adjoint, and in-place suffix preservation. The Python suite
passes 218 tests, including exact small-domain shared-distribution and
rank-count checks.

ASan/UBSan checks pass for the same 21 cases. The research encoder is
instrumented; the linked, unchanged prebuilt BCH objects/library are not.
The four-update one-group bound was independently regenerated at 384-bit
precision and again gives 46.7481695461285 bits.

Measured executable SHA256:
`847007ab88b3155cb12ac0903edcc6234dd09d2d835105b1ad3ef1a74978fdd7`.
Measured `joint.cpp` SHA256:
`b6742112e216a461e0bba93f0be6c1c3fa0cd21175a10e1ded93b04f412ef022`.
Raw logs and receipts remain ignored under `tmp/`; no data is committed.

## Reproduction and Next Step

Set OPENBLAS_NUM_THREADS, OMP_NUM_THREADS, and MKL_NUM_THREADS to 1. From
the repository root:

```text
python -B research/workstreams/permutation_locality/gf16_packets/shared_support.py --output tmp/shared-gf16-q1-screen.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --output tmp/shared-gf16-sparse-screen.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 16 --refined-counts --tilts .0032 .005 .0064 .008 .01 .0128 .016 --max-splits 256 --output tmp/shared-gf16-q16-refined-fine.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 16 --refined-counts --tilts .0032 .005 .0064 .008 .01 .0128 .016 --max-splits 256 --target-bits 40 --proposal-buffer-bits .125 --output tmp/shared-gf16-q16-outward-40.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 --refined-counts --tilts .00016 .00032 .00064 .001 .0016 .0032 .005 .0064 .008 .01 .0128 .016 .024 .032 .048 .064 .096 --max-splits 64 --target-bits 44 --proposal-buffer-bits .125 --output tmp/shared-gf16-r4-d10-q2-16-fill.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 16 17 18 20 24 32 --refined-counts --joint-return-through 3 --lazy-density-through 6 --tilts .0064 .008 .01 .012 .0128 .014 .016 .02 .024 .032 --max-splits 64 --target-bits 44 --proposal-buffer-bits .125 --output tmp/shared-gf16-r4-d10-return-q16-32.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_sparse.py --updates 4 --occupancies 17 18 19 --precision 384 --refined-counts --joint-return-through 3 --lazy-density-through 6 --tilts .008 .01 .0105 .011 .0115 .012 .0125 --max-splits 64 --target-bits 44 --proposal-buffer-bits .125 --output tmp/shared-gf16-r4-d10-q17-19-p384.json
python -B research/workstreams/permutation_locality/gf16_packets/shared_mixture.py --minimum-groups 17 --probe .001 .002 .004 .008 .016 .032 .064 .128 .256 .5 1 --output tmp/shared-gf16-r4-mixture-probe.json
python -B -m unittest discover -s research/workstreams/permutation_locality/gf16_packets -p test_*.py
```

`../run_shared_gf16.sh` builds against the existing BCH objects and library,
then provides `check`, `sanitize`, `confirm`, and `profile` phases. Its third
argument specifies the reference-build directory. These measurements used
`/tmp/spin-shared-gf16-qU2wCQ` and reference `/tmp/spin-joint-9n57TT`.

Prefer **four updates** for the next proof effort. The complete prefix
through eighteen groups has over four bits above the target. Next tighten
the intermediate-support argument, starting at the missing occupancy 19,
and strengthen the shared-group counting interface before attempting a
large dense cover. Retain the certified independent-row encoder until the
replacement has a complete certificate. The current GF16 suite passes
223 tests; the archive helper passes four additional integrity tests.
