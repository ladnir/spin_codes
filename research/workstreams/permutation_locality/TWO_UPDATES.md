# Two updates: a measured proof-oriented alternative

This variant retains the one-column grouped permutation, BCH[256,128],
and the fixed IMT(128,19) expansion and feedback maps. Its recurrence is
y=x+Eq, q'=tau_2(tau_1(q))+Bx, with independent sampled transvections.
There is no flush. The transpose applies transposed updates in reverse order.

## Encoding cost

Peach Ryzen 7950X, GCC 15.2, CPU 15, Release/znver4, K=2^20, 128-bit
elements. Setup, allocation, reference checks, and checksums are excluded.
Each process runs three warmups and 101 in-place calls without input reset.
All measurements are serial under the three shared benchmark locks.
Mode order reverses in repetition two. Entries are medians of three process
medians, not isolated fastest calls.

| Route seed | One update | Two updates |
|---|---:|---:|
| 1 | 5.077606 ms | 5.434293 ms |
| 17 | 5.071726 ms | 5.402463 ms |

The added cost is 6.5--7.0%. Both use mask seed 2; two updates consume twice
as many mask pairs, so these are different sampled codes. The cache-line
routing, BCH kernel, and scratch-page policy are unchanged. The earlier
production control was 9.58--9.65 ms; it was not remeasured in this run.

The research-only LocalMap128S19R2 reuses the existing two-update inner
template with the original unrolled 128-input emission circuit. No production
source changed. Every process checks the inner against explicit dense binary
matrix rows, its adjoint identity, and the full encoder against materialized
routing and BCH. It also checks the unchanged suffix. ASan/UBSan checks pass
at K=2^14 for route seeds 1 and 17. Run check_r2.sh with the remote build
root to repeat the checks and serial comparison.

Executable SHA256:
`299ae14941ee204dd7f4d3d13f47c50c7fe4dbed87d7784c532656eb9139b16f`.
Production library SHA256, unchanged from the earlier comparison:
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.
Raw samples remain outside version control.

## State-envelope transformation

On nonzero states, one transvection has kernel K=(I+U)/2, where U refreshes
uniformly over nonzero states. Since U is a projection,
K^r=2^-r I+(1-2^-r)U. Zero remains zero. Output y is evaluated before mixing.

mixing_rounds.py rescales lazy terms by alpha=2^(1-r) and refresh terms by
beta=2-alpha. Zero-source rows are unchanged. If B bounds a sum of lazy
cancellation L and refresh R, and Rbar separately bounds R, use
alpha B+(beta-alpha)Rbar. Uniform-class coefficients divided by class sizes
provide Rbar. Taking their minimum remains an upper bound.

An empty-step uniform-class self coefficient contains the exact lazy term
exp(-lambda v)/2. Its new bound is beta B minus
(beta-alpha)exp(-lambda v)/2. This subtracts a known exact contribution,
not a loose upper bound. All other entries retain their prior interpretation.

The kernel identity passes exact rational checks for state-space sizes
3 and 7 and r=1,...,4. Direct empty/single-window component checks now use
lazy probability 2^-r and refresh probability 1-2^-r. Each tested operator
passes 330 inequalities using representative states from every expansion-
weight class. These are regression checks; the rescaling argument and prior
envelope bounds justify the other states.

## Complete fourteen-group result

The cover driver also accepts `--occupancies` to replay several occupancies
in one process. It builds the two-update operators once, then independently
constructs and checks each support cover. It does not import one-update
certificates or infer smaller-occupancy results from the fourteen-group result.

To justify this reuse, fix the tilt and shape penalties. Let T_j be the
epoch envelope with j active windows. The region operator for q active
groups is

    R_q = [z^q] (sum_{j=0}^{32} binomial(32,j) z^j T_j)^64 / binomial(2048,q).

Computing more coefficients cannot change a lower coefficient. This holds
for the ordered matrix products; no commutativity of the T_j is assumed.
The implementation also performs the same arithmetic and outward rounding
for each retained coefficient, regardless of the maximum retained degree.
A regression check verifies 28 prefix coefficients using exact rational,
noncommuting matrices.

Each occupancy retains its own group-location factor, support-label
multiplicities, coverage check, and outward replay. The batch reports an
outward sum only if every requested occupancy passes. A batch that omits
any occupancy from 1 through 2048 is not a full-code certificate.

The first one-through-fourteen batch completed with the command below.
Occupancies 2 through 14 passed at 192-bit outward precision. Occupancy 1
did not meet this run's 44-bit search gate, so the batch correctly withheld
an aggregate certificate.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --occupancies 1 2 3 4 5 6 7 8 9 10 11 12 13 14 --mixing-rounds 2 --target-bits 42 --max-splits 200 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --pair-tail --retain-parents --joint-witness --penalties .75 1 --tilts .0004 .001 .002 .0032 .005 .008 .01

| Active groups | Outward margin (bits, rounded) |
|---|---:|
| 2 | 57.087246 |
| 3 | 138.456160 |
| 4 | 112.872420 |
| 5 | 217.467230 |
| 6 | 145.885994 |
| 7 | 304.324749 |
| 8 | 261.432019 |
| 9 | 54.422169 |
| 10 | 321.257578 |
| 11 | 317.700475 |
| 12 | 240.327413 |
| 13 | 181.204948 |
| 14 | 146.593142 |

These are independent two-update replays for all supports, ranks, and group
locations at each listed occupancy. They are not inherited one-update bounds.
Their differing precision and witness grids explain why the fourteen-group
number differs from the earlier replay below.

A finer one-group search using tilts .00016, .00025, .00032, .0004, .0005,
and .00064 passed its 384-bit outward replay with 41.8381878934 bits of
margin. It uses eight retained leaves after 25 splits. Its upper bound
is below 2.543610e-13.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 1 --mixing-rounds 2 --precision 384 --target-bits 39 --max-splits 200 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --pair-tail --retain-parents --joint-witness --penalties 1 --tilts .00016 .00025 .00032 .0004 .0005 .00064

For a conservative union over occupancies 1 through 14, round the one-group
upper to 2.543610e-13 and round every other displayed margin down to an
integer. The resulting outward sum has more than 41.83 bits of margin.
This combines independent two-update results and covers all nonzero
messages supported on at most fourteen groups. Occupancies 15 through
2048 remain outside this aggregate. The goal still requires at least
40 bits for the sum over all occupancies, not for this subset alone.

At exactly fourteen active groups, a 190-leaf cover after 25 splits covers
all supports, ranks, and group locations. At both 192- and 384-bit outward precision its
failure contribution for output weight at most 209715 is below
1.438231e-46, or 152.2843977144 bits. Here N=2^21, so this bounds the event
of relative output weight at most 10%.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 14 --mixing-rounds 2 --max-splits 100 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --pair-tail --retain-parents --joint-witness --penalties .5 .75 1 --tilts .0032 .005 .008 .01

This is not a full-code certificate. Earlier one-through-thirteen results
concern one update and are not transferred. The 384-bit repeat agrees with
the 192-bit result to the displayed precision; add --precision 384 to
the command above to reproduce it.

The optional expansion-window refinement in [WINDOW_HISTOGRAM.md](WINDOW_HISTOGRAM.md)
replays this same occupancy at 192-bit precision with 152.3207736818 bits
of margin. It does not extend the occupancy coverage.

The first 64-group diagnostic tests equal supports 38, 80, 128, 160, 224,
and 256, with tilts .008 and .016 and penalties .75 and 1. Its log2 uppers
are respectively -2189.64, +2008.32, +6641.57, +10123.78, +14411.27,
and +14190.79. Thus this witness grid does not close that range. These
binary64 upper-bound proposals do not exhibit bad words and do not rule
out better tilt choices or a different bound. Higher tilts are the next
diagnostic before attributing the gap to the state envelope itself.

That higher-tilt diagnostic also remains open. See DENSE_DIAGNOSIS.md for
the results, transition sensitivity, and the total-outer-weight refinement.

## Remove the local occupancy ceiling

An epoch has at most 32 active windows, but a region can have 2048 active
groups. placement now accepts maximum_groups separately from local operator
degree. It convolves the available local operators while retaining region
coefficients through the requested global degree. Exact enumeration with
noncommuting rational matrices checks all occupancies on small examples,
including more groups than windows per epoch. The earlier two-group
regression also passes. The cover driver now builds at most 32 local
operators and requests its global degree separately.

This removes an implementation limit, not a proof gap. Large support-box
searches may still be impractical. The next priority is to identify which
remaining ranges need a different bound, not to extend one occupancy at a time.
