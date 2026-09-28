# Diagnose the thirteen-group gap before extending the cover

The twelve-group certificate covers every support and rank. The next full
cover fails, but that alone does not distinguish interval-counting slack
from slack in the state envelope. We therefore tested exact support vectors
using the same weighted outer CDF and the same inner operators.

## Exact support probes

The following binary64 scores include all thirteen-group locations. Each
row fixes all thirteen union supports to u. These are upper-bound proposals,
not outward certificates or evidence that the actual code has bad words.

| u | Log2 contribution upper |
|---|---:|
| 96 | -147.02 |
| 112 | -151.62 |
| 128 | -82.26 |
| 136 | -48.26 |
| 144 | -34.49 |
| 152 | -23.40 |
| 160 | -20.82 |
| 176 | -64.76 |

The mixed vector with ten supports 134, two supports 154, and one support
168 gives -50.05, including its label multiplicity. Finer partitioning alone
cannot make these *tested witnesses* close the difficult equal-support
points. This does not rule out better witnesses for the same envelope.

The exact-support all-one-count screen, which additionally uses the checked
dual-distance premise for OA caps, gives -52.60, -36.07, -25.46, -38.31,
and -256.98 at supports 136, 144, 152, 160, and 176. Its improvement at 176
does not resolve the intermediate range. Neither screen covers all vectors.

Reproduce from the repository root:

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 13 --fresh --window-average --multi-average --joint-witness --penalties .5 .625 .75 1 --tilts .005 .0064 .008 .01 --probe-supports 96 112 128 136 144 152 160 176 --probe-vector 134 134 134 134 134 134 134 134 134 134 154 154 168
    python -B research/workstreams/permutation_locality/occupancy_allones.py --groups 13 --fresh --window-average --multi-average --penalties .5 .625 .75 --tilts .008 .01 --supports 136 144 152 160 176

## Avoid counting the same even subspace at each support

Fix a positive all-one column count J and rank h>=2. As in
ALL_ONE_COLUMNS.md, the even combinations form an (h-1)-dimensional subspace
H. If H has support v, the four-row union has support v+J.

Let A_h(v) be the authenticated upper CDF for rank-h ordered four-tuples.
Let R_h be the number of ordered four-tuples spanning a fixed h-dimensional
space. The number of (h-1)-spaces H with support at most v is bounded by

    L(v) = floor(A_(h-1)(v) / R_(h-1)).

For fixed H, the previous extension argument bounds the number of possible
extension cosets by

    e_J(v) = min(floor(binomial(256-v,J) * 2^D(v) / 2^(h-1)),
                 S(floor(v/2)+J)),

where D(v) bounds the shortened dimension and S is the ordinary nonzero
BCH weight CDF. Set e_J(v)=0 where D(v)<h-1. Each extension gives at most

    t_h = 2^(h-1) product_{i=0}^{h-2}(8-2^i)

spanning four-tuples with the prescribed even subspace. To count tuples
whose union support is at most u, apply the existing CDF folding inequality:

    count_(h,J)(support<=u)
      <= t_h sum_{v=0}^{u-J} (L(v)-L(v-1))
                 max_{v<=w<=u-J} e_J(w),       L(-1)=0.

Intersect this with A_h(u). Sum over h=2,3,4, then add the rank-one count
A_J from the ordinary BCH spectrum. Rank-one tuples with J>0 consist of
four copies of the same nonzero word, so their union support is exactly J.

This bounds each J-prefix directly, instead of adding independent bounds
for every possible exact support. `prefix_flag_caps` intersects these
prefix bounds with the earlier ones before allocating the total count
budget. The optional `--prefix-flags` switch enables the improvement.
Existing certificate commands retain their original default behavior.

Exact enumeration on the three small test codes verifies the improved
weighted CDFs at penalties 1/2, 3/4, and 1. At penalty .625, the saved bits
per group at supports 128, 136, 144, 152, 160, 176, and 192 are respectively
.197, .334, .122, .159, .222, .319, and .494. This is a valid count refinement,
but its size alone does not close the difficult thirteen-group range.

## Identify the useful state-bound target

`occupancy_transition_probe.py` halves selected transition coefficients,
one family at a time, while preserving the actual baseline in a separate
run. These altered coefficients are deliberately **not** upper bounds.
Their scores measure sensitivity only and must never enter a certificate.
The default probe uses thirteen groups, support 152, lambda=.008, and
all-one penalty .625, including the improved prefix counts above.

    python -B research/workstreams/permutation_locality/occupancy_transition_probe.py

At the default point, the baseline log2 score is -25.46. Artificially halving
the mature-density cancellation entry gives -82.31; halving fresh-state
cancellation gives -70.59. Halving mature lazy mass gives -139.98. Density
persistence gives -51.43, uniform-class lazy density gives -58.54, and
fresh-to-density gives -35.47. Halving direct zero-state return barely
changes the score (-25.52). These identify possible proof targets, not
achievable improvements.

## Exact fresh-state cancellation at two and three inputs

A fresh state is dominated by a mixture of distributions p_a: the feedback
of one independent uniform window input of weight a. Fix current column
weights b_1,...,b_j, where j is two or three. Their windows are uniformly
distinct, independently of the fresh input. Let Z_j(b) count ordered
distinct-window assignments and masks of weights b with zero total feedback.
The existing kernel census computes these integers through j=4.

The old window is either distinct from all current windows or equals one
of them. In the first case, cancellation has numerator Z_(j+1)(a,b).
In the second case, suppose it equals current window i. Combine the old and
current masks into a mask D of weight d. For fixed D, the number of pairs
of weights a and b_i with XOR D is

    c(a,b_i,d) = binomial(d,(a-b_i+d)/2)
                * binomial(4-d,(a+b_i-d)/2).

An inadmissible parity or binomial argument gives zero. The case d=0 cannot
cancel: the remaining one or two nonzero windows have injective feedback.
Thus the exact cancellation probability is

    [Z_(j+1)(a,b) + sum_i sum_{d=1}^4 c(a,b_i,d) Z_j(b with b_i=d)]
      / [32 binomial(4,a) (32)_j product_i binomial(4,b_i)].

Here (32)_j is the falling factorial. The denominator allows the old window
to coincide with a current window; imposing distinctness on all j+1 windows
in that denominator would describe a different distribution.

`fresh_collision.py` implements the identity and checks every weight shape
on a small map with injective pairs and nontrivial triple dependencies.
For the actual map, pair injectivity and the existing kernel census provide
the premises. The lazy output factor remains the conservative
exp(-lambda max(0,48-sum_i b_i)); refresh contributes at most that factor
divided by 2(2^19-1). Apply the all-one penalty before maximizing over b and
the fresh mixture component a. Only the fresh-to-zero entry changes.
The optional `--fresh-collision` switch enables this bound.

    python -B research/workstreams/permutation_locality/fresh_collision.py

The initial combined count/collision screen improves log2 scores at
supports 128, 136, 144, 152, 160, and 176 to -91.68, -61.00, -43.78,
-32.44, -31.29, and -78.12. A small additional witness grid improves 144
to -49.35 and 160 to -34.45, but does not close 152. Those are still
selected-point binary64 results, not a thirteen-group certificate.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 13 --fresh --window-average --multi-average --prefix-flags --fresh-collision --joint-witness --penalties .5 .625 .75 1 --tilts .0064 .008 .01 --probe-supports 128 136 144 152 160 176

## Output weights on mature-state cancellation

For mature pointwise density C, a lazy return to zero has incoming state
q=Bx. Its weighted mass is at most

    (C/2) E_x[exp(-lambda wt(EBx+x))].

The average uses the actual distribution of the current input x. Including
Bx=0 terms is conservative because the mature component excludes zero.
`zero_moment.py` enumerates the output-weight histogram for every two- and
three-window weight shape. It enumerates unordered physical window subsets
and all 15 nonzero masks per window, then restores the labeled multiplicity.
For a sorted weight shape, that factor is the product of the factorials of
its repeated-weight counts. Each histogram's total equals
(32)_j product_i binomial(4,b_i).

The program also directly recomputes the first and last physical subsets
using Python integers, evaluating EBx+x from the combined feedback. The
main census uses vectorized XOR and exact integer counts. Only exponential
evaluation uses outward arithmetic. The optional `--zero-moment` refinement
changes the density-to-zero entry at two and three inputs; it does not
assume that the incoming state is uniform or independent of x after fixing
the cancellation event.

    python -B research/workstreams/permutation_locality/zero_moment.py

With all refinements enabled, the selected thirteen-group log2 scores at
supports 144, 152, 160, and 168 become -53.56, -39.44, -41.00, and -65.50.
The witnesses use tilts .0075, .008, and .0085 with all-one penalties .625,
.6875, and .75. These scores do not suffice for a full support union and
have not been replayed as outward certificates. The next substantial proof
target is the mature-mass bound: it forgets expansion-weight information
that the sensitivity probe shows could be valuable. Further small grid
extensions should not substitute for addressing that loss.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 13 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --joint-witness --penalties .625 .6875 .75 --tilts .0075 .008 .0085 --probe-supports 144 152 160 168

The subsequent [mature-tail model](MATURE_TAIL.md) closes thirteen active
groups at 75.22 bits in 192- and 384-bit outward replays. That extension
uses the refinements developed here. Larger occupancies remain open.
Production, the paper, and the measured 1.88--1.90x speedup are unchanged.
