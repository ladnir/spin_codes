# Retaining odd columns and complete column histograms

An even-weight column in four grouped rows contributes zero to their XOR.
If every column is even, the four words are linearly dependent. The earlier
support-only bound did not retain this distinction in the inner analysis.
The following count retains it, but the tested scalar weighting does not
improve the final diagnostic.

## Counting with the XOR word fixed

Fix a binary linear outer code C of length n. Let A_v be a coefficientwise
upper bound on its number of weight-v words, with A_0 = 1. For an ordered
four-tuple of codewords, define f as their XOR and V as its weight. Then V
is exactly the number of odd-weight columns. Write h for the tuple's rank
and u for an upper bound on its union support.

Let M_h be the number of spanning four-tuples in an h-dimensional space.
Let Z_h be the number of spanning three-tuples in that space. Define
Z_4 = 0. Exactly Z_h spanning four-tuples have XOR zero: their first three
words determine the fourth. Each fixed nonzero f is the XOR of

    m_h = (M_h - Z_h)/(2^h - 1)

spanning four-tuples. For h = 1,2,3,4 these multiplicities are 8,56,336,1344.
Thus the XOR-zero contribution is the fraction Z_h/M_h of the corresponding
rank-h support count, and its rank-four contribution is zero.

For nonzero f of weight v, count bases of an h-dimensional subcode D that
start with f. Suppose h >= 2 and D has union support u' <= u. There are

    B_h = product_{i=1}^{h-1} (2^h - 2^i)

ordered extensions (f,b_1,...,b_{h-1}). Each extension position is uniform
over D minus {0,f}. The sum of the weights of all words in D is
2^(h-1)u'. Therefore the average total extension weight is at most

    mu = (h-1)(2^(h-1)u - v)/(2^h - 2).

Every extension has weight at least L = (h-1)d, where d is a valid lower
bound on the outer distance. Choose a threshold T among the possible
extension weights. Let T_next be the next possible weight. If T_next > mu,
at least B_h(T_next-mu)/(T_next-L) extensions have weight at most T.
Round this lower bound upward to a multiple of (h-1)!, since permuting an
extension gives that many ordered extensions. Call the result G_T.

The polynomial (sum_{w>0} A_w z^w)^(h-1) upper-bounds all ordered extension
tuples by total weight, even when those tuples are dependent. If P_T is its
cumulative coefficient through T, the desired four-tuple count is at most

    floor(m_h A_v P_T/G_T).

The implementation minimizes over the eligible thresholds. It also
intersects this result with the existing rank-h support count. For h = 1,
the bound is simply 8 A_v when v <= u.

## A complementary shortening bound

Large v admits another count. Let d_t bound the dimension of C shortened
to any t coordinates. For each f of weight v and each t-set containing its
support, at most [d_t-1 choose h-1]_2 subcodes of rank h contain f.
Here the bracket is the binary Gaussian binomial coefficient.

Every D supported on at most u coordinates appears in at least
binomial(n-u,t-u) such t-sets. Hence, for each t >= u, another bound is

    floor(m_h A_v binomial(n-v,t-v)
          [d_t-1 choose h-1]_2 / binomial(n-u,t-u)).

Use zero when d_t < h. The implementation intersects these bounds over t.
This count complements the basis bound: one helps small XOR weights, the
other helps XOR weights close to the union support.

At union support at most 128, the previous unweighted four-tuple cap has
log2 value 350.9023. The new marginal caps, summed over ranks, are:

| Exact XOR weight V | Log2 count upper bound |
|---|---:|
| 0 | 278.6920 |
| 38 | 315.4878 |
| 48 | 326.3269 |
| 64 | 341.2959 |
| 80 | 350.9023 |
| 96 | 350.9023 |
| 128 | 274.6728 |

These are bounds on distinct subclasses, not replacements for the total
count. Summing all V and intersecting with the old total does not improve
the old total at this support.

## Tested scalar coupling

The research operators now optionally carry eta^V, in addition to the
existing total-weight and all-one-column factors. Each input shape receives
its factor before taking local maxima. The outer sum carries eta^-V and
rho^-J, where J counts all-one columns. A rearrangement bound couples their
two marginal caps. When eta < 1, reflect V to u-V before packing its CDF.
This relaxation does not enforce every joint relation between V and J.

For 64 groups of support 128, with two updates, lambda = 0.032, rho = 0.75,
and exact feedback distributions through six windows, the baseline log2
diagnostic is +3380.7392. The first eta = 1.1 and 1.25 tests, before adding
the shortening refinement, give +3880.7945 and +4563.6071. With both count
refinements, eta = 0.9 and 0.75 give +4353.2700 and +6174.1662.

All tested nontrivial scalar choices worsen this point. Their local
operators are outward bounds, but their matrix-power scores are binary64
diagnostics. They are not certificates and do not exhibit bad codewords.
This experiment does not justify a large scalar-parameter sweep.

The default eta = 1 operators reproduce the earlier 64-group cache within
binary64 rounding. The default 14-group certificate was also replayed at
192-bit precision: all supports pass with upper below 1.438231e-46, or
152.2843977144 bits. This replay covers exactly fourteen active groups,
not every occupancy through fourteen.

## Next route: condition on the full column histogram

The scalar tests still take local worst cases over column-weight shapes.
The actual shared column shuffle supplies more structure. For group i,
let n_i,w count its columns of weight w, for w = 0,1,2,3,4. These counts sum
to 256. The shared uniform column permutation places this multiset uniformly
among the 256 regions. Independent lane shuffles then make each nonzero
four-bit pattern uniform conditional on its weight.

Introduce auxiliary independent categorical draws X_i,r, one for each
group i and region r, with probabilities p_i,w. Choose p_i,w positive
where n_i,w is positive. Let E be the event that every group has its fixed
histogram. Its exact probability is

    Pr[E] = product_i (256! / product_w n_i,w!) product_w p_i,w^n_i,w.

Conditioned on E, the column weights have exactly the required shuffled
distribution. Keep region routing, lane shuffles, and transvections
independent with their original distributions. For the nonnegative output
moment F = exp(-lambda * output_weight),

    E[F | E] <= E[F] / Pr[E].

Before conditioning, different regions are independent. This permits a
region transfer built by averaging categorical input shapes, rather than
taking a separate maximum for every matrix entry. For groups sharing the
same p_i,w, the number of active groups in a region is binomial; their
nonzero weights have probabilities p_i,w/(1-p_i,0). The existing
without-replacement placement recurrence can then average local shapes
with these probabilities. This is a change to the proof, not the encoder.

`categorical_conditioning.py` verifies 256 exact profile/matrix cases on
small examples with noncommuting, nonnegative rational transfer matrices.
It checks the multinomial normalizations, the conditional upper bound,
and reconstruction of the unconditional matrix power. The subsequent
SPIN pilot and a restricted dense certificate are in CATEGORICAL_PILOT.md.

The proposed experiment held a full histogram fixed and tested this
averaged transfer at the difficult 64-group point. Full closure would
additionally need outer counts and coverage for every group histogram,
heterogeneous groups, all occupancies, and an outward total failure sum.

## Reproduction

    python -B research/workstreams/permutation_locality/odd_column_caps.py
    python -B research/workstreams/permutation_locality/odd_column_screen.py --etas 1 1.1 1.25
    python -B research/workstreams/permutation_locality/odd_column_screen.py --etas .9 .75
    python -B research/workstreams/permutation_locality/categorical_conditioning.py

The outer component checks 1,536 exact/inflated-spectrum marginal
inequalities, 29 basis-mean identities, and 110 exact weighted inequalities.
The inner screen checks detailed shapes and 330 direct empty/single-window
component inequalities per tested eta. Its eta = 1 cache comparison uses
the cache-generation command in WEIGHT_AND_REFRESH.md. The current code
includes shortening bounds for both sides of eta = 1; the first two
nontrivial values above record the earlier, weaker count variant.
