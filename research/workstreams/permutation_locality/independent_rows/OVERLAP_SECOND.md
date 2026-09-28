# Constraining cancellation with a second overlap moment

The mean-overlap bound permits a mixture of no cancellation and complete
cancellation. A second moment restricts that extreme mixture. This note
uses the same fixed maps and input distribution as [OVERLAP_MEAN.md](OVERLAP_MEAN.md).
It does not change the encoder or its setup distribution.

Fix a packet-weight shape. Let X be uniform over its D inputs, all of
weight W. For a target state t, write V_t=weight(E(t+BX)) and let H_t
count the common one coordinates of X and E(t+BX). Define d_t=2H_t-W.
Then the output weight is V_t-d_t and |d_t|<=W.

## Exact character calculation

Let a_i give coordinate i of the linear expansion E. Define

    R_il = sum_(x of this shape) x_i x_l (-1)^((a_i+a_l) dot Bx).

Squaring W-2H_t and using x_i^2=x_i gives

    E[d_t^2] = W + (2/D) sum_(i<l) (-1)^((a_i+a_l) dot t) R_il.

Consequently the signed sum gives the exact second moment at t=0.
For every target, independently of t,

    E[d_t^2] <= min(W^2, W + (2/D) sum_(i<l) |R_il|).

This is a second moment around W/2 for H_t, not its variance around
its own mean. No independence between input and feedback is assumed.

Each R_il is a coefficient of the feedback character polynomial with
coordinates i and l forced to one. In one window, use the signed
two-coordinate factor. In distinct windows, multiply two signed
one-coordinate factors. Remove the corresponding unrestricted window
factors by exact truncated polynomial division. All coefficients and
signed sums use unbounded integers.

The complete calculation covers all 8128 coordinate pairs and all 1000
nonempty shapes through ten occupied windows. Ten additional checks use
an independent univariate formula for all-four shapes. For the shape
with ten weight-four packets, W=40 and

    E[d_0^2] = 161272718/4032015,
    E[d_t^2] <= 216797131/5376020 < 40.326698 for every t.

The pointwise bound would be 1600. This improvement concerns the local
input distribution; it does not itself prove a whole-code distance.

## A quadratic upper bound for the output tilt

Fix lambda>0 and a nonnegative rational c. For each nonzero expansion
weight v, choose a nonnegative number a_v satisfying

    a_v >= exp(-lambda*(v-d)) - c*d^2

for every d in {-W,-W+2,...,W}. These are finitely many constraints.
If p_v=Pr[V_t=v] and sigma bounds E[d_t^2 1_(V_t>0)], then

    E[exp(-lambda*(V_t-d_t)) 1_(V_t>0)]
      <= sum_(v>0) p_v*a_v + c*sigma.

For t=0, the exact feedback census supplies p_v and p_0. Since V_0=0
implies H_0=0, the exact restricted second moment is

    sigma = E[d_0^2] - W^2*p_0.

The implementation uses binary64 only to propose c. It rounds c to
an exact nonnegative rational and reconstructs every a_v with outward
arithmetic, checking every discrete overlap value through the formula.
Thus a failed optimizer can weaken the bound but cannot justify an
invalid inequality. The lazy probability alpha=1/4 is applied once.
The all-four penalty is applied before the complete shape maximum.

## Translated-density extension

Let n(s) count feedback values and use the upward counts n'(s) from
[DENSITY_EXTENSION.md](DENSITY_EXTENSION.md), so n(s)<=u*n'(s).
Let G'_v(t) be their translated expansion-class counts. Since a_v>=0,

    E[exp(-lambda*(V_t-d_t)) 1_(V_t>0)]
      <= (u/D) sum_v G'_v(t)*a_v + c*sigma_uniform,

where sigma_uniform bounds the original E[d_t^2] for every target.
The first term uses rounded counts; the second still concerns the
original input law. Taking the maximum over nonzero t and multiplying
by alpha bounds C-to-C. For U_v-to-C, retain only that class, allow it
the entire second-moment budget, and divide by its cardinality.

Here W<=40 and every nonzero v>=48, so intercepts may be bounded by
one. Round them upward to A_v/2^b, where 0<=A_v<=2^b. The positive
integer sum is bounded by D'*2^b; this determines the int64 guard.
Add c*sigma_uniform only after the exact maximum, using outward arithmetic.
Several fixed slopes give alternative valid bounds; take their minimum.

## Replays and remaining scope

At q=80 with support 200 in every active group, tilt .072 and penalty
.9, the following 192-bit outward log2 uppers are rounded upward:

| Local refinement with optimized coupled columns | Log2 upper |
|---|---:|
| Exact ten-window feedback, four-window joint counts, extended density | +158.187720 |
| Mean overlap for zero returns | +133.142758 |
| Mean overlap for zero returns and density | +121.548507 |
| Quadratic overlap for zero returns | +98.248849 |
| Quadratic zero returns plus mean-overlap density | +85.881088 |
| Quadratic zero returns and translated quadratic density | +63.907200 |

All are positive and therefore vacuous as probability bounds. These
are selected support events, not complete occupancy covers. The translated
quadratic-density replay has completed: its 931 density shapes and exact
1000-shape second-moment census give the last row. The positive-count
rounding factor is 15335685049/15335681400. No new occupancy is certified.

```sh
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows/candidates/overlap_second -p 'test_*.py'
python -B research/workstreams/permutation_locality/independent_rows/candidates/overlap_second/replay_second.py --epoch-cache tmp/independent-row-local-families --groups 80 --support 200 --tilt .072 --density-chord --quadratic-density
```

Omit `--quadratic-density` to reproduce the preceding row. Eight
tests in this directory pass. They include exact second-moment checks,
every discrete quadratic constraint, and direct small-state comparisons
for every shape, target and expansion class, with coarse count rounding.
No production encoder, shared-shuffle proof, or certified occupancy
aggregate is changed by these selected-point experiments.
