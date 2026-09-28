# Two active groups: collisions counted, low ranks closed

2026-09-24. The two-column candidate now has an outward bound for every
message with exactly two active groups whose row spans both have rank at
most two. The combined failure contribution is below 9.101189e-20,
or 63.2525 bits, including every pair of group locations.

This extends the single-group coverage in [TWO_COLUMNS.md](TWO_COLUMNS.md).
It does not close higher-rank two-group messages or three or more active
groups. No encoder or production defaults changed. The last matched
performance result remains 5.07--5.10 ms, about 1.89x, not 2x.

## The new obstacle: two windows can share an epoch

Fix two nonzero four-row tuples before setup. Their groups occupy distinct
eight-bit windows in each macroregion. There are 128 epochs, each containing
16 windows. Conditional on both groups contributing nonzero input in that
macroregion, their windows share an epoch with probability 15/2047.

Every single eight-bit window has an injective feedback map. The union of
two windows need not. Exact enumeration of the 120 unordered window pairs
finds 101 maps of rank 16 and 19 maps of rank 15. Each deficient pair has
one nonzero kernel vector; neither half of that vector is zero.

The census retains the two column weights of each group. These sorted pairs
have fourteen nonempty possibilities. For all 196 ordered pairs of shapes,
it enumerates distinct ordered windows and all compatible lane masks.
It records zero feedback, the largest nonzero syndrome atom, and weighted
cancellation counts indexed by the expansion-map weight.

The largest zero-feedback probability is 1/30720, attained when both shapes
are (3,3), conditioned on a common epoch. Including placement gives an upper
bound of 1/(2047*2048) for this particular local cancellation event.
Neither number is a full-code failure probability.

The total number of zero-feedback ordered mask/window choices is 38.
This agrees independently with the rank census: twice the 19 nonzero kernel
vectors. Direct scalar evaluation of the full 128-bit input and feedback map
also reproduces three selected shape-pair censuses.

## Averaging placements without assuming independent epochs

Use the seven-coordinate positive state envelope from the single-group
verifier. Let Z be the empty-epoch operator. Let A and B be the operators
for the two groups, averaged over an individual window. Let C be their joint
operator, averaged over distinct windows in the same epoch.

For E epochs and W windows per epoch, define

    D = sum over 0 <= i < j < E of
        (Z^i A Z^(j-i-1) B Z^(E-1-j)
         + Z^i B Z^(j-i-1) A Z^(E-1-j)),

    S = sum over 0 <= i < E of Z^i C Z^(E-1-i).

The averaged operator is

    R = (W D + (W-1) S) / (E(EW-1)).

Here E=128 and W=16. Distinct epochs contribute W^2 ordered slot choices;
a common epoch contributes W(W-1). The denominator counts all distinct
ordered slots. An independent enumeration with noncommuting small matrices
checks the recurrence for E=1,2,3 and W=2,3.

For the same-epoch operator, output moments use the conservative lower
weight v-a-b when the incoming expansion has weight v and the input weights
are a,b. Returns to zero use the exact cancellation census. Thus rare
feedback cancellation is included rather than discarded.

The outward verifier takes entrywise shape maxima before averaging epochs.
This permits extra choices and gives a valid upper envelope. A separate
diagnostic delays these maxima until after averaging each fixed pair of
shapes; that refinement gives only a small improvement on the current grid.

## Support coefficients and outward replay

For each group, a uniform coordinate shuffle sends a support of size u to a
uniform u-subset of 256 coordinates. The two support subsets are independent.
Let R_ab be the macroregion envelope when the groups occupy a,b columns,
respectively, where a,b are in {0,1,2}.

Introduce auxiliary probabilities p,q and form

    M(p,q) = sum over a,b of
        binom(2,a) p^a (1-p)^(2-a)
        binom(2,b) q^b (1-q)^(2-b) R_ab.

These Bernoulli weights are a device for bounding nonnegative coefficients;
they are not a replacement for the construction's fixed-size supports.
For supports u,v, their conditional exponential moment is at most

    (e_zero M(p,q)^128 1) /
    (Pr[Bin(256,p)=u] Pr[Bin(256,q)=v]).

Multiply by exp(tilt*209715) to obtain a Chernoff bound, then cap at one.
Support buckets have width four. A binomial mass is log-concave in its index,
so its minimum within a bucket occurs at an endpoint. Using those minima
gives a bound for every support in each rectangle, not just grid samples.

The optimizer proposes a tilt from seven fixed choices and p,q rounded to
rationals with denominator 10^9. The verifier recomputes each selected bound
with Arb; optimizer convergence is not part of the proof. Reused positive
upper bounds are rounded upward, and probability denominators use lower
enclosures. Exact small-support enumeration checks the coefficient bound
and the rectangle extension.

For rank one, shell counts are at most 15 times the authenticated BCH shell
caps. For rank two, each bucket uses the authenticated cumulative support
cap at its upper endpoint. Differences of CDF caps are not used as shell
counts. Finally, the verifier includes all binom(2048,2) group pairs and both
rank assignments when the ranks differ.

| Group ranks | Failure upper, including all locations | Margin (bits) |
|---|---:|---:|
| 1,1 | <2.040957e-20 | 65.4093 |
| 1,2 or 2,1 | <3.461718e-20 | 64.6470 |
| 2,2 | <3.598514e-20 | 64.5911 |
| Combined | <9.101189e-20 | 63.2525 |

All three rows passed independent 192-bit and 384-bit outward replays.
The reported bounds agree at both precisions. Setup uses the ideal uniform
coordinate permutations, placements, lane shuffles, and transvections;
these results do not cover the separate heuristic permutation family.

## What remains and where the current bound loses information

The optimized binary64 grid remains vacuous for high-rank pairs. In
particular, the rank-four/rank-four diagnostic has a log2 union bound near
+203, not a small failure bound. Moving shape maxima after region averaging
changes that figure by only about 0.3 bits. Placement averaging alone is not
the main source of this gap.

The current envelope forgets row-weight composition while retaining only
the number of occupied columns. It can repeatedly select unfavorable shapes
without charging their outer multiplicities. The next useful refinement is
to retain total row weight, reuse the joint support/weight counts, and check
whether more feedback-state information is also needed. This is a proposed
way to tighten the proof, not evidence that the missing cases already close.

The eventual certificate must cover all group occupancies. The performance
goal still requires roughly another 5.3% reduction in encoder latency.

## Reproduction

Run from the repository root with NumPy, SciPy, and python-flint installed:

    python -B research/workstreams/permutation_locality/test_two_group.py
    python -B research/workstreams/permutation_locality/two_group_verify.py --precision 192 --ranks 1 1
    python -B research/workstreams/permutation_locality/two_group_verify.py --precision 192 --ranks 1 2
    python -B research/workstreams/permutation_locality/two_group_verify.py --precision 192 --ranks 2 2

Repeat the verifier commands with `--precision 384` for the higher-precision
replay. Each run reconstructs the local census, authenticates 163 BCH
dependencies, and checks 67 exact shortening witnesses. This does not rerun
every historical BCH proof.

The exploratory command is:

    python -B research/workstreams/permutation_locality/two_group_screen.py --step 16 --optimize

Its support grid is not a certificate. The separate verifier covers every
support through the rectangle bounds described above.
