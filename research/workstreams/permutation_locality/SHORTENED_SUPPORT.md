# A bound on delayed activation for 16-row groups

2026-09-23. The uniform 16-row candidate remains approximately 1.22x faster
than the matched production encoder. This iteration improves the proof tools,
not the implementation timing. Production code and the paper are unchanged.

For K=2^20, the new calculation bounds the probability of the following event
by less than 2^-41: some nonzero message confined to one of the 512 groups
keeps the IMT state zero throughout the first 204 of 256 regions. This covers
all messages within each group, not only one active row or sampled messages.
It is not a 10% distance certificate. A message can activate earlier, cancel
later, or span several groups. Those cases remain open.

## Why the ordinary-spectrum bound needed help

The previous report bounds the number of ordered tuples of BCH words by their
span dimension and union support. Applied directly to authenticated BCH caps,
that bound is too loose for high-dimensional tuples. For example, it permits
about 2^981.63 rank-16 tuples with union support at most 96 coordinates.

There is another constraint: all those words lie in a shortened BCH code on
the same coordinates. Any such code has length 96, even weights, and minimum
distance at least 38. A verified polynomial bound limits its dimension to 20.
Counting tuples in these shortened codes improves the preceding cap to about
2^559.96. Both numbers are upper bounds, not estimates of actual counts.

The BCH shell-cap authentication checks 163 dependency files. It reconstructs
the Christoffel caps, checks the retained joint-shell receipts and their exact
reported rational bounds, and checks equality with the existing cap consumer.
This does not rerun every historical LP or the original BCH containment proof.
Those authenticated results remain premises of this calculation.

## Dimension bounds for shortened codes

Fix an even binary linear code D of length u and minimum distance at least 38.
Its Krawtchouk polynomials are

    K_j(w) = sum_i (-1)^i binomial(w,i) binomial(u-w,j-i).

For each j, summing K_j(wt(c)) over c in D gives |D| times the number of
weight-j words in the dual code. In particular, this sum is nonnegative.
Choose nonnegative rational coefficients alpha_j such that

    F(w) = 1 + sum_{j=1}^u alpha_j K_j(w) <= 0

at every even w between 38 and u. Then

    |D| <= sum_{c in D} F(wt(c)) <= F(0).

The first inequality uses the nonnegative dual counts. The second uses the
nonpositive contributions from all nonzero words. Since |D| is a power of two,
this gives an integer upper bound on the dimension of D.

`shortened_bound.py` uses a numerical LP only to propose coefficients.
It converts the proposal to rationals, rescales it if necessary, and verifies
every required inequality exactly. No unverified numerical objective is used.
All 67 proposals at lengths 38 through 104 pass this check. A failed proposal
would contribute no bound. The script also uses the Hamming and Griesmer bounds.

Shortening on u-v coordinates loses at most u-v dimensions. Thus a bound k_v
at a smaller length also gives k_u <= k_v + u-v. Selected resulting bounds are:

| Support size u | Upper dimension k_u |
|---:|---:|
| 80 | 9 |
| 87 | 13 |
| 91 | 16 |
| 96 | 20 |
| 104 | 26 |
| 128 | 50 |
| 160 | 81 |
| 192 | 109 |

These bounds concern any shortened subcode of the fixed outer. They do not
claim such subcodes exist, or use independence of different outer words.

## Counting tuples through shortened codes

Let R(k,g,h) be the number of ordered g-tuples in a k-dimensional binary space
whose span has dimension h. For 1 <= h <= min(k,g),

    R(k,g,h) = product_{j=0}^{h-1}
              ((2^k-2^j)(2^g-2^j)/(2^h-2^j)).

Set R(k,g,h)=0 when h>k. For each fixed u-subset S of the 256 coordinates,
the BCH words supported inside S form a shortened code of dimension at most
k_u. Every tuple with union support at most u lies inside at least one such S.
Therefore the cumulative tuple count satisfies

    T_{g,h}(u) <= binomial(256,u) R(k_u,g,h).

Take the minimum with the spectrum-only bound from the previous report and
the total count R(128,g,h). Caps at larger supports can also cap earlier ones:
the true cumulative count is nondecreasing. The implementation uses this
monotonicity without interpreting differences of caps as shell counts.

For g=16, the combined support lower bounds for ranks 1 through 16 are

    38, 57, 67, 72, 75, 77, 78, 79,
    80, 82, 84, 85, 86, 88, 90, 91.

## Applying the counts to a zero-state prefix

Fix a nonzero message supported within one 16-row group. Its union support
has some size u. A shared uniform coordinate permutation sends that support
to a uniformly sampled u-subset of the 256 regions.

In each active region, the group occupies one aligned 16-position window
within an IMT epoch. Its window is uniform among the eight possibilities,
and its lane permutation is uniform and independent. These choices are
independent across regions. Starting from zero, the state remains zero
through that epoch exactly when Bx=0 for its input x.

The script reads the production weight5_seed0 feedback columns, checks them
against their deterministic construction, and enumerates all subsets in all
eight windows. There are exactly two nonempty zero-feedback choices, both
of weight eight. Hence, for every fixed nonempty group column, the probability
of zero feedback is at most

    q = 2 / (8 binomial(16,8)) = 1/51480.

No assumption about the transvection distribution is needed for this step:
every linear state update maps zero to zero. This argument concerns a single
active group. Other active groups could cancel its feedback within an epoch.

For a prefix of ell complete regions, the probability that this fixed message
keeps the state zero throughout the prefix is at most

    f_ell(u) = sum_j binomial(ell,j) binomial(256-ell,u-j) q^j
                     / binomial(256,u).

The function f_ell is nonincreasing in u. Couple successive supports by adding
one element from a random ordering; the prefix intersection can only grow.
If C_h(u) is the combined cumulative cap for rank h, summation by parts gives
the following union bound over all 512 possible single active groups:

    512 sum_{h=1}^{16} [
        C_h(256) f_ell(256)
        + sum_{u=0}^{255} C_h(u) (f_ell(u)-f_ell(u+1)) ].

Every multiplier of a cumulative cap is nonnegative. All calculations use
integers and rational numbers; logarithms are used only for display.

| Prefix regions ell | log2 of the untruncated union bound |
|---:|---:|
| 180 | 182.316 |
| 192 | 7.659 |
| 204 | -41.060 |
| 216 | -72.095 |
| 224 | -138.406 |
| 240 | -323.522 |
| 255 | -512.366 |

Positive entries are vacuous probability bounds, not predicted failure rates.
At ell=204, the script checks the strict inequality against 2^-41 exactly.
The rank-one contribution dominates there. Without the shortened-code count,
the corresponding bound has log2 about +239.296 and is vacuous.

## Reproduction and next step

Run `python -B research/workstreams/permutation_locality/bch_joint_support.py`
from the repository root, with python-flint and SciPy installed. The script
authenticates the outer caps, repeats the production feedback census, and
rebuilds the exact counting and prefix bounds. It writes no data files.
`shortened_bound.py` separately checks the polynomial recurrence against its
binomial definition and tests bounds on small even codes. The joint-support
tests still cover 16,640 exhaustive tuples. Prefix probabilities are checked
against exhaustive subsets at lengths 1 through 8, including q=0 and q=1.

Next, include the inner's emitted weights and later returns to zero for one
active group. The shortened-support count supplies multiplicities; it does
not replace that inner analysis. Do not transfer the original independent-row
certificate. On implementation, the best confirmed candidate remains 7.86 ms
versus approximately 9.57 ms; the 2x target remains open.
