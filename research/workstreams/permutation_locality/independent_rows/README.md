# Proof exploration: independently shuffled BCH rows

A [first full four-bit certificate](dense_closure/FIRST_CLOSURE.md) now
closes at relative distance greater than 0.5%, with more than 43.744 bits
of setup-failure margin, for the unchanged two-update construction. This
lower-distance result was replayed at 256-bit precision. The original
10% target discussed below remains open beyond the sparse range.

This is a separate ensemble from the shared-shuffle route in the parent
directory. Both remain active; see [PROOF_TRACKS.md](../PROOF_TRACKS.md).
The implementation measures 6.50 ms versus 5.42 ms for the shared route.
This note establishes an outer-counting reduction and reports an initial
inner diagnostic. The later [complete sparse batch](LOW_OCCUPANCIES.md)
verifies every support at each occupancy from 1 through 58. Its combined
margin exceeds 43.74 bits for those messages; unrestricted occupancies
59 through 2048 remain open. [BRIDGE.md](BRIDGE.md) states the reduction
from conditional inner bounds to this averaged outer measure.

[REFINEMENTS.md](REFINEMENTS.md) records a stronger selected-class bound
at occupancy 64. [DENSE_DIAGNOSIS.md](DENSE_DIAGNOSIS.md) identifies the
remaining middle-occupancy gap and the new density-column calculations.

[SHAPE_POTENTIAL.md](SHAPE_POTENTIAL.md) explores a continuation bound
inspired by the companion lifting paper. It preserves whole packet shapes
within a step and derives an exact two-epoch single-packet return moment.
These research bounds do not extend the completed occupancy cover or
change the encoder.

## Why this route may be easier to prove

Fix four BCH words, before setup is sampled. Independently permute the 256
coordinates of each word, then form one four-bit packet in each region.
Each packet contains one coordinate from each row. Subsequent group and
lane shuffles are unchanged.

For fixed row weights w_1,...,w_4, the four active-region sets are independent
uniform subsets of those sizes. Their distribution no longer depends on
the words' original coordinate overlaps or linear relations. This replaces
the old four-word pattern-counting problem by an ordinary weight-spectrum
problem. It does not make positions within a row independent.

For example, four equal weight-38 words share exactly 38 active regions in
the old route. With independent row shuffles, their expected union has
121.381 regions. The probability that their union still has only 38 regions
is binomial(256,38)^(-3), approximately 2^(-453.83). This is a fixed-message
support event, not a distance bound or a union over all messages. Four
all-one words still occupy every packet with weight four.

## Exact expected outer counts

Let B be the row length and A_w the number of outer words of weight w,
including A_0 = 1. For a fixed subset S of v regions, define

    H_A(v) = sum_{w=0}^v A_w binomial(v,w) / binomial(B,w).

This is the expected number of row words whose shuffled support lies in S.
Independence of the four row permutations gives H_A(v)^4 for four-tuples.
Let C_u be the expected number of tuples with union size exactly u. Subset
inversion gives

    C_u = binomial(B,u) sum_{v=0}^u (-1)^(u-v) binomial(u,v) H_A(v)^4.

The expectation is over setup, and the sum is over ordered tuples of outer
words. Linearity of expectation does not require different messages to
have independent setup. All rows use the same fixed BCH code, but their
coordinate permutations are independently sampled.

`support.py` evaluates this identity using a common integer denominator and
exact integer finite differences. No floating subtraction enters the count.
For componentwise shell uppers Abar_w >= A_w, applying the *entire* formula
to the synthetic nonnegative spectrum Abar gives componentwise upper bounds
on C_u. The reason is the original positive sum over independent subsets,
whose weights are products of A_w/binomial(B,w). It is not valid to insert
unrelated containment upper bounds into an alternating sum.

### Using the known total number of words

We also know sum A_w = 2^128. Form a synthetic spectrum Astar by assigning
this total to the lowest weights first, never exceeding the shell caps.
Its weight CDF dominates the true weight CDF. Coupling uniform subsets by
prefixes of a common random ordering shows that their union size is
monotone in each row weight. Consequently, the union CDF computed from
Astar dominates the true union CDF. The all-zero tuple is subtracted once.

This yields an exact rational bound on the expected number of nonzero
tuples with union size at most u. Its successive differences are **not**
shellwise upper bounds. A verifier must use the CDF by summation by parts
or another justified monotone envelope. The final CDF equals 2^512-1.

The checked BCH input has minimum nonzero weight 38. Its shell bounds are
authenticated through the existing 163-file dependency check. Selected
base-2 logarithms of the new CDF upper are:

| Union size at most u | Unweighted count upper, log2 |
|---|---:|
| 80 | 115.864 |
| 96 | 154.330 |
| 128 | 244.431 |
| 144 | 299.618 |
| 160 | 353.326 |
| 176 | 401.005 |
| 192 | 442.005 |
| 256 | 512.000 |

These are averaged tuple counts, not probabilities or security margins.

## Retaining all-one packet penalties

Let J be the number of regions where all four rows are active. The existing
inner envelope can attach a factor rho to every such packet, with
0 < rho <= 1. The outer count must then carry the reciprocal factor
gamma^J, where gamma = 1/rho.

For a fixed v-subset S and a fixed j-subset R of S, the single-row mass of
supports containing R and contained in S is

    H_A(v,j) = sum_{w=j}^v A_w binomial(v-j,w-j) / binomial(B,w).

Expand gamma^J as a sum over subsets of the all-row intersection. The
weighted containment mass is exactly

    F_gamma(v) = sum_{j=0}^v binomial(v,j) (gamma-1)^j H_A(v,j)^4.

Subset inversion of F_gamma gives the weighted union shells. Again,
componentwise spectrum caps are valid through the underlying positive
subset sum. For an additional CDF bound, use J <= u whenever the union has
size at most u, hence multiply the unweighted CDF upper by gamma^u. The
implementation takes the smaller of these two valid weighted CDF bounds.
All calculations are exact rationals. The lowest-weight spectrum alone is
not used to bound the intersection weight: its monotone-union argument
does not justify that stronger claim.

At rho = .75, the log2 weighted CDF bounds at u = 128,144,160 are
246.897, 304.119, and 358.395, respectively. Thus the treatment retains
all-one packets rather than assuming they disappear under shuffling.

## Which inner analysis can be reused

Condition on a group's histogram of packet weights 0,1,2,3,4. Its sequence
of packet weights is uniform over arrangements with that histogram:
independent uniform row subsets are invariant under a common permutation
of region labels. Fresh lane shuffles then make each packet a uniform
four-bit mask of its weight, independently of other lane shuffles.

This is the same conditional input distribution obtained from any fixed
four-row array with that histogram, followed by a shared column shuffle
and independent lane shuffles. Thus inner bounds valid uniformly over such
histograms can be reused. The new ensemble changes the histogram's outer
averaging measure, not the conditional inner experiment. Different groups
remain independent before the unchanged region routing is sampled.

This reduction does not transfer rank-dependent BCH counts, certified
support covers, or previously reported failure margins. Those belong to
the old outer measure. A new full certificate needs its own weighted outer
counts, support cover, and outward sum over every occupancy.

## Initial inner diagnostic

`pilot.py` reads the existing shape-resolved nine-coordinate operator cache
only after checking its source/dependency hash. It takes entrywise maxima
over input shapes after applying rho to each all-one packet. It uses
without-replacement group placement within a region and a positive
coefficient bound to condition each group's union size. The outer factor
uses the new exact weighted CDF, never a shared-route joint-pattern count.

For exactly 64 active groups, tilt .032, rho = .75, bad output weight at
most 209715, and the indicated union size in every active group:

| Union size per group | Log2 contribution upper, binary64 proposal |
|---|---:|
| 80 | -4094.79 |
| 96 | -4214.05 |
| 128 | -2735.70 |
| 144 | -977.86 |
| 160 | +757.91 |
| 176 | +2253.14 |
| 192 | +3455.89 |
| 224 | +4643.09 |
| 256 | +3154.57 |

These are selected homogeneous support vectors, not a cover. They use
binary64 inner arithmetic and have not undergone outward replay. Negative
scores indicate promising proposals; positive scores leave those bounds
open and do not exhibit bad codewords. This nine-coordinate pilot omits
the old track's mature-tail coordinates and later refinements, so its
scores should not be compared directly with the strongest old-track scores.

## Checks and reproduction

The seventeen tests in `test_support.py` include exhaustive weighted-support
enumeration, every triple of row weights for B <= 4, exact alignment
probabilities, shell-cap monotonicity, and total-constrained CDF domination.
They also check a nonsymmetric code where the averaged old/new laws differ
and a permutation-invariant code where they coincide. A further enumeration
checks uniform column-weight orderings conditional on the histogram.
The later tests also check positive total-input-weight tilts on both sides
of one, and reject invalid greedy substitutions for weighted CDFs.
Tests pass using
exact rational equality; no tolerances are used.

```sh
python -B research/workstreams/permutation_locality/independent_rows/test_support.py
python -B research/workstreams/permutation_locality/independent_rows/screen.py --penalty .95 --penalty .75
python -B research/workstreams/permutation_locality/independent_rows/pilot.py --cache tmp/categorical-q64-pilot-checked.npz
```

The optional `screen.py --compare-shared` reports the existing basis-lattice
and shortening-moment shared-route bounds. It deliberately labels them as
excluding the later dual refinements; they are not the strongest old bounds.
Displayed logarithms are rounded numerical summaries of exact count bounds.

The pilot cache can be regenerated by the parent `categorical_screen.py`
using `--cache tmp/categorical-q64-pilot-checked.npz` at its default tilt.
That cache is diagnostic data, not a certificate, and stays untracked.

`verify.py` connects the new counts to the eleven-coordinate inner bounds
and performs separately labeled outward support covers. See
[STRONG_INNER.md](STRONG_INNER.md) for the results and remaining gap.
The separate `--shell-cover` and `--prefix-rank` options combine individual
shell caps with cumulative caps; seven exact and outward tests check their
folding and cover interfaces. [WEIGHT_TILT.md](WEIGHT_TILT.md) describes
the optional `--weight-tilt`, which changes only the bound.

An optional `--operator-cache` stores exact dyadic operator uppers with
source and parameter checks. It is local memoization, not a certificate;
omit it when independently regenerating a proof run. Keep generator sources
unchanged throughout a cached run.

The experimental hill-climb drivers can also reuse source-bound local
families with `--epoch-cache`; see [DENSITY_EXTENSION.md](DENSITY_EXTENSION.md).
[FOUR_BIT_ATTACK.md](FOUR_BIT_ATTACK.md) records the later mixing and
history experiments, including an outward three-update bound for one
96-group support class. It does not extend the complete occupancy cover.
[OVERLAP_MEAN.md](OVERLAP_MEAN.md) records a newer exact overlap calculation
and local cancellation bound. Its 1000-shape census and small exact checks
pass, but it has not yet extended the complete occupancy coverage.

[DENSE_UNIFORM.md](DENSE_UNIFORM.md) gives a separate partial decomposition
bound. [ROW_CONDITIONING.md](ROW_CONDITIONING.md) explores packet-shape
averaging for the remaining intermediate and dense cases. Its outward
eleven-coordinate replay now gives 79.92 bits at q = 64 when every row of
every active group has weight in [64,192], allowing unequal row weights.
Messages with other row types remain open at that occupancy. Neither
partial method extends the unrestricted occupancy range by itself.
Keep improving the shared-route proof in parallel as a
research track; do not replace its driver or treat its certificate as a
certificate for independently shuffled rows.
