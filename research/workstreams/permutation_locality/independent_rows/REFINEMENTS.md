# Refinements of the row-conditioned bound

Three analysis refinements improve the outward bound for one restricted
message class from 79.928 bits to **209.98 bits**, rounded down. They do not
change the encoder or its setup distribution.

The subsequent density-column calculation in
[DENSE_DIAGNOSIS.md](DENSE_DIAGNOSIS.md) raises the same restricted-class
margin to 403.96 bits. The three-refinement checkpoint below is retained
with its original reproduction command.

The class has exactly 64 active groups among 2048 four-row groups. Every
active group has four nonzero BCH rows, and each row has weight in
[64,192]. Row weights may differ. The bound includes all group locations
and all row words in this class; it does **not** cover all messages with
64 active groups.

For K = 2^20 and N = 2^21, the event is that some message in this class
has encoded weight at most H = 209715. Probability is over the ideal
independent-row setup, including row shuffles, region placements, lane
permutations, and the two-transvection inner. The outward first-moment
upper is enclosed by

```text
[6.14137355906991052900563666622175240476309335677339401010e-64 +/- 4.26e-121]
```

Its reported log2 upper is
`-209.984836712759255526963251321972381078990999247676752353`.
In particular, the event probability is less than `6.141373559070e-64`.
This is a probability bound for the stated class, not a distance guarantee
for every sampled code.

## What changes in the analysis

The row-count reduction remains unchanged. With Bernoulli witness p = 1/2,
it dominates the class by a reference experiment with independent row
bits, multiplied by the group-location count and 256 row-density factors.
The reference moment uses tilt lambda = .052. Both parameters are proof
witnesses, not changes to the construction. See
[ROW_CONDITIONING.md](ROW_CONDITIONING.md) for the counting argument.

The inner tracks mature tilted mass M, its low-expansion-weight subset
masses L48 and L56, and a pointwise density coordinate C. Actual measures
satisfy `0 <= L48 <= L56 <= M` and `C <= M`. The latter is a density
inequality; C is not another disjoint state mass.

1. **Class-specific arrivals into the low-weight subsets**
   (`class_tail.py`, `--class-tail`). For each packet shape through eight
   occupied windows, refine the mature contribution to L48 and L56.
   One- and two-window tail censuses retain the source expansion-weight
   class. For larger shapes, the pair bound applies to an arbitrary
   translated state: other packets can change the source class, so the
   class-specific pair bound cannot simply be reused. The resulting
   class bounds are represented by nonnegative coefficients of
   `(M,L56,L48)`. Their cumulative sums are bounded by those of the old
   form. This comparison is valid on the nested-mass inequalities above,
   not by taking an independent minimum of each matrix entry. The code
   preserves exact dyadic coefficient differences and leaves a coupled column
   unchanged when the required representation is unavailable.

2. **Average high-occupancy output moments before maximizing**
   (`averaged_windows.py`, `--averaged-high`). Above eight occupied
   windows, the homogeneous Bernoulli reference supplies independent
   nonzero packet weights. For each exact expansion-window histogram,
   average the four possible packet weights in each window factor first.
   A positive coefficient recurrence then averages over uniformly chosen
   distinct windows. Taking the appropriate maximum or census-weighted
   average gives output-moment bounds for arbitrary mature states, fresh
   states, and uniform expansion classes. The refinement is restricted
   to the unsplit coarse mature-mass fallback. It does not add a feedback
   cancellation claim or assert that an output-tilted state is uniform.

3. **Enforce subset inequalities between regions**
   (`cone_moment.py`, `--clamp`). After each nonnegative region transfer,
   replace the upper for L56 by `min(L56,M)`, then L48 by
   `min(L48,L56)`, and C by `min(C,M)`. Each minimum is still an upper for
   the corresponding true quantity. Outward rounding precedes this
   operation. Repeating it for all 256 regions avoids propagating
   inconsistent auxiliary uppers; the terminal sum excludes the density
   and subset coordinates to avoid double counting.

All three transformations received independent semantic review. The
validation at this checkpoint passed 66 unit tests and 9012 exact
class-tail checks, including coupled-coordinate and direct toy-window
checks. These complement, rather than replace, the domination arguments.

## Reproduction and remaining coverage

Run from the repository root with NumPy, SciPy, and python-flint installed:

```sh
python -B research/workstreams/permutation_locality/independent_rows/row_verify.py --groups 64 --active-rows 4 --intervals 64:192 --tilt .052 --p 1/2 --precision 192 --cut 8 --full-feedback 6 --window-histogram 8 --joint-cancellation --class-tail --averaged-high --clamp --target-bits 52
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows -p 'test_*.py'
python -B research/workstreams/permutation_locality/independent_rows/class_tail.py
```

The replay reconstructs exact censuses and uses outward Arb arithmetic;
binary64 witness screens are not certificate inputs. The bound above uses
all three refinements together and does not assign separate bit gains to
them.

A binary64 sweep at q = 80, 96, 128, and 256, with tilts .052, .08, .12,
.16, and .24, did not close the tested middle-weight classes. This is a
failure of those screened bounds, not proof that the construction fails.
The full goal still requires other active-row patterns, row weights
outside [64,192], mixtures of group types, and all remaining occupancies.
The next useful step is to identify which retained state or shape
information limits those larger-occupancy bounds.
