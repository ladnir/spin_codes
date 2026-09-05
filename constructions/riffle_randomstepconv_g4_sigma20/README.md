# Riffle RandomStepConv g=4 sigma=20

**Status:** `PAUSED_ONE_LAP_BOUNDARY_OBSTRUCTION`

The candidate is paused. At 9% relative output weight, an explicit support
contained in a late suffix nearly saturates the inner upper bound. The active
two-lap successor, **Riffle RandomStepConv-2Lap g=4 sigma=20**, is also paused
because its measured improvement is too small relative to its compute cost.

This candidate replaces four lane-parallel accumulators with one random
packet-level state transition.  Fresh setup randomness at every position
makes the ensemble depend only on whether each input packet and state are
zero.  The inner analysis therefore has two state classes instead of four
accumulator histories or (2^{20}) concrete states.

The immediate objective is an efficient, proof-facing inner bound at the
full target size.  The exact transfer enumerator is not expanded into a
table.  Its rank-one active transition instead gives a positive sum over
terminated live episodes.  The implementation gives each summand its own
coefficient radius and optimizes one output-weight tilt.

The outer plan is recorded but deferred.  The repository contains the full
binary weight spectrum of the extended BCH constituent and the exact MDS
symbol-weight spectrum.  Independent within-block permutations convert each
BCH weight into an exact packet-support polynomial.  The two field parity
relations correlate the BCH values, so the product of the local average
polynomials is not automatically exact.  The intended proof will seek a
codimension-two moment inequality that consumes the full local spectrum.
Replacing every occupied BCH block by a minimum-weight word is only a
fallback.

## Work plan

1. Prove and implement the exact two-state inner transfer enumerator.
2. Derive an outward-rounded tilted coefficient bound for fixed packet
   support (h) and output threshold (D).
3. Validate the bound against exhaustive small instances, then evaluate the
   full (N=524352) packet instance.
4. Compress the inner bound into a certified envelope for every support from
   18 through \(N\).
5. Build the exact local BCH packet-support polynomial from the stored BCH
   spectrum.
6. Prove or refute the codimension-two moment interface for the two field
   parity equations.
7. Combine the outer spectrum and inner bound in a complete first-moment
   calculation.

The active work stops after Step 4 until the inner model has a credible
target-size distance curve.

## Goal 01 status

Goal 01 is complete.  Exhaustive small instances match the exact transfer and
rank-one gap identities.  The target evaluator reports outward-rounded bounds
at 5%, 9%, and 12% relative output weight for packet supports from 18 through
2048.  See `proof/GOAL_01_INNER_TRANSFER_BOUND_REPORT.md` and
`receipts/goal01_inner_transfer_bound.json`.

## Goal 02 status

Goal 02 is complete.  Eighteen support intervals cover every integer support
from 18 through \(N\).  Each interval stores one analytic certificate and two
outward-rounded endpoint evaluations.  Discrete convexity makes the endpoint
maximum a valid cap for the entire interval.  See
`proof/GOAL_02_ALL_SUPPORT_ENVELOPE_REPORT.md` and
`receipts/goal02_all_support_envelope.json`.

## Goal 03 status

Goal 03 is complete as a diagnostic.  Hölder handles one active data block.
A restricted Cauchy inequality handles two or more active data blocks using
only the ordinary BCH spectrum.  All tested occupation shells from 1 through
8, plus 16, are below one at the 5% threshold.  The result does not cover all
occupations.  See `proof/GOAL_03_BCH_MOMENT_GATE_REPORT.md` and
`receipts/goal03_bch_moment_gate.json`.

## Goal 04 status

Goal 04 is complete as a diagnostic.  A closed binomial expression sums all
16384 possible data-block occupations.  A convex secant retains the
support-dependent reciprocal-binomial factor from the inner bound.  The full
5% first moment is (2^{-19.7044}), with support 24--31 and occupation two
dominating.  The evaluation is not yet outward rounded.  See
`proof/GOAL_04_ALL_OCCUPATION_CLOSURE_REPORT.md` and
`receipts/goal04_all_occupation_closure.json`.
