# Reusing directed bounds across support classes

The alternating-pair candidate still lacks a whole-code distance certificate.
This note records a complete local cover at the unchanged 35% threshold.
The cover contains every profile in the declared box, not selected samples.

## One transform, many profiles

Fix the rational witness `(z, alpha, markers)` from the asymmetric checkpoint
in [the selective-pair analysis](HD_SHIFTED_GOOD_PAIRS.md).
The two regional polynomials do not depend on the target profile.
They depend on the witness, the local packet bounds, and the route.
Their coefficients are nonnegative.

[The coverage driver](theory/hd_shifted_coverage.py) computes each regional
polynomial without target-dependent truncation. It takes its directed
inverse transform once and retains all positive residue-alias upper bounds.
For each target box, it selects the required residues and takes entrywise
upper maxima. The stored similarity factors restore the original scaling.
The driver then applies that box's scalar normalization and support count.

This reuse changes neither the construction nor its probability model.
It must retain the same rational markers: changing only a target's marker
divisor would not produce a valid bound for a different witness.
The exact support-count generating function is also cached separately.

## A complete local cover

Write a profile as `(E,O,a,b)`: active even groups, active odd groups,
full pairs in the first matching, and full pairs in the second matching.
The covered box is

    E = 508..516, O = 476..484, a = 446..454, b = 416..424.

All 6,561 profiles in this box are feasible. The driver bisects a box when
its bound exceeds `2^-80`; its two children partition the parent exactly.
It closes this box after 151 evaluations, with 76 accepted leaves and no
unresolved leaves. The sum of their directed bounds is at most

    2^(-39837898990827/500000000000),

or approximately `2^-79.675797981654`. Replays at 256 and 384 Arb bits give
the same rational ceiling. The original radius-two box also closes in one
evaluation, at approximately `2^-120.922956520589`.

Each leaf's support count already includes every indexed support in that
leaf. The final sum is over leaves; it needs no additional profile-count
factor. The bisection tree keeps unfinished boxes explicit if an evaluation
budget is reached. It never labels a partial cover complete.

From the theory directory:

    python -B hd_shifted_coverage.py --radii 2 4 --precision 256
    python -B hd_shifted_coverage.py --radii 2 4 --precision 384
    python -B -m unittest hd_shifted_coverage_tests

Add `--details` to print every leaf. These are local cover receipts, not
whole-code certificates. No native benchmark or encoder change was made
for this result.

## Next reduction

The [orientation-independent bound](HD_SHIFTED_UNORIENTED.md) removes the
even/odd balance from the operator search. It leaves three support
parameters and a two-variable regional transform. That smaller search is
the next candidate for full coverage; the four-parameter machinery remains
available when retaining orientation improves a bound.
