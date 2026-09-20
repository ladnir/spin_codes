# Exact low shells of the hull dual

Updated: 2026-09-05.

The new code-specific results are

\[
d(HQ^\perp)=16,\qquad A_{14}(HQ^\perp)=0,\qquad
A_{16}(HQ^\perp)=16,592.
\]

Adding these two exact shell equations improves the full deterministic M22
margin to 37.56825962340 bits. The resulting cap is
A38(C) <= 227,085,449,947,170. The original 40-bit goal remains open. No
statistical shell bound or parameter change is used.

## Why the previous constraint search was stopped

Two more exact Farkas cuts were generated at successive optima. Their
certificates and mapped warm starts passed exact checks. They changed the
cap only negligibly. The final cut-only cap was 228,870,861,596,719.

A separate feasibility test then imposed the auxiliary requirement
h38 >= 5e12 on the full intersection model, including the reciprocal
constraint and the cuts through round two. This lower bound was used only
to seek a hypothetical witness; it is not claimed as a BCH property.
The solver returned a rational feasible point. Independent checks verified
all original intersection constraints, including the fixed coordinates
eliminated for numerical conditioning.

At that point, 31 h38 times the certified lower coefficient for the true
weight-38 inner tail exceeds 2^-40. Hence that relaxation cannot prove the
target even if its inner-tail calculation is made exact. This is neither
an actual BCH spectrum nor an actual failure-probability lower bound.
The witness predates the new shell equations below and does not establish
insufficiency after adding those equations to the full intersection model.

This evidence redirected the work toward additional code-specific constraints.

## A known supercode gives a completeness check

Retain the notation from BCH_HULL_SPECTRUM_CONSTRAINTS.md. Let A be the
dimension-69 hull of the published dimension-71 code L. Its spectrum is
the weight-0-mod-4 part of L's spectrum. Let D=A-perp. Since A is a subcode
of HQ, the hull dual HQ-perp is a subcode of D.

The exact MacWilliams transform of the known A spectrum gives

\[
A_{14}(D)=65,280,\qquad A_{16}(D)=3,408,432,
\]

and D has no nonzero words below weight 14. We enumerate words in D and
compare the number of distinct words with these known shell cardinalities.
Thus completeness does not depend on a classification of low-weight words.
The published spectrum remains an explicit mathematical input.

Coordinates are indexed by the 256 elements of F256. A binary affine flat
is a coset of a vector subspace over F2. Its indicator is a length-256 word.
Exact basis checks show that A, HP, and HQ are invariant under x -> x+1
and x -> 2x. These permutations generate all maps x -> ax+b with a nonzero.
Their dual codes are invariant as well because coordinate permutations
preserve the binary inner product.

## Weight 14: one complete affine orbit

The script enumerates all 97,155 binary three-dimensional linear subspaces
of F256 through reduced row echelon bases. For each subspace U, it computes
the syndrome of its eight-point indicator against a basis of A. Two
indicators with equal syndrome have a difference in D.

There are exactly 255 matching unordered pairs. Each pair intersects only
at zero, so its symmetric difference has weight 14. Starting from one such
word, the script generates every affine image. The orbit has 65,280 distinct
words. It therefore exhausts the known weight-14 shell of D.

Every orbit word is tested for membership in D and HQ-perp. None belongs
to HQ-perp. Together with HQ-perp subset D, this proves minimum distance
at least 16. The receipt retains a seed, its field support, the orbit size,
and a digest of the complete sorted word set. Replay regenerates the orbit.

## Weight 16: two disjoint enumerated classes

First enumerate all 200,787 four-dimensional linear subspaces. Every
16-point indicator lies in D, as checked by its exact syndrome. Each
subspace has 16 affine cosets. Affine invariance gives 3,212,592 distinct
weight-16 words in D. Distinct affine flats have distinct indicators, and
their direction subspaces are unique. Of the linear four-spaces, 1,037
have indicators in HQ-perp. Their affine cosets contribute 16,592 words.

For the second class, revisit the 255 matching pairs of three-spaces U,V.
Enumerate all 32 cosets of each. The syndrome of each translated indicator
is checked directly. A pair of disjoint cosets gives a weight-16 word in D.
The script retains 195,840 distinct such words; none lies in HQ-perp.

These words are not affine four-flats. The direction spaces U,V intersect
only at zero and span dimension six. Disjointness means the difference of
the two coset representatives lies outside U+V. Their union therefore has
affine span of dimension seven. This proves disjointness of the two classes.

The total count is

\[
3,212,592+195,840=3,408,432=A_{16}(D).
\]

Thus both classes together exhaust D's weight-16 shell, and the subcode
count A16(HQ-perp)=16,592 is exact. Since it is positive, the previously
proved lower distance bound is attained: d(HQ-perp)=16.

## Incorporation into the BCH bound

With the previous half spectra r,s and the length-256 Krawtchouk transform K,
the two new equations are

\[
\bigl(K(r+255s)\bigr)_{14}=0,\qquad
\bigl(K(r+255s)\bigr)_{16}=2^{93}\cdot16,592.
\]

The first new model has 196 variables and 1152 rows. Its exact certificate
gives A38(C) <= 228,602,379,051,008 and a 37.55989408948-bit M22 margin.
Adding the exact weight-16 count gives 1153 rows and the cap stated above.
All primal constraints, dual signs, dual columns, and objective equalities
are checked over rational arithmetic. The bounds use 31 times the floor
of the rational h38 optimum, without orbit-lattice rounding.

The smaller model's exact primal still exceeds the target when multiplied
by the true-tail lower coefficient. Further information is required for
that smaller relaxation. The full intersection model with these new shell
equations has not yet been tested.

Reproduce the checks sequentially:

    python -B code/bch_hull_cut_iteration.py verify 2
    python -B code/bch_hull_cut_iteration.py verify 3
    python -B code/bch_hull_intersection_feasibility.py verify
    python -B code/certify_bch_hull_dual_weight14.py --verify
    python -B code/certify_bch_hull_dual_weight16.py --verify
    python -B code/audit_bch_extended_hull_cap.py prepare_bch_hull_dual16_count_probe oa21_hull_dual16_probe --verify

Next: combine the exact shell equations with the full intersection model,
and attack D's weight-18 shell. Its known size can again certify completeness
of a constructive enumeration. No classification or proposed geometric
construction for weight 18 has yet been proved sufficient.
