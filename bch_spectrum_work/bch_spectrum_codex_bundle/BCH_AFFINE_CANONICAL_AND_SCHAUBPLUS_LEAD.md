# Faster orbit certificates and a stronger-moment lead

Updated: 2026-09-05.

The certified M22 result remains 37.57061695409 bits, with
A38(C) <= 226,659,579,388,193. This note records a faster exact orbit method,
a partial weight-20 result, and a separate provisional literature-input
experiment. The 40-bit goal is not closed.

## Canonical affine orbits without full expansion

Let S be a support of size w >= 2 in F256. For each ordered pair x,y of
distinct elements of S, form

\[
N_{x,y}(S)=\{(z+x)/(y+x):z\in S\}.
\]

Each normalized support contains 0 and 1. Choose its least 256-bit
indicator, using the fixed field-to-coordinate map and integer ordering,
as the canonical representative c(S).

The normalized images are exactly those affine images of S that contain
both 0 and 1. Indeed, an affine map with such an image has unique preimages
x,y of 0,1 and must be z -> (z+x)/(y+x). Thus c(S) is constant on each
affine orbit. Equal canonical representatives imply equal orbits.

Let m(S) count ordered pairs whose normalized image is c(S). The affine
maps taking S to any fixed orbit member form a coset of the stabilizer of
S. Their number is therefore the stabilizer size. Consequently the orbit
has exactly 65,280/m(S) words. This method uses w(w-1) maps per support
instead of expanding all 65,280 maps.

The C++ search produces distinct canonical representatives and claimed
stabilizer sizes. An independent Python verifier reconstructs every
normalized image of each representative, checks its canonicality and
stabilizer, checks membership in the known supercode D, and tests membership
in HQ-perp. The previously verified affine invariance justifies extending
each membership test to its whole orbit. Distinct canonical keys prove
disjointness. Comparing their total orbit size with the known D shell
cardinality gives a completeness check that does not trust the search.

Regression tests compare the method with full affine expansion for supports
having both trivial and nontrivial stabilizers, including a four-element
subfield with stabilizer 12. They also match all 97 canonical keys and orbit
sizes in the independently completed weight-18 shell.

## Weight 20: a certified partial result

The known supercode has A20(D)=966,144,000. The translated-three-flat
candidate search produces 221,336 distinct words. Canonical verification
finds 2,040 disjoint orbits totaling 133,171,200 words, none in HQ-perp.
Therefore

\[
0\le A_{20}(HQ^\perp)\le966,144,000-133,171,200=832,972,800.
\]

This does not prove that the shell is empty. A second construction adds
each affine four-flat to representatives of the complete weight-14 and
weight-18 shells. Its 10,992 weight-20 candidates produce exactly the same
2,040 canonical records, byte for byte. These searches therefore supply
one partial union, not two disjoint contributions. This partial bound has
not been incorporated into the reported M22 certificate.

## Provisional SchaubPlus input

The indexed text of Ponchio and Sala's *A lower bound on the distance of
cyclic codes*, dated February 25, 2003, reports a Table 1 on manuscript
page 16. Its SchaubPlus lower-bound column gives 24 for the dual of the
length-255 dimension-123 BCH code (designed distance 39), and 26 for the
dual of the dimension-131 code (designed distance 37).
[Indexed manuscript](https://citeseerx.ist.psu.edu/document?doi=f4cc68b3184f5a89579bf7375f463e0ae7d4e4ea&repid=rep1&type=pdf).

The full PDF could not be retrieved for page inspection. The table layout
and underlying rank argument have not been independently audited here.
The manuscript is also cited as BCRI Preprint 2003, file BCRI-07.ps, in
reference 11 of [Kaida and Zheng (2015)](https://www.sciencepg.com/article/10.11648/j.pamj.s.2015040201.17).
These facts identify a lead, not a newly certified input.

If the reported bounds apply as stated, the translation-and-puncturing
argument in BCH_OA21_REFINEMENT.md gives d(Q-perp) >= 24 and
d(P-perp) >= 26 for our length-256 extensions. With A(P)=q+255h, the new
moment equations would be

\[
(Kq)_{22}=(Kh)_{22}=0,\qquad (K(q+255h))_{24}=0.
\]

The first pair follows because Q and all its cosets have uniform projections
onto any 23 coordinates. The last equation uses the stronger P-dual bound;
it does not assert strength 25 separately for Q or its cosets.

## Exact arithmetic for the provisional model

Adding those three equations to the audited hull model gives 196 variables
and 1157 normalized rows. Its rational optimum passes every primal and dual
check. Conditional on the new literature inputs, it would give

\[
A_{38}(C)\le90,005,732,934,389
\]

and a full M22 margin of 38.64383138558 bits with the saved other-weight
and occupation bounds. This is not below the required failure threshold.
The weight-38-only true-tail lower contribution at this optimum no longer
exceeds 2^-40. That observation removes one previously demonstrated
obstruction for this stronger model; it does not prove that a complete
first-moment bound can succeed.

The receipt is deliberately named `provisional_audit.json`, omits the
standard certified-cap and aggregate keys, and forces
`original_M22_target_closed` to false. Do not use this folder as the prior
certified envelope in an application audit.

## Replay and next step

Run sequentially:

    python -B code/check_bch_affine_canonical.py --verify
    python -B code/certify_bch_canonical_shell.py 18 hull_weight18_canonical bch256_hull_dual_weight18_canonical --verify
    python -B code/certify_bch_canonical_shell.py 20 hull_weight20_canonical bch256_hull_dual_weight20_partial --verify
    python -B code/certify_bch_canonical_shell.py 20 hull_weight20_flat_canonical bch256_hull_dual_weight20_flat_partial --verify
    python -B code/audit_bch_schaubplus_provisional.py --verify

Priority: verify the new dual-distance inputs by obtaining the primary
manuscript or independently producing a rank-bound certificate. If that
succeeds, test the strengthened full intersection model and a joint weighted
objective before committing to further shell enumeration. The current
insufficiency witness applies to the model without these new equations.
Sharper true-inner tails and improved weight-40/42 bounds may then matter;
none of these prospective improvements is claimed here.
