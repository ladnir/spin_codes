# Eliminating the hull-intersection spectrum

Updated: 2026-09-05.

The first exact elimination cut proves A38(C) <= 228,870,861,644,645.
This improves the preceding cap by only 8,959 words, a relative change below
4e-11. The full M22 margin is 37.55841850571 bits and the 40-bit target remains
open. The principal result is an exact method for extracting small valid
inequalities from the larger, numerically difficult intersection model.

## The intersection and its spectrum

Use HP and HQ for the dimension-85 and dimension-93 hulls from
BCH_HULL_SPECTRUM_CONSTRAINTS.md. Let A be the dimension-69 hull of the
published code L71, whose spectrum a is known. Define J = HP intersect A.
Exact binary nullspaces verify

\[
\dim J=61,\qquad HP+A=HQ.
\]

The punctured cyclic generator quotient for A/J is 0x169. Multiplication
by x modulo that polynomial has one orbit on all 255 nonzero residues.
Consequently, all nonzero A/J cosets have a common spectrum u. If t is the
unknown J spectrum, then a=t+255u. The natural map A/J to HQ/HP is an
isomorphism: its kernel is A intersect HP = J, and both quotients have
dimension eight. The coset inclusions give, coefficientwise,

\[
0\le t\le r,\qquad t\le a,\qquad t+255s\ge a,
\]

where r and s retain their previous meanings. The spectrum t has total mass
2^61, is symmetric under complementation, and vanishes outside multiples of
four. It also vanishes wherever a does. In particular, t0=1 and its positive
weights start at 64 or later.

Taking duals of HP+A=HQ gives HP-perp intersect A-perp = HQ-perp. Both duals
are contained in J-perp. Their union therefore implies

\[
A(J^\perp)\ge A(HP^\perp)+A(A^\perp)-A(HQ^\perp).
\]

Let K denote the length-256 Krawtchouk transform. Substituting the previous
MacWilliams identities and multiplying by 2^93 yields the linear inequality

\[
2^{32}Kt-255Kr+255Ks\ge 2^{93}A(A^\perp).
\]

Self-orthogonality of J additionally gives Kt >= 2^61 t. These constraints
produce a 229-variable, 1329-row model. Exact substitution of the 16 fixed
t coordinates reduces it to 213 variables and 1297 rows. The substitution
preserves the rational feasible set; no constraint was dropped merely for
numerical convenience.

Two bounded exact-solver attempts ended at their CPU limits without a
certificate: 90 seconds for the original model and 120 seconds for the
reduced model, starting at 256 and 512 bits respectively. Neither numerical
infeasibility messages nor timeouts were treated as mathematical evidence.

The hull of the narrow-sense dimension-63 BCH code is not J. Both have
dimension 61, but their intersection has dimension 53. A search over all
invertible cyclic coordinate multipliers found no equivalence. This does
not rule out equivalence under an arbitrary coordinate permutation; no
published spectrum for that other hull was substituted for t.

## Reciprocal coset containment

The exact identity HQ intersect P-perp = HP gives
HQ minus HP contained in Q-perp minus P-perp. Since HQ=r+255s and

\[
A(Q^\perp)-A(P^\perp)=255\,2^{-131}K(q-h),
\]

we can add s <= 2^-131 K(q-h). These 33 inequalities require no new variables.
The resulting 1148-row model solved exactly. Its rational primal and dual
certificates prove that its optimum is unchanged from the anchor model.
Thus this reciprocal containment alone cannot improve the weight-38 bound.

## A certificate that eliminates t

Fix q,h,r,s to the exact optimum returned by the reciprocal model. The
remaining extension problem has only 17 t variables and 150 rows. The solver
reported infeasibility, but its saved status contained no explicit witness.
To obtain independently checkable evidence, we solved the alternative
linear program for nonnegative Farkas multipliers.

Orient every inequality as >= and replace an equality by its two orientations.
Normalize the rows by positive exact rational factors. Let z be the remaining
t coordinates divided by their positive column scales. Write the extension
problem as B z >= d with z >= 0. The saved multiplier vector y satisfies

\[
y\ge0,\qquad y^T B\le0,\qquad y^T d>0.
\]

These three relations are checked with exact fractions. They contradict
feasibility because y^T B z <= 0 for every nonnegative z. The extra constraint
sum y <= 1 only normalizes the search for y; it is not a BCH assumption.

Now unfix q,h,r,s while retaining the exact fixed t coordinates. Replay the
same 16 nonzero row multipliers against the original physical constraints.
All coefficients of the remaining t variables are nonpositive. Eliminating
those terms gives one valid inequality in q,h,r,s alone. Clearing its
denominators yields an integer-coefficient row. Direct substitution of the
old rational primal violates this new row by exactly the previously checked
positive Farkas objective.

This cut is valid for every spectrum satisfying the intersection constraints.
It is not an assumption fitted to the old primal. The old primal only
determines which valid cut the alternative linear program finds.

## Exact result and replay

The smaller model with this cut has 196 variables and 1149 rows. The checker
reconstructs the intersection algebra, Farkas witness, eliminated row, scaled
LP, and rational solution. All primal rows, dual signs, dual columns, and
primal-dual objective equality pass. The resulting cap and M22 margin appear
at the start of this note. The exact feasible primal still exceeds the
2^-40 target when combined with the certified true-tail lower coefficient.

Run these checks sequentially; none starts an optimizer:

    python -B code/audit_bch_hull_intersection_structure.py --verify
    python -B code/audit_bch_extended_hull_cap.py prepare_bch_hull_reciprocal_probe oa21_hull_anchor_probe --verify
    python -B code/audit_bch_extended_hull_cap.py prepare_bch_hull_benders_cut oa21_hull_reciprocal_probe --verify

The receipts are generated/bch256_hull_intersection_structure.json,
generated/oa21_hull_reciprocal_probe/audit.json, and
generated/oa21_hull_benders_cut/audit.json. Models, solver logs, and rational
witnesses are retained in the corresponding generated folders. All runs
were sequential; all solver processes are terminal.

Next: repeat the exact extension test at the new optimum, generate further
elimination cuts if needed, and reuse a correctly mapped simplex basis to
avoid restarting the larger optimization. Stop this line when a feasible
intersection extension proves that its relaxation is exhausted. If that
bound remains above the target, stronger structural information is required.
