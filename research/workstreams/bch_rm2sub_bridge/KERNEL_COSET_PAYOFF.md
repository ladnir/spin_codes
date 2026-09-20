# Limited coset payoff test: park this route

One bounded exact solve gives a gain of only **0.004507 bits**. Its checked
feasible spectrum limits any improvement from this augmented LP, with the
current transfer coefficients fixed, to **less than 0.005635 bits**.
This route is therefore parked. No further rank or codeword search was run.

The new full-score upper bound corresponds to 42.514535 bits at K=2^28,
compared with the frozen 42.510028-bit aggregate certificate. It retains the
same other-shell and higher-occupancy bounds. We did not replace the frozen
aggregate receipt for this negligible improvement.

## What was added

Keep P, Q, and E from `PDUAL_RANK_PROGRESS.md`. Define D=E dual. The exact
generator checks identify D with the even extension of the designed-distance-31
BCH code. Its parameters are [256,139,d>=32], and P is a subcode of dimension
131. The quotient of their punctured generator polynomials is 0x1e7.
Multiplication by X cycles through all 255 nonzero quotient elements.

Thus every nonzero P-coset in D has the same weight spectrum u. Write
p=q+255h for the spectrum of P, using the existing q and h variables.
Then the spectrum of D is p+255u. The new model gives u nonnegative,
complement-symmetric even-weight entries, mass 2^131, and zero entries at
weights below 32.

Let K denote the length-256 Krawtchouk transform. The additional constraints
are the nonnegativity and support constraints for the two spectra

\[
A(E)=2^{-139}K(p+255u),\qquad
A(g+E)=2^{-139}K(p-u),
\]

where g is any word of P dual with nonzero F7. The second spectrum is common
to the 255 nonzero E-cosets in P dual. E has minimum distance at least 32,
while these cosets inherit the current distance-30 lower bound from P dual.
We do **not** set their weight-30 entries to zero.

For completeness, the formulas follow by applying MacWilliams to D and P:
K(p+255u)=2^139 A(E) and Kp=2^131(A(E)+255 A(g+E)). Subtraction gives the
second displayed formula. A small exhaustive-code example checks the same
normalizations with a quotient of size four.

We also include two containment constraints. The hull HP=P intersect P dual
is contained in E, because every word of P has F7=0. Duality then gives
D contained in HP dual. Their spectra were already represented in the old
model, so these inclusions require no new auxiliary spectrum beyond u.
The new model adds 65 nonnegative variables to the old 196-variable LP.

The checker verifies D=E dual by showing that the generator rows of P and
the eight F7 check vectors span D. It replays the previous distance-32 kernel
certificate. The tests separately check F7=0 on every generator row of P.
These are established structure constraints, not reference-spectrum assumptions.

## Exact outcome and its limit

The solver finished in 31.95 seconds, within the 60-second limit. An independent
rational audit checked primal feasibility, dual signs and columns, and exact
primal-dual equality. The objective uses the saved outward Q1 coefficients
for BCH-256 with the selected t64/s20 inner map at K=2^28.

The feasible witness also evaluates the full fixed Q1 objective. Every upper
bound valid over this augmented feasible set must exceed that witness value.
The resulting 0.005635-bit ceiling rules out a meaningful gain from merely
reoptimizing this same model or its objective rounding. It does not rule out
stronger information about the coset, including a future proof that its
weight-30 shell is empty. We are parking that harder question on ROI grounds,
not declaring it false or mathematically impossible.

This is a distance-bound calculation over the same setup randomness as the
existing certificate. It is not a new full SPIN security theorem or evidence
that the heuristic BCH spectrum is accurate.

## Next work and replay

Return to the two main tracks: direct weighted low-shell bounds for BCH-256,
and calibration/evidence for the spectrum growth estimator. Keep E and the
completed rank witnesses available, but do not spend another turn expanding
this particular coset relaxation without a genuinely stronger counting fact.

The implementation is `kernel_coset_payoff.py`; the exact solution and audit
are under `generated/kernel_coset_payoff_v1/`. No solver is invoked by replay:

```text
python -B kernel_coset_payoff.py --verify
python -B -m unittest test_kernel_coset_payoff
```

The CWC presentation distinguishes the checked LP bound, the limit of its
relaxation, and the still-unproved stronger code statement. Frozen source
files and previous certificates were not modified.
