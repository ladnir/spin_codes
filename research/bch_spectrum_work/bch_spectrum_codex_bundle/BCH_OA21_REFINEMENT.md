# Stronger BCH moments: a certified refinement, not M22 closure

Updated: 2026-09-04.

The BCH-derived [256,128] outer admits stronger deterministic moment constraints
than the strength-15 constraints used previously. Strength 21 improves the
complete M22 first-moment bound to approximately 32.90151 bits. The required
40-bit bound remains open. No statistical shell caps enter this refinement.

## Published input and extension argument

Let Q0 be the primitive narrow-sense binary BCH code of length 255 and designed
distance 39. Let Q be its even extension to 256 coordinates. The existing
generator audit establishes dim(Q0)=dim(Q)=123 and Q contained in the fixed
BCH-derived outer C.

Augot and Levy-dit-Vehel report d(Q0 dual) >= 22 in the length-255 table on
page 10 of their author manuscript, *Bounds on the minimum distance of the
duals of BCH codes*. The first column gives the designed distance of the
original BCH code, not its dual. The row for designed distance 39 gives 22
in Schaub's bound column. Pages 8-9 explain that this is an algorithmic lower
bound valid for all nonzero words, rather than a sampling estimate.

Source: [author manuscript](https://www.lix.polytechnique.fr/~augot/IEEE_IT96.pdf),
published in IEEE Transactions on Information Theory 42(4), 1257-1260 (1996),
[DOI](https://doi.org/10.1109/18.508853).
The local source is `paper/augot_levy_dual_bch_1996.pdf`; its hash is in the
receipt. We use this published mathematical result. We have not independently
reproduced the underlying Schaub rank certificate.

The dual of an even extension is not generally the even extension of the
dual. We therefore justify the step to length 256 explicitly.

Index Q by F=GF(256), with the extra coordinate indexed by 0. Define

\[
p_j(c)=\sum_{x\in F}c_x x^j\quad(j\ge1),\qquad
p_0(c)=\sum_{x\in F}c_x.
\]

The code Q consists exactly of the binary vectors satisfying p_j(c)=0 for
0 <= j <= 38. The j=0 equation is even parity; the remaining equations are
the defining BCH root equations. For b in F, translating coordinate labels
by b preserves these equations: expand (x+b)^j by the binomial theorem, and
each term contains a vanishing p_i(c) with i <= j. Thus translations preserve
Q, and consequently preserve Q dual.

Suppose a nonzero v in Q dual had weight below 22. It has a zero coordinate.
Translate that coordinate to 0 and call the resulting vector v'. Then
v' lies in Q dual, has the same nonzero weight, and has v'_0=0. For every
c0 in Q0, orthogonality of v' to the even extension of c0 implies that the
puncturing of v' is orthogonal to c0. This produces a nonzero word in Q0 dual
of weight below 22, contradicting the published bound. Hence

\[
d(Q^\perp)\ge22.
\]

For any set of at most 21 coordinates, projection of Q onto those coordinates
is surjective. Otherwise a nonzero linear functional annihilating the image
would give a word of Q dual supported on that set. All fibers have equal
size, since projection is linear. The same uniform projection property holds
for every affine coset of Q. Thus Q and its cosets are orthogonal arrays of
strength 21.

## Exact LP audit

As in the prior envelope, q_w counts weight-w words in Q, and h_w counts them
in any one nonzero Q-coset inside P. The quotient symmetry gives
A_w(C)=q_w+31h_w. Both half-spectra are complement-symmetric.

For 1 <= j <= 21, the orthogonal-array property gives

\[
\sum_w q_w K_j(w)=\sum_w h_w K_j(w)=0,
\]

where K_j is the binary length-256 Krawtchouk polynomial. Odd degrees vanish
already by complement symmetry. We add degrees 16, 18, and 20 for each
half-spectrum to the prior strength-15 system: six additional equalities.

The checker independently reconstructs these coefficients using the
three-term Krawtchouk recurrence. It re-runs the prior audit of the BCH
containments and all baseline LP rows. Every empirical orbit-search lower
bound is replaced by zero. Each of the three new solutions passes exact
primal feasibility, dual signs, dual feasibility, and primal-dual equality.
Each LP has 402 rows and 130 nonnegative variables.

| Weight | Previous deterministic cap | Strength-21 cap |
| --- | ---: | ---: |
| 38 | 773,397,230,534,890 | 678,661,807,513,927 |
| 40 | 4,210,950,731,605,482 | 3,350,512,214,809,133 |
| 42 | 38,281,539,412,412,374 | 27,886,467,650,792,183 |

At weight 38 we use 31 times the floor of the exact h_38 bound. No
divisibility or orbit-based lattice rounding is used. The generated model
retains its original `provisional` labels; this note supplies the mathematical
applicability argument without modifying the saved model or solutions.

## Consequence for the original random-inner target

Keep the exact RandomStepConv-M22 setup, 8192 outer rows, 2^20 message bits,
2^21 output bits, and cutoff 209716. Reuse the certified full-Arb coefficients
for all 92 possible nonzero shells. Replace the three caps above and retain
the previous deterministic bounds at every other weight.

The exact rational aggregate is

\[
U=Q_1+\frac{128}{127}(Q_2+Q_3+T_{\ge4}).
\]

The factor accounts for removing the old orbit lower bounds, as proved in
RANDOM_INNER_M25_CLOSURE.md. Increasing no cap beyond the audited replacement
envelope preserves that comparison. The new receipt stores U exactly and
reports -log2(U) approximately 32.90151142. The logarithm is diagnostic; exact
rational comparison confirms that U is above 2^-40.

There is also an exact obstruction for this particular relaxation and these
fixed transfer coefficients. The feasible primal spectrum maximizing h_38
has Q1 objective approximately 2^-32.96558, above 2^-40. Therefore no valid
upper bound for that LP objective can be below 2^-40. Jointly optimizing the
same constraints and coefficients cannot establish the target.

This feasible spectrum need not belong to any actual code. The obstruction
does not refute the BCH conjecture or rule out sharper transfer coefficients,
additional algebraic constraints, integrality, or another proof method.

Reproduce the read-only exact checks with:

    python -B code/certify_bch_oa21_probe.py --verify

Receipt: `generated/oa21_closure_probe/audit.json`. The source model,
LP files, and exact solutions are retained beside it. Solves were sequential.

Next, quantify the loss in the one-row Chernoff bound before attempting a
large tail computation. In parallel in the research plan, seek BCH-specific
shell constraints beyond these moment identities; do not merely rerun this
LP hoping a different optimizer will close it. Run experiments sequentially.
The separate certified M25 alternative remains available, but has not been
adopted as a replacement for M22.
