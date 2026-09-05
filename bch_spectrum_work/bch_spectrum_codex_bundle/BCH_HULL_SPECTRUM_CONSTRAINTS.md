# Coupled hull spectra for the remaining weight-38 bound

Updated: 2026-09-05.

The coupled hull model with the derived anchor proves
A38(C) <= 228,870,861,653,604. The resulting full M22 margin is 37.55842 bits,
up from 36.12681 bits. The original 40-bit target is not closed.

## Exact subcodes

Retain P, Q, and the common nonzero Q-coset spectrum h from the existing
model. For an even code E, its hull is R_E = E intersect E-perp. The hull
is the radical of the binary inner product restricted to E.

An independent binary-nullspace calculation constructs explicit bases for
R_P and R_Q. Their dimensions are 85 and 93. Both are doubly even, contain
the all-one word, and have minimum distance at least 40. The calculation
also verifies

\[
R_P\subset R_Q\subset Q,\qquad P\cap R_Q^\perp=Q.
\]

Puncturing the parity coordinate gives cyclic codes. The generator quotient
for R_Q/R_P is 0x169. Multiplication by x modulo this polynomial visits all
255 nonzero residues. Thus the nonzero cosets have a common weight spectrum.
Let r be the spectrum of R_P and s that common coset spectrum. Then

\[
A(R_Q)=r+255s.
\]

Both r and s vanish outside weights divisible by four. They are symmetric
under weight complementation and each has total mass 2^85. The cosets have
no zero word. The audit supplies all bases and punctured defining root sets;
the hulls are not identified with narrow-sense BCH codes of these dimensions.

## Coupled inequalities

For a length-256 spectrum a, write (Ka)_j = sum_w a_w K_j(w), where K_j
is the binary Krawtchouk polynomial. MacWilliams gives

\[
A(R_P^\perp)=2^{-85}Kr,\qquad
A(R_Q^\perp)=2^{-93}K(r+255s).
\]

The containments above imply coefficientwise

\[
r+255s\le q,\quad r\le A(P^\perp),\quad
r+255s\le A(Q^\perp),\quad q\le A(R_Q^\perp).
\]

Also P minus Q is a subset of R_P-perp minus R_Q-perp. Subtracting the
two dual spectra therefore yields

\[
h\le 2^{-93}K(r-s).
\]

This last inequality uses a set difference with the verified intersection
P intersect R_Q-perp = Q; it does not subtract arbitrary containment bounds.

## Signed Fourier constraints

Fix E equal to P or Q, and set f(x)=wt(x)/2 modulo two. Define

\[
F_E(u)=\sum_{x\in E}(-1)^{f(x)+u\cdot x},\qquad
G_j=\sum_{\operatorname{wt}(u)=j}F_E(u).
\]

The previous quadratic-form audit gives sum_E (-1)^f = -2^108. Since f
vanishes on the hull, Fourier cancellation makes F_E vanish outside
R_E-perp. On R_E-perp its values have magnitude M=2^108. One proof squares
the sum and uses the quadratic identity: the square is |E| times the
character sum over the radical, hence |E||R_E| on R_E-perp and zero outside.

For u in E, shifting the summation variable gives
F_E(u)=-M(-1)^f(u). For u in E-perp, F_E(u)=-M. Let e_j, b_j, r_j, and d_j
count weight-j words in E, E-perp, R_E, and R_E-perp, respectively. Therefore

\[
\frac{d_j+G_j/M}{2}\ge
  \mathbf 1_{j\equiv2\ (4)}e_j,
\]

\[
\frac{d_j-G_j/M}{2}\ge
 b_j+\mathbf 1_{j\equiv0\ (4)}e_j-r_j.
\]

The subtraction counts the intersection of E-perp with the weight-0-mod-4
part of E. That intersection is exactly R_E. Finally,
G_j = sum_w (-1)^(w/2) A_w(E) K_j(w), so these are linear spectrum constraints.

The implementation clears every denominator using exact powers of two.
Exhaustive enumeration checks Fourier support, magnitude, and both signed
inclusions on 60 small codes with negative Gauss sums. These toy checks
supplement the algebra; they do not replace its proof for the BCH codes.

## First exact certificate and its limitation

The first model has 196 nonnegative variables and 1017 rows. QSopt_ex
returned a rational optimum. A separate checker reconstructed the model
and verified every primal row, dual sign, dual column, and objective equality.
The reported A38 cap uses 31 times the floor of the rational h38 bound,
without orbit-lattice rounding. All other contributions use the existing
one-shell M22 aggregation.

The exact feasible primal still exceeds the 2^-40 budget when its A38
value is multiplied by the certified true-tail lower bound. Thus even this
stronger relaxation cannot close the target by improving only the inner-tail
calculation. This is not a lower bound on the actual code's failure probability.

Reproduce the completed checks sequentially:

    python -B code/audit_bch_hulls.py --verify
    python -B code/check_bch_hull_fourier.py --verify
    python -B code/audit_bch_hull_cap.py --verify

Receipts: generated/bch256_hull_structure.json,
generated/bch256_hull_fourier_toy.json, generated/oa21_hull_probe/audit.json.

## A derived spectrum anchor

Let L be the published dimension-71 subcode used by the previous LP. Its
polar rank is exactly two, its hull has dimension 69, and its signed sum
is -2^70. Thus L has exactly (2^71-2^70)/2 = 2^69 words whose weights are
divisible by four. Its doubly-even hull already has that size. Consequently,
the hull spectrum is exactly the weight-0-mod-4 part of the published L
spectrum. This obtains a new exact spectrum without codeword enumeration.

Binary elimination verifies R_L subset R_Q. We can therefore add
A(R_L) <= r+255s and A(R_Q-perp) <= A(R_L-perp). The latter known spectrum
has minimum weight 14, improving the root-only lower bound of eight.
The augmented model has 1115 rows. Its inputs retain the same published
spectrum assumption as the original anchor-based proof.

The augmented solver returned an exact optimum. Independent rational checks
passed all 1115 rows, 196 variable nonnegativity checks, dual signs and columns,
and primal-dual objective equality. The improved cap is

\[
A_{38}(C)\le228,870,861,653,604.
\]

Substitution in the saved full M22 aggregation gives 37.55841850566 bits.
This is a deterministic-spectrum result under the same ideal random-inner
model, with the previously stated published anchor and OA21 inputs. It does
not use statistical shell caps. The sufficient A38 <= 10^13 condition remains
unproved. The new exact primal still fails the true-tail lower-bound test,
so this augmented relaxation also remains insufficient by itself.

    python -B code/audit_bch_hull_anchor_cap.py --verify

Receipt: generated/oa21_hull_anchor_probe/audit.json. Both solver runs are
terminal; they ran sequentially. Frozen earlier models and receipts were
not modified.

Next: couple the dimension-61 intersection R_L intersect R_P to both known
R_L and unknown R_P spectra. The rank computation already establishes this
intersection dimension. Its sum with R_P has dimension 93 and equals R_Q.
Its spectrum and coset constraints have not yet been added or certified.
