# Code-specific constraints around a minimum word

Updated: 2026-09-04.

The remaining M22 proof obligation is a bound on A38. The ordinary OA21
relaxation is insufficient even with exact inner tails. We therefore retain
intersection information relative to one explicit minimum word.

This note establishes new exact constraints, not a new bound on A38. The
current deterministic M22 result remains 36.11873 bits. The full 40-bit result
still follows conditionally from A38<=10^13, as stated in
M22_ONE_SHELL_REDUCTION.md.

## A fixed partition with exact local information

Let P and Q be the even extensions of the binary length-255 BCH codes of
designed distances 37 and 39. Their dimensions are 131 and 123. Use the
published Wambach word already represented in the bundle, and extend it by
its parity coordinate. Call this vector v and its 38-coordinate support S.
Let T be the other 218 coordinates.

The new audit verifies v in P but not Q by exact generator-polynomial
division. It constructs all 123 shifted generator rows of Q, including their
parity coordinates. Binary elimination proves that projection onto S has
rank 38. Thus every pattern on S occurs exactly 2^85 times in Q.

Define R as the code on T obtained by requiring zero coordinates on S and
then deleting S. The audit constructs an explicit independent basis of R.
It has dimension 85 and minimum distance at least 40: deleting zero positions
does not change the weight of a word of Q.

The same audit proves

\[
d(R^\perp)\ge7.
\]

Write the generator matrix of R as an 85-by-218 binary matrix. A set of its
columns sums to zero exactly when its indicator is a dual word. The audit
checks that all columns are nonzero and distinct. It then generates and
sorts all 23,653 pair sums and all 1,703,016 triple sums.

- A zero triple or a pair matching a column would give a dependency of size at most three.
- Two distinct equal pair sums would give a dependency of size at most four.
- A triple matching a pair would give a dependency of size at most five.
- Two distinct equal triple sums would give a dependency of size at most six.

None occurs. Conversely, every dependency of size at most six appears in
one of these checks, by partitioning its support into sets of size at most
three. The calculation is exhaustive. Its operations are binary elimination,
unsigned-integer XOR, sorting, and comparisons; no sampling or floating-point
arithmetic enters this certificate.

The complete support, kernel basis, and 218 columns are saved in
`generated/bch256_wambach_shortening.json`. Replay with:

    python -B code/audit_bch_wambach_shortening.py --verify

## Consequences for the split enumerator

For 0<=i<=38 and 0<=j<=218, define

\[
B_{i,j}=\#\{q\in Q:\operatorname{wt}(q|_S)=i,
                         \operatorname{wt}(q|_T)=j\}.
\]

For any fixed pattern a on S, its outside fiber is an affine coset of R.
The dual-distance bound makes R, and every such coset, an orthogonal array
of strength 6. Indeed, a nonsurjective projection onto at most six coordinates
would supply a nonzero dual word supported there. Surjectivity gives equal
fibers because the projection is linear.

Let K_t^(n) be the binary length-n Krawtchouk polynomial. Summing over the
binom(38,i) inside patterns of weight i gives the exact identities

\[
\sum_j B_{i,j}=\binom{38}{i}2^{85},\qquad
\sum_j B_{i,j}K_t^{(218)}(j)=0\quad(1\le t\le6).
\]

These are stronger than the generic strength-21 identities for Q: they remain
valid after fixing all 38 inside bits. They use the chosen support S and do
not assert that every 38-coordinate projection has these properties.

For comparison, dual minimum distance 22 gives only a generic rank bound of
36 on 38-coordinate projections. A dual code supported on 38 positions
cannot have dimension three: its eight words would have pairwise-distance
sum at least 28*22=616, whereas 38 binary coordinates permit at most
38*16=608. Thus the supported dual has dimension at most two. We used this
weaker fact only in the initial exploratory model; the new model uses the
explicit rank-38 certificate instead.

The ordinary spectra are marginals of B:

\[
q_w=\sum_{i+j=w}B_{i,j},\qquad
h_w=\sum_{38-i+j=w}B_{i,j}.
\]

Here h is the spectrum of v+Q. The already audited quotient symmetry makes
this the common spectrum of every nonzero Q-coset in P. Therefore
A_w(C)=q_w+31h_w for the fixed BCH-derived outer.

Support restrictions follow from the two minimum distances. Except for the
zero and all-one words, i+j lies between 40 and 216. The coset weight
38-i+j lies between 38 and 218. Both are even. Complementation gives
B_(i,j)=B_(38-i,218-j).

The models use one variable for each complementary pair. At the single fixed
point (19,109), the stored variable is half of B_(19,109), and its coefficient
is counted twice. This is an exact reparameterization of the continuous LP;
no integrality constraint is imposed on that auxiliary variable.

Finally, the split MacWilliams transform is

\[
B^\perp_{u,t}=2^{-123}\sum_{i,j}B_{i,j}
                  K_u^{(38)}(i)K_t^{(218)}(j).
\]

To derive it, sum the character (-1)^(q dot z) over Q: the sum is 2^123 for
z in Q dual and zero otherwise. Group z by its two support weights. The
character sums over those groups are the two Krawtchouk factors.
Thus every transformed coefficient is nonnegative. Global dual distance 22
sets coefficients of total degree 1 through 21 to zero. The shortening result
also sets every coefficient with outside degree at most 6 to zero, except
the constant coefficient. Full projection rank handles outside degree zero.

An exhaustive RM(1,3) toy check validates the split transform, coset marginal,
and conditional-moment formulas on a small code where every word is known:

    python -B code/check_bch_split38_identities.py

## Optimization status

The initial split model has 1,696 nonnegative variables and 720 rows. Its
50-second and 240-second exact-solver attempts ended at their time limits,
without a certified optimum. Preliminary floating infeasibility messages
were not exact certificates and are not treated as mathematical results.

The code-specific model adds the local moment equalities. It eliminates 40
fixed variables and selects independent equality rows modulo 2^31-1, leaving
1,656 variables and 633 rows. Removing rows can only weaken the relaxation;
we do not claim that modular dependence proves rational equivalence. All
pre-selection rows are retained separately for a later primal check.

The 512-bit, 180-second attempt on this reduced model also reached its time
limit without a solution file. All three bounded solver processes are now
terminal. `generated/bch_split38_attempts.json` records these inconclusive
outcomes; none is an infeasibility certificate.

Sources and models are preserved under `generated/split38_probe` and
`generated/split38_local_probe`. Long LP lines were wrapped without changing
the token stream, to avoid the exact solver's line-length limit.

No A38 cap should be attributed to either model without a checked exact dual
witness. The prepared `code/verify_bch_split38_solution.py` checks primal rows,
dual signs and columns, and exact objective equality when a solution exists.
It also reports whether the primal satisfies the equality rows omitted during
modular selection. The verifier is not evidence of success merely because
its source exists.

Next: obtain a stable exact solve for the code-specific model, or derive
smaller local inequalities from these certified facts. Keep the existing
M22 bounds unchanged until a new certificate passes.
