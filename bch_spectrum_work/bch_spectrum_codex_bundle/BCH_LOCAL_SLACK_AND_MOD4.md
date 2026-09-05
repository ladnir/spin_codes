# Local slack and exact weight classes modulo four

Updated: 2026-09-05.

The new deterministic cap is

\[
A_{38}(C)\le674,685,258,232,216.
\]

This is only a small improvement over 678,661,807,513,927. The full certified
M22 bound improves from 36.11873 to 36.12681 bits. The sufficient condition
A38<=10^13 remains unproved, and the full-closure goal remains active.

## A small formulation with an exact negative result

The existing OA21 dual certificate implies

\[
h_{38}+\sum_w(r^q_w q_w+r^h_w h_w)\le U,
\qquad r^q_w,r^h_w\ge0.
\]

The script reconstructs U and every nonnegative reduced cost r exactly from
the saved rational dual multipliers. If the local moment constraints forced
the second term to be positive, subtracting that amount would improve U.

Use the split enumerator B from BCH_SPLIT38_CONSTRAINTS.md. Convert each
half-spectrum reduced cost into a symmetric full-weight cost: divide it by
two except at weight 128. For a split cell (i,j), add the q cost at i+j to
the h cost at 38-i+j. Denote the resulting nonnegative cost by c_i(j).
Complement symmetry ensures that the global slack equals sum_(i,j) B_i,j c_i(j).

For each inside weight i, the normalized outside distribution has mass one
and six vanishing normalized Krawtchouk moments. The small LP minimizes its
expected c_i(j). A proposed degree-six polynomial below c_i(j) supplies a
lower bound through its constant coefficient. Every pointwise inequality is
checked over exact fractions after the numerical proposal; optimizer output
alone is not a certificate.

All 20 complementary blocks instead admit an exact zero-cost primal witness.
Each saved witness is a nonnegative rational distribution of total mass one,
satisfies all six moment equations, and is supported only where c_i(j)=0.
Therefore the optimum of each small LP is exactly zero. No degree-six local
polynomial can improve this fixed global dual certificate through these
independent block constraints.

The witnesses do not have to share the same global q and h marginals.
Consequently, this result does not prove that the full coupled split LP is
insufficient. It rules out this particular independent-block correction.

The saved rational witnesses can be checked without running an optimizer:

    python -B code/verify_bch_local_slack.py

Receipt: `generated/bch256_local_slack_certificate.json`.

## Exact signed weight sums

For an even binary linear code E, define

\[
f(x)=\operatorname{wt}(x)/2\pmod2,\qquad
S(E)=\sum_{x\in E}(-1)^{f(x)}.
\]

The identity wt(x+y)=wt(x)+wt(y)-2|supp(x) intersect supp(y)| gives

\[
f(x+y)=f(x)+f(y)+x\mathbin\cdot y\pmod2.
\]

Thus f is a quadratic form whose polar form is the binary inner product.
For an even coset representative v, the signed sum on v+E is
(-1)^f(v) times the sum of (-1)^(f(x)+v dot x) over E.

The audit evaluates these sums by binary elimination. If two basis vectors
u,z have u dot z=1, split off their two-dimensional span. Its signed sum is
2(-1)^(g(u)g(z)), where g(x)=f(x)+v dot x. Replace each remaining basis vector
x by x+(x dot z)u+(x dot u)z to make it orthogonal to the pair. Repeat until
the polar form vanishes on the remaining radical.

On that radical, g is linear. Its sum is zero if g is nonzero on any basis
vector; otherwise the sum is the radical's size. This computes the signed
sum exactly without enumerating all codewords. The audit checks the procedure
against 54 exhaustively enumerable code/coset cases.

For the BCH codes P and Q and the fixed Wambach representative v, the result is

\[
S(P)=S(Q)=-2^{108},\qquad S(v+Q)=0.
\]

The polar ranks are 46 for P and 30 for Q. Their radical dimensions are
85 and 93. The common nonzero Q-coset spectrum h therefore has exactly
2^122 words in each of the two even weight classes modulo four. For Q,
the counts in classes zero and two modulo four are respectively
(2^123-2^108)/2 and (2^123+2^108)/2.

The relation S(P)=S(Q)+255 S(v+Q) provides another consistency check using
the already audited quotient symmetry. The signed sum on C is also -2^108,
because its spectrum is q+31h.

Reproduce the signed-sum audit with:

    python -B code/audit_bch_quadratic_sums.py --verify

## Certified small-LP refinement

Adding these two signed identities to the OA21 model yields 404 rows and
130 nonnegative variables. The exact solver found an optimum. A separate
checker verified all primal equations and inequalities, dual signs, every
dual column, and exact primal-dual objective equality.

The resulting A38 cap above uses 31 times the floor of the h38 bound; it
uses no orbit-based lattice rounding. Substituting the smaller cap into the
previous full M22 certificate gives the stated 36.12681-bit bound. No other
shell coefficient or higher-occupation bound changes.

The new exact primal still gives a weight-38 first-moment contribution above
2^-40 when multiplied by the certified true-tail lower bound. Hence OA21
plus these two signed identities remains insufficient even with exact inner
tails. This is a limitation of that relaxation, not a lower bound on the
actual BCH code's bad-setup probability.

    python -B code/audit_bch_mod4_cap.py --verify

Receipt: `generated/oa21_mod4_probe/audit.json`.

## Larger split models remain unresolved

An exact change to square-root-binomial column scales preserved the physical
feasible set but did not yield a certificate within the bounded QSopt_ex run.
The same model was also tested with SoPlex, configured for exact rational
input and solving. Its rational solution was unavailable at the time limit.
SoPlex's process exit code was zero despite that status; exit codes alone
were not treated as proof of success.

A further model removes 70 duplicate or modularly equality-dependent
inequalities, leaving 563 rows. Removing constraints is a safe weakening;
modular dependence is not asserted to prove rational equivalence. A 100-digit
SoPlex attempt on this model also reached its time limit without a rational
solution. Logs and process receipts are retained beside the models.

SoPlex's exact-solving mode uses rational checks and iterative refinement;
see its [official documentation](https://soplex.zib.de/doc-7.0.0/html/EXACT.php).
The local Ubuntu package is version 6.0.4. It was extracted inside this bundle,
not installed globally. These unsuccessful runs establish no new shell cap.

Next: use the exact radical/hull structure or genuinely coupled split
constraints to obtain a stronger counting inequality. Independent OA6 block
corrections and the two signed identities have now been tested and are not
enough. Preserve the original M22 target and the certified one-shell budget.
