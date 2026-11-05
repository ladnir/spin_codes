# Group-weight markers in the positive-moment comparison

The F1 late-alpha control reduces the quiet, all-zero-state component to
+295.013579 log2 bits, but does not close it. Its input-weight floor is
aggregate: eighteen arity-eight groups contribute at least 1296 nonzero
pre-H coordinates. Each group separately contributes at least 72. This
note proves that the positive-moment comparison can retain those group
labels. It does not implement or evaluate the resulting floor calculation.

The reference is `theory/hd_two_sum_late_alpha_source_f1.py`. Its frozen
control uses r=18, arity8, z=691/700, alpha=1/5, and source stratum
F=max(a,b)=1. Restricting regional subsets to f<=1 gives
log2 R_00(18)=-184.866784. The single-line and two-atom controls are
+440.847312 and +295.013579, respectively. These are FLOAT components
of a majorant, not failure lower bounds or whole-code certificates.

## A coefficientwise inequality with labelled groups

Fix a finite group index set G. For each g in G, introduce a formal
variable X_g. Let L be a nonnegative real random variable whose moments
through order |G| are finite, and fix lambda>=0. For I subset G, define

\[
 Z_I(l)=\prod_{g\in I}(1+\lambda l X_g),\qquad
 F_I(\mathbf X)=\mathbb E_L Z_I(L).
\]

Empty products are one. Write P<=_coef Q when every coefficient of Q-P
is nonnegative. For every I,J subset G,

\[
 F_I F_J\ \le_{\rm coef}\ F_{I\cup J}F_{I\cap J}.                 \tag{1}
\]

To prove (1), let K=I intersect J, U=I minus J, and V=J minus I.
Take independent copies L,L' of the same random variable. Direct
expansion gives

\[
\begin{aligned}
 F_{I\cup J}F_{I\cap J}-F_I F_J
 =\frac12\mathbb E\big[&Z_K(L)Z_K(L')\\
 &\cdot(Z_U(L)-Z_U(L'))(Z_V(L)-Z_V(L'))\big].
\end{aligned}                                                     \tag{2}
\]

For a nonempty subset U' of U, the corresponding difference coefficient
is lambda^|U'| times (L^|U'|-L'^|U'|). Its sign is the sign of L-L'.
The analogous coefficient for V has the same sign. Thus their product
has nonnegative coefficients. The remaining factors in (2) also have
nonnegative coefficients. Taking expectations proves (1).

This is a multivariate statement. Setting every X_g equal to X recovers
an aggregate-weight comparison, but is not needed for the proof.

## Concentrating regional incidences without resampling row sets

Fix h groups and n regional coordinates. Each group has a fixed active
row set A_g subset [n] of size a. For coordinate i, put

\[
 I_i=\{g:i\in A_g\}.
\]

Each group occurs in exactly a of the incidence sets I_1,...,I_n. Then

\[
 \prod_{i=1}^n F_{I_i}(\mathbf X)
 \ \le_{\rm coef}\ F_G(\mathbf X)^a.                              \tag{3}
\]

Apply (1) whenever two incidence sets are incomparable under inclusion.
Replace them by their union and intersection. This preserves every
group's incidence count. It strictly increases the bounded integer
sum_i |I_i|^2 by 2|I minus J||J minus I|. Consequently the procedure
terminates with a nested family. In a nested family, every group of
degree a belongs to the same a largest sets. The final family therefore
has a copies of G and n-a empty sets, proving (3).

The row sets A_g were arbitrary and fixed throughout. This comparison
does not sample new row sets or assert that actual regional occupancies
are independent. The independent copies of L implement an auxiliary
positive-moment majorant of the regional curve.

Every coordinate in (3) uses the same law of L. A later good/bad
conditioning step can produce different moment laws across coordinates.
The present lemma does not justify uncrossing those heterogeneous
factors, nor reconstruct a labelled arrangement after an endpoint
coefficient maximum has discarded it.

## The quiet comparison with a marker for each group

For GF(256), take lambda=255 and eta=255^-2. The existing quiet bound
uses, for one group, the positive polynomial

\[
 (1-\eta)+\eta\prod_{i\in A_g}(1+255 X_g).
\]

It admits synthetic weights one and two. It is an upper comparison, not
an exact law for a nonzero quiet vector. Expand this mixture by the set
H of selected groups. Applying (3) to each H gives the labelled bound

\[
 P(\mathbf X)=
 \sum_{H\subseteq G}\eta^{|H|}(1-\eta)^{|G|-|H|}
 F_H(\mathbf X)^a.                                                \tag{4}
\]

For the two-atom law Pr[L=B]=p, Pr[L=0]=1-p, (4) has a compact positive
representation. Let

\[
 \pi_k={a\choose k}p^k(1-p)^{a-k},\qquad
 f_k(X)=(1-\eta)+\eta(1+255 B X)^k.
\]

Exchanging finite positive sums yields

\[
 P(\mathbf X)=\sum_{k=0}^a\pi_k\prod_{g\in G}f_k(X_g).             \tag{5}
\]

The integer k is auxiliary comparison randomness shared by the groups
within this column. It is not a fresh choice of their constituent row
sets. Formula (5) is an algebraic representation, not an authorization
to fractionalize its summands separately.

Positive products and any fixed positive coefficient-selection map
preserve (3) and (4). In particular, a whole-word calculation can retain
the exponent of each X_g before enforcing that group's input-weight
floor. Fractional powers also preserve individual coefficient
inequalities, but merging coefficients and fractionalizing are different
operations. A new implementation must keep the chosen late-alpha
chronology explicit; this lemma alone supplies no evaluated alpha bound.

For sixteen RS[16,8] columns and a=8, the necessary conditions are:

- Every active constituent row has weight at least 9.
- Therefore every group's total pre-H weight is at least 9a=72.
- Therefore the aggregate weight is at least 72r=1296 when r=18.

The implications do not reverse. Requiring every group to be seen once
is weaker still. Group72 floors retain more information than the current
aggregate floor, but do not enforce the eight individual row9 floors.
No improvement in bits is claimed until that distinction is evaluated.

Before implementing group floors, the cheaper order control should keep
the entire row-set sum and whole-word coefficient tail inside the
conditional first moment for one fixed source G. The existing uniform
per-row-set majorant permits this without assuming row-set independence.
For the selected positive component, Jensen can then be applied once;
the exact source count remains outside alpha. The alpha=1 control is
also admissible because the regional bank is unpowered. A subsequent
good/bad route conditioning control is a separate refinement.

Only after those controls should a group-floor implementation be
considered. There is no new provider, witness fit, construction change,
or production DP in this note.

## Exact small checks

`theory/hd_two_sum_group_uncrossing_tests.py` uses rational polynomial
coefficients. It checks (1) and identity (2) for every pair of subsets
on two and three variables, including overlapping sets. It also checks
(3) for all equal-degree three-group/two-of-three incidence arrays, and
checks (5) against explicit selected-group mixtures.

Run from `theory`:

```text
python -B -m unittest hd_two_sum_group_uncrossing_tests
```

These are finite algebra controls, not a numerical proof gate or a
distance certificate. All four tests passed in 0.113 seconds. The
argument above supplies the general lemma.
