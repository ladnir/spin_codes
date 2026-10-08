# Removing the dense nonzero-diagonal counting loss

The alternating five/six-output packet with outer GRS[128,74] has nominal rate
37/44, approximately 84.091%. Its old dense bound fails even at full inverse-input
weight. The failure comes from replacing exact nonzero MDS word counts by a
dimension bound. This note gives a separate refinement; the frozen dense proofs
remain unchanged.

The new [exact helper](theory/gf8_exact_spectrum_dense.py) covers the complete K20
dense interval from density 9/10 with a 184-bit bound. This is one component,
not a whole-code certificate. The lemma below applies to fixed or periodic
puncturing, provided the stated inner-map hypothesis holds.

## Fixed-vector membership and exact support counts

Fix a linear MDS code C of length n and dimension d over GF(q), and put p=n-d.
The outer consists of L rows of C. Independently permute the L row positions in
each of the n coordinate regions. Independently sample a nonzero GF(q) multiplier
for every coordinate. These are the same regional route and input diagonals as
in [the existing membership proof](REGIONAL_GLOBAL_ROUTE.md).

Fix a vector x of length N=nL. Let u_j be its weight in coordinate region j,
alpha_j=u_j/L, W=sum_j u_j, and alpha=W/N. Conditioning on any independently
sampled inner map does not change this remaining route/diagonal distribution.

For a prescribed support of size w, let B_w count codewords with exactly that
support. MDS shortening makes this count independent of the support positions.
Inclusion-exclusion gives the exact integers

\[
B_0=1,\qquad B_w=0\quad(1\le w\le p),
\]
\[
B_w=\sum_{j=0}^{w-p-1}(-1)^j\binom wj
       \left(q^{w-p-j}-1\right)\quad(w>p).
\tag{1}
\]

To obtain (1), fix j additional zero coordinates inside the support. The shortened
subcode has dimension max(w-p-j,0). Subtract its zero word before applying
inclusion-exclusion. The resulting B_w are nonnegative, even though their exact
integer computation contains alternating summands.

Define f_w=B_w/(q-1)^w. Given a row support of size w, f_w is exactly its
membership probability after the independent nonzero diagonal multipliers.
For a positive tilt lambda and independent Bernoulli support indicators with
means alpha_j, the tilted row polynomial is

\[
\Phi(\alpha_1,\ldots,\alpha_n;\lambda)
=\sum_{S\subseteq[n]} f_{|S|}\lambda^{|S|}
 \prod_{j\in S}\alpha_j\prod_{j\notin S}(1-\alpha_j).
\tag{2}
\]

The existing tilted-column conditioning argument, before its old MDS relaxation,
gives

\[
\Pr[x\text{ belongs to the routed, scaled outer}]
\le (L+1)^n\lambda^{-W}\Phi^L.
\tag{3}
\]

In each column, conditioning Bernoulli entries on their prescribed count gives
the required uniform subset. Comparing its tilted count probability to the
binomial mode at alpha_j costs at most L+1. This argument does not require the
regional means to be equal.

## A positive bound uniform over regional profiles

The exact sequence f_w need not be discretely convex after a small tilt. Applying
the old homogenization argument to that sequence without checking convexity would
be invalid. Instead, separate a constant membership baseline from a nonnegative
excess:

\[
b:=q^{-p},\qquad \delta_w:=\max(0,f_w-b).
\]

For every lambda>=1, define

\[
Q(\alpha,\lambda)
=b\bigl(1+(\lambda-1)\alpha\bigr)^n
 +\sum_{w=0}^n\binom nw\delta_w\lambda^w(1-\alpha)^{n-w}.
\tag{4}
\]

Then Phi<=Q for every regional profile with mean alpha.
To prove this, first use f_w<=b+delta_w in (2). The baseline contributes

\[
b\prod_{j=1}^n(1+(\lambda-1)\alpha_j)
\le b(1+(\lambda-1)\alpha)^n
\]

by arithmetic-geometric mean. For the excess of weight w, discard the factors
alpha_j on live coordinates. The remaining support sum is

\[
\sum_{|S|=w}\prod_{j\notin S}(1-\alpha_j)
=e_{n-w}(1-\alpha_1,\ldots,1-\alpha_n)
\le\binom nw(1-\alpha)^{n-w}.
\]

Here e_k is the elementary symmetric polynomial of degree k. Its displayed
bound is Maclaurin's inequality. It also follows by repeatedly replacing two
arguments by their mean: the coefficient of their product is the nonnegative
e_{k-2} of the other arguments. Every term multiplied by this inequality is
nonnegative. Summing proves (4).

No second use of the input diagonal randomness occurs. Exact MDS word counting
replaces an earlier combinatorial overcount; it is not an additional random
constraint. Combining (3) and (4) gives the uniform membership upper bound

\[
(L+1)^n\lambda^{-W}Q(\alpha,\lambda)^L.
\tag{5}
\]

At full density and lambda=1, Q=max(q^-p,f_n). The old bracket instead retained
q^-p(q/(q-1))^n. For the present parameters, this artificial factor was large
enough to prevent the dense proof from closing.

## Puncturing, integer cells, and the output union

Condition on the state diagonals. The retained inner map H has N inputs and M
outputs, has full row rank, and is independent of the input diagonals and route.
For the concrete packet, any fixed subset of its top output rows has an
invertible matching minor in the current-input feedthrough. Those packet minors
give a block-triangular pivot minor for H. This works with unequal retained
output counts in a public chronological period.

Every fixed output has exactly q^(N-M) preimages under H. Apply (5) separately
to each fixed preimage after conditioning on H. The setup-failure union for
GF256 therefore adds 8(N-M) bits and the output-ball exponent 1.27M.
The old exact inequality
(10/9)^100 * 2295^10 < 2^127 establishes that ball exponent at relative distance
1/10. Zero output is included in the ball.

For an integer inverse-input cell [u,v], put a=u/N and c=v/N. The following
exact rational bracket is uniform over the cell:

\[
Q_{[a,c]}(\lambda)
=b(1+(\lambda-1)c)^n
 +\sum_{w=0}^n\binom nw\delta_w\lambda^w(1-a)^{n-w}.
\]

The baseline increases with density and each excess term decreases. Also,
lambda^-W<=lambda^-u. Consequently, a lower bound on the failure margin is

\[
u\log_2\lambda-L\log_2 Q_{[a,c]}(\lambda)
-n\log_2(L+1)-1.27M-8(N-M).
\tag{6}
\]

The helper computes all B_w, delta_w, and endpoint brackets as exact fractions.
A floating optimizer only proposes a rational tilt. Outward Arb logarithms
decide each cell; failing cells are bisected. The cells establish one uniform
fixed-preimage bound. Their count is not another factor in the output union.

## Receipt and independent tests

For GRS[128,74], period (5,6), and requested K=2^20:

| Quantity | Value |
| --- | ---: |
| Outer rows L | 14,208 |
| Padded dimension | 1,051,392 |
| Raw length N | 1,818,624 |
| Retained length M | 1,250,304 |
| Forbidden output threshold | 125,030 |
| Dense weights | 1,636,762 through 1,818,624 |
| Directed cells / precision | 21 / 192 bits |
| Uniform dense margin | 184 bits |

The nominal rate is 37/44. The shortened requested rate is 1024/1221;
these rates should not be conflated. The witness digest is
`61c3a29a382a27cb266b31f1ba857cd50ad45d99e65eab93226a505f82d83136`.
The exact support-count digest is
`df2437289e98ccc055ec12a0ba54c7cf647995194a085894145041dd76002954`.

Eight tests pass. They independently enumerate small prime-field Reed-Solomon
codes, check complete GF256 spectrum cardinalities, and test 1,875 unequal
Bernoulli profiles against (4). They also test 256 profiles against the separate
convex-sequence bound, plus 36 GF256 profile/tilt cases through exact
Poisson-binomial convolution over 128 unequal regional means. Other checks cover
admissible coefficients, endpoint coverage, and removal of the full-support
diagonal overcount. The finite examples test the implementation; the universal
profile claim follows from the argument above.

```text
python -m unittest gf8_exact_spectrum_dense_tests -v
python gf8_exact_spectrum_dense.py --dimension 74 --schedule 5 6 --certify --first 9/10
```

Run from the `theory` directory. To bind an aggregate's canonical metadata,
call `certify_profile(g, first_density=Fraction(9,10))`. Match every construction
field and join the sparse and middle ranges before claiming a whole-code bound.
No native code or benchmark is part of this refinement.
