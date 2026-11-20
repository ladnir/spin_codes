# Exact ingredients for the zero-state second moment

The first-moment obstruction is now established for the actual fixed
construction: at Q=2620, the expected number of bad messages exceeds
2^19600. `ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md` gives the proof and replays.
The remaining question is whether these messages occur in typical setups.
The ingredients below do not yet answer that probability question.

Keep T={c in C:38<=wt(c)<=80}, a=|T|, Q=2620, and L=8192.
Let S contain all messages with exactly Q rows in T and zero rows elsewhere.
Then |S|=binom(L,Q)a^Q. For setup theta, let Z_0(theta) count messages in S
whose every epoch has zero B-syndrome. Such messages have output weight
at most 209600 and are bad for the retained cutoff 209716.

Sample U,V independently from S, independently of setup theta. This is
auxiliary message sampling for the counting identities

\[
\mathbb E_\theta Z_0=|S|\Pr_{U,\theta}[U\text{ stays zero-state}],
\]
\[
\mathbb E_\theta Z_0^2=|S|^2\Pr_{U,V,\theta}
 [U,V\text{ both stay zero-state}]. \tag{1}
\]

Both messages use the **same** row and region permutations. Resampling
those permutations for V would compute a different quantity.
The goal is a verified upper bound on the second probability in (1).
Combined with the existing first-moment lower bound, it would imply

\[
\Pr_\theta[Z_0>0]\ge
 \frac{(\mathbb E_\theta Z_0)^2}{\mathbb E_\theta Z_0^2}. \tag{2}
\]

No useful bound on (2) has yet been certified.
The follow-up `PAIR_FOURIER_AND_PARITY.md` now bounds a finite set of actual
shared-region types and the joint parity event across all regions. It also
introduces an exact-weight-80 subfamily with a certified mixed fourth-moment
bound. Integration over all region types remains open.

## The fixed BCH tail is coordinate-transitive

Index the extended BCH coordinates by x in F256. Translation sends the
coordinate x to x+b. The code P satisfies the even-parity equation and
p_j(c)=sum_x c_x x^j=0 for 1<=j<=36. The binomial expansion shows that
translation preserves these equations. For j=37, its only proper binary
subsets are 0,1,4,5,32,33,36. Their syndromes vanish on P, so translation
also preserves p_37 itself.

Therefore every translation preserves C={c in P:p_37(c) in S_5}, where
S_5={0,...,31} is the fixed five-dimensional field subspace. Here S_5 is
distinct from the message family S defined above. Translations preserve
weight, so they also preserve T and each individual weight shell of C.
The translation group is transitive on all 256 coordinates.

`bch_translation_symmetry.py` checks the eight additive translation
generators on all 131 basis words of P. Each of the 1048 images remains
in P and has the same p_37 syndrome. The exact replay passed.

For a uniform row U_T in T, every coordinate consequently has the same
one-bit probability

\[
p:=\frac{\mathbb E[\operatorname{wt}(U_T)]}{256}.
\]

This does not assert independence between coordinates of U_T.

## One region has no setup-to-setup variation in its candidate fraction

Fix any setup theta and one region r. Under uniform auxiliary U in S,
the Q selected row positions are uniform. Their bits in region r are
independent Bernoulli(p), because the row words are independent and every
coordinate of T has marginal p. Thus the input vector of region r is
exchangeable, with total weight J distributed as Bin(Q,p).

The setup's fixed region permutation does not change this law. If beta_j
is the exact single-region kernel probability from the obstruction note,
the fraction of messages in S passing this region's kernel checks is

\[
\kappa:=\sum_{j=0}^{Q}\binom Qj p^j(1-p)^{Q-j}\beta_j. \tag{3}
\]

Equation (3) holds for every theta, not only on average over setups.
For independent U,V, conditional on the same fixed theta, the probability
that both pass this one region is exactly kappa^2.

It does **not** follow that the full second moment factors over regions.
The 256 coordinates within each BCH row are dependent. An exact toy test
has the same single-region count in every tested setup, but all-region
counts vary among 0,1,2,4. The unresolved dependence is across regions.

## Exact pair spectrum for the inner map

For the fixed map A:F_2^15 -> F_2^128, define

\[
P_{a,b,c}:=\#\{(u,v):\operatorname{wt}(Au)=a,
 \operatorname{wt}(Av)=b,\operatorname{wt}(A(u+v))=c\}.
\]

`inner_pair_spectrum.py` computes all these integers by XOR convolution.
There are 155 nonzero triples, totaling 2^30=1073741824 ordered pairs.
The calculation uses exact integer Walsh transforms on 32768 entries,
not enumeration of all ordered pairs. Integer overflow is excluded by
the bound (2^15)^4<2^63 on every absolute accumulated sum.

The replay checks the total, each marginal, permutation symmetry among
a,b,c, and the diagonal/zero cases. Exhaustive small-code tests check the
transform against direct pair enumeration.

This table gives the complete pair enumerator of ker(B). Introduce four
variables x_00,x_01,x_10,x_11, one for each binary coordinate pair, and define

\[
L_{ab}(x):=\sum_{e,f\in\{0,1\}}(-1)^{ae+bf}x_{ef}.
\]

For a triple (a,b,c), let

\[
n_{10}=(a+c-b)/2,\quad n_{01}=(b+c-a)/2,\quad
n_{11}=(a+b-c)/2,\quad n_{00}=128-n_{10}-n_{01}-n_{11}.
\]

Character orthogonality gives the homogeneous polynomial

\[
K_2(x):=2^{-30}\sum_{a,b,c}P_{a,b,c}
 L_{00}(x)^{n_{00}}L_{01}(x)^{n_{01}}
 L_{10}(x)^{n_{10}}L_{11}(x)^{n_{11}}. \tag{4}
\]

Its coefficient for a type counts ordered pairs in ker(B) having that
type. For coordinate pairs sampled independently with probabilities
p_00,p_01,p_10,p_11, K_2(p) is their exact joint zero-syndrome probability.
`pair_kernel_iid.py` evaluates this identity, and exact toy tests compare
it to direct kernel enumeration. Independent bit marginals give exactly
the product of the two single-kernel probabilities.

The actual region calculation is not iid. Given its four type counts
m=(m_00,m_01,m_10,m_11), summing to 8192, a shared uniform region permutation
has joint zero-syndrome probability

\[
\beta_2(m):=\frac{[x^m]K_2(x)^{64}}
 {8192!/(m_{00}!m_{01}!m_{10}!m_{11}!)}. \tag{5}
\]

The two input vectors must share the same permutation in (5).
Equation (4) is now explicit. The follow-up note certifies (5) at a central
type and a 125-type box. A bound covering the required type distribution,
and its integration over all 256 correlated region types, remain open.

## A certified bound on close BCH tail pairs

Two independent uniform words U_T,V_T in T satisfy

\[
\Pr[\operatorname{wt}(U_T+V_T)\le80]<1/128. \tag{6}
\]

This is an outward-certified statement about the actual BCH tail, not
about a random code. It includes the possibility U_T=V_T.

To prove (6), fix d in C of weight c>0, and split coordinates into its
support and complement. For u in C let i and j be its respective weights
on those two parts. Then wt(u)=i+j and wt(u+d)=c-i+j.

Let K_r^(n)(w) be the binary Krawtchouk polynomial. The audited dual distance
at least 30 makes all products of coordinate polynomials of total degree
at most 28 have the independent-uniform expectation. Consequently the
degree-at-most-14 products

\[
K_r^{(c)}(i)K_s^{(256-c)}(j),\qquad r+s\le14,
\]

are orthogonal under uniform C. Their squared norms are
binom(c,r)binom(256-c,s). The reproducing-kernel argument bounds the number
of C-words at each pair (i,j) by

\[
\frac{2^{128}}{D_c(i,j)},\qquad
D_c(i,j):=\sum_{r+s\le14}
 \frac{K_r^{(c)}(i)^2K_s^{(256-c)}(j)^2}
 {\binom cr\binom{256-c}s}. \tag{7}
\]

Indeed, the kernel polynomial centered at (i,j) has value D_c(i,j) at that
point and squared expectation D_c(i,j). The point's probability is at
most 1/D_c(i,j). This proves (7) uniformly for the fixed support of d.

Sum the rounded point caps over 38<=i+j<=80 and 38<=c-i+j<=80, with both
weights even. Intersect them with the ambient binomial count and the
certified one-weight caps. The resulting I_c bounds |T intersect (T+d)|.
If U_c is the certified C-shell cap and a_0 is the tail lower bound, then

\[
\Pr[\operatorname{wt}(U_T+V_T)\le80]
 \le\frac1{a_0}+\frac1{a_0^2}\sum_{c=38,40,\ldots,80}U_c I_c
 <2^{-7}. \tag{8}
\]

`certify_bch_tail_overlap.py` evaluates (7)--(8) outward at 256 bits.
Its 512-bit replay passed, with diagnostic margin 7.42979623196436 bits.
Exact small-code tests check bivariate orthogonality and every point cap.

Conditional on two occupied-row supports overlapping in r positions,
the r row-pair choices are independent. Thus their number of close pairs
is stochastically bounded by Bin(r,1/128). This controls unusually similar
row pairs, but it is not yet a bound on the complete second moment (1).

## Replay and next step

```powershell
python -B workstreams/bch_rm2sub_bridge/inner_pair_spectrum.py --verify
python -B workstreams/bch_rm2sub_bridge/bch_translation_symmetry.py --verify
python -B workstreams/bch_rm2sub_bridge/certify_bch_tail_overlap.py --verify
python -B workstreams/bch_rm2sub_bridge/test_inner_pair_spectrum.py
python -B workstreams/bch_rm2sub_bridge/test_pair_kernel_iid.py
python -B workstreams/bch_rm2sub_bridge/test_bivariate_christoffel.py
python -B workstreams/bch_rm2sub_bridge/test_single_region_invariance.py
```

Next, bound the shared-permutation pair probability in (5), retaining the
outer row dependence across regions. The exact pair spectrum, translation
symmetry, and close-pair cap are usable inputs. They do not license treating
all region constraints as independent for the two-message experiment.
