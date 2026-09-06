# Shared-permutation bounds for the zero-state family

The actual first-moment obstruction does not determine how often a setup
contains a bad message. We now have bounds on several dependencies needed
for that probability question. None is a full second-moment bound.

Keep the fixed BCH [256,128] code, inner t128_s15, and setup experiment
from `ZERO_STATE_SECOND_MOMENT.md`. Both auxiliary messages use the same
setup. In particular, a region permutation is shared between them.
Write beta_j for one region's kernel probability at input weight j.
Write K_2 for the complete pair enumerator of the 128-bit kernel.
For a region type m=(m_00,m_01,m_10,m_11), with sum 8192, the exact
shared-permutation probability is

\[
\beta_2(m)=\frac{[x^m]K_2(x)^{64}}
 {\binom{8192}{m_{00},m_{01},m_{10},m_{11}}}.
\]

## A discrete Fourier bound on the actual type coefficient

A positive-tilt coefficient bound loses about 15 bits at the central type.
The following finite Fourier sum retains the coefficient's concentration
factor without requiring a full coefficient table.

Fix positive probabilities p=(p_00,p_01,p_10,p_11), with sum one.
Let M_1,M_2,M_3 be positive integers and let zeta_j be a primitive M_j-th
root of unity. Set W(x):=K_2(x)^64. Discrete character orthogonality gives

\[
\frac1{M_1M_2M_3}\sum_{a,b,c}
 W(p_{00},p_{01}\zeta_1^a,p_{10}\zeta_2^b,p_{11}\zeta_3^c)
 \zeta_1^{-am_{01}}\zeta_2^{-bm_{10}}\zeta_3^{-cm_{11}}
 =\sum_{m'\equiv m}[x^{m'}]W(x)p^{m'}.
\]

The sums on the left use 0<=a<M_1, 0<=b<M_2, and 0<=c<M_3.
Congruence on the right concerns the three non-00 coordinates. Homogeneity
determines m'_00. Every coefficient of W is nonnegative, so the right side
is at least its target term. The triangle inequality therefore proves

\[
\beta_2(m)\le
\frac{\displaystyle\frac1{M_1M_2M_3}\sum_{a,b,c}
 |K_2(p_{00},p_{01}\zeta_1^a,p_{10}\zeta_2^b,p_{11}\zeta_3^c)|^{64}}
 {\displaystyle\binom{8192}{m}p^m}. \tag{1}
\]

Aliased coefficients increase this upper bound. There is no degree cutoff,
asymptotic approximation, or assumption that the grid resolves every coefficient.

### Outward numerical evaluation

`certify_pair_type_fourier_tracked.py` evaluates (1) on a 128x128x64 grid.
The tilt probabilities are exact dyadic rationals with denominator 2^40.
Arb encloses every root, then verifies component errors at most 2^-52 for
the binary64 approximations. Producer precision is 256 bits; replay uses
512 bits for the roots and final rational conversion.

Here is the error budget used by `pair_fourier_error.py`. Let u:=2^-53.
Let ell_j be the four computed signed linear forms in the pair MacWilliams
formula. The exact forms are within epsilon_0:=16u=2^-49 in complex norm.
Root error contributes at most 4u in complex L1 norm before weighting.
Probability multiplication and summation of four signed terms fit inside
the remaining budget.

For each j, compute an upward bound R_j on |ell_j|+epsilon_0. The pair
formula has 155 terms with nonnegative exact dyadic coefficients c_d,
sum_d c_d=1, and total degree 128. Define

\[
A_R:=\sum_d c_d\prod_j R_j^{d_j},\qquad
D_R:=\sum_d c_d\prod_j R_j^{d_j}\sum_j d_j/R_j.
\]

The implementation rounds every nonnegative operation in these expressions
upward. It bounds the absolute error of the complex evaluation by

\[
2^{-40}A_R+2^{-49}D_R+2^{-1000}. \tag{2}
\]

The derivative term bounds root and linear-form perturbations along their
line segments. For arithmetic error, each monomial can be expanded into
at most 130 complex multiplication steps. A normal complex multiplication
has error at most 8u times the product of its input magnitudes. The exact
coefficient scaling and summation of 155 terms are covered by 156 factors
of 1+4u. The exact rational inequality

\[
(1+8u)^{130}(1+4u)^{156}-1<2^{-40}
\]

proves the first budget in (2). Propagated subnormal errors are below
2^-1000: intermediate magnitudes are bounded near one, and the expanded
arithmetic graphs contain far fewer than 2^20 operations. This calculation
assumes standard IEEE binary64 round-to-nearest arithmetic and square root.
Complex powers use explicit multiplication, not a library complex power.

After adding (2), the code scales the magnitude by 2^30 and performs six
upward-rounded squarings. This avoids underflow in the relevant 64th powers.
Any value lost at subnormal scale is covered by upward rounding. Each chunk
contains at most 2048 nonnegative terms; its sum is increased by 1+2^-38.
The exact check 1/(1-2047u)<1+2^-38 covers the reduction order. Chunks are
then accumulated as exact rationals, with the 2^1920 scaling removed exactly.

`test_pair_fourier_roundoff_budget.py` checks the rational budgets.
`test_pair_fourier_error.py` compares 48 grid points, including parity peaks,
against 512-bit complex balls. `test_pair_fourier_alias.py` checks the finite
alias identity exactly on a toy polynomial. These tests supplement the
error analysis; sampled agreement is not its justification.

### Certified types and their scope

For m=(6638,736,736,82), both single-message weights are 818. The outward
receipt `generated/pair_type_central_tracked_outward.json` proves

\[
\frac{\beta_2(m)}{\beta_{818}^2}<1.003217.
\]

The stored upper bound is approximately 1.0032168417224794. Its 512-bit
replay passed. The earlier uniform-error certificate, about 1.065659,
remains retained and valid.

The same Fourier numerator bounds other types by an exact change in the
denominator. For a certified bound B(m_0), equation (1) gives

\[
B(m)=B(m_0)\prod_i\frac{m_i!}{m_{0,i}!}p_i^{m_{0,i}-m_i}.
\]

`extend_pair_type_box.py` checks all 125 types with single weights
x,y in {814,816,818,820,822} and overlap k in {80,81,82,83,84}.
For m=(8192-x-y+k,y-k,x-k,k), every ratio B(m)/(beta_x beta_y)
is below 1.101. The largest stored upper bound is about 1.1004688048850855.
The attempted 1.1 threshold failed; the receipt preserves that false flag.
This is an enumerated finite box, not interpolation or coverage of all types.

## Joint parity mixing across all 256 regions

This result concerns only the 256 column parities of each message. It does
not include the remaining kernel constraints. Return temporarily to
T={c in C:38<=wt(c)<=80}, and sample each occupied row independently from T.
Fix any two occupied-row supports of size Q=2620, with intersection size r.
Average over row words and the shared independent row permutations.

For v in F_2^256 define mu_hat(v):=E_{U in T}(-1)^(v dot U).
There are two trivial characters, v=0 and v=1, because T contains even words.
For every other character, exact Krawtchouk calculations and shell caps give

\[
|\mathbb E_\pi\widehat\mu(\pi^{-1}v)|<2/5,\qquad
\mathbb E_\pi\widehat\mu(\pi^{-1}v)^2<3/8. \tag{3}
\]

To obtain the second bound, expand the square using independent U,V in T.
Their close-pair probability Pr[wt(U+V)<=80] is below 1/128. For all even
distances 82..160 and nontrivial character weights, the absolute normalized
Krawtchouk value is at most 23/64. Thus the second expectation is at most
p_close+(1-p_close)23/64<3/8. For the first bound, the exact sum of absolute
Krawtchouk values weighted by certified shell caps gives a bound below 2/5.

For two nontrivial characters on a shared row, Cauchy--Schwarz and (3)
bound their joint character factor by 3/8. Across all occupied rows, its
absolute value is at most

\[
(2/5)^{2(Q-r)}(3/8)^r\le(3/8)^Q.
\]

If exactly one character is nontrivial, the bound is (2/5)^Q. The joint
Fourier inversion has four trivial terms. Consequently the probability
P_even that both messages have even parity in every region satisfies

\[
\left|\frac{P_{\rm even}}{2^{-510}}-1\right|\le\epsilon,
\quad
\epsilon:=(2^{256}-2)(2/5)^Q+
 \frac{(2^{256}-2)^2}{4}(3/8)^Q<2^{-3100}. \tag{4}
\]

`certify_tail_parity_mixing.py` verifies (3)--(4) by exact rational arithmetic.
The bound is uniform over the two occupied supports. It is unweighted:
inserting functions of region weights inside the parity expectation requires
a further argument.

## Fixing the total weight of the candidate family

The first-moment lower bound can use the conserved total input weight.
`zero_state_affine_lower.py` verifies log(beta_j)>=a+bj for every even
j in [64,1536], with b<0. It checks every integer; convexity is not assumed.
On this region-count interval, sum_j J_j<=80Q implies

\[
\prod_{r=1}^{256}\beta_{J_r}\ge\exp(256a+80Qb).
\]

The retained parity and binomial-tail argument gives probability above
2^-256 for this interval with even region counts. The resulting actual
first-moment lower bound for T exceeds 2^20423. Its 512-bit replay passed.

For the next analysis, restrict each occupied row to
T_80:={c in C:wt(c)=80}. The input weight is then exactly 209600.
This is an actual message subfamily of the unchanged encoder.
Subtracting all certified shell upper bounds below 80 from the tail lower
bound proves |T_80|>=a_80, where log_2(a_80) is approximately 98.167565.

`certify_weight80_family.py` proves that this restricted family's actual
zero-state first moment exceeds 2^19593. It also bounds its close-pair
probability by approximately 2^-8.6074. The intersection proof uses the
single Christoffel point (i,j)=(c/2,80-c/2) for each difference weight c.
Both occupied words must have weight exactly 80. Repeated receipt replay
passed. A large expectation remains distinct from a failure-probability bound.

The parity proof also applies directly to T_80. Its largest absolute
nontrivial single-character average is exactly 3/8. Combining this with
its own close-pair bound in (3)--(4) gives relative error below 2^-3300,
uniformly over occupied supports. `verify_weight80_dependence.py` checks
this refinement by exact rational arithmetic. It remains unweighted.

## The first mixed dependence of exact-weight rows

Translation invariance of T_80 gives coordinate marginal p=5/16. Let
U,V be independent uniform T_80 words, and let pi be a shared uniform
coordinate permutation. Set X_i:=U_{pi(i)}-p and Y_i:=V_{pi(i)}-p.
Conditional on pi, X and Y are independent and both have mean zero.
Thus all cross covariances and mixed third moments vanish exactly.

Define v:=p(1-p), and let M be the covariance matrix of U before permutation.
Its diagonal is v, and every row sum is zero because wt(U)=80. Its
permutation average has diagonal v and off-diagonal -v/255. Let D be M
minus this average, and define

\[
\tau:=\frac1{256\cdot255}\sum_{i\ne j}D_{ij}^2.
\]

For i!=j and k!=l, the mixed fourth moment minus its marginal product is
E_pi D_{pi(i),pi(j)}D_{pi(k),pi(l)}. Symmetry and zero row sums give

| Relationship between unordered pairs {i,j}, {k,l} | Mixed fourth correction |
| --- | ---: |
| Equal | tau |
| Share exactly one coordinate | -tau/254 |
| Disjoint | 2 tau/(254*253) |

If either pair is diagonal, the correction is zero. For the shared-vertex
case, sum D_ij D_ik over k!=i,j and use sum_k D_ik=0. Summing over all
ordered pairs then gives the disjoint case. `test_fixed_weight_fourth.py`
checks these identities by exhaustive words and permutations on a transitive
fixed-weight toy family with nonzero tau.

We can bound tau without knowing T_80's full two-coordinate distribution.
Let c:=wt(U+V). Its mean is 110, and

\[
\mathbb E(55-c/2)^2
 =\frac{256^2}{255}v^2+256\cdot255\tau. \tag{5}
\]

For every possible nonzero c in {38,40,...,160}, the single-point
Christoffel bound gives an intersection cap I_c. Together with the shell
cap U_c, it gives Pr[c]<=min(1,U_c I_c/a_80^2); also Pr[c=0]<=1/a_80.
The producer proposes rational numbers alpha,lambda using an LP, then
checks the exact majorant

\[
(55-c/2)^2\le\alpha+\lambda c+h_c,
\quad h_c:=\max(0,(55-c/2)^2-\alpha-\lambda c).
\]

Summing against the actual probabilities, using their mean and upper caps,
certifies the objective independently of the optimizer's numerical accuracy.
`certify_weight80_fourth.py` proves

\[
\mathbb E(55-c/2)^2<225.584,\qquad \tau<0.003274.
\]

For two fixed occupied supports of size Q whose overlap is r, sum their
centered rows. Mixed fourth cumulants add over independent rows, so the
same tensor is multiplied by r. Each coordinate variance is Qv. Therefore
the largest normalized mixed fourth correction is bounded by

\[
\frac{r\tau}{Q^2v^2}\le\frac{\tau}{Qv^2}<\frac1{36000}.
\]

The receipt and exact toy test passed replay. The independent checker
`verify_weight80_dependence.py` reconstructs all 63 difference caps and
the dual majorant using exact fractions, without an optimizer or floating
point. This controls fourth moments, not all higher moments or the
exponentially small joint kernel event.

To replay these last results:

```powershell
python -B workstreams/bch_rm2sub_bridge/certify_weight80_family.py
python -B workstreams/bch_rm2sub_bridge/verify_weight80_dependence.py
python -B workstreams/bch_rm2sub_bridge/test_fixed_weight_fourth.py
python -B -m unittest discover -s workstreams/bch_rm2sub_bridge -p test_pair_fourier*.py
```

## Remaining obligation

We still need an upper bound on the full shared-setup second moment,
including atypical region types and all cross-region dependence. Neither
multiplying the central type ratios nor inserting (4) into a weighted
expectation is justified by these results. The exact-weight-80 family avoids
row-weight fluctuations and supplies a cleaner next experiment. A useful
route must combine its row dependence with the actual coefficients (1),
then normalize against an adequately sharp first-moment lower bound.

Full positive distance coverage remains Q=1..1655, plus the retained dense
partial classes. No new complete occupancy is closed by this note.
