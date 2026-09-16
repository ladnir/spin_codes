# IMT at 10.99% asymptotic relative distance

This proof draft combines the new IMT inequalities with the existing
shared Golay--BA-3 outer argument. Every occupancy regime now has a bound.
The analytic reductions below require review before replacing the paper's
theorem. The numerical certificates do not machine-check these reductions.

## Construction and claim

Let L=128m. Let b be the least positive multiple of 24 satisfying
b >= (39/4) log2(Lb), and put N=Lb. Thus b=(39/4)log2(N)+O(1).
Sample the rate-half Golay--BA-3 outer once at length b, as defined in
[`SINGLE_SAMPLED_BA_RM2SUB_D11.md`](../../paper_architecture/certificates/single_sampled_ba_rm2sub/SINGLE_SAMPLED_BA_RM2SUB_D11.md).
Use that same outer at all L row positions. Independently sample a uniform
coordinate permutation for each row and a uniform permutation of the L bits
in each of the b regions. Concatenate the permuted regions.

The inner has t=128 input bits per epoch and s=19 state bits. Fix the expansion
map A and `weight5_seed0` feedback map B reconstructed by
[`candidate.py`](../asymmetric/bch256/weight5/candidate.py).
The certificates bind their exact columns. All A coordinates and B columns
are nonzero. Each epoch independently samples u from the nonzero s-bit
vectors and v from its orthogonal complement, including zero, then sets
M=I+uv^T. The inner computes

```text
Y_i = X_i + A q_i,
q_(i+1) = M_i q_i + B X_i,     q_0=0.
```

State persists across regions, and the final state is discarded. All setup
randomness is shared by all encoded messages. It is not resampled per message.

**Claim of this draft.** For this ensemble, as m tends to infinity,

```text
Pr[d_min <= floor(0.1099 N)] = o(1).
```

The probability includes the sampled outer, route permutations, and IMT
mixers. The code has rate exactly 1/2. The inner is invertible because its
block-triangular matrix has identity diagonal blocks. Ordinary and transposed
encoding retain O(N) work under the existing outer and routing cost model.
This is a native-length asymptotic claim, not a finite margin certificate.
Its outer is not the fixed BCH-256 code.

## Imported outer event

Write c=39/4, W=[0.104,0.896], and let h denote binary entropy in natural units.
The authenticated outer proof supplies an event G_b of probability 1-o(1).
On G_b all nonzero outer weights lie in bW. If A_b(w) counts outer words of
weight w, then uniformly on these weights,

```text
log A_b(w) <= b a_hat(w/b) + O(log b).
```

Here a_hat is the certified concave, piecewise-affine majorant in
`golay_ba3_concave_majorant.json`. In particular, for ell=0.01281 there exists
epsilon_b=O(log(b)/b), independent of w, such that

```text
A_b(w)/binom(b,w) <= exp[b(-ln(2)/2 + ell + epsilon_b)].       (1)
```

The assembly checker verifies a_hat(x)-h(x)+ln(2)/2 < ell at every
segment endpoint. Convexity on each affine segment proves the whole interval.
The outer event and finite-enumerator remainder are imported mathematical
results, not inferred from a grid. The portable artifact manifest authenticates
their 31 retained dependencies. The outer event is paid once, not L times.

Condition on any outer satisfying G_b for the remainder of the argument.
Let Z_Q count messages with exactly Q nonzero outer rows whose output weight
is at most floor(delta N), where delta=0.1099. Expectations below are over
the remaining route and mixer randomness. We will show sum_Q E Z_Q=o(1).

## Positive occupancy

Set alpha=Q/L. Fix a tuple of active-row relative weights x_1,...,x_Q in W,
and let x be their mean. The existing route argument has this useful form:

```text
moment <= (b+1)^Q (L+1)^b
          exp[N KL(alpha*x || p*y)] * iid_input_moment(p*y).  (2)
```

To obtain (2), first replace row j by iid Bernoulli-x_j bits conditioned on
its weight. The inverse conditioning probability is at most b+1. After a
uniform region permutation, the Poisson-binomial region law is pointwise at
most L+1 times iid Bernoulli-(alpha*x). This follows from the Chernoff upper
bound and binomial type lower bound at each region weight. The original
array has exactly N alpha x ones, so changing the reference probability to
py costs the displayed likelihood ratio. Finally,

```text
KL(alpha*x || p*y) <= KL(alpha || p) + alpha KL(x || y).
```

Concavity of a_hat bounds the number of messages in the fixed weight tuple
by exp[Qb a_hat(x)+O(Q log b)]. There are at most (b+1)^Q weight tuples.
The choice of active positions contributes at most 2^L. All these factors,
including those in (2), are exp(o(N)), uniformly for alpha >= 10^-4.

For each certified box choose its fixed p,y and z in (0,1). Its IMT
Collatz or scalar bound supplies the iid moment exponent log(rho)/128.
Markov's inequality therefore gives the exponent per output bit

```text
alpha*a_hat(x) + KL(alpha||p) + alpha*KL(x||y)
               + log(rho)/128 - delta*log(z).                (3)
```

`DENSE_OUTWARD_D1099_v3.json` covers alpha in [10^-4,1] and x in W with
658 boxes. Exact geometry and 256-bit outward arithmetic bound (3) by
-4.3631360128178675e-7. The 512-bit replay passes every retained bound.
For one affine outer segment, (3) is convex in (alpha,alpha*x); hence four
vertex checks cover its quadrilateral. There are finitely many witnesses,
so their moment prefactors are uniformly bounded. Consequently this entire
range contributes at most L exp(-eta N/2)=o(1), for any fixed positive eta
smaller than the magnitude of the retained maximum.

## Uniform sparse bound

The sparse range needs a different route comparison. Paying (L+1)^b would
not suffice when Q grows slowly. Equation (1) instead dominates the counting
measure on each independently permuted active row by a uniform b-bit row,
with multiplicative factor exp[b(ln(2)/2+ell+epsilon_b)]. This is a bound on
the sum over messages, not a claim that a fixed outer word is uniform.

Under this reference law, each region contains Q uniformly located marked
positions carrying independent fair bits. To generate it, sample iid marks
with probability p=(8/5)alpha and condition their count to equal Q. Give each
mark a fair bit. Without conditioning, the input bits are iid Bernoulli-beta,
where beta=(4/5)alpha. The regions are independent under the reference law.

For integer Q>=1 and alpha=Q/L<1, Q is a mode of Bin(L,alpha). Chebyshev's
inequality puts probability at least 3/4 in fewer than 4 sqrt(Q)+1 integer
positions around Q. Therefore its modal mass is at least 1/(8 sqrt(Q)).
The exact binomial likelihood ratio then gives

```text
Pr[Bin(L,p)=Q] >= exp[-L KL(alpha||p)] / (8 sqrt(Q)).          (4)
```

Set z=1-(8/5)alpha. The exact polynomial certificate proves

```text
T(beta,z) w(alpha) <= (1-96 alpha) w(alpha),
                         0 < alpha <= 10^-4.                (5)
```

T is the complete seven-state occupation envelope, not an entrywise mixture
of transfer families. The witness has zero-state coordinate 1; its other
coordinates are affine and at least 1/2048 throughout the interval.
The certificate uses x=10000 alpha in [0,1]. It constructs rational
polynomial upper bounds for all row actions. After removing each residual's
zero at x=0, an exact coefficient bound proves strict negativity. The
certificate also proves every maximum replacement used in those polynomials.
Thus (5) controls the whole interval, not only its first-order expansion.
It implies an iid-input moment bound of 2048(1-96 alpha)^(N/128).

Combining (1), (4), (5), and Markov's inequality yields

```text
E Z_Q <= 2048 binom(L,Q)
  * exp[Qb(ln(2)/2+ell+epsilon_b)] * (8 sqrt(Q))^b
  * exp[N KL(alpha||p)] * (1-96 alpha)^(N/128) * z^(-delta N).
                                                               (6)
```

Use log binom(L,Q) <= Q log L+Q and b>=c log2(N). Also,

```text
KL(alpha||(8/5)alpha)/alpha <= ln(5/8) + (3/5)/(1-0.00016),
-log(1-(8/5)alpha)/alpha <= (8/5)/(1-0.00016),
log(1-96 alpha)/(128 alpha) <= -96/128.
```

The resulting uniform coefficient, even at the larger delta=0.11, is

```text
C = ln(2)/2 + ln(5/8)
    + (3/5 + 0.11*8/5)/(1-0.00016) - 96/128
    + ell + ln(2)/c
  < -0.013403840578015145.                                    (7)
```

Dividing the logarithm of (6) by Qb adds only
epsilon_b + 1/b + ln(8 sqrt(Q))/Q + ln(2048)/(Qb) to C.
For Q>=4096, the third term decreases and is below 0.002. For sufficiently
large b, the other terms together are below 0.004, uniformly in Q.
Hence E Z_Q <= exp(-0.006 Qb) whenever
4096<=Q<=10^-4 L. Summing this geometric bound gives o(1).
No lower bound on the rate at which Q tends to infinity is required.

## Fixed occupancy: reduction to two states

We now handle the remaining finite set Q<4096. The following reduction is
valid for each fixed Q. Its constants may depend on Q; we do not apply it
to growing Q.

Write M=2^19-1 and r=1/M. On the nonzero states, an empty epoch has transition
P=(I+Pi)/2, where Pi replaces the state by a uniform nonzero state.
Let f(q)=wt(Aq). Since each A coordinate is nonzero, its stationary mean is
mu=128 p0, where p0=2^18/M. For an arbitrary initial nonzero state,

```text
E[f(q_j)|q_i] = mu + 2^(-(j-i)) (f(q_i)-mu),   j>=i.
```

For S_g=sum_{i=0}^{g-1} f(q_i), this identity gives

```text
|E S_g - mu*g| <= 2t,
Var(S_g) <= 3t^2*g,
E|S_g-mu*g| <= t sqrt(3g)+2t.                                (8)
```

For the variance bound, use
Cov(f(q_i),f(q_j))=2^(-(j-i)) Var(f(q_i)) and Var(f(q_i))<=t^2.
These statements do not assume a stationary initial state.

Fix theta>0 and use z_L=exp(-theta/L). Let W_0(z_L) be the full-state
weighted empty-epoch kernel. For each pair of live states q,q', the
Lipschitz bound for exp(-x), (8), and P^g=2^-g I+(1-2^-g)Pi give

```text
|W_0(z_L)^g(q,q') - exp(-theta*mu*g/L)/M|
    <= theta*t*(sqrt(3g)+2)/L + 2^-g.                        (9)
```

The zero state remains zero exactly. Set gamma=p0 theta and
D_gamma(u)=diag(1,exp(-gamma*u)). Define J to lift the classes zero/live
to all states, and R to project them to delta_0/uniform-live. Thus RJ=I_2.
For gaps g>=ceil(3 log2 L) and g<=L/t, (9) approximates the empty kernel by
J D_gamma(tg/L) R with uniform entry error O(L^-1/2).

Let W_1(z_L) be the kernel for one input bit at a uniformly selected local
position. Its weight differs from W_1(1) by O(theta*t/L) in row norm.
Invertibility of every mixer preserves the uniform-live distribution before
adding the nonzero B column. It follows exactly that

```text
R W_1(1) J = P0 = [[0,1],[r,1-r]].                           (10)
```

Consider a region with a fixed number a<=Q of uniformly placed one-bits.
The probability of a shared epoch is O_Q(1/L). The probability of a boundary
or inter-impulse gap shorter than ceil(3 log2 L) epochs is O_Q(log(L)/L).
Outside these events, replace every empty gap using (9) and each impulse
using (10). Products telescope with bounded row norms, since the number of
factors is fixed. The ordered uniform epoch sites converge by Riemann sums
to uniform spacings. The limiting region kernel is J K_a(theta) R, where

```text
K_a(theta) = a! integral_{u_i>=0, sum_{i=0}^a u_i=1}
  D_gamma(u_0) P0 D_gamma(u_1) ... P0 D_gamma(u_a) du.          (11)
```

The integral uses a-dimensional simplex coordinates; K_0=D_gamma(1).
The additive error in every entry is
O_{Q,M,theta}(L^-1/2+log(L)/L), uniformly over entering and exiting states.
Short-gap and collision events contribute at most their probability because
the actual weighted kernels have row sums at most one.

This also gives a componentwise multiplicative upper comparison. For a=0,
the zero/live cross entries vanish exactly in both kernels. For a=1, the
zero-to-zero entry vanishes exactly because B has no zero column. Every
other coarse entry is positive; for a>=2 all entries are positive.
Consequently, for each fixed Q and theta,

```text
F_{a,L}(z_L) <= (1+epsilon_L) J K_a(theta) R,
epsilon_L=O_{Q,M,theta}(L^-1/2+log(L)/L),    0<=a<=Q.          (12)
```

Here F is the actual full-state region kernel. Since b=O(log L), the
multiplicative loss across b regions is exp(o(1)). This comparison is
uniform over all regional input patterns and arbitrary boundary states.
It therefore remains valid under coefficient extraction enforcing row weights.

For convenience replace P0 by the nonnegative upper matrix
P_plus=[[0,1],[r,1]]. Every term of (11) increases entrywise.
Write K_a^+ for these larger matrices. They are not stochastic kernels.

## The fixed-Q inequalities

For Q=1,2 use the fair-row counting comparison (1). Each region has a
binomial number of one-bits, so its continuum matrix is 2^-Q T_Q, where
T_Q=sum_{a=0}^Q binom(Q,a) K_a^+. With gamma=Q ln(2), theta=gamma/p0, set

```text
a0=exp(-gamma), f=(1-a0)/gamma,
g=2*(gamma-1+a0)/gamma^2, e=2*(1-(1+gamma)*a0)/gamma^2.
K_0^+ = diag(1,a0),
K_1^+ = [[0,f],[r*f,a0]],
K_2^+ = [[r*g,e],[r*e,r*e+a0]].
```

The exponent per outer bit, including the choice of active rows, is

```text
log rho(T_Q) + delta*theta - Q ln(2)/2 + Q ell + Q ln(2)/c.   (13)
```

`ONE_TWO_OUTWARD.json` checks (13) at delta=0.11 with 256-bit outward
arithmetic and 512-bit replay. Dividing by Q, the upper bounds are below
-0.10918 for Q=1 and -0.10914 for Q=2. Thus E Z_Q tends to zero in both cases.

For fixed Q>=3 retain each row's weight through a fugacity u_j>0.
The continuum generating matrix is

```text
T_Q(u,theta) = sum_{e in {0,1}^Q} (product_j u_j^e_j) K_|e|^+(theta).
```

Coefficient extraction bounds the moment for row weights w_j by
product_j[binom(b,w_j)^-1 u_j^-w_j] times e_0 T_Q^b 1,
with the loss in (12). To bound T_Q, set

```text
sigma=127/250, tau=133/125, v=(1,3/1600), theta=Q*tau/p0,
G_Q = sup_{0<=y<=1} [-tau*y + (1+1/Q)ln(1-y+y/sigma)],
m_v(u) = max(1+u*sigma*v_1, u*r/v_1+sigma*(1+u)).
```

Here is why the old product bound applies even though P_plus is not
stochastic. Insert all Q potential impulse sites, using I at an unused site
and P_plus at a used site. On a fixed class path, let l be the number of live
spacings. Their total duration Y has distribution Beta(l,Q+1-l), with
deterministic endpoint interpretations. The elementary beta-integral identity

```text
E[(1-Y+Y/sigma)^(-(Q+1))] = sigma^l
```

and the definition of G_Q bound its weight by exp(QG_Q) sigma^l. The identity
follows by substituting u=y/[sigma+(1-sigma)y] into the beta density.
This pathwise argument uses nonnegative transition coefficients only.
Summing paths, averaging the Q! orders, and using the weighted row norm
||H||_v=max_i (Hv)_i/v_i gives

```text
||T_Q(u,Q*tau/p0)||_v <= exp(QG_Q) product_j m_v(u_j).
```

Since sigma<1, G_Q<=G_3 for Q>=3. For a row weight x=w/b the exponent
per active-row bit is therefore at most

```text
a_hat(x)-h(x)-x ln u + ln m_v(u) + G_3 + delta*tau/p0 + ln(2)/c.
                                                               (14)
```

Choose a fixed rational u separately for each affine outer segment.
Expression (14) is convex in x on that segment. `FIXED_LIMIT.json` checks
both endpoints of all 31 segments using P_plus, not P0. At delta=0.11 its
maximum upper bound is -0.0008679786206340028. The uniform outer remainder,
the (b+1)^Q weight tuples, and (12) contribute o_Q(b). Hence for each fixed
Q>=3, E Z_Q <= exp[-0.0008679 Qb+o_Q(b)].

## Summing the regimes

Use the fixed-Q result only for the finite set 1<=Q<4096. Every term tends
to zero, so their sum does too. Equation (6) and the fixed cutoff handle
4096<=Q<=10^-4 L by a uniform geometric sum. The positive-occupancy bound
handles the rest. Therefore sum_{Q=1}^L E Z_Q=o(1), uniformly for outers
satisfying G_b. Markov's inequality and Pr[G_b^c]=o(1) prove the claim.

The explicit cutoff matters: separately proving every fixed Q and every
growing sequence would not by itself justify summing all occupations.

## Evidence and review scope

`verify_imt_asymptotic.py` authenticates the retained inputs, checks coverage,
replays the small outward inequalities, and reconstructs the sparse polynomial
certificate. It binds this draft and the tests to its assembly record.
It does not formally verify the outer argument, transfer derivation, continuum
lemma, or asymptotic union. Those are mathematical review obligations.

The local receipts retain their original conditional scope and are not
rewritten to claim a theorem. No paper statement or implementation default
changes here. Recovering exactly 11% is a separate tightening task near dense
occupancy; it is not needed for the 10.99% draft.
