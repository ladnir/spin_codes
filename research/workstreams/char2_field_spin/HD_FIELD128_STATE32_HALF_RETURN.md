# Half-aware return counts for state 32

This note sharpens a local return cap for the conditional state-32
construction. It is not a distance certificate and does not change the
state-30 certificate. No coefficient computation or numerical validation
of this refinement is claimed. The older
[half-rank note](HD_FIELD128_HALF_RANK_RETURN.md) concerns different states.

## Setting and the exact count

Let E=GF(256), F=GF(2^128), q=|F|. Let A be the 64-by-32 Lagrange map
from nodes128..159 to nodes0..63, and C=A^T. Order output coordinates as
H0=0..31 and H1=32..63. The [Gram identity](HD_FIELD128_STATE32_GRAM.md) is

    P=A C=[I B; B^T I],       B B^T=B^T B=I,

where B is a nonzero scalar multiple of a Cauchy matrix over E. Every
square minor of B is nonsingular. Fix gamma of degree16 over E satisfying
the principal-nonsingularity property G(gamma) in the
[random-global-gamma note](HD_FIELD128_RANDOM_GAMMA32.md). Put

    lambda=1+gamma,       L=I+gamma P=[lambda I, gamma B;
                                    gamma B^T, lambda I].

Both gamma and lambda are nonzero. G(gamma) is an explicit hypothesis,
not a consequence of the degree condition at state32. If gamma is sampled
uniformly from F\GF(2^64), its failure probability is at most
16(2^62-1)/(q-2^64), charged once per setup, not once per packet or support.

Fix a potential input support S, with a_i=|S intersect H_i| and
j=a0+a1. Here 0<=a_i<=32. Let X be uniform on F^S and zero elsewhere.
This is the local counting experiment, not a claim that an actual outer
codeword has independent coordinates. Write R=C X. The incoming state
that returns to zero is Q=gamma R; for R!=0 its output is Y=L X.

Choose active zero sets D_i subset S intersect H_i and inactive zero
sets E_i subset H_i\S, of sizes c_i and e_i. Write c=c0+c1, e=e0+e1,
Z=D0 union D1 union E0 union E1, and u=min(32,j-c). For
k_*=rank_F L[Z,S], exact linear counting gives

    Pr[R!=0 and Y_Z=0] = (q^(j-k_*)-q^(j-c-u))/q^j.                 (1)

Indeed, Y_Z=0 has q^(j-k_*) solutions. When R=0, Y=X, so precisely the
coordinates D must vanish and C_(S\D) X_(S\D)=0. This kernel has
dimension j-c-u by the MDS property of A. It is contained in the first
kernel, hence k_*<=c+u. This subtraction is necessary: a lower bound
k=c+u proves the nonzero-syndrome event empty.

For e>=32 the event is empty because A_(E0 union E1) has full column
rank. It is also empty for c=j by G(gamma). Below take e<32 and c<j.
The [joint return count](HD_FIELD128_GAMMA_RETURN.md), conditional on G,
already implies the rank lower bound

    k_old=c+u-min(u,32-e).                                      (2)

## Two half-aware rank bounds

The inactive rows E0 constrain only X1, and E1 only X0. Their rank and
the two residual dimensions are exactly

    r_E=min(e1,a0)+min(e0,a1),
    d0=a0-min(e1,a0),       d1=a1-min(e0,a1).

Choose E-bases for these two inactive kernels. On the residual variables,
the coordinate projections D_i have ranks t_i=min(c_i,d_i). The cross
projections D0 on the half1 kernel and D1 on the half0 kernel have ranks
v0=min(c0,d1), v1=min(c1,d0). To see these assertions, stack coordinate
rows with Cauchy rows, or stack the distinct Cauchy rows E_i and D_i;
subtract the rank of the inactive constraints. Cauchy minors give the
stated residual ranks, independently of the particular sets.

Define f(p,q)=max(p+min(q,15),q+min(p,15)). Then

    k_* >= r_E+f(t0,t1),       k_* >= r_E+f(v0,v1).               (3)

For the first bound select p<=t0 and q<=t1 independent coordinate
projections and residual columns, normalized over E. The resulting
square minor has form [lambda I_p, gamma K; gamma H, lambda I_q],
with H and K over E; they need not be transposes. For p>=q its determinant
is

    lambda^(p-q) det(lambda^2 I_q+gamma^2 H K).

As a polynomial in gamma, the second factor is the square of an
E-polynomial of degree at most q, with constant term1. If q<=15 it
cannot vanish at a degree16 element. The case q>=p is symmetric.
Maximizing p+q subject to min(p,q)<=15 proves the first bound.

For the second bound normalize the independent cross projections and
reverse the two column groups. The minor is
[gamma I_p, lambda K; lambda H, gamma I_q]. Its determinant is
gamma^|p-q| times the square of a polynomial of degree at most min(p,q).
That polynomial is nonzero: at gamma=1 its square equals1. The same
degree argument and maximization prove the second bound. These two
arguments require no additional bad-gamma event.

A useful third bound restricts the residual variables to one half:

    k_* >= r_E+max(min(d0,c),min(d1,c)).                          (4)

For half0, stack coordinate rows D0 with Cauchy rows E1 union D1.
Eliminating c0 coordinate pivots gives rank
c0+min(a0-c0,e1+c1)=min(a0,c+e1). Subtracting min(a0,e1) leaves
min(d0,c). Half1 is symmetric. These restriction bounds are combined
by a maximum, not added.

Take k to be the maximum of (2), the two expressions in (3), and (4),
and set delta=c+u-k. Each is a proved rank lower bound, so delta>=0.
Equation (1) gives the uniform half-count cap

    p_bar=q^(-k)(1-q^(-delta)).                                 (5)

In particular delta=0 means exactly zero, not a clipped approximation.
For a0=a1=16, c0=e0=16, c1=e1=15, equation (2) gives k_old=31.
But r_E=31, (d0,d1)=(1,0), and (3) gives k=32=c+u. This deletes the
old positive q^-31(1-q^-1) contribution for every such choice of sets.

## Positive moment and transfer use

For 0<z<1, apply z^wt(Y)=product_h[z+(1-z)1_{Y_h=0}]. For fixed
(a0,a1), sum (5) over the allowed half-counts to obtain

    U_(a0,a1)(z) = sum_{c_i,e_i; c<j,e<32}
      product_i binom(a_i,c_i) binom(32-a_i,e_i)
      z^(64-c-e)(1-z)^(c+e) p_bar.

Here 0<=c_i<=a_i and 0<=e_i<=32-a_i; omitted cases are exactly zero.
All summands are nonnegative. Retain the minimum with every existing
valid basic/joint return cap. This bound is uniform over supports with
the stated half sizes; it does not average over independent half draws.

If the incoming weighted-state measure is independent of the comparison
input X, has nonzero mass at most W, and bounds each nonzero atom by
W/(q^32-1), its return contribution is at most
W U_(a0,a1)/(q^32-1). Thus only the return entry J/(q^32-1) changes.
The quiet, birth and live caps, and the refresh-density factor, are
unchanged. One must not subtract a return upper bound from a live cap.

The smallest occupancy-only interface takes the minimum of the old J_j
and max_a U_(a,j-a), where max(0,j-32)<=a<=min(32,j). This keeps the
existing 65 local occupancy matrices, but may lose much of the gain.
Retaining half sizes requires the actual paired-route multiplicities:
each half has16 pair slots, and fair swaps do not cross halves. For
h=(s0,f0,s1,f1), let a_i=s_i+2f_i and

    m(h)=product_i binom(16,s_i) binom(16-s_i,f_i),  s_i+f_i<=16.

Let P_h be the resulting local 2-by-2 transfer cap and 0<alpha<=1.
In a chronological two-packet coefficient, m(h)m(h') multiplies
(P_h P_h')^circ(alpha) outside the power. Singleton/full marker exponents
are the totals over both halves and packets. Do not average raw P_h
over halves before alpha. Regional normalization, original-r support
counts, projective prefactor, state continuity and terminal treatment
remain unchanged; fresh coefficients and coverage would still be needed.

First test the optimistic replacement J=0 in the existing proof pipeline.
It is not a valid return upper bound or certificate. If even that
replacement fails, sharpening only J cannot close that fixed-witness
calculation; inspect live/quiet caps, composition or coefficient loss
instead. A successful ablation only motivates this refinement. No claim
here establishes a 45% distance bound or an encoder performance gain.
