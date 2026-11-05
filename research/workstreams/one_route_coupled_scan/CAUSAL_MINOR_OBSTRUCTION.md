# A causal-minor obstruction for the additive-coset inner

The single-observation Cauchy lemma in
[SINGLE_DIAGONAL_RANK.md](SINGLE_DIAGONAL_RANK.md) does not extend to
arbitrary observation times. Both batch sizes 64 and 128 have a nonzero
2-by-2 impulse submatrix of rank one. Its four coordinates occur in the
actual full-Cantor K20 schedule, and none of its entries is a causal zero.

This is a counterexample to a proposed rank argument, not a low-distance
codeword or a failure of the whole code. The outer, random diagonal,
and the other inverse constraints are not imposed in this example.

## The matrix being tested

Work in the native field J=GF(2^32), with its AES-byte subfield E=GF(256).
Let Delta(x)=x^2+x. A full Cantor chain satisfies gamma_0=1 and
Delta(gamma_i)=gamma_(i-1). Define

    V_h=span_F2(gamma_0,...,gamma_(h-1))=ker Delta^h,
    b=2^h,     L_h(X)=Delta^h(X).

The monic vanishing polynomial of V_h is L_h, and its derivative is one.
Each output batch uses one coset a_t+V_h. Its moment matrix is
C_t[i,j]=(a_t+beta_j)^i, with beta_j the ordered points of V_h.
For a fixed output sequence, the inverse syndrome gives

    V_t=C_t^(-1) sum_(s<t) C_s Y_s.

For an observation point x in batch t and a source point z in batch s<t,
the impulse coefficient is

    K(x,z)=L_h(x+z)/(x+z).                                (1)

The denominator is nonzero because distinct batches have disjoint cosets.
For one fixed observation batch, the numerator in (1) is constant down
each source column. That matrix is a scaled Cauchy matrix. With several
observation batches, the numerator can create additional dependencies.

## A two-dimensional eigenspace produces a rank-one minor

Suppose U is a two-dimensional F2-subspace of J and lambda!=0 satisfies

    L_h(v)=lambda*v for every v in U.                     (2)

Then U intersect V_h={0}. For any a in J, the four points in a+U lie
in four distinct V_h cosets. Suppose all four cosets occur in the schedule.
Order these points by their batch times. Use the first two as source
columns and the last two as observation rows.

Every observation follows both sources. Every row-column difference is
a nonzero element of U, so (1)--(2) give the matrix

    [ lambda  lambda ]
    [ lambda  lambda ].                                  (3)

Its rank is one. The pattern of causal nonzeros would permit rank two,
and the state dimension b is larger than two. Therefore neither causal
zeros nor the usual rank-b cut through the state explains this defect.
Reordering these same four cosets does not remove the example: select
the first two and last two again after reordering.

## Batch 64: a GF(4) line

For h=6,

    L_6(X)=X^64+X^16+X^4+X.

This polynomial is GF(4)-linear. Choose any d in E outside V_6 and set
U=d*GF(4). Such d exists because |E|=256>|V_6|=64. For each v in U,

    L_6(v)=(L_6(d)/d)*v.

The scalar lambda=L_6(d)/d is nonzero. This proves (2).

The native tower uses u^2+u=delta with delta=0x20 in E and
Tr_(E/F2)(delta)=1. Consequently

    Delta^8(u)=1,     Delta^9(u)=0.

Choose a=u. Since U is contained in E=V_8 and u is outside V_8,
all four points lie in V_9 and none lies in V_6. This places the
counterexample well inside the full-Cantor K20 prefix.

## Batch 128: an explicit odd-level example

For h=7, all binary binomial coefficients of 7 are odd, so

    L_7(X)=sum_(i=0)^7 X^(2^i).

On E this is the absolute trace. Let zeta=0x02 in the AES representation
and mu=zeta^(-1)=0x8d. The minimal polynomial of zeta is

    X^8+X^4+X^3+X+1.

Thus mu has reciprocal minimal polynomial

    X^8+X^7+X^5+X^4+1.

Its degree-seven and degree-six coefficients give, respectively,

    sum_i mu^(2^i)=1,
    sum_(0<=j<i<=7) mu^(2^i+2^j)=0.                      (4)

Choose w in GF(2^16) satisfying w^2+w=mu. Such roots exist in that
quadratic extension and not in E, because the first trace in (4) is one.
Set

    d=mu,     e=mu*w,     U=span_F2(d,e).

The elements d,e are F2-independent. The first identity in (4) gives
L_7(d)=1=zeta*d. Iterating w^2=w+mu gives

    w^(2^i)=w+sum_(j=0)^(i-1) mu^(2^j).

Substitution into L_7(mu*w), followed by both identities in (4), yields

    L_7(mu*w)
      =w sum_i mu^(2^i)
        +sum_(0<=j<i<=7) mu^(2^i+2^j)
      =w=zeta*e.

Therefore (2) holds with lambda=zeta, including for d+e by F2-linearity.
The odd value h=7 does not prevent the obstruction.

This construction also has a direct native-tower realization. The native
field uses byte order [1,u,v,uv], with

    u^2+u=delta,     v^2+v=delta*u,     delta=0x20.

Choose the smaller byte t satisfying t^2+t=mu+delta, and take w=u+t.
The byte equation is solvable because its right side has trace zero.
It gives w^2+w=mu in exactly the native AES tower representation.

Now choose a=v. The first quadratic subfield is GF(2^16)=V_16. The
second quadratic has conjugation v -> v+1, so

    Delta^16(v)=v^(2^16)+v=1,     Delta^17(v)=0.

Hence a is in V_17 outside V_16, whereas U is contained in V_16.
All four points of a+U lie in V_17, are outside V_7, and occupy distinct
V_7 cosets. This proves actual-prefix inclusion without assigning an
arbitrary isomorphic representation to the byte constants.

## Why these cosets occur at K20

At K=2^20 and rate one half, N=2^21. In the full-Cantor schedule, the
quotient directions are gamma_h,gamma_(h+1),..., and batch t uses the
quotient coordinates g(t+1), where g(m)=m xor (m>>1).

The first 2^(21-h)-1 positive Gray indices enumerate every nonzero
element of V_21/V_h. All these indices occur among the n=N/b batches.
The examples above lie in V_9 and V_17, respectively. Neither uses the
omitted predecessor coset V_h. Thus all four required cosets occur in
the actual prefix, before any coordinates outside V_21 are needed.

This inclusion is stated for `SPIN_OR_FULL_CANTOR=1`. The algebraic
eigenspace examples are field identities, but another finite coset
schedule needs its own inclusion check.

## Consequences for the proof search

The valid single-observation lemma and the selected-time block-triangular
rank bound remain unchanged. In (3), the two observation rows belong to
different batches. It is therefore not a counterexample to either lemma.

A stronger global rank theorem must allow algebraic defects beyond the
rank-b state cut. Treating the causal impulse matrix as a generic Cauchy
matrix or a generic network realization is unsound for these points.
One possible next interface is a lower bound with an explicit defect
term for repeated eigenspace constraints, combined with the actual outer
parity checks. No such bound is established here. Merely changing the
order of a schedule containing the four cosets cannot repair (3).

## Bounded exact regression

`theory/test_causal_minor.py` uses the native AES quadratic tower. It
reconstructs the full Cantor chain, checks the trace and subtrace
identities, and constructs both examples. It converts each point to
Cantor coordinates and computes its actual Gray batch and lane index.
Finally it checks distinct cosets, strict source-before-observation
ordering, and the four equal nonzero impulse entries.

The script writes no files and runs serially with `python -B`.
The root agent's exact replay passed both cases. In chronological order,
the verified zero-based batch/lane coordinates are

| Batch size | First source | Second source | First observation | Second observation |
|---|---|---|---|---|
|64|(3,19)|(4,27)|(5,56)|(6,48)|
|128|(555,82)|(556,55)|(861,81)|(862,52)|

The replay checks all four equal nonzero entries in each case, not only
the zero determinant. This is an algebra regression, not a whole-code
distance or native-kernel performance test.
