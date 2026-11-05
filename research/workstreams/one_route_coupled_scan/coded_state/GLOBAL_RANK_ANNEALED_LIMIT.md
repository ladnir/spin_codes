# What a rank-weighted flow reformulation does and does not gain

2026-11-01. A uniform flow tail is stronger than the geometry statement
needed by a particular codeword. One can avoid that uniform event by
retaining each support's cut score inside the field average. The exact
algebra below shows that this returns the original joint support/rank
enumerator. It does not create a gain merely by changing the proof
network's temporal capacity. A second identity isolates the genuinely
geometric quantity that a different global argument must control.

The encoder and all setup distributions remain fixed. Use
`N=2097152`, `K=N/2`, `W=943718`, `q=2^32`, and 8192 outer rows of
length 256 and dimension 128. The flow and field conventions are in
[QUARTET_GENERAL_FLOW_COVER.md](QUARTET_GENERAL_FLOW_COVER.md).

## The exact support-specific cancellation

Fix the routing and quartet partitions and an actual state support
`S`. Write `L,D,a_T,F,R` for its volume, internal deaths, terminal
population, footprint volume, and capped outer-row rank. If its
eligible set is `H`, with size `h`, then

```
R=128h+|V outside H|.
```

For network capacities `(a,b,c)`, put `p=1-a/c` and `beta=b/c`.
The normalized score of this particular cut is

```
g_H(S)=(1-p)L-beta*D-|V outside H|.
```

The flow deficit is the maximum of these scores over state supports.
Before taking that maximum, however, the rank relation is an identity:

```
R=(1-p)L-beta*D+128h-g_H(S).                         (1)
```

At `x=q^-p`, the fixed-support field bound, including output
cancellations and the free terminal, therefore becomes

```
gamma/(q-1) * q^(pW-128h)
 * (1+q^(p-1))^(L-D-a_T)
 * q^(g_H(S)-(p-beta)D).                            (2)
```

Its remaining score satisfies

```
g_H(S)-(p-beta)D
 =(1-p)L-pD-|V outside H|.                          (3)
```

In particular, the temporal-capacity parameter `beta` cancels exactly.
Summing (2) over supports while preserving their own cut scores is
the original joint field/rank sum, expressed differently. Changing
`beta` helps only if it enables a stronger classwise maximum bound,
count, or other inequality before this cancellation. It supplies no
independent annealed gain by itself.

This is not a prohibition on averaging rank and amplitudes jointly.
Such an average may be sharper than a uniform all-`H` event. It still
needs a theorem controlling the joint support, rank, and output
weights. Moving the row-set union into that same expression does not
prove the needed theorem.

## Eliminating the row set gives an exact but degenerate scalar cut

The capped rank is already the minimum over row-auxiliary choices:

```
R(V)=min_H [128|H|+|V outside H|].                   (4)
```

A network with source-to-state capacity `a`, forward death arcs of
capacity `b`, state/quartet footprint implications, capacity-`c`
footprint-to-row arcs, and row-to-sink capacity `128c` has deficit

```
Delta=max_S [aL-bD-cR(V)].                          (5)
```

Thus the row-set minimization can be eliminated exactly, without an
explicit union over `H`. This does not give a positive rank margin.
Every row has degree at most 256, so deterministically

```
min(128,degree)>=degree/2,
R(V)>=F/2>=L/2.                                    (6)
```

If `a/c<=1/2` and `b>=0`, (5) is identically zero for every setup:
all scores are nonpositive and the empty state set has score zero.
If `a/c>1/2`, the full state set has `L=N`, `D=0`, and `R=K`, giving

```
Delta>=aN-cK.
```

These are actual extremal cuts. They are not a numerical failure to
find a useful potential. Linear volume penalties expose the convex
lower hull through `(L,R)=(0,0)` and `(N,K)`. They cannot by themselves
certify an intermediate-volume rank profile strictly above the chord
`R=L/2`. Adding a nonnegative death penalty leaves both endpoint cuts
unchanged.

A constrained or nonlinear cut argument remains possible. The missing
information is the volume/death class or another statistic separating
nontrivial cuts from both endpoint phases, rather than the row-set
union alone.

## Exact rank gain as two geometric boundary costs

Let `d_r` be the footprint degree of outer row `r`. The elementary
identity

```
2*min(128,d_r)-d_r=min(d_r,256-d_r)
```

gives

```
2R-L=(F-L)+sum_r min(d_r,256-d_r).                   (7)
```

The first term counts footprint positions not occupied by the state.
The second is the number of footprint-label disagreements with the
best row-constant binary labeling, choosing each row independently.
Equivalently, if `H` is a minimizing row set in (4), then

```
2R-L=|S symmetric_difference H_slots|
     +2*|(V outside S) outside H|.                  (8)
```

Here `H_slots` contains the 256 routed positions of each row in `H`.
Equation (8) uses `S subset V`; it is an equality, not an independence
assumption. At a degree-128 tie either row choice gives the same value.

Consequently a small value of `2R-L` requires both limited quartet
footprint expansion and limited row-minority boundary. These two
requirements concern the same state support and must be controlled
jointly. Bounding their worst cases separately can select incompatible
supports and lose the desired rank gain.

## A precise nonlinear proof target

For the sampled routing and partitions, define the profile

```
I(l,d)=min_(S:L=l,D=d)
       [(F-L)+sum_r min(d_r,256-d_r)].               (9)
```

Then every support in that class satisfies
`R>=(l+I(l,d))/2`. A quantitative lower bound on (9), together with a
count or transfer for the same support classes, can be inserted into
the full field/output enumerator. It need not pass through a uniform
all-row-set deficit event.

This is a boundary problem on the state, quartet, and outer-row factor
graph. It retains the quartet union and the capacity-128 row cost.
A graph-edge relaxation that charges one covered input once for each
of its possible predecessors can lose a substantial constant factor;
an ordinary edge-expansion estimate is not automatically a bound of
the needed strength for (9).

The empty and full states have zero boundary in (7), so a useful
statement must depend on the population class or distinguish both
phases. No quantitative nonlinear expansion bound sufficient for the
remaining codeword classes is proved here. The exact identities
identify what a new global argument must retain, and explain why
another unconstrained scalar cut or a capacity reparameterization
cannot replace that argument.

## The fixed-row energy and saturation of cut witnesses

For an arbitrary designated row set `H`, write
`R_H=128h+|V outside H|`; this need not equal the actual capped rank.
Its energy is exactly

```
E_H=2R_H-L
   =|S symmetric_difference H_slots|
    +2*|(V outside S) outside H|.
```

In the half-row unit network `(a,b,c)=(1,1,2)`,

```
delta_H=256h-min_S(E_H+D),
tau_H=delta_H/2.
```

The empty and full states have energies `256h` and `N-256h`,
respectively, and both have zero internal deaths. At half rows both
energies are `N/2`. Thus a deficit measures an energy gain against two
constant phases. The spill term is a quartet-union cost, not just a
nearest-neighbor interaction. This interpretation does not make
different deviations independent or establish a contour probability
bound.

There is also an exact compression of unrestricted cut witnesses.
For a canonical footprint `V=V(S)`, define `A_t(V)` to contain the
positions in `V_t` whose next quartet is entirely in `V_(t+1)`;
at the free endpoint use `A_T(V)=V_T`. Then

```
S subset A(V),             V(A(V))=V.
```

The second equality follows in both directions: footprint monotonicity
gives one inclusion, and the defining next-quartet condition gives
the other. Consider the maximal coordinate runs of `A(V)`. Filling an
occupied run increases volume and cannot increase deaths. Filling a
previously empty closed run of length `ell` changes
`aL-bD-c|V outside H|` by `a*ell-b`; filling an empty open run changes
it by `a*ell`. Hence, whenever `a>=b`, every unrestricted cut maximum
has a representative

```
S=A(V).
```

This preserves the footprint and its actual eligible set exactly.
It removes the choices of which allowed runs are occupied from this
cut problem. It does not turn an arbitrary designated `H` into the
actual eligible set: that requires enlarging or changing the footprint.
It also need not preserve an upper volume bound, a death cap, or the
sufficient dimension-obstruction cutoff `R+D+1<=W`.

## Projecting temporal constants does not create row discrepancy

Let `B` be the `N`-by-8192 routed-row embedding: a row value is copied
to its 256 input positions. For every routing realization,

```
B^transpose B=256 I.
```

Let `P_time` be the orthogonal projection onto vectors constant within
each of the 4096 updates. Its rank is at most 4096. Therefore
`ker(P_time B)` has dimension at least 4096. Every vector in that
kernel is centered, since the sums in all updates vanish. On this
subspace the residual embedding preserves its full norm, giving

```
||(I-P_time)B||_2=16.
```

The same conclusion holds if only region-constant modes are removed.
Thus temporal-mode removal by itself cannot yield a small uniform
spectral norm from the row repetition count. A useful estimate must
include quartet/time energy or additional restrictions on the binary
state and row labels. This observation does not rule out such a joint
spectral or nonlinear boundary theorem; no adequate constant for one
has been proved here.

## The next discriminating test

The immediate next test is the existing adaptive comparison with its
exact physical population action. It should retain the population
through all sixteen updates before regional quota bounds are applied.
This removes the exponential-envelope or bin projection whose loss
dominated the completed controls. It tests the current majorant before
introducing another proof representation.

The exact action has a positive polynomial evaluation. Let `T` and `S`
be the hole and selected homogeneous quartics, and put `m=128`. For a
physical vector `f`, multiply each coefficient of `T^k` with old-color
degrees `(i,j,4k-i-j)` by `f(4k,i+j)`; denote this operation by
`T^k odot f_k`. Define

```
Tpow_0=1,       H_0=f(0,0),
Tpow_k=T*Tpow_(k-1),
H_k=S*H_(k-1)+binom(m,k)*(Tpow_k odot f_k).
```

Induction gives

```
H_m=sum_(k=0)^m binom(m,k)*(T^k odot f_k)*S^(m-k).
```

Extract the old-population coefficient and divide by its multinomial
once. The result is the exact action on `f`; the recurrence uses only
two positive quartic multiplications per stage. It imposes no
exponential shape on `f` and requires no dense population matrix.
The identity has passed an independent algebraic audit. Its full-size
regional bound remains an implementation and certification task.

If that test still fails, a precise remaining question is whether the
cut comparison can be restricted to canonical footprints whose
eligible set is the designated `H`, with the already certified
footprint/rank caps and a lower state-volume condition. The saturation
lemma preserves those footprint restrictions and lower volume. It
does not preserve a death cap or upper state volume. No probability
bound for this restricted family is asserted here.

The cancellation, cut identities, saturation, and row-isometry
statements are deterministic. Independent review of the complete note
passed, including those statements, their scopes, and the Horner
identity. No numerical jobs, helper implementation, or construction
changes were made for this note. The full-code certificate remains
open.
