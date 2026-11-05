# A limitation of the coarse active-union gate

Improving regional route moments alone cannot close the complete bound
in [ACTIVE_UNION_GATE.md](ACTIVE_UNION_GATE.md). One small activity class
already makes its sum at least one, even with exact route probabilities.
This is a failure of that counting envelope, not a low-distance codeword.

## A guarded layout exists for every route

Use the current parameters

    K=2^20, N=2^21, B=64, k=32, p=4, b=128,
    M=N/B=32768, L=M/b=256,
    e=floor(0.35N)=734003.

Fix one source group of p outer rows and any 33 distinct regions. Consider
an arbitrary realization of the packet route. In each selected region,
the group's packet occupies one inner batch because p divides b.
Choose a two-batch interval inside that region containing this batch.
Such an interval always exists: use the next batch unless the packet
is in the last batch, in which case use the preceding batch.

These 33 intervals are disjoint, though intervals in adjacent regions
may abut. Regard them as a marked layout with one excursion per interval.
Closed excursions have length two; the last can be treated as free in
the upper bound. The layout parameters are

    c=33,       v=66,       u=bv=8448,       a=b(v-c)=4224.

Let U be its whole-batch union. Each of the four selected source rows
has a coordinate in U in all 33 selected regions. The MDS shortening
identity therefore gives

    D_out(U)>=p*(33-k)=4.                                 (1)

Thus every route has at least one layout in this fixed (v,c) class with
positive shortened dimension. The chosen layout may depend on the route;
the complete union over layouts includes it.

## The coarse conditional amplitude cap is already one

The active-union gate defines

    C=2^a q^((b-a)/2) V(u,e,c)/(q-1),
    V(u,e,c)=sum_(w=w0)^min(e,u) binom(u,w) q^(w/2),
    w0=(b+1)(c-1)+1.

Here w0=4129 and u=8448<e. The term w=u is included, so

    V(u,e,c)>=q^(u/2),
    C>=2^a q^((b+u-a)/2)/(q-1)
      =2^a q^(b(c+1)/2)/(q-1)>1.                         (2)

For every route, the gate's conditional bound for a fixed layout is

    1[D_out(U)>0] min(1,C q^D_out(U)).

By (2), this is exactly 1[D_out(U)>0]. Consequently the strongest result
obtainable from that conditional bound, even before introducing any
moment approximation, is the route probability of its guard.

Let L_(66,33) denote all marked layouts in this class. Equation (1) gives
the pointwise inequality

    sum_(U in L_(66,33)) 1[D_out(U)>0] >= 1.

Taking expectation proves

    sum_(U in L_(66,33)) Pr_route[D_out(U)>0] >= 1.        (3)

The sum may count several marked layouts with the same union; this only
increases its value. Exact regional occupancies, exact guarded moments,
or a sharper count of these layouts cannot reduce the sum in (3) below
one. Fractionalizing the same conditional cap cannot improve it either.

## What this rules out, and what remains useful

The construction does not supply an output realizing the marked layout.
In particular, D_out(U)>0 does not imply that the inner can keep a word
inside U while meeting all required resets. The coarse amplitude bound
discarded too much of that compatibility information.

This obstruction differs from [BAD_ROUTE_OBSTRUCTION.md](BAD_ROUTE_OBSTRUCTION.md).
That note constructs actual low-weight words for specially aligned routes.
Here every route supplies an outer shortening witness, but no compatible
inner word is asserted.

The coarse gate can still certify particular ranges of v,c. It cannot
serve as the sole complete gate. Do not continue a regional-only operator
implementation on the assumption that sharper route moments repair it.

A remaining interface retains outer compatibility before aggregating masks.
For exact output support W and inverse-randomness support A, put
S=A union W, with S contained in U. The overlap-rank bound in
[AFFINE_COMPATIBILITY.md](AFFINE_COMPATIBILITY.md) gives

    expected eligible outputs
      <= q^(|W|-R_overlap-rank P_S)
       = q^(D_out(S)-|A minus W|-R_overlap)
      <= q^(D_out(U)-|A minus W|-R_overlap).               (4)

Equation (4) still shares the same D_out(U) across masks, but retains
deterministic ranks that the coarse volume discarded. A complete sum of
its weighted masks remains open; no closure is claimed here. The next
step is to bound that sum or preserve equivalent joint constraints,
not to refine the guard in isolation.

This is an analytic obstruction and requires no numerical search or
benchmark. No execution or whole-code distance claim accompanies it.
