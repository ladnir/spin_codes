# General row-count flow, adaptive cover, and field bridge

2026-11-01. The half-row flow interface extends to every eligible-row
count using the same two population statistics. This note proves that
extension and specifies its endpoints and probability normalization.
It provides an interface for the remaining classes, not their required
flow-tail estimates. The fixed encoder and its terminal convention do
not change.

Use `n=512`, `T=4096`, `N=nT=2097152`, `M=8192`, `K=N/2`,
`W=943718`, and `q=2^32`. There are 256 regions of 16 updates each.
For a state support, use volume `L`, internal deaths `D`, terminal
population `a_T`, footprint `V`, and rank envelope
`R=sum_rows min(128,degree_V)`. Its eligible row set consists of the
rows of degree at least 129. If this set is `H`, with `|H|=h`, then

```
R=128h+|V outside H|.                                  (1)
```

The local input labels `H_t` are the positions routed to these rows;
write `J_t=H_t^c`. Each region has exactly `h` labels in `H` and
`M-h` labels in `J`. All field statements use the accepted-field
density loss `gamma<=1024/989`.

## A network with an exact short matching

Fix integers

```
a>=1,       0<=b<=a,       c>=a+b.
```

Give each source-to-state edge capacity `a`, each same-coordinate
forward temporal edge capacity `b`, and each non-`H` input-to-sink
edge capacity `c`. A state has infinite-capacity edges to its own input
and to the four inputs in its next quartet. There are no temporal
edges after time `T`. Its integer deficit is exactly

```
delta_H=aN-maxflow
       =max_S [aL-bD-c|V(S) outside H|].               (2)
```

Indeed a finite cut with state set `S` has capacity
`a(N-L)+bD+c|V(S) outside H|`. Input vertices outside the forced
footprint can be removed from the cut without increasing its capacity.

This maximum flow equals a matching with `a` unit copies of each
source: `a-b` copies may use their own input or the next quartet;
`b` copies may additionally use the quartet two updates ahead.
Each non-`H` input has capacity `c`; neighbors outside the time range
are omitted.

To prove equality, choose an integral maximum flow minimizing total
temporal flow. If a state receives positive temporal flow while its
own source edge is unsaturated, replace an incoming path prefix by
unused own-source capacity. This preserves the flow value and strictly
reduces temporal flow. Therefore any state with incoming temporal flow
uses all `a` units of its own source edge. Since its outgoing temporal
flow is at most `b<=a`, pair all outgoing temporal units with its own
source units, and pair all incoming temporal units with input exits.
Every resulting path uses at most one temporal edge. A carried path
exiting at the next own input already has a direct next-quartet route;
the other carried paths use the two-ahead quartet. At most `b` own
units carry, so they can be assigned to the long copies. Conversely,
every short matching embeds in the network: only long copies use a
temporal edge, with at most `b` copies sharing that edge.

This proof uses integral capacities and `b<=a`. When an intended
anchor has temporal capacity larger than the source capacity, one may
restrict temporal capacity to `b=a`. The resulting deficit upper-bounds
the unrestricted deficit and still gives (2), but the restriction may
make the probability bound weaker.

## Rank margin and the low-death field theorem

Put

```
p=1-a/c,       tau_H=delta_H/c.
```

Here `p` is a chosen network anchor; it need not equal `h/M`.
Equations (1)-(2) imply, for every support with eligible set `H`,

```
R >= (a/c)L-(b/c)D+128h-tau_H.                       (3)
```

Suppose the geometry satisfies, simultaneously for every row set of
size `h`,

```
tau_H <= 128h-pW-s.                                 (4)
```

The fixed-support death-cover field bound, including output
cancellations and the free terminal state, is

```
gamma/(q-1) * q^-R * x^-W * (1+qx)^L
 * [qx^2/(1+qx)]^D * [qx/(1+qx)]^a_T.
```

Set `x=q^-p` and `A=1+q^(p-1)`. Substitution of (3)-(4) gives

```
gamma/(q-1) * q^-s * A^(L-D-a_T)
 * q^((a+b-c)D/c)
 <= gamma/(q-1) * q^-s * A^N.                        (5)
```

The last inequality uses `c>=a+b`. It requires no condition
`L+D<=W`. Therefore, for `D<=d`, the conditional field failure is at
most

```
gamma/(q-1) * q^-s * A^N
 * sum_(i=0)^min(N,2d+n) binom(N,i).                 (6)
```

The support census uses its `2D+a_T` binary changes from the empty
initial state. There is no further row-set union inside (6), since
(4) already holds for every eligible set of this size. See
[QUARTET_FLOW_FIELD_BRIDGE.md](QUARTET_FLOW_FIELD_BRIDGE.md) for the
underlying fixed-support bound and the exact half-row constants.

For `p<=1/2`, the common estimate `A^N<2^47` applies. It does not
apply when `p>1/2`; in that range retain the actual factor `A^N`.
The extra death factor in (5) can also be retained. With
`eta=(c-a-b)/c>=0`, its exact unrestricted support sum is

```
[sum_d binom(T+1,2d+1) q^(-eta*d)]^n,
```

and a low-death restriction can be applied to the corresponding
polynomial coefficient sum. This is an optional sharper census,
not a claimed closure.

For comparison, the sufficient dimension-obstruction conditions
`L>R` and `R+D+1<=W` force

```
delta_H >= 128ch-(c-a)W+c.                           (7)
```

To see this, use `L>=R+1`, `c-a>=b`, and `R+D<=W-1` in (2).
Excluding (7) alone does not prove distance. The stronger margin (4)
is what supplies the complete low-death amplitude bound (6).

## Adaptive normalization for arbitrary capacities and row counts

For the short matching, the dual input score is
`c|C|-|N^-(C)|`, where source copies count separately. The left
deficit equals its maximum plus `aN-c|J|`. A maximum can be chosen
with `C_1=J_1` and, for each later update,

```
C_t=U_t intersect J_t,
U_t is a union of current quartets meeting J_t.
```

Adding a missing `J` input in an already selected quartet gains
capacity `c` and introduces at most `a` own source copies; all earlier
neighbors are already present. Thus saturation cannot decrease the
score. The same argument justifies `C_1=J_1`.

Write `V_t=U_t^c`. In reverse time, the known future is

```
A_t=V_(t+1),       W_t=V_(t+1) intersect V_(t+2).
```

Both formal future hole sets are full. For a current quartet define

```
u=|A_t intersect Q|,        w=|W_t intersect Q|,
j=|J_t intersect Q|,
r_A=|J_t intersect A_t intersect Q|,
r_W=|J_t intersect W_t intersect Q|,
B=cj-(a-b)r_A-b r_W,
A_loss=max(0,4a+bu-B).                              (8)
```

The immediate gain from selection is exactly `B>=0`. With every
remaining earlier choice fixed, selecting the quartet adds at most
`4(a-b)+bu` predecessor copies at time `t-1` and `4b` at time `t-2`.
Consequently the difference of the two complete scores, and also of
their maxima over the remaining common choice space, belongs to
`[B-4a-bu,B]`, hence to `[-A_loss,B]`. Missing early time layers only
improve the lower endpoint.

For `y>1` and positive `A_loss,B`, set

```
w0=(1-y^-A_loss)/(1-y^(-A_loss-B)),
w1=(1-y^-B)/(1-y^(-A_loss-B)),
g=y^-B*w0,       p_select=w1.                        (9)
```

The weighted-maximum inequality gives an upper bound on each branch
maximum, and `g+p_select=1`. If `A_loss=0`, force selection. If
`B=0`, force a hole. A quartet with `j=0` is always a hole. These
degenerate rules are valid because the respective branch difference
has a fixed sign. The induction over reverse-time decisions is the
one in
[QUARTET_ADAPTIVE_HALL_NORMALIZATION.md](QUARTET_ADAPTIVE_HALL_NORMALIZATION.md).

The baseline current source coverage is
`(a-b)|A_t^c|+b|W_t^c|`. Factoring the selected score out of each
weighted branch leaves the conditional probabilities (9). Its common
exponent at time `t` is

```
c|J_t|-an+(a-b)|H_t intersect A_t|+b|H_t intersect W_t|.
```

The forced first-layer witness has the same expression. The dual
offset `aN-c|J|` cancels the sum of its first two terms, for every
row count. Thus, for each fixed routing and partition realization,

```
y^delta_H <= E_adaptive coins[
 y^(sum_t ((a-b)|H_t intersect A_t|+b|H_t intersect W_t|))].  (10)
```

This is a positive moment majorant. A sampled cost need not be a
feasible flow or an upper bound on the actual deficit.

## Population kernel, endpoints, and one row-set union

The same population state suffices:
`(ell,omega)=(|A_t|,|W_t|)`, with `ell` divisible by four and
`0<=omega<=ell<=n`. At width 512 there are 33,153 safe states.
Let `x` mark non-`H` inputs. For an old quartet type with `w` positions
in `W`, `e` in `A outside W`, and `z` outside `A`, sum over
`r,s,t` non-`H` labels in those three parts. Put
`u=w+e`, `j=r+s+t`, `r_A=r+s`, and `r_W=r` in (8)-(9).
Its coefficient has multiplicity

```
binom(w,r) binom(e,s) binom(z,t)
```

and weight

```
x^(r+s+t) y^((a-b)(u-r-s)+b(w-r))
 * [p_select*(new_Z)^4 + g*(new_W)^u*(new_E)^z].      (11)
```

Multiply (11) by
`multinomial(4;w,e,z) old_W^w old_E^e old_Z^z`, sum the 15 old
quartet types, raise to `n/4`, and divide extracted coefficients by
the old population multinomial. This gives the exact positive
partition-averaged operator `K_x`. Its physical initial state is
`(n,n)`. The first-layer endpoint is

```
q_x(ell,omega)=(1+x)^(n-ell)
 *(y^(a-b)+x)^(ell-omega)*(y^a+x)^omega.              (12)
```

There is no additional factor `y^(an)`: the full-hole start already
charges the final source layer.

For fixed `h`, put `k=M-h` and `C_h=binom(M,h)`. Define

```
B_h=[x^k]K_x^16/C_h,
b_h=[x^k]K_x^15 q_x/C_h.
```

Independent regional routing and fresh partitions give the fixed-row
set moment bound `e_(n,n) B_h^255 b_h`. Markov's inequality and one
union over row sets give, for every threshold `d_star`,

```
Pr_geometry[exists H, |H|=h, delta_H>=d_star]
 <= C_h*y^-d_star * e_(n,n) B_h^255 b_h.             (13)
```

For (4), one may use
`d_star=floor(c*(128h-pW-s))+1`; any smaller positive threshold
is conservative. The disjoint-bin comparison and the order of
regional Cauchy minima remain those in
[QUARTET_ADAPTIVE_POPULATION_COMPARISON.md](QUARTET_ADAPTIVE_POPULATION_COMPARISON.md).
For fixed `(a,b,c)`, the evaluated population actions do not depend on
`h`: only the quota, normalization, and tail threshold change. Changing
the capacities changes the adaptive probabilities and requires new
actions.

## Remaining scope

The half-row choice `(a,b,c)=(1,1,2)` recovers the current operator
and the deficit target 60518 for the `D<=65536` field bridge. No
required global flow-tail estimate has been proved, including that
primary target. Other anchors are necessary to address all row counts;
for `p>1/2`, taking `b=a` is a restricted network whose sharpness must
be established separately.

The total input capacity gives a useful necessary check before any
tail calculation:

```
tau_H >= max(0,256h-pN).
```

For the half anchor and the simple census margin `s=22170`, this and
(4) require `3860<=h<=4332`. Even allowing arbitrary anchors, the
largest margin permitted by this capacity check is
`h*(K-W)/M`, attained at `p=h/M`. The same unweighted census therefore
requires `h>=1733`, even with an ideal flow. Smaller row counts need a
smaller death census, the retained weighted-death factor, or additional
support compatibility. These are limitations of this margin/census
interface, not obstructions to the construction.

The existing high-rank class covers `R>=1009256`, hence eligible
counts `h>=7885`; the large-footprint class has its separate bound.
These facts do not cover all larger-death histories or all remaining
`h<7885`. Completing the fixed-construction goal requires an explicit
accounting of those remaining row/death classes, their geometry-tail
allowances, and their conditional field bounds. The present note
provides that common interface without claiming the required estimates.

Next work: audit the bounded primary population comparison, distinguish
its bin/tilt losses from the underlying adaptive moment, and obtain a
rigorous primary tail before extending the same interface across the
remaining row and death ranges. No code or numerical jobs were added
for this generalization.

Independent full-note review passed the integral path compression,
arbitrary-row-count dual offset and adaptive induction, polynomial and
endpoint, field residual factor, structural threshold, and capacity
range calculations. The displayed threshold 60518 is conservative for
the exact cap 60518, whose strict failure event starts at 60519. This
review establishes the interface, not any missing probability estimate.
