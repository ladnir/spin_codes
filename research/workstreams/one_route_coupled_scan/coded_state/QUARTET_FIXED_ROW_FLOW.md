# A fixed-row flow certificate for every state mask

2026-10-31. Counting support histories has left large positive bounds.
A flow certificate can instead control every state mask simultaneously
for one fixed set of outer rows. The remaining probabilistic question is
whether the random routing and quartet partitions admit enough flow.
This note establishes that reduction and its window version. It does not
assert a flow-deficit tail bound or a whole-code certificate.

Use `M=8192`, `K=128M=1048576`, `n=512`, `T=4096`, `N=nT`, and
`W=943718`. Fix the route, all quartet partitions, and a row set `H` of
size `h`. Put `p=h/M`. A coordinate is an `H` input when its routed outer
row belongs to `H`. Each region contains exactly `h` such inputs, so the
whole encoder has exactly `pN` of them.

## The network and its exact cut identity

Create a vertex for every state coordinate `(t,i)` and a separate vertex
for every input coordinate. Add a source and sink, with these arcs:

* source to each state vertex, capacity `1-p`;
* each state vertex to its same-position input vertex, infinite capacity;
* state `(t,i)` to all four input vertices in its quartet at update
  `t+1`, infinite capacity, for `t<T`;
* state `(t,i)` to state `(t+1,i)`, capacity `p`, for `t<T`;
* each non-`H` input vertex to the sink, capacity one.

An `H` input has no sink capacity. There are no next-update or temporal
arcs after `T`. Total source capacity and total input-to-sink capacity
are both `(1-p)N`. A finite implementation may replace infinity by a
capacity greater than `(1-p)N`.

For a source-side state set `S`, use the standard statistics

```
L(S)=sum_t |S_t|,
D(S)=sum_(t=1)^T |S_(t-1)\S_t|,       S_0=empty,
V_t(S)=S_t union Gamma_t(S_(t-1)).
```

The infinite arcs force every input in `V(S)` onto the source side.
Choosing exactly those inputs minimizes the remaining input cost; extra
`H` inputs have no effect. The temporal arcs crossing the cut contribute
`pD(S)`. Hence the minimum cut has value

```
(1-p)N + min_S [|V(S)\H|+pD(S)-(1-p)L(S)].           (1)
```

Let `Flow_H` be the maximum flow and define its deficit

```
tau_H=(1-p)N-Flow_H.
```

The empty and full state masks both give zero in the bracket of (1).
Max-flow/min-cut therefore proves the exact identity

```
tau_H=max_S [(1-p)L(S)-pD(S)-|V(S)\H|] >= 0.          (2)
```

Thus any constructed feasible flow of value `(1-p)N-tau` proves, for
every state mask at once,

```
|V(S)\H|+pD(S) >= (1-p)L(S)-tau.                    (3)
```

This is a certificate about footprint geometry. It uses no field-map
amplitudes and does not change the encoder.

## A structural obstruction forces a large flow deficit

For any footprint, define its eligible row set by
`H={row i:|V intersect row i|>=129}`. Then its rank satisfies

```
R(V)=128h+|V\H|=pK+|V\H|.
```

Suppose a mask satisfies the dimension/Singleton obstruction

```
L>=R+1,                 R+D<=W-1.
```

For this eligible `H`, substitution in (2) gives

```
tau_H >= (1-p)L-pD-(R-pK)
      >= p(K-R-D)+(1-p)
      >= p(K-W)+1.                                  (4)
```

Consequently, for each fixed `h`, the simultaneous inequalities

```
tau_H < p(K-W)+1       for every H of size h           (5)
```

exclude every structural obstruction with that eligible-row count.
There is no separate union over `R,L,D,f,F` or state histories. A large
deficit does not conversely establish an obstruction: (2) is a scalar
relaxation of the two structural inequalities.

At `h=4096`, `p=1/2`, so the threshold in (4) is `52430`. Multiplying
all capacities by two gives integer capacities:

```
source demand per state:          1,
non-H input sink capacity:        2,
temporal capacity:                1.
```

If the integer maximum flow is `V_flow`, the integer deficit is
`N-V_flow=2*tau_H`. A structural obstruction therefore requires integer
deficit at least `104860`.

This fixed-`H` formulation does not repeat the unconstrained cut problem
whose empty/full convex hull defeated scalarization. The row set is fixed
before the flow is constructed. Empty and full masks give deficit zero;
a bad mask must improve the cut by the positive amount in (4).

## Independent windows and their boundary costs

Choose a window length `ell` divisible by 16 and dividing `T`. Partition
the updates into `B=T/ell` consecutive windows. In each window keep only
the network vertices and arcs entirely inside it. In particular, omit
next-input and temporal arcs leaving its last update. Retain the
same-position input arcs at that last update.

Each window contains whole routing regions. Its source and sink totals
are exactly `(1-p)n*ell`. Let its maximum-flow deficit be `tau_j`.
The union of the window flows is a feasible global flow, so

```
tau_H <= sum_(j=1)^B tau_j.                          (6)
```

Boundary losses are already included in `tau_j`; do not add them a
second time. For example, the last-layer sources at `H` inputs cannot
reach a sink inside the window. Therefore

```
tau_j >= (1-p)*(number of H inputs in its last layer).
```

For fixed `H`, that unavoidable component has expectation
`p(1-p)n`. At `p=1/2` it is 128, and its worst possible value is 256.
This is only the terminal-source component, not an upper bound on the
entire boundary or interior deficit.

The primary per-window thresholds obtained by equally dividing `52430`
are

| Window length | Number of windows | Deficit threshold per window |
| ---: | ---: | ---: |
| 16 | 256 | 204.8046875 |
| 32 | 128 | 409.609375 |
| 64 | 64 | 819.21875 |

These thresholds are in the original, unscaled capacities. They describe
the required average deficit; a proof need not bound every window by its
individual threshold. Larger windows reduce how often a boundary is
introduced and permit more temporal transport.

Condition on one fixed row set `H` before exposing the setup. Its `H`
labels are an independent uniform `h`-subset of the slots in each region.
Quartet partitions are also independent between windows. Hence the
restricted window deficits are independent and identically distributed.
For any `a>1`, (6) gives

```
Pr[tau_H >= p(K-W)+1]
 <= a^(-p(K-W)-1) * (E a^(tau_window))^B.             (7)
```

The expectation in (7) uses exact regional quotas. Replacing those quotas
by independent slot labels would require a separate comparison.
Symmetry and a single row-set union give

```
Pr[some obstruction with eligible count h]
 <= binom(M,h) * Pr_fixed_H[tau_H >= p(K-W)+1].        (8)
```

A whole theorem must also cover the different `h` values. In the primary
integer network, a window deficit `delta=2*tau_window` would supply the
certificate objective

```
lambda*104860 - B*log2(E 2^(lambda*delta))
                 - log2 binom(8192,4096),            (9)
```

for `lambda>0`. A sufficiently large positive value in (9) supplies the
desired failure exponent for that whole eligible-row class. No support
history count appears.

## An exact short matching for the half-row case

At `p=1/2`, a short matching represents the full flow exactly. Give each
state source one unit of demand and each non-`H` input capacity two.
Allow that source to use its own input, its next quartet's inputs, or the
quartet inputs two updates later. The last choice follows its temporal
arc once and then a next-quartet arc. No source uses two temporal arcs.

Every temporal arc is then used only by the source at its tail, so its
unit capacity is respected. Any fractional bipartite matching using
these neighbors is a feasible flow in the full network. Away from window
boundaries, both sides of this matching graph have degree at most nine.
Its deficiency equals the deficit of the full network, as follows.

Choose an integral maximum flow that minimizes the total flow on temporal
arcs. Suppose a state vertex has one unit on both its incoming and
outgoing temporal arcs, but no flow on its own source arc. Follow the
incoming unit back to its source. Cancel that path prefix and inject the
unit through the unused source arc at the current vertex. This preserves
the flow value and strictly reduces temporal flow, a contradiction.

Thus every vertex carrying one unit both in and out also uses its own
source arc. Conservation then gives one unit exiting directly to an
input. At each such vertex, pair the incoming temporal unit with that
input exit, and pair its own source unit with the outgoing temporal unit.
At a vertex without outgoing temporal flow, every incoming temporal unit
already exits to an input. These pairings decompose the flow into paths
that use at most one temporal arc.

A path with one temporal arc from `(t,i)` exits either to its own input
at `t+1` or to its quartet inputs at `t+2`. The former input already lies
in its direct quartet at `t+1`. Hence every path is a permitted matching
edge, and every used source and input keeps its load. The matching flow
therefore attains the full maximum flow. This proof needs unit source
and temporal capacities and integral input capacities. It does not need
randomness, regional quotas, or a boundary correction.

For future eligible-row classes with `p<=1/2`, the same path argument
has a fractional form. Write `q=1-p>=p`. A flow minimizing temporal use
must saturate a vertex's own source capacity `q` whenever it receives
positive temporal flow; otherwise an incoming prefix can be replaced.
Since outgoing temporal flow is at most `p<=q`, it can be supplied
entirely by that own source. All incoming temporal flow can exit to
inputs. Thus one temporal edge still suffices. The corresponding
matching must retain a per-source capacity `p` on its two-update choices,
in addition to total source capacity `q` and input capacity one. An
auxiliary vertex can enforce this constraint. It is redundant in the
half-row case after integer scaling. No probability bound for the
additional row counts is claimed here.

Hall deficits in this graph can be decomposed over connected components
of the source-neighborhood overlap graph. That graph has maximum degree
at most 72. This suggests counting small deficient components rather
than arbitrary state masks. Large components still require a separate
argument using the regional balance and temporal geometry.

For scale only, an interior source has nine distinct candidate inputs.
For fixed partitions, the probability that all nine are `H` inputs is at
most `p^9`, by sampling without replacement separately in each region.
At `p=1/2`, the expected number of such interior isolated sources is at most
`N/512=4096`. This does not upper-bound the total matching deficiency;
larger deficient components have not been controlled.

## A feasible backward algorithm and a failed moment comparison

For the half-row case, a simple backward algorithm constructs a feasible
integer flow without using temporal arcs. After serving the future
sources, each next-layer input has remaining capacity zero, one, or two.
Zero means that it is an `H` input. A fresh quartet has total remaining
capacity `C`. If its current sources include `r` at `H` inputs, serve
those from the future capacity first. They have no usable own input.
The algorithm loses `(r-C)_+` sources. It then serves

```
j=min(4-r,(C-r)_+)
```

non-`H` current sources from the future capacity. Their own inputs retain
capacity two. The other non-`H` sources use one unit of their own inputs,
leaving capacity one. Excess future capacity is discarded.

This allocation is feasible because every current source in the quartet
can use every one of its next-layer inputs. Those input quartets are
disjoint. At the terminal, every non-`H` source uses its own input and
every `H` source is lost. The total algorithmic loss upper-bounds the
optimal integer deficit.

Under independent half-density `H` labels, the next-layer capacities
depend on later labels and partitions, but not on the fresh partition
used for the current allocation. The current `H` labels are independent
as well. This permits a moment analysis before conditioning on the
regional quotas.

Let `a,b` count zero and two capacities in one next-layer quartet, so
`C=4-a+b`. For `t>=1` and positive `alpha,beta`, use potential
`alpha^(number of twos)*beta^(number of zeros)`. Its 15 local ratios are

```
A_(a,b) = sum_(r=0)^4 binom(4,r)/16
             * t^((r-C)_+) * alpha^min(4-r,(C-r)_+) * beta^r
             / (alpha^b*beta^a).                    (10)
```

Writing `Lambda=max A_(a,b)`, the resulting iid-label loss moment is at
most

```
Lambda^(128*(T-1))
 * {[(1+t*beta)/2]*max(1,alpha^-1,beta^-1)}^512.       (11)
```

The terminal factor counts the terminal `H` losses; the maximum bounds
the unrestricted remaining capacity vector at the other end. Conditioning
on exact half-region quotas and taking the row-set union would multiply
the Chernoff bound by

```
binom(M,M/2) * [binom(M,M/2)*2^(-M)]^(-256).
```

The pointwise maximum in (11) cannot certify the target for any marker
choice. Mix the input compositions `(a,b)=(4,0),(0,0),(0,4)` with weights
one half, one quarter, and one quarter. Their average old zero and two
counts are two and one. Their average new counts are also two and one,
so both potential drifts cancel. Their mean local loss is one: only the
all-zero-capacity quartet loses sources, with conditional mean two.
Jensen's inequality therefore gives

```
max_(a,b) log2 A_(a,b) >= log2 t.
```

The endpoint factor in braces in (11) is at least one. Since
`128*(T-1)=524160>104860`, the product of (11) with the target Chernoff
factor `t^-104860` cannot decay. No pointwise-marker search is warranted.
This proves a limitation of that comparison, not of partition averaging
or of the greedy algorithm itself.

## Verification and the first window outcome

[quartet_fixed_row_flow.py](quartet_fixed_row_flow.py) implements the
bounded integer network and the backward algorithm. Source review passed
the next-partition indexing, free endpoint, unit-source Dinic routine,
exact regional sampling, and independent reconstruction of the canonical
state-cut identity. The parent ran the exact checks serially: 18,432 state
histories, 72 minimum cuts, six terminal checks, five quota checks, and
76 comparisons with the feasible greedy flow passed.

For one width-512, 32-update window with seed zero and exact regional
quotas, the parent obtained

```
Maximum integer flow:           15926
Integer deficit:                  458
Physical deficit tau:             229
Allocated physical window share:  409.609375
Greedy integer loss:             1657
Greedy physical upper bound:      828.5
```

All cut identities passed. The returned state cut had `L=2814`, `D=32`,
`|V\H|=1162`, and `F=3053`, giving
`2814-32-2*1162=458`. This is a window cut, not a global dimension
obstruction. The optimal deficit has useful margin in this instance;
the simple greedy flow does not meet the allocated window share.

The parent then ran one 64-update window, also with seed zero and exact
regional quotas. It returned maximum integer flow `32301`, integer
deficit `467`, and physical deficit `233.5`, below the allocated share
`819.21875`. The independent cut identity passed with `L=18351`, `D=42`,
`F=18605`, and `|V\H|=8921`: `18351-42-2*8921=467`. The greedy loss was
`3021`, or `1510.5` in physical units, again above the allocated share.
The bounded run took 9.36 seconds.

The parent also solved the short matching on that same 64-update
instance. It returned exactly the full flow value `32301`, with deficit
`467`, in 8.846 seconds. Its assignments used 10,724 own-input edges,
15,701 next-quartet edges, and 5,876 two-update edges. Independent route
and input-load checks passed. The equality is explained by the theorem
above, rather than inferred from this diagnostic.

Read-only review of
[quartet_short_matching.py](quartet_short_matching.py) passed the two
partition indices, the Hall neighborhood
`N_t=S_t union Gamma_t(S_(t-1) union S_(t-2))`, the canonical cut, and
the conflict-free temporal embedding. That Hall cut has value
`n*T-|S|+2|N(S)\\H|`; it has no death term.

The two window instances use different draw sizes. Their identical seed
does not make them the same boundary experiment, and the similar optimal
deficits are not a measured decomposition into boundary and interior
costs. Both are structural diagnostics, not probability estimates.

No additional diagnostic is proposed here. The eventual proof needs a
certified window tail or a stronger constructive flow, such as a
matching argument with controlled short augmenting paths, followed by
the single row-set union (8). The pointwise comparison (11) is not a
proposed next search.
