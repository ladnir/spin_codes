# Deterministic pair mixing before the input diagonal

A cheap pair mixer can spread an isolated raw input across two inner-input
coordinates. This note gives a safe proof interface for that change.
The first conservative screen does not improve the three-quarter-rate
bound, so no implementation or certificate is promoted.

## Order of operations

Keep the existing GF256 MDS outer and regional route. Inside each
eight-symbol packet, apply four disjoint copies of a fixed 2-by-2 MDS
matrix `F`. Every entry and the determinant of `F` are nonzero. Apply the
existing independent nonzero input diagonal **after** these pair mixers.
Then use the unchanged repaired cubic eight-input, four-state inner,
including its fresh independent state diagonal.

For example, `F(x0,x1)=(x0+x1,x0+c*x1)` is MDS when `c` is neither zero
nor one. This mixer needs no new randomness. It does not increase the
source dimension or justify treating its outputs as independent inputs.

For an actual raw pair, the support transforms as follows:

| Raw support size | Mixed support size |
|---:|---:|
| 0 | 0 |
| 1 | exactly 2 |
| 2 | 1 or 2 |

The last row cannot be assigned an independent cancellation probability.
For example, raw source forms `(f,f)` are both live but their mixed first
coordinate is identically zero. The subsequent input diagonal cannot
remove that cancellation.

## Cover the mixed zero flats first

Fix the regional route, exact raw input support, past randomness, and the
common span of previously imposed equations. Work in the resulting source
quotient. Every prescribed active raw coordinate must be nonzero on the
witness. The current input and state diagonals remain unexposed.

Let `a` raw pairs have one active coordinate, and let `d` have two. The raw
packet weight is `w=a+2d`. A feasible witness cannot make either output of
a singleton pair zero. Nor can it make both outputs of a double pair zero.

For a witness, let `Z` be its exact set of zero mixed coordinates. Then
`Z` contains at most one coordinate from each double pair. Write `z=|Z|`
and let `k` be the rank of these mixed forms in the current quotient.
Thus `0<=k<=z<=d`. On this branch, the mixed support has size

`u=2a+2d-z`.

The span of the forms indexed by `Z` determines `Z`: every mixed form in
that span also vanishes on the witness. Choose a canonical independent
basis of this span. A rank-`k` basis contains at most one coordinate from
each double pair. There are at most `binom(d,k) 2^k` possible bases.
This also allows one rank-zero flat containing deterministic cancellations.

Impose the basis equations before revealing the current diagonals. They
increase equation rank by `k`, giving the potential factor `tau^-k`.
Require all other active mixed coordinates to remain live on the witness.
The current diagonals are still fresh and independent on this quotient.

The existing local matrix `T_u(y,tau)` applies uniformly to these mixed
forms, even when they depend on each other or on the incoming state.
If a flat kills a required raw form or the required nonzero incoming
state, that branch has no witness; retaining it only enlarges the bound.

Consequently, a valid two-state matrix for fixed pair counts is

`U_(a,d) = sum_(k=0)^d binom(d,k) 2^k tau^-k
            max_(k<=z<=d) T_(2a+2d-z)`.

The maximum is entrywise. The formula is a positive upper cover, not an
exact transition distribution. It counts each deterministic zero-flat
rank before applying the fresh-diagonal bound. No marginal probabilities
are multiplied as though those events were independent.

## Preserve the raw support markers

The number of raw supports with pair counts `(a,d)` is

`binom(4,d) binom(4-d,a) 2^a`.

The packet generating matrix therefore uses

`sum_(a+d<=4) binom(4,d) binom(4-d,a) 2^a x^(a+2d) U_(a,d)`.

The exponent of `x` is the **raw** weight. The outer shortened-space
dimension, active-row count, and regional placement normalization remain
those of the original raw source. In particular, doubling singleton
support does not double the outer source dimension.

The helper `weighted_matrices` returns a sum over placements at each raw
weight. It is not a pointwise bound for every support of that weight.
This sum can be used directly in the regional coefficient product, whose
positive expansion already sums those packet placements.

## First numerical gate

The screen uses the GF256 `[128,96]` outer, requested `K=2^20`, one regional
route, 10% target distance, and the repaired cubic 8/4 inner. Padding gives
`L=10944`, `N=1400832`, and requested rate `128/171`.

| Sampled class | Pair-premix margin |
|---|---:|
| 64 active outer rows | -2,322.50 bits |
| 128 active outer rows | -4,690.48 bits |
| Inverse raw-input density .10 | -78,465.44 bits |
| Inverse raw-input density .14 | -113,510.17 bits |

All four optimizations reported convergence. Without the premixer, the
same middle densities gave about -62,763 and -89,761 bits. Thus this
conservative cover is worse, not merely short of 40 bits. These are failed
upper-bound screens, not codeword counterexamples or global optimizer
optimality claims.

Singleton spreading improves birth output weight from one to two. However,
the reused return bound changes from the especially strong weight-one
minimum eight to the weight-two minimum six. The mixed-flat cover and its
weight maximization also add slack.

## A specific refinement, if revisited

The two-state cover forgets a useful distinction: a singleton raw input
creates two live mixed coordinates but only one source-form dimension.
A rank-cap transfer should therefore send a singleton birth to cap one,
not cap two.

More generally, after imposing a mixed zero flat of rank `k`, the mixed
input span has rank at most `w-k`. Suppose the surviving incoming state
forms have rank at most `h`, and `i` selected output-zero equations act on
mixed active coordinates. A safe successor cap can combine the previous
bound `h+u-i` with the source-span bound `h+w-k` and the state dimension.
This could retain singleton spreading without charging it as two source
dimensions. It needs its own continuation-cone audit; the present script
does not claim that improvement.

## Validation and scope

```powershell
python research/workstreams/char2_field_spin/theory/gf8_pair_premix_tests.py
python research/workstreams/char2_field_spin/theory/gf8_pair_premix_screen.py
python research/workstreams/char2_field_spin/theory/gf8_pair_premix_screen.py --mode sparse
```

Five tests pass. They check raw-support multiplicities, exact singleton
behavior, 768 GF4 mixed-flat profiles, deterministic rank-zero cancellation,
and agreement between the floating packet and an exact rational sum.
An independent audit checked another 256 four-pair GF4 profiles of rank
three, covering 615 exact flats and 1024 matrix-entry comparisons.
The audit also checked the quotient chronology and flat-count argument.

No native implementation, whole-code certificate, or frozen result was
changed. The bounded screen is complete; the rank-cap refinement is the
only identified next step for this particular proposal.
