# Occupation one at the paper's outer-growth schedule

> **Later result:** The [packet-window obstruction](PACKET_WINDOW_OBSTRUCTION.md)
> proves that the full `>11%` theorem is false for four-bit packets at
> this growth constant. The occupation-one certificate here remains
> valid for its stated class.

## Result and probability space

Use the growing Turbo member in [TURBO_FAMILY_BRIDGE.md](TURBO_FAMILY_BRIDGE.md):
one shared rate-half Golay--BAA outer of length `b`, `R=b/4` packet
regions, `L` outer positions, uniform packet and regional shuffles, and
the fixed s22 inner with independent uniform `GL(22,2)` updates. Let
`L` be a multiple of 16, `N=bL`, and take `b` to be a multiple of 24
with `b >= (39/4) log2 N`. Let `Z_{d,1}` count nonzero messages active in
exactly one outer position whose output weight is at most
`d=floor(0.11N)`.

The [packet-floor proof](PACKET_FLOOR_TAIL.md) gives an event `G_b`,
determined only by the two sampled BAA permutations, with
`Pr[G_b]=1-o(1)`. On `G_b`, every nonzero outer word has at least
`ceil(0.19R)` active four-bit packets. The certificates here establish

```
E[ Z_{d,1} 1[G_b] ] = o(1).
```

The expectation includes the shared BAA outer, route permutations,
and inner GL updates. Thus
`Pr[Z_{d,1}>0] <= Pr[not G_b]+E[Z_{d,1}1[G_b]]=o(1)`.
This is a proof for occupation one under the **ideal uniform setup law**.
It does not cover two or more active outer positions.

## Common transfer bound

Fix one nonzero message in one outer position. Let `v` be its active
packet count and `rho=v/R`. Packet labels do not affect the zero/live
state transition: every physical packet slot has rank-four feedback.
Write `M=2^22-1` and `p_0=2^21/M`. For `gamma>0`, let
`E=e^(-gamma)` and `J=(1-E)/gamma`. The continuum region matrices are

```
K_0 = [[1,0],[0,E]],
K_1 = [[0,J],[J/M,(1-1/M)E]].
```

For `y>0`, put `S(y,gamma)=K_0+yK_1` and let `lambda_S` be its Perron
root. A coefficient bound and the packet-subset denominator give the
inner contribution per outer bit

```
I(rho;y,gamma)
 = (ln lambda_S - rho ln y - h(rho)
    + gamma * 0.11/p_0)/4.                         (1)
```

The finite-step upper transfer in
[Q1_PACKET_TRANSFER.md](Q1_PACKET_TRANSFER.md) uses the exact s22
output-weight spectrum. It differs from the continuum matrices by
`o(1)` per region when `b=O(log L)`; its exceptional-step factor is
`exp(O(b/L))`. All witnesses below are drawn from finite sets, so the
correction is uniform and changes the normalized exponent by `o(1)`.
Positive 2x2 matrices have finite Perron prefactors for each fixed
witness. Binomial coefficients have uniform polynomial Stirling
remainders over the stated feasible ranges.

## Sparse Golay source messages

Split source messages by `j`, the number of active Golay blocks, with
`eta=1/50000`. For `1<=j<=eta b`, their number is at most
`C(b/24,j)4095^j`. Uniformly over this range, its normalized log is at
most

```
sigma = eta ln(4095e/(24eta)).
```

On `G_b`, `rho>=0.19`. Choose the fixed rational witnesses `y=6` and
`gamma=2`. For these witnesses, `I(rho;y,gamma)` is convex in `rho`, so
its maximum on `[0.19,1]` occurs at an endpoint. The outward-rounded
check `verify_packet_q1_sparse.py` gives

```
ln(2)/(39/4) + sigma + max_{rho in [0.19,1]} I(rho;6,2)
 < -0.02522.
```

The factor `L` for the active outer position costs at most
`exp(b ln(2)/(39/4)+o(b))`. The sparse-source first moment therefore
decays exponentially in `b`, up to polynomial factors. This argument
uses only the packet floor, not the detailed BAA packet spectrum.

## Dense Golay source messages

For `j>eta b`, let `alpha=a/b` be the Golay output weight fraction and
`beta=c/b` the first accumulator output weight fraction. Then
`alpha>=8eta=1/6250` and
`beta in [alpha/2,1-alpha/2]`. The paper's Golay and first-accumulator
exponents are `g(alpha)` and `pi(alpha,beta)`.

The exact final-accumulator packet matrix is

```
T(x,u) = [[1+u(6x^2+x^4), 4ux(1+x^2)],
          [x+u(3x+4x^3), u(1+6x^2+x^4)]].
```

For positive `x,u`, its coefficient bound contributes

```
O(beta,rho;x,u)
 = (ln lambda_T - rho ln u)/4 - beta ln x - h(beta),
```

where `lambda_T` is the Perron root. The dense-source exponent to
upper-bound over the feasible three-dimensional domain is

```
ln(2)/(39/4) + g(alpha) + pi(alpha,beta)
  + O(beta,rho;x,u) + I(rho;y,gamma).               (2)
```

`verify_packet_q1_interval.py` covers
`alpha in [1/6250,1]`, `beta in [0,1]`, and `rho in [0.19,1]`,
discarding only boxes disjoint from accumulator feasibility. Each
accepted box chooses fixed positive **rational** `x,u,y,gamma` and
upper-bounds (2) using outward-rounded interval arithmetic. Golay,
accumulator, and entropy bounds reuse the paper's interval verifier.
The completed run processes 1,663 boxes, accepts 1,664 leaves, reaches
depth 27, and has maximum accepted upper below `-9.10e-6` nats per
outer bit. The strict uniform gap absorbs polynomial coefficient
factors and `o(b)` finite-step corrections as `b` grows.

The binary64 saddle and finite-length screens were useful for choosing
the witnesses, but neither is used as the certificate. Source messages
with `j<=eta b` and `j>eta b` are both covered above. Selection of the
shared outer is charged once, before multiplying by `L`.

Reproduce the four proof gates from the repository root, in order:

```
python research/workstreams/turbo_spin_migration/verify_packet_tail_sparse.py
python research/workstreams/turbo_spin_migration/verify_packet_tail_interval.py --rho 19/100 --eta 1/50000
python research/workstreams/turbo_spin_migration/verify_packet_q1_sparse.py
python research/workstreams/turbo_spin_migration/verify_packet_q1_interval.py
```

## Remaining theorem work

The growing Turbo member now has an occupation-one asymptotic bound at
distance 11% and the existing `39/4` schedule. The full distance theorem
still needs occupation two, fixed higher occupation, growing sparse
occupation, and positive relative occupation. Those classes must retain
the shared-outer dependence and account for packet collisions and
feedback cancellations between different active positions.
