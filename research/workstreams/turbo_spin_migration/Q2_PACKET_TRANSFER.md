# Two active outer positions: packet-region transfer

> **Later result:** The [packet-window obstruction](PACKET_WINDOW_OBSTRUCTION.md)
> rules out the full `>11%` theorem for this four-bit growing family at
> the paper's `39/4` schedule, despite the closed Q2 and Q3 bounds below.

## Why this class differs

The [occupation-one proof](Q1_SCHEDULE_CERTIFICATE.md) uses one active
four-bit packet at most in each region. Two active outer positions can
send two packets to one region. They usually occupy distinct physical
steps; if they share a step, their feedback syndromes can cancel. The
same shared BAA outer also supplies both local words, so its expected
packet spectrum cannot simply be squared.

## Limiting two-packet kernel

Keep `M=2^22-1`, `r=1/M`, `E=e^(-gamma)`, and
`J=(1-E)/gamma`. Let `P=[[0,1],[r,1-r]]` be the zero/live transition at
one isolated active packet, and `D(t)=diag(1,e^(-gamma t))` the Laplace
factor for an empty fraction `t` of a region. Conditional on two
packets reaching a region, their distinct uniform slots become two
independent uniform points as `L->infinity`. In order-statistic form,

```
K_2 = 2 integral_{0<u<v<1} D(u) P D(v-u) P D(1-v) du dv.
```

Define

```
H = 2(gamma-1+E)/gamma^2,
G = 2(1-(1+gamma)E)/gamma^2.
```

Direct integration gives

```
K_2 = [[r H,       (1-r) G],
       [r(1-r) G,  r G+(1-r)^2 E]].
```

At `gamma=0` by continuity, `H=G=1` and `K_2=P^2`.
Together with the one-packet matrices `K_0,K_1`, this determines the
continuum two-block transfer. If the two outer words have packet counts
`v_1,v_2`, their independently shuffled packet supports give

```
F_{R,v_1,v_2}(gamma)
 = [x^v_1 y^v_2] (1,0)
   (K_0+(x+y)K_1+xy K_2)^R (1,1)^T
   / (C(R,v_1) C(R,v_2)).
```

This formula conditions on the two realized packet counts. Packet
labels do not enter its limiting transition because every isolated
nonzero packet has nonzero feedback in every physical slot.

## Finite-step collision comparison

Write `L=16n` and use the exact s22 output spectrum to define the
empty-step Laplace factor `q(theta)` in the occupation-one note. Take
`theta=gamma/(4Lp_0)`, `D_n=diag(1,q(theta))`, and
`P=[[0,1],[1/M,1-1/M]]`. An active packet at physical step `i`, counted
from zero, has the proxy transfer `D_n^i (D_n P) D_n^(n-i-1)`.
Conditional on two packets reaching a region, the common regional
permutation puts them in two distinct uniform slots among the `L`
packet slots. Their step pair is sampled without replacement. For
every nonnegative transfer entry, its expectation under this law is
at most `L/(L-1)` times the expectation under two independent uniform
slots. Under the latter law the physical steps are independent uniform
elements of `[n]`. If the steps differ, the proxy
transfer is the ordered product with two `D_n P` factors. If they
coincide, their combined feedback may be zero, including complete
cancellation. Upper-bound that step's zero/live transition by the
all-ones matrix `U`; the colliding proxy transfer is
`D_n^i (D_n U) D_n^(n-i-1)`.

Let `K_{2,n}^+` be the average of these `n^2` transfers, with the
factor `L/(L-1)` for the actual slot law. The collision
part has coefficient exactly `1/n=16/L`. It is entrywise at most
`1/n` times an all-ones matrix because every `D_n` entry is at most
one. For a fixed `gamma>0`, `q(theta)=1-gamma/n+O(n^-2)` from the exact
finite output spectrum. The noncollision Riemann sum therefore differs
entrywise from the integral defining `K_2` by `O_gamma(1/n)`. All four
entries of `K_2` are positive for fixed `gamma>0` and fixed `M`; hence

```
K_{2,n}^+ <= (1+C_gamma,M/n) K_2       entrywise.
```

The same comparison holds for `K_{0,n}` and `K_{1,n}`, preserving their
exact zero entries. Mark every active physical step and its immediate
successor. There are at most four marked steps per region, including a
collision. Every unmarked live pre-state is uniform nonzero after the
intervening uniform `GL(22,2)` update, conditional on its zero/live
history. Replacing the true weighted output on marked steps by the
proxy weight changes the upper moment by at most `exp(64 theta)` per
marked step. For two outer words with packet counts `v_1,v_2`, the full
finite moment is consequently bounded by the coefficient transfer
with `K_{a,n}^+`, multiplied by
`exp(128 theta(v_1+v_2))`. Since `R=O(log L)`, both the matrix comparison
and exceptional-step factor contribute `exp(O(b/L))` for each fixed
set of positive witnesses. This gives the required continuum exponent
in the `L->infinity` limit. The collision event is included in the
upper bound; no noncancellation assumption is made there.

## Shared-outer selection and first screen

For a realized shared outer, let `B_O(v)` count its nonzero words with
`v` active packets. Besides the proved 19% support floor, select the
event

```
B_O(v) <= b^2 E[B_O(v)]  for every v.
```

Markov and a union bound over `R+1=O(b)` packet counts charge this event
once with probability `O(1/b)`. On it, products of two realized packet
spectra can be bounded by `b^4 E[B_O(v_1)]E[B_O(v_2)]`. The factor is
polynomial and does not affect a strict asymptotic exponent. This step
preserves the fact that one outer realization is reused in all positions.

`q2_inner_exponent_screen.py` minimizes a Perron coefficient bound for
the continuum kernel. At `v_1/R=v_2/R=0.19`, distance `0.11`, and the
paper's `39/4` schedule, **two position choices plus the inner** have
exponent about `-0.05636` nats per outer bit, before paying for the
outer words. A second screen folds the exact finite BAA packet spectrum
into both outer words. On sampled packet-count pairs, the worst
normalized exponent is about `-0.0932` at `b=432` and `-0.0836` at
`b=864`, both at the packet-floor pair. The latter screen is in
`q2_spectrum_screen.py`. These are binary64 diagnostics, not proofs for
all packet-count pairs. The next gate is a uniform
outer-spectrum/inner exponent certificate over both packet fractions,
followed by the finite-step collision lemma.

An initial six-variable outward-rounded branch verifier is in
`verify_packet_q2_interval.py`. It combines two copies of the Golay,
first-accumulator, and packet-spectrum exponent with the two-packet
inner transfer. After 10,000 boxes its largest unresolved upper is
`+1.37` nats per outer bit, despite the negative sampled objective.
The broad boxes lose cancellations between each outer's weight and
packet fraction. This is a **bound-resolution problem**, not evidence
of a bad code. The next approach should first certify a one-dimensional
packet-spectrum envelope for the selected BAA outer and then verify
the remaining two packet fractions jointly. The direct verifier is
preserved as a control and diagnostic.

## Closed occupation-two exponent

`verify_packet_spectrum_envelope.py` partitions the allowed packet
fraction `[0.19,1]` into 32 rational bands. For each band it branches
over the Golay source-weight and first-accumulator-weight fractions.
Every leaf evaluates `g(alpha)+pi(alpha,beta)` and the final packet
Chernoff transfer with fixed positive rational witnesses, using the
paper's outward-rounded interval bounds. Taking the largest live leaf
upper after 1,000 splits is valid regardless of whether the refinement
budget has converged. The generated `packet_envelope_32.json` records
the 32 certified upper bounds on
`b^-1 ln E[B_O(v)]` for every `v/R` in each band. Polynomial coefficient
and summation factors vanish after division by `b`.

`verify_packet_q2_envelope.py` then checks every unordered pair of
bands, 528 pairs in total. It uses a positive rational witness for
each pair in the two-packet Perron coefficient bound, with outward
interval arithmetic over both packet fractions. Every pair closes
without subdivision. The largest certified upper for the complete
Q2 exponent, including two active-position choices, both outer word
counts, and the 11% low-output Chernoff charge, is
`-0.0320290481` nats per outer bit. Thus the selected ideal setup law
satisfies `E[Z_{floor(0.11N),2} 1[G_b]]=o(1)` under the same
`b >= (39/4) log2 N` schedule as Q1. This is an asymptotic result;
the one-time selected-outer failure remains `o(1)`.

The first derivation of `K_2` had `rH` in its live-to-live entry. Direct
integration shows the correct term is `rG`. An independent numerical
quadrature of all four entries at `gamma=2` agrees with the corrected
matrix within `3e-23`; the certificate was rerun after this correction.

Reproduce the two exponent gates from the repository root:

```
python research/workstreams/turbo_spin_migration/verify_packet_spectrum_envelope.py --bands 32 --max-boxes-per-band 1000 --output research/workstreams/turbo_spin_migration/packet_envelope_32.json
python research/workstreams/turbo_spin_migration/verify_packet_q2_envelope.py research/workstreams/turbo_spin_migration/packet_envelope_32.json
```

Occupation two is no longer the exponent obstruction. The full
minimum-distance theorem still requires fixed `Q>=3`, growing sparse
occupation, and positive occupation. The finite-step comparison above
must also receive a separate proof audit before paper migration.

## Occupation three

The same continuum method closes `Q=3`. For three packets in one
region, let `k` be the number of live spacings among the four ordered
spacings. Their total live-time fraction has the `Beta(k,4-k)` law,
with endpoint interpretations for `k=0,4`. Enumerating the four-state
paths yields a positive two-by-two `K_3`; `verify_packet_q3_envelope.py`
evaluates its closed Laplace integrals with outward interval arithmetic.
For packet fugacities `u_1,u_2,u_3`, the region polynomial is

```
T_3 = K_0 + (u_1+u_2+u_3)K_1
          + (u_1u_2+u_1u_3+u_2u_3)K_2 + u_1u_2u_3 K_3.
```

The checker combines this matrix with the same 32-band selected-outer
packet envelope and enumerates all 5,984 unordered triples of bands.
Its largest outward upper for the complete three-position exponent is
`-0.0480973770` nats per outer bit. The fixed-Q collision comparison
extends to three positions: sampling distinct physical slots differs
from independent slots by a factor `L^3/(L)_3=1+O(1/L)`, and step
collisions have probability `O(1/L)`. The `K_3` entries are positive,
so an all-ones collision transfer incurs an entrywise `1+O(1/L)`
factor. Together with the marked-step correction this is `exp(o(1))`
over `R=O(log L)` regions. Thus occupation three has the same
asymptotic first-moment conclusion as occupation two.

```
python research/workstreams/turbo_spin_migration/verify_packet_q3_envelope.py research/workstreams/turbo_spin_migration/packet_envelope_32.json
```

The earlier product-norm attempt for all fixed `Q>=3` does not close:
it replaces the exact multi-packet kernels by a coarse per-row factor.
`fixed_q_product_screen.py` records a positive worst-band diagnostic
even though the exact `Q=3` transfer closes. It should not be used as
evidence against the construction. Fixed `Q>=4` and both growing
occupation regimes remain open.
