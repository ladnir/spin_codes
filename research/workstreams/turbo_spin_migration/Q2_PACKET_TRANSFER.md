# Two active outer positions: packet-region transfer

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
packets reaching a region, their independent uniform slots become two
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
       [r(1-r) G,  r H+(1-r)^2 E]].
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

## Finite-step collisions

Within a region, two independently shuffled packet slots fall in the
same 64-bit step with probability `16/L`. On that event their combined
feedback may vanish, so the isolated-packet transition `P^2` cannot be
used. A finite-step upper transfer can charge the collision with a
matrix whose entries are at most one. Every entry of `K_2` is strictly
positive for a fixed positive `gamma`, and its smallest entries are
proportional to `1/M`. Thus, for fixed witnesses, the collision changes
the two-packet region matrix by a relative `O(M/L)` entrywise factor.
The packet-step output and one-step state transients admit the same
bounded-weight correction as in occupation one. Over
`R=O(log L)` regions the logarithmic correction is `o(1)`.

This is the route to a finite-step lemma, not yet a written full
certificate: the collision matrix, boundary-state handling, and uniform
entrywise constants must be made explicit before relying on it.

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
