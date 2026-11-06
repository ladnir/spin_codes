# One active block in a growing Turbo family

## Question and construction

The finite Turbo encoder uses a 256-bit RS8 outer, a four-bit packet route,
and a width-64, 22-bit-state inner. A fixed outer length cannot give the
claimed asymptotic distance under this route and inner
([pair-return obstruction](PAIR_RETURN_OBSTRUCTION.md)). The proposed growing
member keeps the route and inner and replaces the outer with the paper's
Golay--BAA constituent of length `b`, where `24 | b` and `4 | b`. This note
tests the first proof class: messages active in exactly one outer block.

Write `R=b/4`, and let `L` be the number of outer blocks. The selected BAA
outer is sampled once and reused in all `L` positions. For each encoded
block, a uniform packet shuffle sends its `R` four-bit packets to distinct
regions; each region independently shuffles the `L` block slots. The inner
processes the resulting 64-bit steps with its fixed maps and independent
uniform `GL(22,2)` state updates. The map check in
`verify_pair_return_map.py` establishes that every nonzero label in each of
the sixteen physical packet slots has nonzero feedback. Therefore a single
active packet always turns a zero state on, whatever its label.

## Exact outer packet enumerator

Let `B_O(v)` count nonzero messages whose output under the *realized* BAA
outer `O` has exactly `v` nonzero four-bit packets. We need packet count,
not just binary weight: four nonzero bits may occupy one packet or four.

The Golay direct sum has weight counts `G_b(a)=[z^a]G(z)^(b/24)`, with
`G(z)=1+759z^8+2576z^12+759z^16+z^24`. Let `P_b(a,c)` be the paper's
exact accumulator weight transition. After the second BAA permutation, an
input of weight `c` is uniform over the `C(b,c)` binary words of that
weight. The final accumulator output has `c` transitions from initial state
zero. Define `D_b(c,v)` as the number of length-`b` binary output words
with exactly `c` transitions and `v` nonzero four-bit packets.

`D_b` has a two-state packet transfer. Its state is the preceding output
bit `p`. For a four-bit output packet `q=(q_1,...,q_4)`, set

```
t(p,q) = (q_1 xor p) + (q_2 xor q_1) + (q_3 xor q_2) + (q_4 xor q_3),
a(q)   = 1[q != 0000],
next   = q_4.
```

Let `T(x,u)_{p,p'}` sum `x^t(p,q) u^a(q)` over all packets `q` ending in
`p'`. Then

```
D_b(c,v) = [x^c u^v] (1,0) T(x,u)^R (1,1)^T,
E[B_O(v)] = sum_{a,c} G_b(a) P_b(a,c) D_b(c,v)/C(b,c).
```

All counts in `D_b` are exact integers. In particular,
`sum_v D_b(c,v)=C(b,c)`. The implementation in
`packet_baa_q1_screen.py` checks this identity at every `c` and directly
enumerates all 8-bit inputs for a small independent check. This enumerator
uses the actual BAA distribution; it adds no mixer to the encoder.

## One-block inner transfer

Consider `L -> infinity` with `b` fixed. The active block's packet slots
in distinct regions become isolated 64-bit steps. A live state is uniform
over the `M=2^22-1` nonzero states after an intervening empty step. At an
active packet, a live state returns to zero with probability exactly `1/M`.
The mean bit density emitted by a live state at an empty step is
`p_0=2^21/M`; every row of the chosen expansion map is nonzero. Packet
steps themselves occupy `O(b)` bits and vanish from the relative output
weight as `L` grows. The nonuniform transient immediately after a packet
also vanishes on this scale.

For a Laplace parameter `gamma>0`, put `E=e^(-gamma)` and
`J=(1-E)/gamma`. The rows and columns of the following matrices are
`(zero,live)`:

```
K_0 = [[1, 0], [0, E]],
K_1 = [[0, J], [J/M, (1-1/M)E]].
```

`K_0` describes an empty region. `K_1` averages the active packet's
uniform position within its region. Since the block's `v` active packets
occupy a uniform `v`-subset of the `R` regions, its limiting Laplace
moment is

```
F_{R,v}(gamma) = [u^v] (1,0)(K_0+u K_1)^R(1,1)^T / C(R,v).
```

This is an exact *continuum expression* for the ideal uniform setup law,
conditional on packet count `v`. It is not yet a finite-`L` bound. To use
it in a theorem we need a uniform approximation or direct finite-step
transfer, including packets that share or abut physical steps.

For relative output distance `delta`, the corresponding one-block
Chernoff first-moment expression is

```
L exp(gamma delta R/p_0) sum_{v>=1} B_O(v) F_{R,v}(gamma).
```

It is an upper bound on the probability of a low-weight one-block word
only after the finite-step relation to the limiting moment is justified.

## Shared-outer selection and numerical screen

An unconditional expectation is misleading here. Rare bad **shared outer
setups** can have extremely low-packet-weight words; their expectation is
multiplied by `L`, even though the paper already removes those setups with
probability `1-o(1)`. Its good-outer event gives binary weight at least
`0.104b` for every nonzero word. Therefore every such word has packet
support at least `ceil(ceil(0.104b)/4)`. On that event, replacing the
restricted expectation by the unconditional `E[B_O(v)]` only above this
packet cutoff gives a valid, but potentially loose, upper envelope.
The same event also gives binary weight at most `0.896b`; the present
packet-only envelope does not exploit that upper bound. For multiple
active blocks, one can additionally select a realized per-packet spectrum
majorant `B_O(v) <= b^2 E[B_O(v)]` at a one-time `O(1/b)` probability cost.

The following values use the exact integer packet census, then binary64
floating-point BAA folding and continuum optimization. They fix the
*comparison scale* `K=2^20` input bits (`L=2K/b`), even where this is not
an integer native size. They are diagnostic values, **not certificates**.

| `b` | active-packet floor | `log2` one-block expression | `log2` expression per outer position | dominant packet counts |
| ---: | ---: | ---: | ---: | --- |
| 216 | 6 | -1.23 | -14.48 | 6, 7, 50, 51 |
| 432 | 12 | -13.61 | -25.85 | 12, 13, 101, 100 |
| 648 | 17 | -20.22 | -31.88 | 17, 18, 152, 151 |
| 864 | 23 | -31.03 | -42.28 | 23 and nearby; run script for full ordering |

The decreasing per-position charge is encouraging: growing the outer
attacks the fixed-width obstruction while leaving the packet route and
inner intact. However, the paper's present schedule
`b=(39/4) log2 N+O(1)` would make `log2 L` grow at rate about `b/9.75`.
These finite screens show only about `0.05--0.07` bits of per-position
decay per additional outer bit, while that schedule requires more than
`1/9.75` bits. The current **bound** would need a larger growth constant,
roughly 15--21 in these finite screens. This is neither a proved
asymptotic threshold nor evidence that the actual construction fails at
`39/4`: the coarse good-outer cutoff admits bad binary-weight words with
packet support above the cutoff, and the high-packet-count terms remain
large. A joint binary-weight/packet-support envelope is the next useful
refinement.

## Proof gates before paper migration

1. Bound the finite-`L` packet-route/inner transfer uniformly in growing
   `b`, retaining packet collisions and the actual s22 output spectrum.
2. Select a joint binary-weight/packet-support envelope for the shared BAA
   outer. Determine the true one-block decay and a sufficient schedule
   `b=c log2 N` with a certified margin over `log2 L`.
3. Extend from one block to fixed, sparse growing, and positive block
   occupation. Reuse the selected outer across positions; never replace
   the expectation of a power by a power of the mean.
4. Certify all numerical inequalities with outward interval bounds and
   retain the finite RS8 certificates independently.

The current result identifies an exact outer statistic and an explicit
inner transfer for the same Turbo architecture. It does not yet replace
the paper's asymptotic theorem.
