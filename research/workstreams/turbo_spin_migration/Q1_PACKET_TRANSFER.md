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

### A finite-step upper transfer

The direct transfer is available for this one-block class. Let `L=16n`
and set `theta=gamma/(4Lp_0)`. Enumerating the chosen s22 expansion map
over all `2^22-1` nonzero states gives its exact 64-bit output-weight
spectrum:

| Weight | 16 | 24 | 28 | 32 | 36 | 40 | 48 | 64 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Count | 2,604 | 291,648 | 888,832 | 1,828,134 | 888,832 | 291,648 | 2,604 | 1 |

The script `verify_pair_return_map.py --spectrum` reproduces the count
from `productionMap.h`. Define

```
q(theta) = sum_{a != 0} exp(-theta wt(Aa))/(2^22-1),
Q = q(theta)^n,
J_n = (1/n) sum_{h=0}^{n-1} q(theta)^h.
```

If every step with a live *pre-step* state were charged its mean Laplace
factor `q`, the empty- and active-region transfers would be exactly

```
K_0^(n) = [[1, 0], [0, Q]],
K_1^(n) = [[0, J_n], [q J_n/M, (1-1/M)Q]].
```

The active packet occupies a uniform one of the `n` physical steps in
its region. The `q` in the live-to-zero entry charges its active step;
the zero-to-live entry does not, because that step starts in state zero.
GL updates give the zero/live transition probabilities exactly even
when the incoming live state is not uniform.

To relate this transfer to actual output weight, mark every active step
and its immediate successor. There are at most `2v` marked steps. At
every unmarked live step the pre-state is independent uniform nonzero:
an empty update after a nonzero state refreshes it, and a later packet's
return event is independent of that earlier state. Conditional on the
zero/live trajectory, unmarked live-step output weights therefore have
factor `q`; other unmarked steps have zero output. Each marked step has
at most 64 output bits. Replacing its true factor by the `q` factor
changes the moment by at most `exp(64 theta)` in the upper direction.
Consequently, for fixed packet labels and `v` active packets,

```
E[exp(-theta W) | v]
  <= exp(128 theta v)
     [u^v] (1,0)(K_0^(n)+u K_1^(n))^R(1,1)^T / C(R,v).
```

The expectation includes packet and region shuffles and uniform GL
updates. The inequality is for the ideal sampled setup law, not the
seeded implementation. It needs a written lemma-level audit before use
in a paper; the numerical implementation uses binary64, not outward
rounding. As `L -> infinity` with `b=O(log L)`, its entrywise transfer
converges to `K_0,K_1` above and the correction is `exp(O(b/L))`.

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
inner intact. With *only* the paper's binary-weight selection, its present schedule
`b=(39/4) log2 N+O(1)` would make `log2 L` grow at rate about `b/9.75`.
These finite screens show only about `0.05--0.07` bits of per-position
decay per additional outer bit, while that schedule requires more than
`1/9.75` bits. The current **bound** would need a larger growth constant,
roughly 15--21 in these finite screens. This is neither a proved
asymptotic threshold nor evidence that the actual construction fails at
`39/4`: the coarse good-outer cutoff admits bad binary-weight words with
packet support above the cutoff. An exact joint binary-weight/packet-
support count improves the screen, but the lowest allowed packet count
still dominates. For example, at `b=432`, counting binary weight exactly
through `v=20` changes `log2 Q1` from `-13.61` to `-18.29`, while `v=12`
still contributes most of the bound.

## Selecting packet support once

A more effective selection condition is that every nonzero word of the
one shared BAA outer occupy at least `ceil(rho R)` packets. Markov bounds
the one-time setup failure probability by the exact expected number
`sum_{1<=v<ceil(rho R)} E[B_O(v)]`. Unlike the one-block first moment,
this event is **not multiplied by `L`**. The following binary64 screens
use `rho=0.16` and the unconditioned exact packet enumerator:

| `b` | packet floor | expected bad outer words | `log2` conditional continuum Q1 | indicative `b/log2 L` threshold |
| ---: | ---: | ---: | ---: | ---: |
| 216 | 9 | 0.003759 | -14.04 | 7.92 |
| 432 | 18 | 0.001552 | -41.55 | 8.03 |
| 648 | 26 | 0.000589 | -63.23 | 8.65 |
| 864 | 35 | 0.000338 | -89.15 | 8.61 |

The threshold column divides `b` by the magnitude of the computed
per-position log2 charge. It is a diagnostic finite-size ratio, not a
proved asymptotic constant. At native `L=4864` and `b=432`, the
finite-step upper formula gives conditional `log2 Q1=-40.57`; its
separate packet-floor setup-failure expectation is still `0.001552`.
Neither number certifies the complete code distance.

Splitting the `rho=0.16` bad-word expectation by the number of active
Golay source blocks identifies the asymptotic proof split. At `b=432`,
one active Golay block contributes `0.001441` of the total `0.001552`;
two contribute `0.0000998`, three `0.0000101`, and the rest less. At
`b=864`, one block contributes `0.000327` of `0.000338`. A dense-only
large-deviation exponent cannot establish this selection event: it must
also control the sparse source classes, as the existing BAA binary-tail
proof does. The numerical decrease is consistent with a vanishing
failure probability but is not itself a proof of `o(1)`.
The [packet-floor tail note](PACKET_FLOOR_TAIL.md) proves the sparse-source
contribution is `O(b^(-3/2))` and certifies a strictly negative
dense-source exponent with outward-rounded interval arithmetic. Thus
the additional packet-floor setup event has probability `1-o(1)`.

This packet selection makes it plausible to retain the paper's
`39/4` schedule with Turbo's route and inner. It changes only which
sampled outer realizations the proof accepts; it does not alter the
encoder or add setup work. The one-block schedule inequality and all
`Q>=2` message classes remain open.

## Proof gates before paper migration

1. Audit and formalize the finite-step one-block transfer above. It
   uses the exact s22 output spectrum and a factor for exceptional steps.
2. Certify the one-block schedule inequality using the selected 16%
   packet floor and a realized packet-spectrum majorant.
3. Extend from one block to fixed, sparse growing, and positive block
   occupation, including packet collisions between blocks. Reuse the
   selected outer across positions; never replace
   the expectation of a power by a power of the mean.
4. Certify all numerical inequalities with outward interval bounds and
   retain the finite RS8 certificates independently.

The current result identifies an exact outer statistic and an explicit
inner transfer for the same Turbo architecture. It does not yet replace
the paper's asymptotic theorem.
