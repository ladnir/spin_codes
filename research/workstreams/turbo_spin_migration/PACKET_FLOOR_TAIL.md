# A packet-support floor for the shared Golay--BAA outer

## Result and scope

Let `b` grow through multiples of 24 and let `R=b/4`. Sample one
Golay--BAA constituent `O_b=Acc_b o tau_2 o Acc_b o tau_1 o H_{24,b}`
using independent uniform permutations. Let `B_O(v)` count its nonzero
words with `v` active four-bit output packets. We prove

```
E[ sum_{1 <= v < (4/25)R} B_O(v) ] = o(1).             (1)
```

Rounding the strict inequality to the integer packet floor is understood.
Markov then shows that one sampled, shared outer has packet support
at least `ceil((4/25)R)` for every nonzero word with probability `1-o(1)`.
This condition improves Turbo's one-active-block bound without changing
the encoder. The exact finite packet enumerator is in
[Q1_PACKET_TRANSFER.md](Q1_PACKET_TRANSFER.md).

We split messages by `j`, the number of active 24-bit Golay source
blocks, as the paper does for its binary-weight tail. The sparse part
below closes by adapting that proof. The dense part uses an
outward-rounded interval certificate for the packet-aware exponent.
This result concerns the selected outer only; it does not prove Turbo's
global distance.

## Sparse source blocks: an analytic bound

Fix `rho=4/25=0.16`, `eta=1/5000`, and `1<=j<=eta b`. A word with fewer
than `rho R` active packets has binary weight less than `rho b`. Thus the
paper's one-accumulator lower-tail inequality applies to the final
accumulator. For a uniformly interleaved input of weight `c`, it gives

```
Pr[final packet support < rho R | c]
 <= Pr[final binary weight < rho b | c]
 <= (16/21)^floor(c/2)
 <= zeta^(-1) zeta^c,
zeta = sqrt(16/21).
```

The initial Golay output has even weight `a=2l`, where
`4j<=l<=12j`. Let `H_1` be the weight after the first accumulator.
The paper's sparse-moment calculation, with this new `zeta`, gives

```
E[zeta^H_1 | a=2l]
 <= (kappa l/(b-2l+1))^l,
kappa = 4 zeta/(1-zeta).
```

The right-hand side decreases with `l` throughout `l/b<=12eta`:
the logarithmic derivative is at most
`ln(kappa x/(1-2x))+1+2x/(1-2x)` at `x=12eta`, which is less than
`-1.6` using `kappa<28`. There are at most
`(4095 e b/(24j))^j` messages with `j` active Golay blocks. Hence their
total bad-packet expectation is at most

```
zeta^(-1) sum_{1<=j<=eta b} [C_* (j/b)^3]^j,
C_* = (4095 e/24) (4 kappa/(1-8eta))^4.
```

Here `zeta<7/8` follows by squaring rational numbers, so `kappa<28`.
Using `e<68/25` and `eta=1/5000` gives `C_* eta^3<0.59`;
`verify_packet_tail_sparse.py` checks these rational inequalities.
Split the sum at `sqrt(b)`. For smaller `j` it is
bounded by a geometric series with ratio `O(b^(-3/2))`; for larger `j`
it is bounded by a geometric tail with ratio below `0.59`. Therefore
the sparse contribution to (1) is `O(b^(-3/2))`. This argument needs
no new BAA assumption and no numerical optimization.

## Dense source blocks: packet-aware exponent

For `j>eta b`, write `alpha=a/b`, `beta=c/b`, where `a` is the Golay
output weight and `c` the first accumulator output weight. Then
`alpha>=8eta=0.0016`. The paper supplies the Golay exponent `g(alpha)`
and first-accumulator exponent `pi(alpha,beta)`.

The final accumulator requires a packet statistic. Its output bit
immediately before a packet is `p in {0,1}`. A packet `q` contributes
`t(p,q)` transitions and activity `1[q!=0000]`, as defined in the
one-block note. Summing its sixteen possibilities gives the positive
matrix

```
T(x,u) = [ 1+u(6x^2+x^4)       4ux(1+x^2)       ]
         [ x+u(3x+4x^3)       u(1+6x^2+x^4)  ].
```

For `x>0` and `0<u<1`, a coefficient bound and the packet-tail
Chernoff factor give the following normalized exponent for the expected
number of dense-source bad words:

```
g(alpha) + pi(alpha,beta)
  + (1/4) ln lambda(T(x,u)) - beta ln x
  - (rho/4) ln u - h(beta),                          (2)
```

where `lambda` is the Perron root and `h` is binary entropy in nats.
For each feasible `(alpha,beta)`, we may minimize (2) over `x,u`.
A strictly negative *uniform* supremum over
`alpha in [0.0016,1]` and
`beta in [alpha/2,1-alpha/2]` would make the dense contribution
`exp(-Omega(b))`, up to polynomial coefficient factors. The paper's
uniform coefficient and Stirling bounds handle the endpoints; each
accepted box uses a fixed positive matrix with a finite Perron prefactor.

`verify_packet_tail_interval.py` certifies the uniform strict inequality.
It covers the full rectangle `alpha in [1/625,1]`, `beta in [0,1]`,
discarding only boxes disjoint from the accumulator feasibility region.
For each remaining box it chooses fixed positive rational `x,u` and
upper-bounds (2) with outward-rounded interval arithmetic. The verifier
reuses the paper's certified Golay and first-accumulator interval bounds.
Its completed run processes 995 boxes, accepts 996 leaves, reaches depth
24, and has maximum accepted upper bound below `-3.26e-7` nats per
outer bit. Thus the dense contribution is `exp(-Omega(b))`, since the
coefficient and Stirling remainders are polynomial in `b` as in the
paper's outer-tail proof. `packet_tail_exponent.py` remains a faster
binary64 diagnostic; it is not the certificate.

## Finite diagnostic and next gate

The exact packet enumerator gives expected numbers of words below the
16% packet floor of `0.003759`, `0.001552`, `0.000589`, and `0.000338`
at `b=216,432,648,864`. At `b=432`, one active Golay source block
contributes `0.001441` of `0.001552`; at `b=864`, it contributes
`0.000327` of `0.000338`. These numbers accord with the sparse bound.
They are binary64 evaluations, not outward-rounded certificates.

The sparse and dense bounds together establish (1), hence the shared
outer's 16% packet-support floor with probability `1-o(1)`. Multi-block
packet collisions and inner-state cancellations still require proofs
before replacing Structured SPIN's asymptotic theorem.

## Stronger selected floors

The same argument certifies packet-support floors of 17%, 18%, and 19%.
For each row, `verify_packet_tail_sparse.py` checks the rational sparse
constants and `verify_packet_tail_interval.py` checks (2) over the
entire feasible dense domain with outward-rounded interval arithmetic.

| `rho` | sparse split `eta` | dense accepted leaves | maximum dense upper (nats per bit) |
| ---: | ---: | ---: | ---: |
| 0.16 | 1/5,000 | 996 | below -3.26e-7 |
| 0.17 | 1/10,000 | 1,332 | below -1.36e-6 |
| 0.18 | 1/20,000 | 1,734 | below -3.34e-7 |
| 0.19 | 1/50,000 | 2,296 | below -2.27e-7 |

Thus the strongest currently certified statement replaces `4/25` in
(1) by `19/100`: the one sampled Golay--BAA outer has packet support
at least `ceil(0.19R)` for every nonzero word with probability `1-o(1)`.
The proof accepts a more restrictive outer event; setup and encoding
algorithms are unchanged. This improves the one-block inner bound, but
does not by itself prove the full Turbo distance theorem.
