# Goal 02 report: full-size support-33 gate

## Status

The support-33 gate is unresolved. The current upper bound misses the
\(2^{-40}\) target, while the strongest explicit lower family remains below
the target.

This result is evidence against using one parallel accumulator as the deployed
constituent. It is not a proof or refutation of the finite construction.

## Full-size upper bound

Fix one authenticated support-33 outer word. Its labeled active packets undergo
the uniform global packet permutation. Goal 01 gives the exact placement
generating function after fixing the active-value order.

For \(0<z,\tau<1\), the Chernoff and coefficient bounds give

\[
\Pr[W\le d]
\le
z^{-d}\tau^{-(n-h)}
\frac{1}{\binom nh(1-\tau)}
\mathbb E_{\mathrm{ord}}
\left[
\prod_{j=1}^{h}
\frac{z^{a_j}}{1-\tau z^{a_j}}
\right].
\tag{1}
\]

Here \(n=524{,}352\), \(h=33\), and \(d=188{,}766\). The value \(a_j\) is
the bit weight of the prefix XOR after active packet \(j\).

A mixed-radix dynamic program evaluates the expectation in (1) over every
labeled active-value order. The state contains the used multiplicity of each
nonzero packet value. Those multiplicities determine the current XOR state.

The 26 outer words reduce to 14 distinct packet-value multisets. Numerical
optimization gives

\[
\sum_{x\in\mathcal X_{33}}
\Pr[W_x\le d]
\le 2^{-20.623793}
\]

under the numerical evaluation of (1). The bound is far above \(2^{-40}\).
It therefore cannot certify this shell.

The worst multiset belongs to

```text
terminal-w3-s33-one-i0-t6b16d09000000000
```

and has value counts

\[
\{1^3,5^9,7^3,10^9,11^3,12^3,14^3\}.
\]

Its 33 packets have total XOR zero. Its individual optimized bound is

\[
\Pr[W\le d]\le 2^{-21.830205}.
\tag{2}
\]

An independent top-down recursion reproduces (2) within
\(2.7\cdot10^{-9}\) bits at the recorded parameters.

## Explicit lower families

The matching refutation probe restricts every active-order prefix to state
weight at most \(c\). It then limits the total gap length following nonzero
prefixes. These two conditions imply output weight at most \(d\).

For the zero-final-state multiset above, the resulting lower bounds are:

| Maximum prefix weight \(c\) | Lower-family contribution |
|---:|---:|
| 2 | \(2^{-63.519846}\) |
| 3 | \(2^{-56.305469}\) |
| 4 | \(2^{-63.212782}\) |

The \(c=3\) family is strongest. It remains 16.305 bits below the refutation
threshold. Therefore it does not refute the finite construction.

Combining the current results gives the diagnostic bracket

\[
2^{-56.305469}
\le
\Pr[W\le d]
\le
2^{-21.830205}
\]

for the worst authenticated word. The bracket is too wide for a construction
decision.

## Asymptotic obstruction

A single parallel accumulator has a general late-window event. Fix any input
with \(h\) active packets. If all active packets occupy the final \(\rho n\)
packet positions, then

\[
W\le 4\rho n.
\]

For fixed \(h\), the probability of this event approaches \(\rho^h\) as
\(n\) grows. This probability does not vanish with \(n\).

The zero-final-state word has an additional translated-window event. If all
active packets lie in any interval of length \(\rho n\), the output vanishes
outside that interval. The same constant-support phenomenon underlies the
known weakness of one-accumulator serial constructions.

Consequently, one parallel accumulator cannot yield a standard
high-probability linear-distance theorem when every length has even one outer
word of bounded packet support. The late-window event gives that word a fixed
positive failure probability. This statement concerns an asymptotic family.
It does not decide the fixed \(N=2^{21}\) and 40-bit target.

## Decision

Riffle ParallelAcc g=4 succeeds as an analytical stress test. Its exact
constituent interface is substantially simpler than the two-lap BCH interface.

The construction does not yet look competitive as a replacement constituent.
Its first full-size sparse shell already defeats the natural moment proof, and
the late-window event prevents the desired clean asymptotic theorem.

The next bounded finite-size experiment should use importance sampling tilted
by the optimizer in (1). Apply it only to the zero-final-state multiset. That
experiment should estimate whether its true probability lies above or below
\(2^{-40}\) before any wider ParallelAcc proof work continues.
