# Goal 02: full-size support-33 moment gate

## Question

Does the parallel accumulator suppress the 26 authenticated support-33 outer
words strongly enough to justify further proof work?

## Bound

Fix one outer word with active values \(v_1,\ldots,v_h\), where \(h=33\).
The uniform packet permutation independently determines their labeled order and
their support positions.

For parameters \(0<z,\tau<1\), combine the Chernoff inequality with the
coefficient bound from Goal 01. If \(W\) is the binary output weight, then

\[
\Pr[W\le d]
\le
z^{-d}\tau^{-(n-h)}
\frac{1}{\binom nh(1-\tau)}
\mathbb E_{\mathrm{ord}}
\left[
\prod_{j=1}^{h}
\frac{z^{a_j}}{1-\tau z^{a_j}}
\right],
\]

where \(a_j\) is the bit weight of the prefix XOR after ordered packet \(j\).

The expectation is evaluated over every labeled ordering. Its dynamic-program
state contains only the remaining multiplicities and current four-bit XOR.

## Gate

Evaluate the optimized numerical bound for all 26 words and aggregate their
contributions. A numerical result below \(2^{-40}\) permits a later rigorous
hardening. A result near or above \(2^{-40}\) rejects this proof route unless a
materially sharper inequality is identified.

This goal does not cover other outer supports or prove linear distance.
