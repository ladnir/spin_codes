# A window obstruction for direct 2K encoding

The unpunctured GRS[64,32] construction with four state elements cannot
provide 40% distance with a 40-bit setup margin at K=2^20. The obstruction below gives a lower bound
on actual failure probability. It is not a failed upper-bound certificate.

## Construction and window

Fix a field F and R separate message rows, each encoded by a [b,k] MDS
code over F. Collect each outer coordinate into one region of R symbols.
Independently sample a uniform permutation of each region.

A causal linear inner consumes packets of w symbols and emits w symbols.
It carries s field elements of state, initially zero. Zero input and zero
state produce zero output and zero next state. Assume the unpunctured
encoder is injective; the current full-Cauchy inner satisfies this because
its input-to-output block is invertible at every step.

Choose W divisible by w, and retain the first W input positions of each
region as a window. For outer row i, let t_i count its coordinates routed
into these windows. The space of outer messages supported entirely inside
the windows has dimension

\[
S:=\sum_{i=1}^{R}\max\{0,t_i-(b-k)\}.
\]

This follows from MDS shortening, separately in each message row. For any
fixed routing and inner coefficients, require the state to be zero at the
end of each window. These are at most bs homogeneous linear constraints.
If S>bs, a nonzero message satisfies every constraint. Its inner output
vanishes after each window and before the next. Injectivity makes its
output nonzero, so the output has weight at most bW.

This argument holds for every realization of the local maps and diagonals.
It uses no randomness from them.

## Failure probability

Write p:=W/R. For each fixed row, independent regional permutations give
t_i distributed as Binomial(b,p). Rows need not be independent. Therefore

\[
\mathbb E[S]
=R\sum_{t=b-k+1}^{b}(t-b+k)\binom bt p^t(1-p)^{b-t}.
\]

Exactly bW outer coordinates enter the windows. Since
\((t-(b-k))_+\le kt/b\), we have \(S\le kW\). Thus, when kW>bs,

\[
\Pr_{\text{regional permutations}}[S>bs]
\ge
\max\left\{0,\frac{\mathbb E[S]-bs}{kW-bs}\right\}.
\]

Choose W/R strictly below the target relative distance. The displayed
probability is then a lower bound on the probability of violating that
distance, regardless of the field size or the inner coefficients.

## Concrete direct-2K instance

For b=64, k=32, R=32768, w=s=4, the encoder has K=2^20 and N=2K.
There is no shortening or padding in this instance. Exact rational
evaluation gives the following conservative rounded bounds.

| Target distance | W | E[S], approximately | Failure probability | Setup margin cannot exceed |
|---|---:|---:|---:|---:|
| 35% | 11468 | 282.915 | 7.33e-5 | 13.74 bits |
| 40% | 13104 | 2835.266 | 6.15e-3 | 7.35 bits |

The probabilities in the table are lower bounds, not estimates of the
actual failure rate. In particular, increasing the diagonal alphabet
cannot repair this direct-2K geometry.

The argument does not rule out the larger-outer, punctured 4K candidate.
A zero lower bound from this test also does not establish good distance.

Reproduction: `theory/high_distance_window_obstruction.py` computes the
bounds using exact fractions. Its tests check the expectation and lower
bound against exhaustive small routing examples. The independent narrow
inner tests verify the full-Cauchy ranks and agreement with the width-four
transfer formula.
