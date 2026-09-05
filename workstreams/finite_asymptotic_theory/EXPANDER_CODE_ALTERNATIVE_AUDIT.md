# Two-sided regular Expand--Convolute alternative at \(k=2^{20}\)

## Correct use in the present design

The primary opportunity is to use a small expander-based code as the one
sampled \([512,256]\) outer constituent, then reuse that constituent in every
Structured SPIN row. The end-to-end expander construction below is a secondary
alternative, not the proposed replacement for BA.

`BLOCK_EXPAND_CONSTITUENT_ROUTE.md` defines an exact-size degree-14,
memory-15 block EC comparison and an iterated Block Expand--\(t\) family. The
current local target is Block Expand--5: the degree-14 regional map followed
by five independently interleaved accumulators. Conditional on a factor-two
shell-variance bound, its outward mean, cap, and full Structured SPIN transfer
give 42.6254 bits. The sole open obligation is the length-512 factor-two
variance theorem for one reused sample.

## Global alternative

The expander-code project contains a direct alternative to the repeated-block
SPIN construction. A strengthened invocation of its existing outward verifier
closes the present finite target with more than 50 bits of code-sampling
margin. The code is linear-time, but it is not the structured SPIN code and
has not received an implementation benchmark.

The verified parent parameters are

\[
 k_+=1{,}048{,}585,\quad n_+=2{,}097{,}170,\quad
 (d_L,d_R)=(10,5),\quad m=15,
\]

with bad-event cutoff \(L_+=228{,}607\). The in-memory certificate generator
used exact regional transfers for message supports 1 through 32. Its
192-bit Arb verifier returned

\[
 \Pr[d_{\min}(G_+)\le228{,}607]
 \le 7.703111446150614\mathbin\cdot10^{-16}
 <2^{-50.2054082202}. \tag{1}
\]

The prior checked-in certificate used exact transfers only through support
16. At the new cutoff, its coefficient relaxation starts at support 17 and
limits the proof to 34.7518 bits. Extending the exact range removes that
proof discontinuity. It does not change the sampled code.

## Parent construction

Let \(\ell=k_+/5=209{,}717\). The binary expander matrix
\(B\in\mathbb F_2^{k_+\times n_+}\) has ten regions of \(\ell\) columns. In
each region, setup samples a permutation of the \(k_+\) message positions,
splits the permuted list into groups of five, and connects each group to one
column. Every message coordinate has degree ten, and every intermediate
coordinate has degree five.

The wrapped convolution \(C\in\mathbb F_2^{n_+\times n_+}\) has memory 15.
For input \(u\), it initializes \(y_j=0\) for \(j\le0\). At each position
\(t\), setup samples independent bits \(b_{t,1},\ldots,b_{t,14}\), and the
encoder computes

\[
 y_t=u_t+y_{t-15}+\sum_{j=1}^{14}b_{t,j}y_{t-j}. \tag{2}
\]

The coefficient of \(u_t\) is one, so \(C\) is invertible. The parent
generator is \(G_+=BC\). The ten expander permutations and all convolution
bits are sampled once and are shared by every encoded message.

Equation (1) takes probability over these expander and convolution choices.
It is a first-moment distance proof for the complete sampled code. It does
not require concentration of a repeated local constituent.

## Exact power-of-two wrapper

The two-sided regular construction requires an odd message length divisible
by five. The smallest admissible length above \(2^{20}\) is \(k_+\). Fix nine
parent message coordinates to zero. This restriction leaves exactly
\(k=2^{20}\) message coordinates.

Next, puncture any fixed set of eighteen parent output coordinates. The
resulting length is

\[
 N=n_+-18=2^{21}.
\]

If the parent satisfies \(d_{\min}(G_+)\ge228{,}608\), shortening cannot
decrease its distance and puncturing decreases it by at most eighteen.
Therefore the wrapped code satisfies

\[
 d_{\min}\ge228{,}608-18=228{,}590,
 \qquad \frac{d_{\min}}N\ge0.1090002059. \tag{3}
\]

The same inequality also proves that puncturing preserves dimension. Thus,
subject to the parent event in (1), the wrapper produces a binary
\([2^{21},2^{20}]\) code with the target distance. The wrapper is
deterministic and adds no failure probability.

## Verification evidence

The source implementation is in
`C:/Users/peter/.codex/worktrees/30ff/permute_conv/expander_codes/`.
The decisive calculation invoked
`binary_biregular_ec_certificate.verify_candidate` with

```text
k = 1048585
cutoff = 228607
left_degree = 10
right_degree = 5
memory = 15
exact_limit = 32
target_bits = 40
precision_bits = 192
```

The outward contribution ledger was

| Support regime | Probability upper bound |
|---|---:|
| exact supports 1 through 32 | \(7.703111446150614\cdot10^{-16}\) |
| positive-coefficient outer blocks | \(7.964771199433166\cdot10^{-40}\) |
| central local-limit blocks | below \(10^{-1945}\) |
| full support | below \(10^{-344439}\) |

The weight-one term is the largest exact contribution. The verifier treats
the optimizer-selected decimal markers as fixed inputs and evaluates every
bound with outward Arb arithmetic.

## Scope and remaining work

This result changes the construction. It replaces the BCH constituent,
structured routing, and final SPIN inner map with one global two-sided
regular expander followed by a wrapped random convolution. It is a compelling
alternative when the application needs only rate, distance, and linear-time
encoding. It is not a proof about the spectrum of BCH--BA constituents.

The parent markers are now frozen in
`expander_ec_d10_m15_k20_parent_candidate.json`. The independent structural
checker accepts its complete support partition. The outward verifier
reproduces (1), and `certify_expander_ec_k20_wrapper.py` binds that result to
the deterministic wrapper. The combined receipt and artifact manifest are
`expander_ec_d10_m15_k20_wrapper_outward.json` and
`FINITE_K20_EXPANDER_EC_D10_M15_MANIFEST.json`.

The remaining task is implementation comparison. Ordinary and transposed
encoding both have linear operation counts: the transpose uses the reverse
recurrence for \(C^{\mathsf T}\) and traverses the same sparse edges for
\(B^{\mathsf T}\). However, ten sparse edge updates per input and the
memory-15 recurrence
may exceed the current performance target. Benchmark this construction
against the optimized BCH and Golay--BA encoders before changing the primary
design.
