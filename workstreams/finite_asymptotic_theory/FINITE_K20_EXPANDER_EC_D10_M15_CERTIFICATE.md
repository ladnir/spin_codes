# Finite \(k=2^{20}\) certificate for two-sided regular Expand--Convolute

## Theorem

There is an efficiently sampled binary linear-code ensemble with message
length \(k=2^{20}\), output length \(N=2^{21}\), and

\[
 \Pr[d_{\min}<228{,}590]
 <2^{-50.2054082202}<2^{-40}. \tag{1}
\]

Consequently, a sampled code has relative minimum distance at least

\[
 \frac{228{,}590}{2^{21}}=0.10900020599365234
\]

except with the probability in (1). Ordinary and transposed encoding take
\(O(k)\) word operations when the expander degrees and convolution memory are
fixed.

## Parent ensemble

Set

\[
 k_+=1{,}048{,}585,\qquad n_+=2{,}097{,}170,
 \qquad \ell=k_+/5=209{,}717.
\]

The parent generator is \(G_+=BC\). The binary matrix
\(B\in\mathbb F_2^{k_+\times n_+}\) is sampled from the two-sided regular
ensemble with left degree 10 and right degree 5. Its columns are partitioned
into ten consecutive regions of length \(\ell\). Independently in each
region, setup samples a permutation of the \(k_+\) message positions, divides
the permuted list into \(\ell\) groups of five, and connects each group to one
column.

Independently sample a wrapped convolution \(C\) of memory 15. For an
intermediate word \(u\in\mathbb F_2^{n_+}\), initialize \(y_j=0\) for
\(j\le0\). At each position \(t\), sample independent bits
\(b_{t,1},\ldots,b_{t,14}\) and define

\[
 y_t:=u_t+y_{t-15}+\sum_{j=1}^{14}b_{t,j}y_{t-j}. \tag{2}
\]

The coefficient of \(u_t\) is one. Thus \(C\) is invertible. Setup samples
the ten regional permutations and all coefficients in (2) once. Every
message uses the resulting fixed generator.

## Parent distance bound

For \(L_+=228{,}607\), the exact first-moment verifier partitions all nonzero
message supports as follows:

- exact regional transfers cover supports 1 through 32;
- 108 positive-coefficient blocks cover the two outer ranges;
- 47 complemented local-limit blocks cover the central range; and
- a separate transfer covers the all-one support.

The independent structural checker verifies that these ranges are disjoint
and cover every support from 1 through \(k_+\). For fixed decimal markers,
192-bit Arb arithmetic gives

\[
\begin{aligned}
 U_{\rm exact}&<7.703111446150614\mathbin\cdot10^{-16},\\
 U_{\rm outer}&<7.964771199433167\mathbin\cdot10^{-40},\\
 U_{\rm central}&<6.026\mathbin\cdot10^{-1946},\\
 U_{\rm full}&<8.730\mathbin\cdot10^{-344440}.
\end{aligned}
\]

Their outward-rounded sum is below
\(7.703111446150614\mathbin\cdot10^{-16}\). Markov's inequality for the
number of low-weight nonzero messages therefore proves

\[
 \Pr[d_{\min}(G_+)\le228{,}607]
 <2^{-50.2054082202}. \tag{3}
\]

The weight-one message term dominates (3). Extending the exact transfer from
support 16 to support 32 is essential only to the proof. It leaves the parent
ensemble unchanged.

## Power-of-two wrapper

Fix nine parent message coordinates to zero. Since the good event in (3)
implies that \(G_+\) is injective, this shortening leaves a message space of
dimension \(k_+-9=2^{20}\).

Puncture any fixed set of eighteen output coordinates. The output length is
\(n_+-18=2^{21}\). On the good parent event,

\[
 d_{\min}\ge (228{,}607+1)-18=228{,}590.
\]

This distance also exceeds zero, so puncturing preserves the shortened
dimension. The wrapper is deterministic and adds no failure probability.
Equation (1) follows from (3).

## Scope

This theorem concerns a global Expand--Convolute code. It does not certify
the BCH--BA constituent or the structured SPIN routing. It gives a direct
linear-time alternative when minimum distance is the required property.

The theorem is an ensemble statement and provides no decoding algorithm.
The operation count has not been benchmarked. The transposed convolution is
the reverse-scan adjoint of (2), and the transposed expander traverses the same
edges in reverse. In particular, ten sparse edge updates per message
coordinate and the memory-15 recurrence may have a larger constant than the
optimized BCH-based encoders.

## Artifacts

The frozen parent markers and support partition are in
`expander_ec_d10_m15_k20_parent_candidate.json`. The independent checker and
outward parent verifier are in the expander-code worktree named by
`expander_ec_d10_m15_k20_wrapper_outward.json`; that receipt records their
absolute paths and SHA-256 digests.

`generate_expander_ec_k20_candidate.py` reproduces the parent candidate.
`certify_expander_ec_k20_wrapper.py` reruns the structural and outward checks,
verifies the wrapper arithmetic, and writes the combined receipt.
