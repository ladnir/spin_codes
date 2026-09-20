# Riffle SpectrumPerm-TransposeBitShuffle-RandomStepConv g

Fix a binary linear code

\[
 C:\mathbb F_2^{B/2}\longrightarrow\mathbb F_2^B
\]

with weight spectrum \(A_0,\ldots,A_B\). Thus \(A_w\) counts messages whose
codeword has weight \(w\). Injectivity gives \(A_0=1\) and
\(\sum_w A_w=2^{B/2}\).

Let the complete message contain \(n\) bits. Define

\[
 L:=\frac{2n}{B},
 \qquad
 P:=\frac{L}{g}.
\]

Setup samples independent coordinate permutations

\[
 \pi_j\gets S_B
 \qquad (j=1,\ldots,L).
\]

The permutations remain fixed for all encodings. Outer block \(j\) computes

\[
 c_j:=C(m_j),
 \qquad
 \widetilde c_j:=\pi_j(c_j).
\]

Place the words \(\widetilde c_j\) in an \(L\)-by-\(B\) matrix and transpose
the matrix. Setup independently permutes all \(L\) bits in each transposed
row. The encoder partitions every permuted row into consecutive \(g\)-bit
inputs.

The one-lap inner encoder has state \(s_i\in\mathbb F_2^\sigma\). Setup
samples an independent unrestricted linear map

\[
 M_i:\mathbb F_2^{g+\sigma}\longrightarrow
     \mathbb F_2^{g+\sigma}
\]

at every inner position. Starting with \(s_1=0\), the encoder computes

\[
 (y_i,s_{i+1}):=M_i(x_i,s_i).
\]

The codeword is the concatenation of the \(y_i\). The encoder discards the
terminal state.

## Why the spectrum is sufficient for a first moment

Fix a message difference. Suppose its nonzero outer words have weights
\(w_1,\ldots,w_a\). For each active block \(j\), the support of
\(\pi_j(c_j)\) is a uniform \(w_j\)-subset of the \(B\) rows. These supports
are independent across outer blocks.

The distribution for this fixed message therefore depends on each outer
word only through its weight. The spectrum supplies the number of messages
with each weight. Correlations between different messages do not enter a
first-moment union bound.

The coordinate permutations must be independent across outer blocks. A
shared permutation would preserve relative supports and require more than
the individual weight spectrum.

The rows are not independent. A word of weight \(w\) activates exactly \(w\)
rows. The proof must retain this sampling-without-replacement law or replace
it with a proved upper bound.

## Claim scope

This construction has no distance certificate yet. The first experiment
keeps the random inner ensemble and determines which outer spectra suffice.
The ideal random-code spectrum is only a benchmark. It is not the spectrum
of a specified fixed code.
