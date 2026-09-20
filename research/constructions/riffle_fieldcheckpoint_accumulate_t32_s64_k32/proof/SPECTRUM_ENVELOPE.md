# Spectrum-density envelope

This lemma compresses the regular outer-spectrum sum to one active-block
count.

Let \(C\) be a binary linear \([B,k]\) code with spectrum \((A_w)\). For a
fixed nonzero codeword of weight \(w\), an independent coordinate
permutation makes its support uniform among the \({B\choose w}\) subsets of
weight \(w\). Therefore, the expected number of permuted nonzero codewords
with support \(S\) is

\[
 \mu(S):=\frac{A_{|S|}}{{B\choose |S|}}.
\]

Fix a collection \(\mathcal R\) of regular supports and define

\[
 \eta:=
 \max_{S\in\mathcal R}
 \frac{2^B\mu(S)}{2^k-1}.
\]

For every nonnegative function \(f\) on support tuples and every \(a\ge1\),

\[
 \sum_{(S_1,\ldots,S_a)\in\mathcal R^a}
 \left(\prod_{j=1}^a\mu(S_j)\right)
 f(S_1,\ldots,S_a)
 \le
 \eta^a
 \left(\frac{2^k-1}{2^B}\right)^a
 \sum_{(S_1,\ldots,S_a)\in(2^{[B]})^a}
 f(S_1,\ldots,S_a).
\]

The inequality follows by applying the pointwise bound on \(\mu(S_j)\) to
every regular support tuple. The right side is the support sum induced by
\(a\) independent uniform \(B\)-bit outer words, multiplied by
\(\eta^a(2^k-1)^a\).

For the modeled \([256,128,38]\)-shaped spectrum, let \(\mathcal R\) contain
the even weights from 38 through 218. Its multiplicities are proportional to
\({256\choose w}\). Consequently,

\[
 \eta=2
\]

up to the negligible difference between \(2^{128}\) and
\(2^{128}-1\). The implemented logarithmic calculation gives
\(\log_2\eta=1.00000000000011\).

The support of weight 256 is excluded. Its multiplicity is one, whereas its
point mass is much larger than the uniform-support density. The full
certificate must treat blocks carrying that all-one word as a separate
occupation type.
