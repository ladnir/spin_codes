# Riffle S-Stripe RandomStepConv g=4

This candidate applies the S-stripe permutation to four-bit packets. It is
distinct from the packet-width-eight family.

Let the message contain \(n=2^{20}\) bits. Setup independently samples each
outer constituent uniformly from the linear injections

\[
C_i:\mathbb F_2^{B/2}\longrightarrow\mathbb F_2^B.
\]

Define

\[
L:=\frac{2n}{B},
\qquad
q:=\frac{B}{4}.
\]

Write the encoded packets as an \(L\)-by-\(q\) matrix. In the full-stripe
instance, setup independently permutes each packet-coordinate column across
the \(L\) outer blocks. The encoder emits the \(q\) columns consecutively.

The one-lap inner encoder has state \(s_i\in\mathbb F_2^\sigma\). Setup
independently samples an unrestricted binary matrix \(M_i\) of dimension
\(4+\sigma\) at every packet position. The encoder initializes \(s_1=0\) and
computes

\[
(y_i,s_{i+1}):=M_i(x_i,s_i).
\]

The codeword is the concatenation of the four-bit packets \(y_i\). The final
state is discarded.

The outer injections, column permutations, and inner matrices remain fixed
for all encodings. The current random outer ensemble is a proof model rather
than an implementation recommendation.
