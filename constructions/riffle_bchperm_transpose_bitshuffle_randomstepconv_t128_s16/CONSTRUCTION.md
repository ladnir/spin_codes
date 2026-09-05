# Riffle BCHPerm-TransposeBitShuffle-RandomStepConv t=128 s=16

This benchmark keeps the modeled binary `[256,128]` outer code, the sampled
coordinate permutations, the bit transpose, and the sampled permutation of
each transposed region.

The transposed word contains 16,384 epochs. An epoch receives

\[
X_i\in\mathbb F_2^{128}
\]

and a state

\[
Q_i\in\mathbb F_2^{16}.
\]

Setup independently samples

\[
M_i\gets\mathbb F_2^{144\times144}
\]

for every epoch. The encoder initializes \(Q_0:=0\) and computes

\[
(Y_i,Q_{i+1}):=M_i(X_i,Q_i),
\qquad
Y_i\in\mathbb F_2^{128}.
\]

The encoder emits \(Y_0\|\cdots\|Y_{16383}\) and discards
\(Q_{16384}\). The matrices remain fixed across encoding operations.

For every fixed nonzero \((X_i,Q_i)\), the pair
\((Y_i,Q_{i+1})\) is uniform in \(\mathbb F_2^{144}\) over the setup
randomness. In particular, conditioning on \(Q_{i+1}=0\) leaves \(Y_i\)
uniform. The zero-state identity class of SplitState-PreAddMul is therefore
absent.

This benchmark samples unrestricted matrices. Requiring invertibility or
reusing one matrix at every epoch defines a different benchmark.
