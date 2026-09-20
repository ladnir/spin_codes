# Construction

Fix binary linear maps

\[
A:\mathbb F_2^{32}\to\mathbb F_2^{128},
\qquad
B:\mathbb F_2^{128}\to\mathbb F_2^{32},
\qquad BA=0.
\]

The fixed map `B` has systematic form \([H\mid I_{32}]\). The 96 columns of
\(H\) are distinct weight-three words. The 32 remaining columns are the
standard basis. The fixed map `A` is systematic on its first 32 output
coordinates. Its middle 64 coordinates are obtained from weight-seven
columns followed by two fixed accumulations with an intervening fixed
permutation. Its last 32 coordinates are the unique parity coordinates that
enforce \(BA=0\). The receipt files `fixed_b32_columns.txt` and
`fixed_a32_generators.txt` specify both maps exactly.

Setup samples an independent scalar
\(\alpha_i\gets\operatorname{GF}(2^{32})^*\) for each epoch. An epoch
receives input \(X_i\) and state \(Q_i\). It emits

\[
Y_i:=X_i+A(Q_i)
\]

and updates the state by

\[
Q_{i+1}:=\alpha_i Q_i+B(X_i).
\]

The first epoch starts with \(Q_0=0\).

For fixed setup randomness, the epoch map is binary linear and invertible.
Since \(BA=0\), the output determines \(B(X_i)=B(Y_i)\). Therefore

\[
Q_i=\alpha_i^{-1}(Q_{i+1}+B(Y_i)),
\qquad
X_i=Y_i+A(Q_i).
\]

The parent SplitState update is
\(Q_{i+1}=\alpha_i(Q_i+B(X_i))\). PreAddMul uses the same field
multiplication and XOR. Only their order changes.

The outer encoder, per-block coordinate permutations, bit transpose, and
independent permutation of each region are unchanged.
