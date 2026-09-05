# Riffle SplitState-PreAddMul-RM2Sub t=128 s=20

This folder records the fixed-state-size-20 experiment for the transposed
PCG evaluator. The inner step is

\[
Y_i=X_i+A(Q_i),\qquad
Q_{i+1}=\alpha_iQ_i+B(X_i),
\]

where (A:\mathbb F_2^{20}\to\mathbb F_2^{128}),
(B:\mathbb F_2^{128}\to\mathbb F_2^{20}), and (BA=0). Each fresh
nonzero scalar \(\alpha_i\) acts on the 20-bit state field.

The selected RM(2,7) subcode has seed 2719617393. The selector exactly
enumerated 1,024 candidates and minimized their low-weight image spectra.
The selected image has minimum distance 48 and spectrum

\[
1+10416z^{48}+220416z^{56}+586910z^{64}
+220416z^{72}+10416z^{80}+z^{128}.
\]

The exact MacWilliams transform gives a minimum kernel distance of six for
\(B\). In particular, its kernel contains no weight-four input.

For output length (N=2^{21}), the regular-occupation calculation closes at

\[
D=228023,
\qquad
D/N=0.108729839\ldots.
\]

The modeled margins at this distance are 61.6379 bits for occupation one,
108.2055 bits for occupations 2 through 100, and 42.7641 bits for
occupations 101 through 8192. The last range is dominated by occupation
8192. At (D=228024), that class has only 39.7251 bits of margin.
Therefore (D=228023) is the 40-bit endpoint of the current regular-class
calculation.

The result does not reach 11 percent. At 11 percent, the all-active class
fails the current bound by about 8,018 bits. Thus increasing the state from
16 to 20 removes the earlier persistent-zero-state obstruction, but a dense
live-state class becomes limiting near the random-convolution benchmark.

No larger state closes 11 percent under the current outer-spectrum bound.
For the all-active occupation, the Bernoulli reference makes each fresh
epoch input uniform. Adding (A(Q_i)) therefore leaves each epoch output
uniform, independently of the state size and constituent. The modeled even
outer spectrum costs one extra envelope bit per active outer block. Across
8,192 blocks, this produces an 8,192-bit loss. An ideal random-coset inner
at (s=64) still has (-8017.5650) bits of margin.

Consequently, the next useful change is to preserve the outer block's even
parity in the dense analysis, or to use a non-even outer constituent. A
larger inner state cannot remove the present dense-class loss.

The follow-up `Riffle ParityShear-12 SplitState s=19` variant implements the
first option with a sparse invertible outer shear. Its regular calculation
closes at 11 percent with 46.3294 bits.

This is not yet an end-to-end theorem. The calculation assumes the modeled
even-floor outer spectrum, uses nearest binary64 arithmetic, and does not
rerun the special all-one and mixed all-one classes.

See `PROOF_STATUS.md` for the receipt map and exact scope.
