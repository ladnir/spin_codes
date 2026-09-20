# Proof status

Fix output length \(N=2^{21}\) and target distance
\(d=\lfloor0.09N\rfloor\).  The current calculation uses the modeled
complement-symmetric even spectrum for a binary \([256,128,38]\) outer code.
Arithmetic uses nearest binary64 values in the log semiring.

The floating-point certificate is complete.  It bounds the expected number
of nonzero messages of output weight at most 188743 by

\[
2^{-55.0716220969}.
\]

The bound covers every regular group-rank histogram, every zero/live state
trajectory, and every number and placement of exceptional all-one outer
words.

## Exact epoch reduction

After the spectrum-density reduction, an epoch contains \(A\) independent
fair candidate bits at arbitrary coordinates.  Define

\[
D:=2^{256}-1,
\qquad
p_A:=2^{-A},
\qquad
q_A(z):=\left(\frac{1+z}{2}\right)^A,
\qquad
S(z):=(1+z)^{256}.
\]

The exact transfer between the zero and uniform-nonzero state classes is

\[
F_A(z)=
\begin{pmatrix}
p_A & q_A(z)-p_A\\
(1-p_A)/D & (S(z)-q_A(z)-1+p_A)/D
\end{pmatrix}.
\]

This formula depends on \(A\), not on the candidate coordinates.  Therefore,
the 256-bit state removes packet-lane alignment from the inner analysis.  The
remaining packet dependence is the distribution of the candidate counts
among the 32 epochs of each region.

## Certificate reduction

The checker replaces each nonzero packet by one output bit.  A live epoch is
bounded by a normalized uniform 256-bit renewal transition and the global
factor \((2^{256}/(2^{256}-1))^{8192}\).  Its base-two logarithm is below
\(1.1\mathbin{\cdot}10^{-73}\).

Conditioned on \(H\) nonzero packets in one region, the checker computes the
exact shuffled region matrix \(R_H\).  It then takes the entrywise suffix
envelope \(\max_{K\ge H}R_K\).  The proved packing coupling is needed only for
the scalar count \(H\), so the small counterexamples to direct matrix packing
do not apply.

The regular class has 55.0716221062 bits of margin.  The all-one class has
82.2820153825 bits.  Their aggregate has 55.0716220969 bits.

## Earlier diagnostics

The 64-bit-state parent has a genuine medium-density obstruction.  Activating
all 2048 blocks in one packet lane leaves approximately 132,900 bits of
uncovered modeled multiplicity at 9% distance.

The earlier maximally packed exact-transfer diagnostic gave 55.9514 bits of
margin.  It was not a certificate because direct matrix packing is false.
The suffix-envelope certificate loses only 0.8798 bits and removes that gap.

The marked partial-group recurrence agrees with exhaustive enumeration in a
six-slot instance to error below \(6.2\times10^{-16}\).  The exact epoch
matrices are stochastic at \(z=1\) to displayed binary64 precision.

## Invalid shortcuts

Maximal group packing does not always maximize the final moment.  Exhaustive
toy instances contain small counterexamples.  Consequently, the packed-group
receipt is a diagnostic and not a certificate.

An adversarial epoch-composition bound discards the packet-permutation
probabilities and fails by more than 200,000 bits in the middle range.  A
zero/live bound that treats every return to zero as catastrophic also fails.
Neither relaxation is a viable certificate route.

## Remaining obligations

An outward-rounded implementation remains to turn the floating-point
certificate into a machine-checked numerical theorem.  An explicit outer-code
spectrum theorem remains to replace the modeled \([256,128,38]\) spectrum.

## Artifacts

- `receipts/packed_groups_full_regular_delta09.json`
- `receipts/all_profiles_regular_support_relaxation_delta09.json`
- `receipts/all_one_completion_delta09.json`
- `proof/FLOATING_POINT_CERTIFICATE.md`
- `../../scripts/analyze_riffle_packet4_fieldcheckpoint_packed_groups_s256.py`
- `../../scripts/analyze_riffle_packet4_fieldcheckpoint_packed_residue.py`
- `../../scripts/certify_riffle_packet4_fieldcheckpoint_s256_support.py`
- `../../scripts/certify_riffle_packet4_fieldcheckpoint_s256_allone.py`
- `../../scripts/verify_riffle_packet4_fieldcheckpoint_s256_certificate.py`
