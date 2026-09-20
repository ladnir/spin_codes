# Power-of-two repeated EBCH128 with RandomStepConv

## Current result

The unpadded power-of-two instance now has a complete outward certificate at
relative distance 10.9%. Its parameters are

\[
 L=16384,\qquad k=64L=2^{20},\qquad N=128L=2^{21},
 \qquad D=\lceil0.109N\rceil=228590.
\]

For RandomStepConv-M20, the verifier proves

\[
 \Pr[d_{\min}<D]<2^{-7.58436084118}<2^{-5}.
\]

Occupation one is dominant. The proof uses the exact BCH spectrum for that
occupation. It uses the existing parity-pivot transfer through occupation
16319 and an exact-complement recurrence for the final 65 occupations.

The theorem, receipt, and proof details are in
`FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_CERTIFICATE.md` and
`ebch128_randomstepconv_g1_s20_pow2_outward_d109.json`.
`BITTRANSPOSE_RANDOMSTEP_CONV_PROOF_TEMPLATE.md` freezes this proof
architecture for later outer, routing, and inner substitutions.

At M19, the complete binary64 diagnostic has 0.759 bits of margin. No outward
M19 claim has been made. Lowering the distance target may strengthen that
sparse margin, but M20 is the current certified parameter.

## Historical 11% target

The current outer has no zero padding. Its parameters are

\[
 L=16384,\qquad k=64L=2^{20},\qquad N=128L=2^{21},
 \qquad D=\lceil0.11N\rceil=230687.
\]

The earlier 26.19-bit certificate used \(L=16560\). It does not certify this
power-of-two instance. A fresh certificate is required.

The RandomStepConv-M30 diagnostic does not close with the existing
parity-pivot relaxation. Its aggregate logarithmic upper bound is

\[
 \log_2 U=24670.227066088886.
\]

Occupation \(Q=L=16384\) is dominant. The first nonclosing occupation is
\(Q=16287\). Every occupation through \(Q=16286\) has a negative pointwise
bound in this diagnostic.

This failure belongs to the proof relaxation. It is not a counterexample to
the code or to RandomStepConv.

## Why removing 176 rows changes the proof

At exact rate one half, relative distance 11% is close to the random-linear
Gilbert--Varshamov crossing. For a uniform random binary linear code with
the same \((k,N,D)\), the elementary first moment has approximately
188.379 bits of margin.

The current dense proof replaces each permuted nontrivial BCH word by a
uniform even row. The worst-shell density factor costs

\[
 0.0674346406456\ldots
\]

bits per active row. At \(Q=L\), this pointwise comparison alone costs
approximately 1104.849 bits. That loss exceeds the complete random-code
margin before the coefficient relaxation is considered. Consequently,
changing the Chernoff grid or increasing the inner memory cannot repair the
dense proof.

This loss is highly pessimistic. After excluding the zero and all-one
endpoints, the exact BCH body has only 0.721 aggregate KL bits and 1.443
aggregate order-two Renyi bits relative to the uniform-even body at full
occupation. These divergence values do not prove the needed event bound,
but they identify exact shell averaging as the appropriate sharpening.

The 176 zero rows had supplied rate slack. They were not needed to define
the repeated BCH outer, but they made the pointwise spectrum comparison
affordable.

## Memory diagnostics at occupation one

Occupation one remains the sparse bottleneck for inner memory. The exact
BCH-shell diagnostic gives the following nearest-binary64 results.

| Memory | Diagnostic margin |
|---:|---:|
| 19 | 0.177 bits |
| 22 | 16.886 bits |
| 23 | 20.382 bits |
| 24 | 22.601 bits |
| 25 | 23.436 bits |
| 27 | 25.449 bits |
| 30 | 26.327 bits |

The memory-23 value uses a refined Chernoff grid around its optimum. These
values are diagnostics, not outward certificates. They establish a useful
proof target: memory 23 is the smallest tested value whose current
occupation-one analysis exceeds 20 bits. The memory-22 result does not prove
that memory 22 is impossible; it only fails the present 20-bit bound.

## Exact-complement and spectrum-mixture progress

The near-complement transfer computes each region coefficient by its number
of noncandidate positions. At full occupation and M30, this change reduces
the old pointwise bound from +24670.23 bits to +944.86 bits before replacing
the worst-shell spectrum factor. Replacing that factor by its baseline mass
would give -159.99 bits.

The positive-defect decomposition writes the BCH counting measure below a
uniform-even baseline of mass \(2^{64}\) plus a positive residual. The
residual-to-baseline mass ratio is

\[
 0.003906250000018967\ldots.
\]

For M23, the diagnostic now covers all original occupations with at most 128
inactive rows and all residual counts through 512. Their partial sum is

\[
 \log_2 U_{\rm boundary}=-66.9515958938.
\]

The dominant covered term has 63 residual rows and no original inactive
rows. The term with 512 residual rows is below \(2^{-966.9}\). Combining
this boundary sum with the refined occupation-one value and the old middle
transfer gives the partial all-message bound

\[
 \log_2 U_{\rm partial}=-20.3756340391.
\]

This is not yet a certificate. The positive-defect expansion still needs a
bound for more than 512 residual rows. Deleting every residual row becomes
too weak in that remote tail. Each residual row has support at least 24, so
the planned tail bound retains a uniformly selected subset of 24 support
positions instead of deleting the complete row.

M23 is therefore a plausible minimum, but its current diagnostic headroom is
only 0.376 bits after combining the covered ranges. M24 has a 22.601-bit
occupation-one diagnostic and is the safer first outward-certificate target.
After M24 closes, the proof can return to M23. M22 remains a later sharpening
target because its present occupation-one analysis reaches only 16.886 bits.

## Evidence

- M20 10.9% outward certificate:
  `ebch128_randomstepconv_g1_s20_pow2_outward_d109.json`.
- Hash-bound certificate manifest:
  `FINITE_K20_EBCH128_POW2_RANDOMSTEP_CONV_D109_MANIFEST.json`.
- M20 all-occupation and exact-complement witness receipts:
  `ebch128_randomstepconv_g1_s20_pow2_parity_pivot_d109_diagnostic.json` and
  `ebch128_randomstepconv_s20_pow2_dense_complement_r64_d109_u074.json`.
- M19 10.9% diagnostics:
  `ebch128_randomstepconv_g1_s19_pow2_q1_exact_d109_diagnostic.json`,
  `ebch128_randomstepconv_g1_s19_pow2_parity_pivot_d109_diagnostic.json`, and
  `ebch128_randomstepconv_s19_pow2_dense_complement_r64_d109_u074.json`.
- M30 all-occupation diagnostic:
  `ebch128_randomstepconv_g1_s30_pow2_parity_pivot_d11_diagnostic.json`.
- Random-code and exact-spectrum budget:
  `ebch128_pow2_spectrum_budget.json`, generated by
  `analyze_ebch128_pow2_spectrum_budget.py`.
- M19 occupation-one diagnostic:
  `ebch128_randomstepconv_g1_s19_pow2_q1_exact_d11_diagnostic.json`.
- M22 occupation-one diagnostic:
  `ebch128_randomstepconv_g1_s22_pow2_q1_exact_d11_diagnostic.json`.
- Refined M23 occupation-one diagnostic:
  `ebch128_randomstepconv_g1_s23_pow2_q1_exact_d11_diagnostic.json`.
- M25 and M27 occupation-one diagnostics:
  `ebch128_randomstepconv_g1_s25_pow2_q1_exact_d11_diagnostic.json` and
  `ebch128_randomstepconv_g1_s27_pow2_q1_exact_d11_diagnostic.json`.
- M24 occupation-one diagnostic:
  `ebch128_randomstepconv_g1_s24_pow2_q1_exact_d11_diagnostic.json`.
- M23 full old-transfer diagnostic:
  `ebch128_randomstepconv_g1_s23_pow2_parity_pivot_d11_diagnostic.json`.
- M23 boundary mixture through 128 original inactive rows and 512 residual
  rows:
  `ebch128_randomstepconv_s23_pow2_spectrum_mixture_r128_h512.json`.
- Exact-complement and spectrum-mixture evaluators:
  `evaluate_ebch128_randomstepconv_dense_complement.py` and
  `evaluate_ebch128_randomstepconv_spectrum_mixture.py`.
- Optimized M30 full-occupation complement diagnostic:
  `ebch128_randomstepconv_s30_pow2_dense_complement_qfull_u074.json`.

The evaluators accept exact rational distance arguments. The outward verifier
accepts the same length, memory, and distance parameters. The 10.9% receipt
above supersedes the earlier statement that no power-of-two receipt existed.
