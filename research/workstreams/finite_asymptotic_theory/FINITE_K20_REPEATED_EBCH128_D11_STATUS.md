# Repeated EBCH128 at finite \(k=2^{20}\): 11% status

## Decision

The repeated-EBCH128 candidate does not yet have a complete 11% distance
certificate. The current evidence closes the low-occupation endpoint only in
nearest-binary64 arithmetic. A separate high-occupation lemma remains
necessary. This gap is mathematical; outward rounding cannot repair it.

The earlier description of the 26.71-bit occupation-one result as a complete
certificate was incorrect. The number applies only to messages supported on
one outer block.

## Candidate and target event

Fix one full-support binary linear \([128,64,22]\) extended BCH code
\(C\). Use the same code \(C\) in every outer position. Set

\[
  k=2^{20},\qquad L=k/64=16384,\qquad N=128L=2^{21},
\]

and

\[
  D=\lfloor0.11N\rfloor=230686.
\]

For each outer position, sample an independent uniform permutation of the
128 code coordinates. For each transposed region, sample an independent
uniform permutation of its \(L\) positions. For each RM2Sub-S19 epoch, sample
an independent uniform multiplier in \(\mathbb F_{2^{19}}^*\). The outer code
is fixed; it is not resampled at different positions.

Let \(\mathcal C(\omega)\) be the resulting rate-one-half code, where
\(\omega\) denotes only the sampled permutations and multipliers. The desired
finite statement is

\[
 \Pr_\omega\!\left[
   \exists m\ne0:\operatorname{wt}(\mathcal C(\omega)m)\le D
 \right]\le 2^{-\lambda}.
 \tag{1}
\]

Any positive \(\lambda\) proves that some setup has minimum distance at least
\(D+1=230687>0.11N\). A cryptographic setup statement should retain a
material margin. The present low-occupation calculation suggests
\(\lambda\approx26\), but it does not establish (1).

## Imported exact outer data

The exact spectrum is stored in `scripts/EBCH128_64.wd`. Its coefficients
sum to \(2^{64}\), its first nonzero weight is 22, and it contains the unique
all-one word. The spectrum is complement symmetric.

For every nonzero, non-all-one shell, define

\[
 C_{\rm body}:=
 \max_{0<w<128}
 \frac{A_w2^{127}}{\binom{128}{w}}.
 \tag{2}
\]

The exact integer ratios in (2) give

\[
 \log_2 C_{\rm body}\approx64.0674346406.
 \tag{3}
\]

The maximizing shells are 24 and 104. Equation (2) pointwise dominates the
permuted BCH counting measure by \(C_{\rm body}\) times the uniform-even row
law.

## Low occupation

The three-state occupation-one calculation sums every exact BCH shell. Its
nearest-binary64 result is

\[
 \mathbb E[Z_{D,1}]\lesssim2^{-26.7140676240}.
 \tag{4}
\]

The dominant shell is weight 22. The receipt is
`ebch128_direct_rm2sub_s19_k20_d11_one_active_three_state_diagnostic.json`.

For occupations 2 through 100, an exact uniform-shell region recurrence and
the pointwise exact-spectrum envelope give

\[
 \sum_{Q=2}^{100}\mathbb E[Z_{D,Q}]
 \lesssim2^{-28.4502917457}.
 \tag{5}
\]

The dominant occupation in (5) is \(Q=2\). The receipt is
`ebch128_direct_rm2sub_s19_k20_d11_q2_100_diagnostic.json`.

Combining (4) and (5) gives the diagnostic low-occupation margin

\[
 -\log_2\!\left(
  2^{-26.7140676240}+2^{-28.4502917457}
 \right)
 \approx26.3353848579.
 \tag{6}
\]

Equations (4)--(6) are not outward certificates.

## Why the even-row reduction is insufficient

A uniform even row has 127 independent fair coordinates; its final
coordinate is their parity. Omitting the parity-determined transposed region
leaves 127 independent regions. Since a Chernoff variable \(z\) lies in
\((0,1)\), omitting output coordinates gives a valid upper bound.

This reduction removes the former one-bit conditioning charge per active row.
It works well through most of the occupation range. Sampled binary64 margins
are positive from \(Q=100\) through \(Q=16100\). The same calculation has
margin about \(-4853\) bits at \(Q=16200\) and about \(-14545\) bits at
\(Q=L\). See `ebch128_even_body_dense_high_k20_d11_diagnostic.json`.

This endpoint failure is expected. At rate one half and relative distance
11%, the ideal-random-code first-moment margin at this exact length is only
about 178 bits. Dropping one region discards 16384 output coordinates. The
resulting loss is much larger than any rounding allowance.

## Rejected high-occupation shortcut

At \(Q=L\), a uniformly selected nonzero word of any full-support linear
\([128,64]\) code has one-coordinate probability

\[
 p=\frac{2^{63}}{2^{64}-1}.
\]

One may compute the exact uniform-shell transfer for one region, apply a
common positive weighted norm, and use Hölder's inequality across the 128
dependent regions. This route uses the exact one-coordinate marginal and
does not charge (3). Nevertheless, the Hölder moment is dominated by rare
column weights. The binary64 witness at \(z=0.11\) fails by approximately
1.77 million bits. The receipt is
`ebch128_rm2sub_s19_k20_d11_qL_column_holder_diagnostic.json`.

The column-Hölder route is therefore rejected. Witness refinement and
outward rounding cannot recover its loss.

## Exact remaining proof obligation

A full proof of (1) requires a nonnegative 128-region transfer that preserves
the dependence introduced by the BCH row law. It must satisfy all of the
following conditions.

1. It uses the exact repeated-code row distribution, or a positive majorant
   whose total excess is small enough for the approximately 178-bit ideal
   endpoint budget.
2. It retains the parity-determined region. A proof that simply discards one
   region cannot close the all-active endpoint.
3. It covers every occupation \(101\le Q\le16384\), including mixtures of
   ordinary BCH words and the unique all-one word.
4. Its RM2Sub transfer is valid for correlated region inputs and arbitrary
   entering-state distributions. The fresh nonzero multiplier may be used,
   but independence between the parity region and the entering state may not
   be assumed without proof.
5. A final verifier evaluates the complete occupation sum with outward
   arithmetic and binds the BCH spectrum and frozen RM2Sub receipts by hash.

Two plausible implementations of this obligation remain.

- Track the exact BCH shell generating polynomial through all 128 regions,
  using a finite state that preserves the row-parity dependence.
- Prove a parity-region decoupling lemma for RM2Sub-S19, then append a sharp
  final-region transfer to the existing 127-region body calculation.

Neither lemma is currently proved.

## Proof status

- **Exact imported fact:** the repeated outer is the published
  \([128,64,22]\) extended BCH code with the supplied full spectrum.
- **Algebraically valid reduction:** pointwise domination by (2), and omission
  of one parity-determined region for \(z\in(0,1)\).
- **Diagnostic:** the 26.71-bit occupation-one margin, the 28.45-bit aggregate
  margin for occupations 2 through 100, and the sampled dense margins.
- **Rejected:** a common-norm column-Hölder proof at high occupation.
- **Open:** the complete 11% finite certificate, outward verification, the
  all-one mixed tail, and the exact high-occupation dependence lemma.

## Recommendation

Do not describe this candidate as having a certified 11% distance. Preserve
the approximately 26-bit low-occupation result as a witness for a future
exact high-occupation transfer. Resume only if the repeated EBCH128
implementation remains competitive enough to justify that proof work.
