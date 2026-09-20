# Proof status

## Current claim

Fix the selected maps (A) and (B) recorded in this folder. Assume the
modeled even-floor spectrum for a binary [256,128] outer code. Setup samples
the outer coordinate permutations, the region permutations, and the epoch
scalars independently.

For (N=2^{21}) and (D=228023), the current numerical recurrence bounds
the expected number of nonzero regular messages with output weight at most
(D) by approximately (2^{-42.7641}). Consequently, subject to the model
and numerical limitations below, some fixed setup has no such regular
message.

The next integer distance does not meet the 40-bit target. At (D=228024),
the all-active regular occupation has 39.7251 bits of margin.

## No finite-state closure at 11 percent

Consider occupation 8192 under the current Bernoulli outer-spectrum
envelope. Conditioned on the history before epoch (i), the epoch input
(X_i) is uniform in (mathbb F_2^{128}). The state (Q_i) and therefore
(A(Q_i)) are fixed under this conditioning. Hence

\[
Y_i=X_i+A(Q_i)
\]

is uniform in (mathbb F_2^{128}). This conclusion does not depend on the
state size, the maps (A,B), or the state update. Repeating the argument
over the epochs shows that increasing (s) cannot improve the all-active
output law used by this bound.

The modeled even outer spectrum has a Bernoulli-envelope mass of exactly
(2^{129}) per active block. An unconstrained [256,128] reference would
contribute (2^{128}). The current analysis therefore pays one extra bit
for each of the 8192 active blocks.

At 11 percent, an ideal random-coset inner gives the following margins.

- At (s=20), the all-active margin is (-8017.565021) bits.
- At (s=64), the all-active margin is (-8017.565020) bits.

Thus no finite (s) reaches 40 bits at 11 percent under the current outer
envelope. The obstruction is the lost even-parity information, not the
inner state.

## Fixed constituent audit

The selection receipt records seed 2719617393. Exact enumeration proves the
following properties.

- The image of (A) has mass (2^{20}) and minimum distance 48.
- All 128 coordinate forms are nonzero and distinct.
- No three coordinate forms are dependent.
- The first three factorial moments match the coordinate identities.
- The degree argument proves (BA=0).

The exact MacWilliams transform of the image spectrum gives the complete
kernel spectrum of (B). The kernel has minimum distance six and contains
no word of weight four.

## Numerical coverage

The regular calculation is split only to refine the Chernoff grids.

- Occupation one has aggregate margin 61.637943729 bits.
- Occupations 2 through 100 have aggregate margin 108.205500 bits.
- Occupations 101 through 8192 have aggregate margin 42.764071 bits.
- Occupation 8192 is the dominant class in the last range.

The union of these three ranges remains dominated by the final value to the
reported precision.

## Scope

The claim is conditional on the modeled outer spectrum. The calculation
uses nearest binary64 arithmetic rather than outward-rounded intervals. It
does not cover the outer all-one word or mixtures containing that word.
Those classes require the existing separate tail analyses to be rerun for
this state size and distance.

## Receipts

- `receipts/s20_rm2sub_spectrum_search_1024_selection.json`
- `receipts/s20_rm2sub_spectrum_search_1024_a_histogram.csv`
- `receipts/s20_rm2sub_spectrum_search_1024_a_spectrum_audit.json`
- `receipts/s20_rm2sub_spectrum_search_1024_b_kernel_spectrum.json`
- `receipts/s20_rm2sub_best_exact_one_active_D228023.json`
- `receipts/s20_rm2sub_best_exact_regular_2_100_D228023.json`
- `receipts/s20_rm2sub_best_exact_regular_101_8192_D228023_direct.json`
- `receipts/s20_rm2sub_best_all_active_D228024.json`
- `receipts/ideal_coset_s20_all_active_d11.json`
- `receipts/ideal_coset_s64_all_active_d11.json`
