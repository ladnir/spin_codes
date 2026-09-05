# Proof status

This variant is retained as the state-size threshold and endpoint-obstruction
record. The current deployment candidate is the distinct
`Riffle ParityFanout-31x33 s=19` construction in the adjacent folder. Its odd
source and target sets remove the exceptional all-one endpoint, and its full
floating-point ledger reaches 52.716450 bits at relative distance 0.11.

## Minimum-state conclusion

Within the current RM2Sub inner family and the recorded ParityShear-12
outer model, the smallest state size that passes the regular-class
calculation at relative distance 0.11 is

\[
s=19.
\]

The state-size-19 calculation has 46.329381 bits of regular-class margin.
State size 18 is separated from the threshold by a concrete failure: at
occupation 3100, its pointwise margin is -7567.8584 bits. This is not a
rounding-scale miss and does not justify searching smaller state sizes in
the same bound.

This minimum does not apply to the unmodified outer spectrum. That
spectrum loses one envelope bit per active outer block. Even an
ideal-random inner at state sizes 20 and 64 then misses the 0.11 target by
about 8017.565 bits, so increasing the state size alone cannot repair the
recorded bound.

The conclusion is a minimum for the current modeled family and numerical
certificate, not yet an unconditional end-to-end theorem. The remaining
gates listed below still apply.

## Regular-class claim

Assume the modeled even-floor spectrum for the source [256,128] outer code.
Sample each ParityShear-12 map independently across outer-block instances.
Also sample the coordinate permutations, region permutations, and epoch
scalars independently.

For output length (N=2^{21}) and distance

\[
D=\lfloor0.11N\rfloor=230686,
\]

the current recurrence bounds the expected number of nonzero regular
messages of output weight at most (D) by approximately (2^{-46.3294}).
Consequently, subject to the model and numerical limitations below, some
fixed setup has no such regular message.

## Outer shear

For each block, setup samples (j\gets\{1,\ldots,256\}) and then samples a
12-subset (S\) of the remaining coordinates. The map adds the parity on
(S) to coordinate (j).

For a source word of weight (w), the probability of changing its weight
is an exact hypergeometric parity probability. Averaging these transitions
over the source spectrum gives the recorded transformed spectrum. Its
Bernoulli-envelope exponent is 128.015626884 bits per active block. The
unmodified even spectrum requires 129 bits per active block.

The ten-input shear has exponent 128.034617921. Across 8192 active blocks,
that residual loss causes a margin of (-109.2452) bits. The twelve-input
shear reduces the same loss enough to give 46.3294 bits.

## Fixed inner constituent

The (s=19) selector exactly enumerated 256 candidate image spectra and
selected seed 2719678504.

- The image of (A) has mass (2^{19}) and minimum distance 48.
- All coordinate forms are nonzero and distinct.
- No three coordinate forms are dependent.
- The first three factorial moments match the coordinate identities.
- The degree argument proves (BA=0).
- The kernel of (B) has minimum distance six and no weight-four word.

## Numerical coverage

- Occupation one: 59.485911 bits.
- Occupations 2 through 100: 108.806146 bits aggregate.
- Occupations 101 through 8192: 46.329381 bits aggregate.
- Pure all-one checkpoints: at least 2388.957424 bits.

State size 18 is refuted by occupation 3100, whose pointwise margin is
(-7567.8584) bits on the recorded grid. State size 19 gives thousands of
bits on the corresponding medium-occupation probe.

## Remaining gates

1. Replace the modeled outer spectrum with an explicit constituent or a
   proved spectrum bound.
2. Cover general mixtures containing at least two all-one outer words.
3. Replace nearest binary64 evaluation by outward-rounded arithmetic.
4. Benchmark and optimize the integrated (s=19) transposed kernel.

## Receipts

The receipts currently live in the adjacent state-search folder.

- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/parity_break_outputsubset12_outer256_d38_expected_spectrum.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/s19_rm2sub_selection.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/s19_rm2sub_a_spectrum_audit.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/s19_rm2sub_b_kernel_spectrum.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/parity_break_k12_s19_one_active_d11.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/parity_break_k12_s19_regular_2_100_d11.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/parity_break_k12_s19_regular_101_8192_d11.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/parity_break_s18_medium_probe_d11.json`
- `../riffle_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s20/receipts/min_state/parity_break_k10_s19_all_active_d11.json`
