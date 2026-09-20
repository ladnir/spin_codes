# Riffle BCHPerm-TransposeBitShuffle-SplitState

**Status: diagnostic coset exploration.** The established
FieldCheckpoint (K=8) inner is both cheaper and further along in the proof.
The SplitState calculation reaches the 9% target only after assuming an
averaged bound on nontrivial affine cosets. That bound is not yet proved for
the fixed constituent.

## Construction

The outer encoder applies one binary ([256,128]) block code to each message
block. The current spectrum model has minimum distance 38. The encoder
independently permutes the 256 output coordinates of every outer block.

The encoder transposes the outer output into 256 regions. Each region contains
one coordinate from every outer block. The encoder independently permutes all
bits in each region and divides the result into 256-bit epochs.

The inner state has 64 bits. One epoch uses the fixed maps

\[
A:\mathbb F_2^{64}\to\mathbb F_2^{256},
\qquad
B:\mathbb F_2^{256}\to\mathbb F_2^{64},
\qquad BA=0.
\]

For input (X) and state (Q), the epoch outputs (X+A(Q)). It compresses
the output with (B), combines the result with (Q), and randomizes every
nonzero next state by a fixed nonzero field scalar for that epoch.

The recorded pair has a sparse nested implementation. Its local audit gives
4.2422 XORs per output bit for (A), 2.25 for (B), and 1 for combining
the input with (A(Q)). The 64-bit state randomizer raises the current logical
proxy to 10.6147 XORs per output bit. This proxy is not a wall-clock result.

## Comparison at the 9% target

The comparison uses output length (2^{21}) and the same modeled outer
spectrum. A negative pointwise exponent is a successful union bound.

| Inner model | Regular margin | Exponent at 1,024 active blocks | Logical XOR proxy |
|---|---:|---:|---:|
| FieldCheckpoint (K=32), 1,024-bit epoch | not closed | +16,038.90 bits | 1.7806 |
| FieldCheckpoint (K=8), 256-bit epoch | 55.9507 bits | -30,347.96 bits | 4.1226 |
| SplitState, fixed coset bound | not closed | +5,549.82 bits | 10.6147 |
| SplitState, combined coset envelope | 79.3952 bits | -67,523.98 bits | 10.6147 |

The fixed SplitState analysis improves on FieldCheckpoint (K=32), but that
is the wrong performance frontier. FieldCheckpoint (K=8) already closes
the same occupation and every other regular occupation. It also uses about
39% of the SplitState logical XOR proxy.

The outer density envelope is loose for one active block. The final regular
diagnostic therefore uses the exact modeled spectrum at occupation one. It
uses the density envelope from occupation two onward. Occupations two and
three have 126.20 and 94.88 bits of pointwise margin. The sum over all regular
occupations then has 79.3952 bits of floating-point margin.

Both regular calculations exclude configurations containing the special
all-one outer word. The SplitState result additionally depends on the
combined coset envelope described below.

## SplitState parameter choice

The current landscape supports (B=256) as the outer block size.

- At (B=128), the audited (t=256,s=64) inner has only 39.12 bits of
  one-active margin. Aggregating the remaining cases can only reduce it.
- At (B=256), the one-active margin is 79.40 bits. The combined envelope
  also closes all regular higher occupations.
- At (B=512), the one-active margin rises to about 149.87 bits, but the
  larger outer block is unnecessary for the 40-bit target.

Among the explored SplitState geometries, (t=128,s=32) has the best logical cost
proxy and 80.50 bits of one-active margin at (B=256). That point is only an
extrapolation. No matching activation table or fixed nested constituent has
been audited. Neither SplitState point currently displaces FieldCheckpoint
(K=8).

## Missing coset statement

The fixed proof treats every nonzero epoch input as an adversarial shift of
the code (A(\mathbb F_2^{64})). The resulting exceptional term is about
(2^{-64}) per epoch. This term caps the dense-regime exponent and prevents
an end-to-end proof, regardless of the outer block size.

The bit shuffle gives a stronger probability space. Conditioned on a regular
outer occupation, the candidate positions in each epoch form a uniform
subset, and their values are independent fair bits under the density
envelope. The required lemma must average each affine coset over this input
law. For every Chernoff parameter used by the certificate, the lemma should
upper-bound the shifted live-state moment by the smaller of these bounds:

1. the existing distance-40 and three-moment bound; and
2. the corresponding random-like affine-coset moment.

The saved combined-envelope calculation assumes this statement. It does not
assume a fresh random matrix during encoding.

## Evidence

- `../riffle_bchperm_transpose_bitshuffle_structuredstepconv/receipts/inner_global_parameter_landscape.json`
- `../riffle_bchperm_transpose_bitshuffle_structuredstepconv/receipts/splitstate_regular_bulk_comparison.json`
- `../riffle_bchperm_transpose_bitshuffle_structuredstepconv/receipts/splitstate_combined_B256_low_occupations.json`
- `../riffle_bchperm_transpose_bitshuffle_structuredstepconv/receipts/splitstate_combined_B256_checkpoints.json`
- `../riffle_bchperm_transpose_bitshuffle_structuredstepconv/receipts/splitstate_hybridcoset_B256_full.json`
- `../riffle_ldpcsplitstate_g4_t256_s64/receipts/zero_state_activation_table.json`

All reported distance values are nearest-binary64 diagnostics. The modeled
outer spectrum, the averaged coset statement, the all-one cases, and outward
rounding remain outside the current certificate.
