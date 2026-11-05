# Smaller outer groups: a bounded B128 screen

The tested fixed witness does not close the middle-weight bound for the
smaller outer. At the new active-group weight 2048, the FLOAT log2 bound
is **+4845.925682130852**. The sparse screens at weights one and two are
below -40 bits. These are numerical proposal screens, not directed
certificates. In particular, the positive middle-weight bound is neither
a codeword counterexample nor an impossibility result for this construction.

No native implementation, benchmark, witness fit, or full coverage search
was performed. The existing B256 construction and proof defaults are
unchanged. The recommended next step is to finish the B256 certificate;
if B128 is revisited, improve this single middle-weight gate before
starting native work or a broad proof search.

## Construction and measure

Fix the byte field of order q=256, source dimension K=2^20, and output
length N=2^21. Each proposed group contains four RS[32,16] constituents,
so it maps 64 input symbols to 128 output symbols. Independent uniformly
sampled invertible Toeplitz4 matrices mix the four symbols in each output
column. There are 16384 groups and 128 regions. The route uses alternating
adjacent-pair matchings with an independent fair swap for each routed
pair in each region. The 64-input, 16-state, h4 inner, initial zero state,
and discarded terminal state are unchanged.

For a fixed nonzero target group, let h be its number of nonzero columns.
If h<17, the MDS distance makes membership impossible. Otherwise, condition
on the mixers at h-16 nonzero target columns. These columns and the 32-h
zero columns specify 16 positions of each constituent RS word, hence at
most one source word. At each remaining nonzero target column, an
independent mixer hits its required value with probability at most

\[
\rho_T=\frac{q^{-4}}{1-1/q}.
\]

Thus pointwise membership is at most \(\rho_T^{16}\). Relative to the
uniform distribution on 128 output coordinates, the group density is

\[
\beta_{128}=q^{128}\rho_T^{16}
            =q^{64}(1-1/q)^{-16}
            =256^{80}/255^{16}.
\]

Its log2 value is 512.0903450102583, half the sharper B256 Toeplitz
density. Injectivity separates the zero input from nonzero groups;
independent group setup gives the usual product domination. The scalar
projective factor remains 1/255. Any extension-field transfer relies on
the underlying GF256-linear map, not transitivity on arbitrary vectors
over the extension field.

Each new region has 8192 routed pairs and 256 inner packets. Let T be
the unchanged positive two-packet polynomial, including the good-good
composition before the fractional power. The regional polynomial is
now T^128, rather than T^64. If R0 and R1 are the resulting normalized
regional matrices for the ordered matchings, their product is
(R0 R1)^64, rather than (R0 R1)^128. The total number of packets remains
32768. A common-potential bound includes its initial-state factor once,
not once per region.

At matched active-coordinate density, new weight 2048 corresponds to old
weight 1024. The outer prefix is identical with the sharper Toeplitz
density, but the support count grows: the difference between
log2 binomial(16384,2048) and log2 binomial(8192,1024) is about 4452.38
bits. The regional matrices also change, so this entropy difference
cannot simply be added to an old final bound.

## The bounded screens

The isolated helper is [hd_shifted_b128_screen.py](theory/hd_shifted_b128_screen.py).
It uses a distinct construction identifier ending in
`swapped-pairs-b128-k20-rate-half` and records its geometry in every result.

For sparse weights, positive polynomial multiplication is truncated at
input weight two. This coefficient truncation is exact algebraically;
its numerical evaluation is FLOAT. For an even cycle of G groups, the
joint matching counts are G for the sole weight-one class; at weight
two they are binomial(G,2)-G for (a,b)=(0,0), and G/2 each for (1,0) and
(0,1). Both ordered classes are retained because the initial zero state
does not justify swapping noncommuting regional operators.

| Active groups | z (exact saved literal) | alpha | FLOAT log2 bound |
|---:|---|---:|---:|
| 1 | 2250852543970203/2251799813685248 | 1 | -47.57606805706871 |
| 2 | 8999622690980997/9007199254740992 | 1 | -96.99029274927628 |

These predetermined z values halve the old singleton log(z) for new
weight one and retain it for new weight two. Neither witness was fitted.

The sole medium screen reuses the strong frozen B256 weight-1024 bank:
22 exact marker pairs, grid 128, and the same z, alpha, and common state
potential. It rebuilds the regional tables with exponent 128, normalizes
each of the 1025 full-pair histograms using 8192 pairs, and takes
entrywise minima only after normalization. The cyclic-run sum includes
every matching profile and every even/odd split at new weight 2048.

The resulting support-sum logarithm is -179294.2070066937 and the outer
prefix is 184140.13268882455, giving +4845.925682130852. For comparison,
the old weight-1024 bank gives -165.86870150166214 with the inherited
diagonal density, or -222.03569227065245 after the sharper Toeplitz
prefix. The new bound is worse by about 5067.96 bits than the latter
control. The FFT screen uses a normalized numerical cushion; it is not
a directed error estimate. No nearby weights are covered by this test.

## Cost and reproduction

A source-level operation count for the current two-coset, fused-last
schedule estimates 245760 fewer RS GFNI products, a reduction from
1064960 to 819200. Mixer work stays at 1179648 GFNI products, and the
prepared hot keys remain 4.5 MiB. The combined outer count therefore
falls by about 10.95%. Payload traffic and inner work are unchanged;
halving the region count doubles each region's length and stride.
These are operation counts, not a speed prediction.

The inspected g8/RS[32,16] experiments retain 256-symbol groups, while
older B128 direct-GRS studies use different outers, routes, or inners.
Neither is a native test of this exact g4 byte-field candidate.

From `theory`, the following commands create new artifacts exclusively;
use unused output names:

```text
python -B hd_shifted_b128_screen.py --output NEW_b128_sparse.json
python -B hd_shifted_b128_screen.py --medium-source hd_shifted_potential_band_1024_frontier.json --output NEW_b128_medium2048.json
python -B -m unittest hd_shifted_b128_screen_tests hd_shifted_swap_sparse_tests hd_shifted_run_count_screen_tests
```

The last command passed 26 tests, including nine B128 tests. They check
small-cycle exhaustive counts, unchanged old-dimension normalization,
alpha-one paired/unpaired equality, ordered initial-state composition,
the two regional exponents, and the one-weight/no-fit medium guard.
Optional artifact checks skip when receipts are absent; unit tests do
not regenerate the medium bank.

Frozen artifact SHA-256 hashes:

- `theory/hd_shifted_b128_sparse_screen.json`:
  `0569f81147c9a88eb5dc10e55faf2b68fca75408b31fe6cc8a954b31d4a37796`.
- `theory/hd_shifted_b128_medium2048_screen.json`:
  `74525a7bc1094f023ef0c8ef0f1a4d165ef78c60c5994bbf797284873f559ae8`.
- Its source, `theory/hd_shifted_potential_band_1024_frontier.json`:
  `23aee6ac737c2d2a2d5c3d4205775f27bcd1da91b537bdff03748f56a1949d99`.
