# Finite 11% route with one repeated BCH subcode

The optimized outward certificate is complete for the proof-model ensemble.
See `FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_CERTIFICATE.md`. The earlier
256-layer certificate remains valid but is no longer the preferred instance.

## Result

The construction has

\[
 (B,K_B,L,N,D)=(250,124,8576,2{,}144{,}000,235{,}840).
\]

It repeats one fixed \([250,124,\ge38]\) BCH subcode. Each row receives an
independent composition of 56 ParityFanout-31x33 layers and an independent
local coordinate permutation. The construction retains the independent
region permutations and RM2Sub-S19 multipliers.

Fixing 14,848 parent input coordinates gives exactly \(2^{20}\) message bits.
The certified distance-failure probability is below \(2^{-40}\). The
displayed receipt endpoints combine to about 41.51 bits.

## Proof decomposition

Let \(G\) be the event that no wrapped row code contains a nonzero word of
weight at most 12 or at least 238. The outward receipts prove:

| contribution | strict bound | displayed margin |
|---|---:|---:|
| \(G^c\) | \(2^{-45}\) | 45.238916 bits |
| occupation \(Q=1\) | \(2^{-53}\) | 53.947385 bits |
| occupations \(2\le Q\le31\) | \(2^{-52}\) | 53.044345 bits |
| occupations \(32\le Q\le8576\), on \(G\) | \(2^{-41}\) | 41.620995 bits |

On central weights 13 through 237, the expected per-vector multiplicity is
strictly below \(2^{389/1250}\) times the uniform-injection multiplicity.
This shell-sensitive split replaces the former global 0.001-bit comparison.

The coarse strict bounds satisfy

\[
 2^{-45}+2^{-53}+2^{-52}+2^{-41}<2^{-40}.
\]

## Linear-time and operation count

The proof-model encoder remains linear-time. The unfused fanout proxy is

\[
 56(31+33-1)=3528
\]

scalar XORs per row, or 30,256,128 over all rows. The 256-layer proof used
138,313,728 scalar XORs. The optimized proof removes 78.125% of that proxy,
a factor of approximately 4.57.

The 52-layer candidate fails the present tail-plus-central pointwise proof.
This does not prove that 56 layers are necessary. A multi-band argument could
possibly reduce the count.

## Reproducible artifacts

- `bch250_124_parityfanout31x33_l56_cutoff12_spectrum_outward.json`
- `bch250_124_parityfanout31x33_l56_q1_outward_d11.json`
- `bch250_124_parityfanout31x33_l56_q2_31_outward_d11.json`
- `bch250_124_parityfanout31x33_l56_cutoff12_dense_q32_8576_outward_d11.json`
- `FINITE_K20_BCH250_124_ROWLOCAL_FANOUT56_MANIFEST.json`

## Optimized BCH256 performance proxy

The active proxy now preserves the frozen 256-by-8192 packed, tiled, and fused
transpose encoder. It changes only the fanout schedule multiplicity and the
number of layers applied in the hot outer loop.

On one pinned Ryzen 9 7950X logical processor, the zero-layer path takes
10.650813 ms. The 56-layer path takes 20.454284 ms. Fanout-56 therefore adds
9.803471 ms by the difference of medians and multiplies total time by
1.920443. Both paths pass staged-output checks outside the fused method.

This measurement is not the exact BCH250 encoder. The proxy has 8,192 rows,
whereas the proof model has 8,576 rows. Linear row-count scaling estimates a
10.263009-ms fanout increment for the exact row count, but this estimate is
not a benchmark. See `FINITE_K20_BCH256_FANOUT56_PROXY_BENCHMARK.md`.

## Remaining implementation obligations

1. Replace the benchmark schedule generator with a proof-compatible setup
   sampler or authenticate a table produced by such a sampler.
2. Implement the ordinary forward encoder.
3. Prove complete ordinary, transposed, and proof-model equivalence.
4. Implement and measure the exact BCH250 path on the target host if the
   proxy remains performance-competitive.

The one-bit-plane mask kernel remains documented in
`FINITE_K20_BCH250_FANOUT56_IMPLEMENTATION_BENCHMARK.md`. On one pinned
i7-13700H logical processor, it adds 2.171294 ms over the post-BCH zero-layer
copy control. This microbenchmark is not the production 128-bit-block cost.

This is a finite \(k=2^{20}\) theorem. It does not give an asymptotic theorem
with fixed outer length 250.
