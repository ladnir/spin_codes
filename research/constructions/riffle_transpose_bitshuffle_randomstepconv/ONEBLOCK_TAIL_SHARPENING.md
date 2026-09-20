# Sharpening the one-active-block tail

The complete bit-transpose certificates use a Chernoff bound on the inner
output weight. Every certified boundary is dominated by one active outer
block. This note measures the loss in that Chernoff bound.

Fix one active outer block. Let \(P=L/g\) be the number of convolution
positions in each transposed row. Let \(M_0(z)\) be the tilted transition at
an unmarked position. Let \(M_1(z)\) be the transition at the unique position
that contains the candidate bit. The candidate bit is uniform under the
random-outer relaxation.

The exact row transfer is

\[
 R(z):=\frac1P\sum_{j=0}^{P-1}M_0(z)^jM_1(z)M_0(z)^{P-1-j}.
\]

The probability generating function for the complete inner output weight
\(W\) is

\[
 F(z):=e_0^{\mathsf T}R(z)^B\mathbf 1.
\]

The script differentiates this expression to locate the saddle exactly in
binary64 arithmetic. It then estimates \(\Pr[W\le d]\) in two ways. The
lattice saddle includes the Gaussian and geometric tail prefactors. The
tilted Fourier method numerically inverts \(F\) after centering the weight at
\(d\).

## Validation

Three scaled instances admit an exact weight dynamic program. On those
instances, the lattice estimate exceeds the exact log probability by only
0.017 to 0.031 bits. The tilted Fourier method reproduces the exact \(g=8\)
tail to \(1.4\times10^{-14}\) bits.

At full size, transforms of lengths \(2^{20}\) and \(2^{21}\) give margins
that differ by at most \(6.4\times10^{-12}\) bits. The first cyclic alias lies
at least \(2^{20}\) weights beyond the target. Negative numerical mass is at
most \(2.2\times10^{-11}\) in the tilted distribution.

These checks establish strong numerical evidence. They are not an
outward-rounded certificate for the FFT computation.

## Revised numerical boundary

The following table compares the optimized Chernoff margin with the tilted
Fourier margin. Each margin includes the one-active-block outer multiplicity.

| \(B\) | \(g\) | \(\sigma\) | Chernoff \(\lambda\) | Fourier \(\lambda\) | All \(a\ge2\) terms |
|---:|---:|---:|---:|---:|---:|
| 256 | 4 | 17 | 37.2787 | 42.4010 | at most \(2^{-79.37}\) |
| 256 | 8 | 16 | 37.2820 | 42.4043 | at most \(2^{-79.38}\) |
| 512 | 4 | 14 | 34.9674 | 40.2973 | at most \(2^{-129.14}\) |
| 512 | 8 | 13 | 34.9959 | 40.3259 | at most \(2^{-129.22}\) |

Replacing the one-block term by the Fourier value changes the complete sum
by less than the displayed precision. The two-block term is the runner-up in
all four cases.

The next lower states do not approach 40 bits under the lattice estimate:

| \(B\) | \(g\) | \(\sigma\) | Lattice \(\lambda\) |
|---:|---:|---:|---:|
| 256 | 4 | 16 | 30.7564 |
| 256 | 8 | 15 | 30.7613 |
| 512 | 4 | 13 | -2.8793 |
| 512 | 8 | 12 | -2.8378 |
| 1024 | 4 | 12 | 3.3954 |
| 1024 | 8 | 11 | 3.5608 |

The numerical frontier therefore moves downward by exactly one state bit at
\(B=256\) and \(B=512\). It does not move at \(B=1024\). The one-state
advantage of \(g=8\) over \(g=4\) remains unchanged.

## Revised cost proxy

Use the dense-product proxy \(nB+(2n/g)(g+\sigma)^2\).

| \(B\) | \(g\) | \(\sigma\) | Total products |
|---:|---:|---:|---:|
| 1024 | 4 | 13 | 1,225,261,056 |
| 1024 | 8 | 12 | 1,178,599,424 |
| 512 | 4 | 14 | 706,740,224 |
| 512 | 8 | 13 | 652,476,416 |
| 256 | 4 | 17 | 499,646,464 |
| 256 | 8 | 16 | 419,430,400 |

The lowest numerical point remains \(B=256,g=8\). Sharpening the tail reduces
its proxy from 432,275,456 to 419,430,400 products.

## Scope and next proof step

The exact generating function and the full occupation calculations use the
same random-ensemble relaxation as the existing certificates. The new
limitation is numerical: the complex FFT has no outward-rounded error bound.

The cleanest certification target is \(B=256,g=8,\sigma=16\), which has 2.40
bits of numerical room. A rigorous tilted Fourier upper bound would combine
with the existing outward-rounded checker for occupations \(a\ge2\). The
\(B=512,g=8,\sigma=13\) point has only 0.326 bits of room and should follow
only after the wider-margin target closes.

Artifacts:

- `scripts/analyze_riffle_transpose_bitshuffle_oneblock_saddle.py`;
- `receipts/oneblock_fourier_g8_lower_cells.json`;
- `receipts/oneblock_fourier_g8_lower_cells_m2.json`;
- `receipts/oneblock_fourier_g4_lower_cells.json`;
- full-spectrum exploratory receipts for the four displayed cells.

Recommended next goal: build an outward-rounded upper bound for the tilted
Fourier inversion at \(B=256,g=8,\sigma=16\), then certify the sum of all
occupations \(a\ge2\) with the existing checker.
