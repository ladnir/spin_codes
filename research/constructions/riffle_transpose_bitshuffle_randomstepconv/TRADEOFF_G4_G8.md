# Certified bit-transpose tradeoff for g=4 and g=8

Fix message length \(n=2^{20}\), output length \(2^{21}\), and distance
threshold 188743. The outer code is sampled independently for every rate-half
block. Every transposed bit row is independently shuffled before consecutive
bits are formed into packets.

| Outer length \(B\) | \(g=4\) state \(\sigma\) | Certified \(\lambda\) | \(g=8\) state \(\sigma\) | Certified \(\lambda\) |
|---:|---:|---:|---:|---:|
| 1024 | 13 | 83.087891407589 | 12 | 83.201070356493 |
| 512 | 15 | 57.641932041422 | 14 | 57.659313166668 |
| 256 | 18 | 45.400358851794 | 17 | 45.402590000096 |

The \(g=4\) columns inherit the complete packed-support certificates: packing
candidate bits into fixed packets maximizes the low-output probability, so
the same upper bound applies to a fresh bit shuffle. The exact bit-shuffle
enumerator independently reproduces all three margins. The \(g=8\) columns
use direct outward-rounded coefficients for every active-block occupation.

In all six cells, one active outer block is dominant. Consequently the two
packet widths have almost identical margin at their boundary. Increasing
\(g\) from four to eight saves exactly one state bit at each tested outer
length.

## Sharper numerical frontier

`ONEBLOCK_TAIL_SHARPENING.md` replaces the dominant Chernoff term with a
tilted Fourier inversion of the exact one-block generating function. This
numerical calculation supports one smaller state at \(B=256\) and \(B=512\):

| \(B\) | \(g=4\) numerical \(\sigma,\lambda\) | \(g=8\) numerical \(\sigma,\lambda\) |
|---:|---:|---:|
| 1024 | \(13,83.09\), certified above | \(12,83.20\), certified above |
| 512 | \(14,40.2973\) | \(13,40.3259\) |
| 256 | \(17,42.4010\) | \(16,42.4043\) |

These four new cells are not yet certificates because the complex FFT lacks
outward-rounded error control. The occupations with at least two active
blocks remain far below the one-block term.

## Why packet-transpose g=8 behaves differently

Fixed packet groups allow many active outer blocks to be packed into the same
groups in every transposed row. At \(B=1024,g=8,\sigma=12\), this produces a
bulk obstruction near 424 active outer blocks and an exploratory exponent of
about -49417 bits. The first certified packet-transpose point at this outer
length is only \(\sigma=102\), with margin 140.2893 bits. Its dominant event
is still the 424-block bulk event. At \(\sigma=103\), dominance finally moves
to one active block.

Bit transpose redraws those collisions independently in every row. That
removes the persistent packed-group event: already at \(B=1024,\sigma=12\),
the complete bit-transpose certificate is positive and one active block is
dominant.

## Dense-operation proxy

Count each dense binary matrix-vector bit product as one unit. The outer cost
is \(nB\). The inner has \(2n/g\) steps and costs \((g+\sigma)^2\) products per
step.

| \(B\) | \(g\) | \(\sigma\) | Outer products | Inner products | Total proxy |
|---:|---:|---:|---:|---:|---:|
| 1024 | 4 | 13 | 1,073,741,824 | 151,519,232 | 1,225,261,056 |
| 1024 | 8 | 12 | 1,073,741,824 | 104,857,600 | 1,178,599,424 |
| 512 | 4 | 15 | 536,870,912 | 189,267,968 | 726,138,880 |
| 512 | 8 | 14 | 536,870,912 | 126,877,696 | 663,748,608 |
| 256 | 4 | 18 | 268,435,456 | 253,755,392 | 522,190,848 |
| 256 | 8 | 17 | 268,435,456 | 163,840,000 | 432,275,456 |

Under this proxy, bit-transpose \(g=8\) wins at all three outer lengths. The
advantage grows from about 3.8% at \(B=1024\) to about 17.2% at \(B=256\).
The proxy excludes setup expansion, memory traffic, and bit-permutation cost.
Bit transpose shuffles \(2n\) individual items in total, independent of
\(g\); packet transpose shuffles only \(2n/g\) packet items.

## Receipts and checkers

- `receipts/g8_b1024_sigma12_full_interval.json`;
- `receipts/g8_b512_sigma14_full_interval.json`;
- `receipts/g8_b256_sigma17_full_interval.json`;
- `scripts/analyze_riffle_transpose_bitshuffle.py`;
- `scripts/certify_riffle_transpose_bitshuffle_full.py`.

Recommended next goal: certify the tilted Fourier tail at
\(g=8,B=256,\sigma=16\). Then benchmark the permutation and inner kernels for
\(g=4\) and \(g=8\). The proof proxy favors \(g=8\), but the extra bit-shuffle
traffic could change the implementation winner.
