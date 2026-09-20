# Goal 01 results

Fix \(n=2^{20}\), \(B=256\), and output distance threshold 188743. Goal 01
evaluates messages that activate exactly one outer block.

The row-weight coefficient calculation is exact in floating-point reference
arithmetic. Exhaustive small cases agree within \(2.1\times10^{-16}\). The
final Chernoff grid is a numerical diagnostic, not an outward-rounded
certificate.

## Ideal random-code spectrum

The first benchmark uses the expected spectrum of a uniform random
\([256,128]\) linear injection.

| \(g\) | \(\sigma\) | One-block \(\lambda\) | Dominant weight |
|---:|---:|---:|---:|
| 4 | 18 | 46.6215 | 26 |
| 4 | 17 | 38.1852 | 35 |
| 8 | 17 | 46.6239 | 26 |
| 8 | 16 | 38.1886 | 35 |

At the higher states, the expected multiplicity at weight 26 is below one:

\[
 \log_2\overline A_{26}=-10.2775.
\]

This term represents rare bad codes in the random ensemble. A fixed code
with minimum distance above 26 removes it completely.

The lower-state cases remain about 1.81 bits below the 40-bit target under
the Chernoff calculation. Their dominant weight is 35, where the ideal
expected multiplicity is about \(2^{15.56}\).

## Modeled distance-38 spectrum

The second benchmark uses a complement-symmetric even spectrum with minimum
distance 38. Its interior multiplicities are proportional to
\({256\choose w}\), then normalized to contain \(2^{128}\) words. The model
has one word at weights zero and 256.

This spectrum is real-valued. No claim is made that it is the spectrum of an
explicit code.

| \(g\) | \(\sigma\) | One-block \(\lambda\) | Dominant weight |
|---:|---:|---:|---:|
| 4 | 18 | 51.9477 | 38 |
| 4 | 17 | 40.5495 | 38 |
| 8 | 17 | 51.9530 | 38 |
| 8 | 16 | 40.5539 | 38 |

The higher-state points have about 12 bits of one-block room. If all other
modeled weights remain fixed, the multiplicity at weight 38 may increase by
about 12.25 bits before the one-block sum reaches \(2^{-40}\).

The lower-state points are fragile. Under the same comparison, weight 38 may
increase by only 0.68 bits. An actual spectrum must therefore be known quite
accurately before using \(g=4,\sigma=17\) or \(g=8,\sigma=16\).

## Interpretation

Per-block coordinate permutations remove every support property except
weight for a fixed message. They also remove rare low-weight members of the
random-matrix ensemble once a strong fixed code is selected.

The experiment supports the proposed outer replacement. It does not yet
close the construction for three reasons:

1. The distance-38 spectrum is a model rather than an explicit code.
2. Goal 01 covers one active outer block only.
3. The inner encoder still uses independent random matrices.

The next proof step should condition on two outer weights \((w_1,w_2)\). The
independent coordinate permutations choose two uniform row subsets. Their
intersection size gives a one-dimensional compression of the joint support
law. Bit transpose then averages each row using only whether zero, one, or
two nonzero bits are present.

Artifacts:

- `scripts/analyze_riffle_spectrumperm_bitshuffle_oneblock.py`;
- receipts `goal01_random_spectrum_b256_*.json`;
- receipts `goal01_d38model_b256_*.json`.

Recommended next goal: implement the exact two-active-block evaluator from
\((w_1,w_2)\) and their random intersection size. Use the distance-38 model
first, then expose an exact-spectrum input for a selected constituent code.
