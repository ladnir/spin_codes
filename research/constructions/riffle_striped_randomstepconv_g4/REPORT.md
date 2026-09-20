# Initial g=4 full-stripe exploration

## Comparison point

Under one global packet permutation, the first observed packet-width-four
contour was

\[
B=4096,
\qquad
\sigma=20,
\qquad
\lambda\approx47.25.
\]

The preceding outer length \(B=2048\) had a high-memory diagnostic near 32.62
bits. Additional inner memory could not move that global-permutation bound to
40 bits.

## Full-stripe results

The structured evaluator computes every regional candidate coefficient
exactly. It also sums every active-outer-block occupation. The output-weight
tail uses a finite-grid Chernoff bound.

| Outer length \(B\) | Full stripes \(S=B/4\) | State bits \(\sigma\) | Margin \(\lambda\) | Dominant occupation |
|---:|---:|---:|---:|---:|
| 512 | 128 | 20 | +118.52 | 1 |
| 512 | 128 | 40 | +140.48 | 1 |
| 1024 | 256 | 14 | +35.37 | 1 |
| 1024 | 256 | 15 | +110.03 | 1 |
| 2048 | 512 | 12 | -470.79 | 15 |
| 2048 | 512 | 13 | +84.32 | 1 |
| 2048 | 512 | 40 | +634.91 | 1 |

Thus the first observed 40-bit contours are:

\[
(B,\sigma)=(1024,15)
\quad\text{and}\quad
(B,\sigma)=(2048,13).
\]

The point \((B,\sigma)=(512,20)\) also clears the target. Its exact memory
boundary has not been located.

## Interpretation

At packet width four, the outer message entropy per packet is half the value
at packet width eight. Concentrating packets near stripe boundaries is
therefore more expensive relative to the number of outer messages.

Full striping removes the old high-memory floors at outer lengths 512, 1024,
and 2048 in the current bound. As memory increases, the dominant profile
quickly collapses from several active outer blocks to one active block.

The most balanced point found so far is \(B=1024,\sigma=15\). Relative to the
global-permutation point \(B=4096,\sigma=20\), it uses a four-times smaller
random outer constituent and five fewer state bits.

## Scope

The positive values are provisional first-moment upper bounds. They retain
the correct proof direction, but floating-point operations require outward
rounding before any minimum-distance claim is closed.

Negative values only show failure of the displayed upper bound. They do not
refute the corresponding construction.
