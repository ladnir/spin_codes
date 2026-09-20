# Certified g=4 outer-memory tradeoff

Fix \(n=2^{20}\), packet width \(g=4\), and distance threshold 188743 in an
output of length \(2^{21}\). The following boundary points have complete
outward-rounded occupation certificates.

| Outer length \(B\) | Inner state \(\sigma\) | Certified \(\lambda\) |
|---:|---:|---:|
| 1024 | 13 | 83.087891407589 |
| 512 | 15 | 57.641932041422 |
| 256 | 18 | 45.400358851794 |

Every row clears the target \(\lambda=40\). In all three certificates, one
active outer block is the dominant occupation.

The adjacent lower-memory cells do not reach 40 bits in the current
one-block bound:

| Outer length \(B\) | Inner state \(\sigma\) | One-block \(\lambda\) |
|---:|---:|---:|
| 1024 | 12 | -2.29 |
| 512 | 14 | 34.97 |
| 256 | 17 | 37.28 |

For \(B=128\), the one-block bound approaches only about 21.38 bits as
\(\sigma\) grows. Thus the current proof method has no 40-bit point at that
outer length.

## Operation-count proxy

This comparison treats every dense binary matrix-vector bit product as one
unit. It is not a measured implementation cost.

There are \(L=2n/B\) rate-half outer blocks. A dense outer multiplication
uses \(B^2/2\) bit products, so all outer blocks use \(nB\) products. The
inner has \(2n/g\) steps. A dense random step map uses
\((g+\sigma)^2\) products.

| \(B\) | \(\sigma\) | Outer products | Inner products | Total proxy |
|---:|---:|---:|---:|---:|
| 1024 | 13 | 1,073,741,824 | 151,519,232 | 1,225,261,056 |
| 512 | 15 | 536,870,912 | 189,267,968 | 726,138,880 |
| 256 | 18 | 268,435,456 | 253,755,392 | 522,190,848 |

Under this proxy, \(B=256,\sigma=18\) is cheapest despite its larger inner
state. Relative to \(B=1024,\sigma=13\), it cuts the estimated total by about
57%. The \(B=512,\sigma=15\) point is the middle choice.

All three variants permute the same total number \(2n/g\) of packet items.
Smaller \(B\) produces fewer rows with larger permutations. A Fisher-Yates
implementation performs approximately the same number of swaps, although
the random indices require more bits.

These counts omit word-level batching, memory traffic, setup expansion, and
the cost of generating the random maps. They are useful for choosing which
instances to benchmark, not for predicting wall-clock time.

## Receipts

- `receipts/g4_b1024_sigma13_full_interval.json`;
- `receipts/g4_b512_sigma15_full_interval.json`;
- `receipts/g4_b256_sigma18_full_interval.json`;
- `receipts/g4_sparse_tradeoff.json`.

Recommended next goal: benchmark the three certified instances with the same
low-level kernels. If simultaneous reduction of \(B\) and \(\sigma\) is
required, add padded intermediate outer lengths between 512 and 1024 and
repeat the certificate calculation.

