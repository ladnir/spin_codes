# Certified packet-transpose g=8 tradeoff

Fix message length \(n=2^{20}\), output length \(2^{21}\), packet width
\(g=8\), and distance threshold 188743. Complete outward-rounded occupation
sums certify these power-of-two outer points.

| Outer length \(B\) | Inner state \(\sigma\) | Certified \(\lambda\) | Dominant active blocks |
|---:|---:|---:|---:|
| 1024 | 102 | 140.289349191542 | 424 |
| 2048 | 54 | 601.228272652283 | 1 |
| 4096 | 30 | 1011.906471547732 | 1 |
| 8192 | 17 | 2238.062870150003 | 1 |
| 16384 | 11 | 3626.465383813988 | 1 |
| 32768 | 8 | 5191.959819043644 | 1 |

The immediately lower states fail in the exploratory full-spectrum
calculation:

| \(B\) | Failing \(\sigma\) | Exploratory \(\lambda\) | Dominant active blocks |
|---:|---:|---:|---:|
| 1024 | 101 | -371.70 | 424 |
| 2048 | 53 | -197.61 | 208 |
| 4096 | 29 | -366.69 | 104 |
| 8192 | 16 | -1812.88 | 48 |
| 16384 | 10 | -4565 | 24 |
| 32768 | 7 | -3423 | 16 |

The transition is sharp. At the failing cells, many active outer blocks are
packed into a small number of fixed groups. Those collisions persist in all
\(B\) transposed rows. One additional state bit suppresses that bulk event;
except at \(B=1024,\sigma=102\), dominance then jumps to one active outer
block.

This behavior is not an intrinsic limitation of width-eight convolution.
`../riffle_transpose_bitshuffle_randomstepconv/TRADEOFF_G4_G8.md` shows that
reshuffling individual bits before packet formation removes the persistent
collision event and certifies \(B=1024,g=8,\sigma=12\) at 83.2011 bits.

The full receipts are named `g8_b<B>_sigma<sigma>_full_interval.json` under
`receipts/`. The adjacent exploratory receipts use the suffix
`_boundary.json`.

Recommended next goal: benchmark packet transpose against bit transpose at
the same output size. The distance result strongly favors bit transpose at
\(g=8\), while packet transpose shuffles eight times fewer items.
