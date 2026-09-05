# Proof status

Fix message length (n=2^{20}), output length (N=2^{21}), and threshold
(d=\lfloor0.09N\rfloor=188743). The current calculation uses the modeled
complement-symmetric spectrum for a binary ([256,128,38]) outer code.

For each fixed nonzero message, probability is over the coordinate
permutations, region permutations, and checkpoint multipliers sampled during
setup. The first-moment calculation bounds the expected number of nonzero
messages whose encoded output has weight at most (d).

## Regular outer words

The spectrum-density envelope replaces each regular outer word by a uniform
256-bit word at a factor smaller than two. This reduces the regular outer sum
to the number (a) of active outer blocks.

The log-domain recurrence computes the exact tilted two-state transfer for
every region occupation. The following ranges cover every value of (a).

| Regular active blocks | Method | Aggregate margin |
|---:|---|---:|
| 1–511 | exact log-domain recurrence | 55.9507 bits |
| 512–4096 | exact log-domain recurrence | 19399.3836 bits |
| 4097–7241 | exact log-domain recurrence | 28938.7956 bits |
| 7242–8192 | invertibility and Hamming-ball volume | 168.1408 bits |

The low range dominates the sum. Within the exact middle range, the weakest
point occurs near (a=3240) and retains approximately 19399 bits. The
middle-density obstruction observed for K=32 is absent.

These margins use nearest binary64 arithmetic and optimized Chernoff tilts.
They are not outward-rounded certificate values.

## Remaining cases

The density envelope excludes the unique all-one outer word. The analytic
dense bound sums every all-one count when at least 7249 regular blocks are
active. Configurations with fewer regular blocks and at least one all-one
block remain open.

A full construction claim also requires an explicit outer code with the
stated spectrum envelope. The final checker must use fixed tilts and
outward-rounded arithmetic.
