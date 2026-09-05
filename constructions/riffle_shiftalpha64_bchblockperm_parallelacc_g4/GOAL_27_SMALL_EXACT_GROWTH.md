# Goal 27: measure growth with the current exact enumerator

## Result

The exact proof-gym enumerator shows very slow distance growth. The
expected-count-one crossing equals nine through data-block count eight. It
then equals ten through the largest legal instance, data-block count fifteen.

The calculation does not support linear growth in this finite family. It is
also not a refutation of the target construction. The proof gym uses
`GF(16)` and extended BCH `[8,4,4]`, and its shifted coefficient schedule ends
at fifteen data blocks.

## Experiment

For each data-block count $B$ from 1 through 15, the enumerator includes two
field parity symbols. The binary output length is

$$
n=8(B+2).
$$

The calculation averages the independent permutations inside the BCH blocks
and the global packet permutation. It sums every nonzero message in
`GF(16)^B`. Accumulator paths are clipped after output weight 16. Since
accumulator output weight never decreases, every coefficient through weight
16 remains exact.

Let $D_*(B)$ be the least $D$ for which the expected number of nonzero words
of output weight at most $D$ is at least one.

## Exact ladder

| $B$ | Binary length $n$ | $D_*(B)$ | $D_*(B)/n$ |
|---:|---:|---:|---:|
| 1 | 24 | 9 | 37.50% |
| 2 | 32 | 9 | 28.13% |
| 3 | 40 | 9 | 22.50% |
| 4 | 48 | 9 | 18.75% |
| 5 | 56 | 9 | 16.07% |
| 6 | 64 | 9 | 14.06% |
| 7 | 72 | 9 | 12.50% |
| 8 | 80 | 9 | 11.25% |
| 9 | 88 | 10 | 11.36% |
| 10 | 96 | 10 | 10.42% |
| 11 | 104 | 10 | 9.62% |
| 12 | 112 | 10 | 8.93% |
| 13 | 120 | 10 | 8.33% |
| 14 | 128 | 10 | 7.81% |
| 15 | 136 | 10 | 7.35% |

Every instance has minimum nonzero output weight six. The expected count at
weight six stays near $2^{-4.8}$ throughout the ladder.

The cumulative expected-count logs near the crossing are:

| $B$ | $n$ | through 8 | through 9 | through 10 |
|---:|---:|---:|---:|---:|
| 1 | 24 | -0.262 | 0.890 | 1.879 |
| 4 | 48 | -0.690 | 0.648 | 2.332 |
| 8 | 80 | -1.080 | 0.064 | 1.751 |
| 9 | 88 | -1.133 | -0.027 | 1.658 |
| 12 | 112 | -1.239 | -0.229 | 1.450 |
| 15 | 136 | -1.294 | -0.361 | 1.311 |

Thus the weight-nine tail slowly falls below one, while the weight-ten tail
remains above one. The movement is real but small.

## Interpretation

The absolute crossing gains one unit while the output length grows from 24
to 136. Consequently, the observed relative crossing decreases rather than
stabilizing.

A descriptive fit over $B=6,ldots,15$ makes the slow scale visible. The
base-two log of the weight-ten expected count falls by approximately 0.51
bits per doubling of $B$. Extending that fitted line would place its crossing
near $B=86$. Such an extension is not part of the construction: the present
`GF(16)` coefficient schedule is invalid beyond $B=15$.

The current data are consistent with logarithmic or otherwise sublinear
growth. They provide no evidence for linear distance. The next useful
diagnostic is not a larger proof bound. It is a decomposition of the exact
$B=15$ tail by outer occupation and binary input weight. That decomposition
would identify which shell keeps the crossing at ten.

## Reproduction

Run

```powershell
python scripts/analyze_riffle_small_exact_growth.py
```

The script writes
`receipts/goal27_small_exact_growth_b1_through_b15.json`.
