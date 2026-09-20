# Goal 01: entropy ledger and sigma boundary

## Objective

Minimize true setup entropy subject to an end-to-end first moment at most
\(2^{-40}\). Do not count pseudorandom expansion as information-theoretic
randomness.

## Lossless map-family reductions

The outer proof needs every fixed nonzero 256-bit input to produce a uniform
nonzero 512-bit output. This output law needs at least
\(\log_2(2^{512}-1)\) setup bits per independently sampled block. The
FieldPair512 map attains this lower bound.

The inner proof needs every live 22-bit input-state value to produce a uniform
22-bit output-state value. This transition needs at least 22 fresh setup bits
per step. Scalar multiplication in \(\mathbb F_{2^{22}}\) attains this lower
bound.

Thus both map families are entropy-optimal under the current exact laws.

## Entropy ledger

| Setup object | Entropy in bits |
|---|---:|
| 4096 FieldPair512 maps | 2,097,152 |
| 524352 scalar inner maps | 11,535,744 |
| Uniform global packet permutation | 9,206,311.134 |
| Two uniform 128-bit parity permutations | 1,432.323 |
| Total | 22,840,639.457 |

The total is approximately 2.723 MiB. The dense-matrix baseline used
848,105,407.457 bits, so the reduction factor is 37.13.

## State-size boundary

At \(\sigma=18\), three exact-saddle intervals give

```text
supports 50--399:     -53.130967
supports 400--1400:   -43.600759
supports 1401--3200:  -56.600958
combined:             -43.598633
```

Supports below 50 and above 3200 are negligible in the sampled Chernoff scan.
They still require a formal cover.

At \(\sigma=17\), the central interval from 400 through 1400 has exponent
`+2.273900`. Therefore the current first-moment method cannot certify
\(\sigma=17\).

## Next entropy target

The uniform global permutation now consumes 40.3% of the setup entropy. The
next reduction should replace it with a structured interleaver or a
low-entropy permutation family. Such a replacement changes the placement law
and needs a new proof.
