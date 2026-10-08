# Sparse certificate for fixed seven-of-eight output puncturing

The candidate uses a GF(256) MDS [128,84] outer, one regional route, and the full-Cauchy width8/state8 inner.
Each packet retains output coordinates 0 through 6 and discards coordinate 7.
The retained code has nominal rate 3/4.
This note certifies its sparse active-row classes only.
The whole-code claim additionally needs the middle and dense components and their checked joins.

The shared construction identity is in `theory/gf8_punctured_profile.py`.
At requested K=2^20, padding gives 1,053,696 input records and 12,544 outer rows.
The route consumes 1,605,632 records; the punctured encoder emits 1,404,928 records.
The distance event includes output weights through 140,492, strictly below one tenth of the emitted length.
The normalized F,A,C matrices and eta byte have SHA-256 digest
`6e3cf3e5d6b5876838afc5b2d7230f931974f363f9f67b402843f677521450a1`.

## Retained-output local bound

Restricting the original full-Cauchy packet to seven selected output rows preserves the selected-row rank and Schur-minor arguments.
The transfer sums selected zero sets among those seven rows only.
The birth and return minimum output weight for input weight w is max(0,8-w).
In particular, a weight-eight input may produce zero retained output; the bound includes this possibility.
It does not assume that the punctured inner is globally invertible.

`theory/gf8_punctured_screen.py` supplies exact rational local coefficients when its two marks are rational.
The sparse replay reuses the existing positive width-eight matrix-polynomial multiplication and regional binomial-ratio sum.
For A active outer rows, the source parity penalty is tau^(-44A-1).
The route uses 128 regions of length 12,544.
The output-weight mark uses the punctured output length and threshold, not the route length.

## Directed result

`theory/gf8_punctured_sparse_certificate.py` covers every integer A from 1 through 713.
It uses rational marks, 4096-bit outward Arb arithmetic, and adaptive positive interval bounds.
The replay accepts 49 disjoint intervals and gives at least **116 bits** for their combined failure upper bound.
The singleton class alone has at least 232 bits.
The partition digest is
`c1ef785aeb8db864918cdfaf3ab487e20a52c07df44309cdcba9065277a7f280`.

Each interval uses the entrywise monotonicity of the positive regional matrix in A.
The source factor is maximized independently over the interval.
The replay multiplies by the number of integer classes and adds every interval bound.
The returned receipt contains the exact dyadic sum and every rational mark; no optimizer participates in acceptance.

The outer's minimum row weight is 45.
Therefore every source with A at least 714 has raw weight at least 32,130.
A middle component beginning at raw weight 32,113 overlaps the sparse component and leaves no gap.
This arithmetic does not establish the middle component by itself.

Reproduce from the repository root:

```text
python -B research/workstreams/char2_field_spin/theory/gf8_punctured_sparse_certificate.py
```

The three dedicated tests check shared geometry, exact local coefficients, a directed singleton, and the complete short range A=1 through 16.
All pass.
