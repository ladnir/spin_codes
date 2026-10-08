# Three-quarter-rate GF256 code: complete K20 certificate

The fixed seven-of-eight punctured construction closes the proof target.
At nominal rate 3/4, its minimum retained-symbol distance is at least 10%, except with setup-failure probability at most **2^-116**.
This is a mathematical certificate under the stated setup distribution.
It is not a native timing result or a claim that an unverified implementation realizes the certified map.

## Certified code and probability space

The outer is a fixed GF(256) MDS [128,84] code, repeated over 12,544 rows.
For each of its 128 output coordinates, setup samples an independent uniform permutation of those rows.
The resulting single regional routing stage feeds the normalized full-Cauchy width8/state8 inner.
Each packet independently samples nonzero GF(256) input and state diagonal entries.
The state starts at zero; its terminal value is discarded.
The encoder retains packet output coordinates 0 through 6 and discards coordinate 7.

The full16 Cauchy matrix uses row nodes 0 through 15 and column nodes 128 through 143 in the AES field representation.
Its bottom-state block is invertibly normalized to identity.
The checked nonzero row/column scalings produce the F,A,C matrices identified by
`6e3cf3e5d6b5876838afc5b2d7230f931974f363f9f67b402843f677521450a1`.
The digest covers row-major F,A,C bytes followed by eta=2.
`theory/gf8_punctured_profile.py` is the shared construction identity.

For requested K=1,048,576, the certified padded code has:

| Quantity | Value |
|---|---:|
| Dimension | 1,053,696 |
| Coordinates before puncturing | 1,605,632 |
| Retained output length | 1,404,928 |
| Nominal rate | 3/4 |
| Excluded retained-output weights | 0 through 140,492 |
| Certified minimum distance | 140,493 |
| Setup-failure upper bound | 2^-116 |

Because 140,492 is strictly below one tenth of the retained length, the distance claim is at least 10%.
Restricting to the requested-dimensional subcode preserves the minimum-distance guarantee.
The corresponding requested rate is 256/343.

Sixteen byte lanes can use the same sampled scalar code in parallel.
If a vector-valued input is nonzero, at least one scalar lane is nonzero and has the certified support size.
Thus the same setup event covers the block-symbol distance of all lanes; no sixteen-fold failure union is needed.

## Three components cover every nonzero source

Let A count nonzero outer rows, and let W be the raw outer-output weight before the inner input diagonals.
The MDS outer has minimum row weight 45, so W is at least 45A.
The independently replayed components are:

| Component | Covered classes | Failure margin |
|---|---|---:|
| Sparse | A=1 through 713 | 116 bits |
| Middle | W=32,113 through 722,534 | 2,782 bits |
| Dense | W=722,535 through 1,605,632 | 2,539 bits |

A source outside the sparse range has A at least 714 and therefore W at least 32,130.
That weight lies inside the middle or dense range.
The middle range ends immediately before the dense range begins.
The overlap at the sparse/middle join is harmless; no nonzero class is omitted.

The aggregate adds the actual sparse dyadic upper bound, the directed middle decimal upper bound as an exact rational, and the dense dyadic upper bound.
Their sum is at most 2^-116.
It does not merely compare rounded bit counts.

## Audit of the dense preimage factor

Discarding outputs invalidates a direct appeal to invertibility of the retained inner.
The proof instead uses the full, unpunctured causal inner.
Its current-input feedthrough F is invertible, so the complete causal map is invertible for every fixed state-diagonal realization.

Fix a retained low-weight output and a completion of its discarded coordinates.
After conditioning on the state diagonals, the inverse full output is fixed independently of the input diagonals and regional permutations.
The existing fixed-vector outer-membership bound therefore applies.
Every retained vector has exactly 256^200704 full-output completions.
The dense component charges this entire factor, whose base-two logarithm is 1,605,632.
It also charges the retained-output Hamming-ball bound.

The inverse weight W depends on the fixed full output and conditioned state diagonals, not on the remaining random input diagonals or route.
The dense replay is uniform over all eligible W.
Its interval partition therefore needs no additional union over cells after the output/completion union.

## Replay and validation

Run from the repository root:

```text
python -B research/workstreams/char2_field_spin/theory/gf8_punctured_certificate.py
```

The helper checks the exact Cauchy factorization, state normalization, nonzero row/column scalings, frozen matrix digest, and causal feedthrough invertibility.
It replays all three components and checks every shared metadata field, both coverage joins, the strict distance threshold, and the dense completion factor.
The sparse arithmetic uses rational witnesses and 4096-bit outward Arb operations.
The middle arithmetic uses rational witnesses and 80-digit directed intervals.
The dense arithmetic uses exact rational endpoint brackets and 192-bit outward Arb operations.

The independent audit reran all ten middle/dense tests and four aggregate tests successfully.
The aggregate tests reject changed geometry, matrix identity, missing completion factors, invalid probability claims, and gaps between covers.
The complete replay returned `complete=true`, `rigorous=true`, `whole_code=true`, and `margin_bits_lower=116`.

The component partition digests are:

- Sparse: `c1ef785aeb8db864918cdfaf3ab487e20a52c07df44309cdcba9065277a7f280`.
- Middle: `48815bda6f5eeb8b785e8f1b38d6603a5f121132545675d0c5e63cb3a9ad5822`.
- Dense: `7bf2bfa264e77f44e00779fcb46dc71f20f809e727c10d92e17742c8c6d031bb`.

The separate [native report](GF8_PUNCTURED_NATIVE.md) now records successful
scalar, release, sanitizer, and large-size checks, with serial whole-call
measurements. Those checks bind the implementation to this map; they do
not change the probability space or certify the benchmark's fixed seed.
