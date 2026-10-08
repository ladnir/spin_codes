# Directed middle and dense certificates for the punctured code

Keeping seven outputs from each eight-input packet permits a stronger outer at the same final rate.
This checkpoint certifies the middle and dense input classes for the resulting GF(256) code.
It does not report an implementation time.

## Fixed construction

The outer applies a fixed MDS `[128,84]` code to each message row.
Each outer coordinate forms a region, and each region receives one independent uniform permutation.
The inner starts from zero state and processes eight input symbols per packet with eight persistent state symbols.
At each packet it samples independent nonzero input and state diagonals, then computes

\[
X=D_{\rm input}x,\qquad R=D_{\rm state}r,
\qquad Y=FX+2AR,\qquad r'=R+CX.
\]

The maps are `fullstate_matrices(normalized=True)` from `gf8_mds_packet.py`.
They come from a full 16-by-16 Cauchy matrix, with an invertible change of basis on all state rows.
The encoder retains output coordinates 0 through 6 and discards coordinate 7 from every packet.
It discards the final state.
There is no second permutation or random puncturing.
Setup is sampled once and fixed for all messages.

The shared identifier is
`gf256-full16-cauchy-w8-s8-keep7-grs128-84-regional-nonzero-diagonals`.
The ordered `F,A,C,eta` bytes have SHA-256
`6e3cf3e5d6b5876838afc5b2d7230f931974f363f9f67b402843f677521450a1`.
All component receipts use the same metadata helper.

At requested `K=2^20`, natural padding gives:

| Quantity | Value |
|---|---:|
| Outer rows | 12,544 |
| Padded message dimension | 1,053,696 |
| Inner input length `N_in` | 1,605,632 |
| Retained output length `N_out` | 1,404,928 |
| Nominal rate | 3/4 |
| Requested-dimension rate | 256/343 |
| Largest forbidden output weight | 140,492 |

Thus excluding every nonzero message with output weight at most 140,492 proves relative symbol distance at least 10%.
Sixteen parallel GF(256) lanes can use the same sampled code.
A nonzero lane already supplies the scalar distance guarantee, so no additional lane union is required.

## Local transfer after puncturing

For `w` live input coordinates, `h` live state coordinates, and `z` selected retained-output zeros, the numerical ranks are

\[
\min(z,w+h),\qquad \min(8+z,w+h),
\]

for the selected outputs and for the selected outputs together with all eight state-return equations, respectively.
Keeping seven fixed rows preserves these statements.
They do not require the state-normalized packet matrix itself to be fully superregular.

Eliminating `min(w,z)` input pivots leaves `max(0,z-w)` superregular state equations.
The audited diagonal-state moment bound handles those equations before the remaining input diagonals are exposed.
The rectangular birth and return maps give retained output weight at least `max(0,8-w)`.
The zero lower bound at `w=8` is intentional: puncturing can remove the only nonzero output of a local packet.
The global proof must cover those events rather than exclude them by assumption.

`gf8_punctured_screen.local_matrices` implements the resulting positive two-state cover.
It counts pure returns and selected-output/return equations using the established canonical incoming-zero-flat cover.
The directed replay does not require the newer designated-live return refinement.

## Middle certificate

`gf8_punctured_middle_certificate.py` covers every integer input weight from 32,113 through 722,534.
This is the half-open raw-input density interval `[.02,.45)`.
The local matrices use exact rational marks; all matrix powers and source bounds use outward rounding at 80 decimal digits.

The outer parity is 44, so every nonzero row has at least 45 live coordinates.
For source marks `tau,lambda`, put `beta=tau*lambda`.
When `beta>=2`, the existing Bernoulli comparison permits the source polynomial

\[
(1-\alpha)^{128}
+\tau^{-44}\sum_{h=45}^{128}\binom{128}{h}
(\beta\alpha)^h(1-\alpha)^{128-h},
\qquad \alpha=W/N_{\rm in}.
\]

For a density cell, the empty term is maximized at its lower endpoint.
Each occupied term is maximized at `h/128`, clipped to the cell.
The script sums these positive termwise upper bounds and raises the result to the number of rows.
If `beta<2`, it falls back to the established relaxed source polynomial.
There is no numerical assumption about where the source sum attains its maximum.

The complete replay uses 430 cells and gives

\[
\Pr[\text{bad message in the middle classes}]
\le
2.3956738651723222382231874260523429016350255300753346620351500878570554852452462
\cdot10^{-838}
<2^{-2782}.
\]

The cell-partition SHA-256 is
`48815bda6f5eeb8b785e8f1b38d6603a5f121132545675d0c5e63cb3a9ad5822`.

## Dense certificate with omitted-output completions

The old dense membership theorem applies to a fixed **full** inner output.
Puncturing requires an additional union over its omitted coordinates.

First fix the state diagonals and a full output vector of length `N_in`.
Because `F` is invertible, the causal recursion determines the diagonal-scaled inputs `X` uniquely.
This inverse computation does not use the input diagonals or regional permutations.
Their remaining randomness is therefore available to the fixed-output membership theorem.

Each retained output vector has exactly `256^(N_in-N_out)` full-output completions.
The radius-10% retained-output ball has fewer than `2^(1.27*N_out)` vectors.
The combined union factor is at most

\[
2^{\,1.27N_{\rm out}+8(N_{\rm in}-N_{\rm out})}.
\]

Here `N_in-N_out=200704`, so the completion factor contributes 1,605,632 bits to the exponent.
It must not be omitted merely because the kept output has smaller length.

`gf8_punctured_dense.py` applies the unchanged fixed-full-output membership bound to each completion.
Arb outward rounding at 192 bits certifies every input weight from 722,535 through 1,605,632 using 112 cells.
The uniform dense failure bound is at most `2^-2539` after the output and completion union.
There is no further union over density cells: conditioned on the state diagonals, each full output has one inverse-input weight.
The helper takes the worst verified cell bound before performing the single output union.

Its witness SHA-256 is
`7bf2bfa264e77f44e00779fcb46dc71f20f809e727c10d92e17742c8c6d031bb`.

## Joining the components

Every nonzero outer row has weight at least 45.
Consequently every input of weight below 32,113 has at most 713 active rows.
The separate sparse receipt covers active counts 1 through 713.
The middle endpoint 722,534 is immediately followed by the dense start 722,535.
These interfaces give complete coverage once the three receipts are checked and their failure bounds are added.
Some sparse messages are also covered by the other components; overlap is harmless for this upper bound.

The two new component suites each have five passing tests.
They check exact geometry, local-transfer agreement, density-cell domination, full-output completion counts, directed replays, and rejected floating witnesses.
The preimage argument and the directed middle implementation also received an independent read-through audit.
No frozen certificate or production encoder was changed.

```text
python research/workstreams/char2_field_spin/theory/gf8_punctured_middle_tests.py
python research/workstreams/char2_field_spin/theory/gf8_punctured_dense_tests.py
python research/workstreams/char2_field_spin/theory/gf8_punctured_middle_certificate.py
python research/workstreams/char2_field_spin/theory/gf8_punctured_dense.py
```

The next step is the independent whole-code receipt, followed by a correctness-checked native implementation and a serialized performance comparison.
