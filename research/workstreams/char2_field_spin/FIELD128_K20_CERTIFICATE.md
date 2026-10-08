# GRS[128,87] at K20: complete distance certificate

The ideal regional ensemble defined below has relative distance at least 10% at K20, except with probability less than `2^-67`. Its requested dimension is `2^20`, its output length is `1548288`, and its effective rate is `128/189 > 2/3`.

The [aggregate replay](theory/field128_whole_certificate.py) recomputes and joins the sparse, middle, and dense certificates. It completed successfully with component margins 67, 408, and 13690 bits. This is a whole-code result for the specified K20 ensemble, not a claim about the naturally repadded K16 or K18 profiles. Encoding performance remains a separate implementation question.

## Exact construction and setup distribution

Let `F=GF(256)`, represented by AES-field bytes with modulus `X^8+X^4+X^3+X+1`. Addition of byte labels below is field addition, hence XOR.

The padded source has `L=12096` rows, each containing 87 symbols. For source row `r`, interpolate the unique polynomial `f_r` of degree below 87 from its values at byte labels `0,...,86`. Its outer word is `(f_r(0),...,f_r(127))`, the fixed systematic GRS[128,87] codeword. Each outer coordinate forms one region containing its `L` row values.

Setup samples the following mutually independent randomness:

- For each of the 128 regions, a uniform permutation of its `L` rows. Together these form one regional routing stage.
- One uniform element of `F*` for each of the `N=128L` routed input coordinates.
- Eight uniform elements of `F*` for each sixteen-coordinate packet, forming that packet's state diagonal.

All these choices are sampled once at setup and then fixed. Repeated encodings use the same linear code. Independence across packet indices is a setup assumption, not a requirement to resample the code on each call.

The deterministic packet maps are those checked in [state_capacity_reference.py](theory/state_capacity_reference.py). To specify them explicitly, let

\[
 W=(0,1,92,93,224,225,188,189),\qquad V=16+W,
 \qquad (p_0,p_1,p_2,p_3)=(0,2,6,4).
\]

Here `W` is an ordered eight-element additive subgroup, not a GF8 subfield. Phase `a` has the ordered input labels

\[
 U_a=(p_a+W,\ p_a+80+W).
\]

For `u` in `U_a` and `v` in `V`, define

\[
 (A_a)_{u,v}=(u+v)^7+(u+v)^3+(u+v)+1,
 \qquad C_a=A_a^T.
\]

The reference proves that both maps are superregular and `C_a A_a=0`. These are the symmetric maps; no asymmetric injection or alternative state basis is substituted.

Process the routed vector in consecutive packets, using phase `a=t mod 4` at packet index `t`. Let `X_t` be that packet after its input diagonal, and let `D_t` be its eight-coordinate state diagonal. With `Q_0=0`, compute

\[
 R_t=D_tQ_t,\qquad
 Y_t=X_t+2A_aR_t,\qquad
 Q_{t+1}=R_t+C_aX_t.
\]

Concatenate the `Y_t` values and discard the final state. The state is not reset at regional boundaries. Identity feedthrough makes the map from `X` to `Y` invertible: given the earlier inputs, each packet's input is recovered from its output and known incoming state. The input diagonal is also invertible.

The probability bound concerns this ideal independent setup distribution. It does not prove that a particular seeded setup is good, or that a pseudorandom setup implementation has the same statistical law.

## Geometry and bad event

Natural padding to a multiple of `87*64` gives

\[
 K_{\rm req}=2^{20}=1048576,\quad
 K_{\rm pad}=87L=1052352,\quad
 N=128L=1548288.
\]

The design rate is `87/128`; the requested rate after padding is

\[
 \frac{K_{\rm req}}N=\frac{128}{189}.
\]

For a sampled setup, let `G` be its padded encoder over `F`. The bad event is

\[
 \exists x\in F^{K_{\rm pad}}\setminus\{0\}:
 \operatorname{wt}(Gx)\le e,
 \qquad e=\lceil N/10\rceil-1=154828.
\]

Weight counts nonzero field coordinates. Outside this event, the padded code has minimum distance at least `154829`, and therefore relative distance at least 10%.

## Complete class coverage

Let `A` be the number of nonzero outer rows, and let `W` be the support size entering the inner map. Routing and nonzero input diagonals preserve that support size. Each nonzero outer row has weight at least `128-87+1=42`, so `W>=42A`.

| Class | Complete integer range | Failure margin | Replay |
|---|---|---:|---|
| Sparse | `1 <= A <= 1105`, with no restriction on `W` | 67 bits | [Sparse certificate](theory/field128_sparse_certificate.py) |
| Middle | `46449 <= W <= 371589` | 408 bits | [Middle certificate](FIELD128_DENSITY_CERTIFICATE.md) |
| Dense | `371590 <= W <= 1548288` | 13690 bits | [Dense certificate](FIELD128_DENSE_CERTIFICATE.md) |

Every nonzero source has `A>=1`. If it is outside the sparse class, then

\[
 W\ge42\cdot1106=46452>46449.
\]

It therefore lies in the middle or dense class. Their weight intervals meet exactly at `371589/371590`. Overlap with the sparse class is harmless for a union bound; no independence between these three bad events is assumed.

The sparse replay covers every active-row count with 129 outward-rounded intervals. Its partition hash is `ec0b0f8a52d326e2b115519a832fef5ff8f9f7a193a57dc868bbc5d53f411d5b`. The middle replay covers all 325141 weights in 210 directed cells, with hash `ea5af80dd068950da774d2daf8ed3cd5577588120652a17b7dc817952204fe83`. The dense replay uses 152 exact rational density cells and one uniform membership bound over their union.

The aggregate checks that all receipts use the same outer, field, packet width, state size, padding, output length, and forbidden threshold. It adds the sparse exact dyadic bound, the middle outward decimal endpoint interpreted as an exact fraction, and the dense bound `2^-13690`. It does not substitute the minimum component margin or sum rounded bit margins.

The resulting upward-rounded bound is

\[
 \Pr[\mathrm{bad}]
 \le 3.8940860675637015699105930599988818654381707827758864323492950113200209500005810
 \cdot10^{-21}<2^{-67}<2^{-40}.
\]

Thus the whole-code claim follows from an exact sum of matching finite bounds, not from interpolation between successful floating samples.

## Source order and shortening

The prototype uses the fixed four-row source order

\[
 \operatorname{sourceIndex}(r,j)
 =348\lfloor r/4\rfloor+4j+(r\bmod4),
 \qquad 0\le j<87.
\]

This is a coordinate permutation of the padded source. It does not change the padded code. Keeping only the first `K_req` coordinates in that order and setting the rest to zero restricts the code to a subspace. A subcode cannot have smaller minimum distance. This argument permits a partially occupied final tile; it does not require padding to remove whole rows.

The theorem therefore applies to the requested source size and that source order. It does not certify naturally repadded K16 or K18 instances, whose regional lengths and output lengths differ. Shortening the fixed K20 code to fewer source coordinates would retain this K20 output length and is a different statement.

Separate [K16/K18 certificates](FIELD128_SIZE_CERTIFICATES.md) now cover those naturally repadded instances at 10% distance with 63/79 bits. They do not modify this K20 replay.

## Record widths and scalar extension

Fix a setup outside the base-field bad event. Let `E` be any extension of `F`, and let `beta_1,...,beta_d` be an `F`-basis of `E`. Every nonzero source vector over `E` has an expansion

\[
 x=\sum_{j=1}^d\beta_j x_j,
 \qquad x_j\in F^{K_{\rm req}},
\]

with at least one nonzero `x_j`. Since every encoder coefficient lies in `F`,

\[
 Gx=\sum_{j=1}^d\beta_jGx_j.
\]

The support of `Gx` contains the support of each component `Gx_j`. A nonzero component therefore supplies the same distance lower bound `154829`. A good base-field setup works simultaneously for all compatible scalar extensions; no extra failure factor is needed.

In particular, the supported 8-, 32-, 64-, and 128-bit records correspond to `GF(2^m)` with `m` equal to those bit widths. Their extension degrees over GF256 are `1,4,8,16`. Bytewise GF256 arithmetic also realizes the same componentwise linear map on those records, without general multiplication in the larger field. The distance counts nonzero records, not individual nonzero bits.

The scalar-extension argument applies to the complete base-field theorem. Individual sparse, middle, or dense classes need not be preserved when taking a component, so those partial certificates are not lifted separately.

## Replay and review

```text
python -B research/workstreams/char2_field_spin/theory/field128_whole_certificate.py
python -B research/workstreams/char2_field_spin/theory/field128_whole_certificate_tests.py
```

The full aggregate completed successfully. All five aggregate tests pass; they exercise exact summation, geometry mismatches, missing classes, nonrigorous receipts, and invalid probability bounds. The sparse, middle, and dense component suites have 4, 8, and 4 passing tests respectively. The sparse and middle implementations also received independent mathematical audits and separate full replays.

The next research gate is implementation competitiveness for this certified geometry. Any change to the setup distribution, packet maps, state coordinates, routing structure, or finite geometry requires a matching proof audit or replay; an arithmetic refactor preserving the exact encoder does not change the theorem.
