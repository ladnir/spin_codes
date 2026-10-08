# GRS[128,87] at K16 and K18: complete distance certificates

The naturally padded K16 and K18 instances of the certified GRS[128,87] ensemble have relative distance at least 10%, except with probabilities below `2^-63` and `2^-79`. Both requested rates are exactly `2/3`. These are complete finite certificates, not extrapolations from the separate [K20 certificate](FIELD128_K20_CERTIFICATE.md).

The [new replay](theory/field128_size_certificate.py) imports the frozen packet, source, and interval helpers without changing their constants. Geometry is passed explicitly. Two proof refinements suffice: an exact modal-conditioning factor for the middle class, and a stronger exact dense membership exponent.

## Ensemble and finite geometry

The field, fixed GRS outer, packet maps, and setup distribution are unchanged from the [K20 construction](FIELD128_K20_CERTIFICATE.md#exact-construction-and-setup-distribution). In particular, the field is AES-represented `F=GF(256)`. Each source row contains 87 values at labels `0,...,86`; interpolation extends them to labels `0,...,127`.

For each regional length `L`, setup samples 128 independent uniform permutations of `L` rows. Their product is one regional routing stage. Setup also samples independent uniform nonzero input diagonals and eight independent uniform nonzero state-diagonal entries per sixteen-coordinate packet. All setup choices are mutually independent and remain fixed across encodings.

The deterministic maps are the symmetric phase maps `A_a,C_a=A_a^T` in [state_capacity_reference.py](theory/state_capacity_reference.py). The phase order is `(0,2,6,4)`, with packet index taken modulo four. These maps satisfy superregularity and `C_a A_a=0`. With input-diagonal-adjusted packet `X_t`, state diagonal `D_t`, and `Q_0=0`, the recurrence is

\[
 R_t=D_tQ_t,\qquad Y_t=X_t+2A_aR_t,\qquad
 Q_{t+1}=R_t+C_aX_t.
\]

The state persists across regions; the final state is discarded. Identity feedthrough makes the inner map invertible. No packet map, diagonal distribution, or routing structure changes in these size replays.

Let `K_req=2^k`, let `L` be the least multiple of 64 with `87L>=K_req`, and set `K_pad=87L`, `N=128L`. The forbidden output threshold is `e=ceil(N/10)-1`.

| Size | `L` | `K_pad` | `N` | `e` | Certified minimum distance | Requested rate |
|---|---:|---:|---:|---:|---:|---:|
| K16 | 768 | 66816 | 98304 | 9830 | 9831 | `2/3` |
| K18 | 3072 | 267264 | 393216 | 39321 | 39322 | `2/3` |

For each size, the bad event is the existence of a nonzero padded source whose output has weight at most `e`. Probability is over the ideal setup just specified. The theorem does not identify a particular good seed or prove statistical equivalence of a pseudorandom setup implementation.

## A sharper conditioning factor

The previous middle comparison charged `(L+1)^128` for conditioning regional counts. The following bound improves that factor without assuming equal regional densities.

For an even positive integer `L`, define

\[
 m_L(h)=\binom Lh(h/L)^h(1-h/L)^{L-h},\qquad
 \kappa_L=\frac{2^L}{\binom L{L/2}},\qquad 0\le h\le L,
\]

with `0^0=1`. Then

\[
 m_L(h)\ge m_L(L/2)=\kappa_L^{-1}.
\]

To prove this, define `g(0)=0` and `g(n)=n log n-log(n!)` for positive integers `n`. Expanding the binomial probability gives

\[
 \log m_L(h)=\log(L!)-L\log L+g(h)+g(L-h).
\]

For `n>=1`, the discrete increment is `g(n+1)-g(n)=n log(1+1/n)`; the increment at zero is zero. These increments increase because `log(1+x)>x/(1+x)` for `x>0`. Thus `g` is discretely convex. The symmetric sum `g(h)+g(L-h)` is minimized at `h=L/2`, proving the claim.

We next check where this lemma applies. Fix an inner-input support with regional counts `h_1,...,h_128` and total weight `W`. Inverse regional routing produces independent uniform subsets of `h_j` row labels in region `j`. Instead, sample independent Bernoulli columns with their respective parameters `alpha_j=h_j/L`. Let `E` be the event that every column has its prescribed count. Conditioned on `E`, the Bernoulli law is exactly the routed subset law, and

\[
 \Pr[E]=\prod_{j=1}^{128}m_L(h_j)\ge\kappa_L^{-128}>0.
\]

For any nonnegative function `G` of the support matrix and any `lambda>0`,

\[
 \mathbb E_{\rm route}G
 =\lambda^{-W}
   \frac{\mathbb E_{\rm Bern}[\lambda^{\mathrm{wt}}G\,1_E]}{\Pr[E]}
 \le \kappa_L^{128}\lambda^{-W}
   \mathbb E_{\rm Bern}[\lambda^{\mathrm{wt}}G].
\]

Every conditioned count is modal under its own Bernoulli parameter. The lemma does not lower-bound an arbitrary count probability under a fixed common parameter.

For the source comparison, let `r` be a row's support size and fix `1<=tau<=255`. Its feasible rank factor is one at `r=0`, zero for `1<=r<=41`, and `tau^(r-41)` for `r>41`. The product of these factors over rows is nonnegative, so the conditioning inequality applies. Relaxing each row polynomial gives

\[
 \prod_j(1-\alpha_j)
 +\tau^{-41}\prod_j(1-\alpha_j+\tau\lambda\alpha_j).
\]

Only now apply arithmetic--geometric mean separately to the two products. With `alpha=W/N`, the result is bounded by

\[
 (1-\alpha)^{128}
 +\tau^{-41}(1-\alpha+\tau\lambda\alpha)^{128}.
\]

Thus the existing middle formula remains valid with `(L+1)^128` replaced by `kappa_L^128`. No regional profile is omitted, and no union over profiles is added. The packet transfer remains the audited joint-return transfer from the [K20 middle certificate](FIELD128_DENSITY_CERTIFICATE.md).

The replay computes `kappa_L` as an exact fraction. It reuses eight frozen rational packet witnesses and the existing exact source-tilt retargeting rule. Each adaptive weight interval uses the directed source bound at its density endpoints, the correct endpoint for the extraction power, and the number of integer weights in the interval. Successful point samples are not used as coverage evidence.

## Exact dense exponent

The dense replay keeps the original conditioning factor `(L+1)^128`. It strengthens the per-coordinate membership exponent from `1.28` to `1.30` using all 152 existing rational tilts.

For cell index `j=48,...,199`, let `B_j` be the rational upper bracket returned by `interval_base(128,87,j,lambda_j)`. The replay verifies the exact integer-rational inequality

\[
 B_j^{200}2^{260\cdot128}\le\lambda_j^{128j}.
\]

This proves the required exponent `13/10` throughout every density cell in `[.24,1]`. The output ball contributes at most `2^(1.27N)`, as certified by

\[
 (10/9)^{100}2295^{10}<2^{127}.
\]

Consequently the dense failure bound is at most

\[
 (L+1)^{128}2^{-.03N}
 \le 2^{-\left(\lfloor .03N\rfloor-128\,\operatorname{bitlength}(L)\right)}.
\]

The resulting dense margins are 1669 bits at K16 and 10260 bits at K18. This step does not rely on the new modal-conditioning lemma.

## Complete coverage and exact sum

Let `A` be the number of nonzero outer rows and `W` the inner-input weight. The GRS minimum distance gives `W>=42A`. The sparse replay covers every weight for its stated active-row range.

| Size | Sparse active rows | Middle input weights | Dense input weights | Component margins | Whole margin |
|---|---|---|---|---|---:|
| K16 | `1..122` | `5161..23592` | `23593..98304` | `67,63,1669` bits | 63 bits |
| K18 | `1..280` | `11797..94371` | `94372..393216` | `99,79,10260` bits | 79 bits |

Outside the K16 sparse class, `W>=42*123=5166>=5161`. Outside the K18 sparse class, `W>=42*281=11802>=11797`. In each case, the middle and dense intervals are adjacent. These inequalities cover every nonzero source; overlap with the sparse class is harmless.

The aggregate checks the outer, field, packet width, state size, construction identifier, source sizes, regional length, output length, and forbidden threshold of every component. It also checks all class endpoints. It adds the sparse dyadic bound, the middle directed decimal endpoint interpreted as an exact fraction, and the dense dyadic bound. It does not add rounded bit margins.

The exact sums have the following upward-rounded decimal bounds:

\[
 \Pr[\mathrm{bad}_{16}]\le6.988111\cdot10^{-20}<2^{-63},\qquad
 \Pr[\mathrm{bad}_{18}]\le1.370907\cdot10^{-24}<2^{-79}.
\]

Both bounds are below `2^-40`. No independence between the sparse, middle, and dense bad events is assumed.

## Shortening and record widths

The padded theorem permits the fixed four-row source order

\[
 \operatorname{sourceIndex}(r,j)
 =348\lfloor r/4\rfloor+4j+(r\bmod4),\qquad 0\le j<87.
\]

This is a source-coordinate permutation. Setting coordinates after `K_req` to zero selects a subcode, so distance cannot decrease. The argument allows a partially occupied final tile; it does not require whole-row padding.

Fix a setup outside the complete base-field bad event. Over any extension of GF256, expand a nonzero source in a GF256 basis. At least one component is nonzero, and the encoder acts componentwise because all coefficients lie in GF256. The output support contains that component's support. Hence the same minimum distance holds over every compatible extension, with no additional failure probability.

In particular, 8-, 32-, 64-, and 128-bit records use extension degrees `1,4,8,16`. Bytewise GF256 multiplication realizes the same linear map without general multiplication in the larger field. Distance counts nonzero records. Scalar extension applies to the complete theorem, not separately to the three support classes.

## Replay and review

```text
python -B research/workstreams/char2_field_spin/theory/field128_size_certificate.py
python -B research/workstreams/char2_field_spin/theory/field128_size_certificate_tests.py
```

Both full size replays completed successfully and were independently repeated. All ten tests pass. They include both full replays, exact small modal and source-profile checks, all dense rational witnesses, a packet-prefactor reference, and rejection of corrupted geometry, construction, coverage, or probability bounds.

The sparse partitions contain 28 and 48 intervals; their hashes are `507d1ec6401404b4506e11436d316a0d308fd4acc81eefece4a0d4ed2fe431ab` and `1a69bcf89fcd694c299ded95210f9fc9a68cebfc73dcd2d49bad65b8af1e1bd8`. The middle partitions contain 191 and 220 cells; their hashes are `ca13f18f6dfd644eb71ca7b0e46274874f0498cf127f44f88a6c85137ad71f77` and `83ff5fee53441de52a1d0461fe1e5023017a6c429604ce54ddec631a2aa6c3dc`. The script recomputes the bounds; these hashes are not substitutes for the arithmetic replay.

The frozen K20 scripts remain unchanged. The next task is implementation competitiveness for this now-certified three-size family, without changing its ensemble. These proof replays make no new runtime claim.
