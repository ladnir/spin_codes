# Smaller outers with the state-eight joint return cover

The bounded K20 screen identifies GRS[128,87] as a proof candidate: every tested active-row class has more than 208 bits of floating headroom. GRS[128,88] has a narrow trough near 240 active rows, where the evaluated bound gives only 1.27 bits. GRS[64,44] has a much larger unresolved range. The stronger joint return operator does not improve these particular sparse samples: its pure-return alternative is selected at their optimized marks.

These are floating evaluations of finite upper-bound formulas, not distance certificates. The sampled classes do not cover all active-row counts or all input weights. An inadequate upper bound does not establish that a construction fails. No native implementation or runtime experiment was performed in this screen.

The subsequent [GRS[128,87] middle replay](FIELD128_DENSITY_CERTIFICATE.md) closes every integer input weight in `[.03N,.24N)` at K20 with a directed bound below `2^-408`. It does not turn the sparse samples below into a complete sparse certificate.

## Construction and finite formula

The principal screen keeps the certified inner construction unchanged: width 16, eight state coordinates, GF256, symmetric Cauchy maps, eta 2, and independent nonzero input and state diagonals. It changes only the systematic GRS outer and its natural padding. There is still one regional route.

Write `b,k` for the outer length and dimension, `p=b-k`, and `A` for the number of active outer rows. At K20, the regional length is `L=64 ceil(2^20/(64k))`, the output length is `N=bL`, and the forbidden output weight is `e=ceil(N/10)-1`. For marks `0<y<=1` and `1<=tau<=255`, let `T_w(y,tau)` be the two-state packet cover for `w` active input coordinates. Its return entry is the minimum of the complete pure-return and split selected-zero covers in [state8_joint_return.py](theory/state8_joint_return.py). The other entries are unchanged.

Define the positive matrix polynomial and the regional matrix by

\[
 F_h=[x^h]\left(\sum_{w=0}^{16}\binom{16}{w}x^wT_w(y,\tau)\right)^{L/16},
 \qquad
 R_A=\sum_{h=0}^{A}\frac{\binom Ah}{\binom Lh}\tau^h F_h.
\]

The existing active-row rank cover gives the finite class bound

\[
 B_A=\binom LA\tau^{-pA-1}y^{-e}(1,0)R_A^b\mathbf1,
\]

where `1` is the all-ones column vector. The implementation retains the without-replacement factor `C(A,h)/C(L,h)`. It truncates positive matrix polynomials only above degree `A`, which cannot affect the retained coefficients. This is the same regional/source comparison as the [directed sparse replay](theory/state8_sparse_certificate.py); it is not an ordinary word first moment. The reported margin is `-log2(B_A)` evaluated in floating point. The search parametrizes `log(y)=-exp(logtilt) A/N` and `tau=exp(logtau)`.

## Executed state-eight samples

All columns use the joint local operator and K20. The 5/8-rate column is a backup comparison, not the preferred higher-rate target.

| Active rows A | GRS[64,44] | GRS[128,88] | GRS[128,87] | GRS[64,40] |
|---:|---:|---:|---:|---:|
| 1 | 92.215 | 200.729 | 208.708 | 124.056 |
| 2 | 182.568 | 400.558 | 416.517 | 246.251 |
| 4 | 365.337 | 802.283 | 834.200 | 492.702 |
| 8 | 547.377 | 1225.302 | 1286.237 | 797.677 |
| 16 | 331.478 | 878.935 | 947.001 | 612.764 |
| 32 | 153.083 | 694.976 | 799.029 | 584.128 |
| 64 | -136.794 | 439.314 | 621.861 | 622.771 |
| 128 | -583.747 | 154.168 | 500.140 | 854.163 |
| 256 | -1247.125 | 4.589 | 685.486 | 1563.247 |
| 512 | -2102.943 | 755.212 | 2137.274 | 3510.739 |

The corresponding `(L,N,e)` values are `(23872,1527808,152780)`, `(11968,1531904,153190)`, `(12096,1548288,154828)`, and `(26240,1679360,167935)`.

The additional GRS[128,88] trough checks were:

| A | 160 | 192 | 224 | 240 | 288 | 320 | 384 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| Margin | 75.449 | 26.548 | 3.680 | 1.268 | 27.812 | 72.336 | 222.539 |

The `A=240` search also used an independently distant initial point. This supports stopping a complete replay of the present bound for GRS[128,88]; it does not prove that no better marks or bound exist.

For GRS[128,87], additional checks at `A=96,160,192` gave `533.790,504.183,538.692` bits. All thirteen selected marks are frozen in `WITNESSES_128_87` in [joint_outer_screen.py](theory/joint_outer_screen.py). At every selected mark in these tables, the joint improvement over the improved pure-return cover was zero to the reported floating precision. The successful change is one additional outer parity coordinate, not a new sparse benefit from selected output zeros.

## Distinct state-six coarse comparison

The state-six comparison uses the generic operator in [superregular_return_screen.py](theory/superregular_return_screen.py), not the state-eight joint operator. Its proposed maps have six distinct residue nodes outside all 64 packet labels. With `Q(X)` their degree-six vanishing polynomial, set `A[u,j]=Q(u)/(u+v_j)` and `C=A^T`. Each column is a degree-five polynomial. Products of columns have degree at most ten, so their sum on each affine sixteen-label space is zero and `CA=0`. The maps are row-scaled Cauchy and therefore superregular.

The generic return cover uses minimum returned weight `max(1,11-w)` and six diagonal rank charges, with the weighted-flat factor. No state-eight singleton or joint-output assertion is imported.

| A | State 6, GRS[128,87] | State 6, GRS[128,86] |
|---:|---:|---:|
| 1 | 208.708 | not screened |
| 4 | 833.042 | not screened |
| 16 | 569.992 | not screened |
| 64 | 203.511 | 387.072 |
| 128 | 0.998 | 356.478 |
| 256 | -46.004 | 663.524 |
| 512 | not screened | 2100.149 |

State-six GRS[128,87] does not clear the 40-bit target at the tested trough. GRS[128,86], at rate `43/64`, clears these four tested classes, but neither its unsampled sparse classes nor its full input-weight range has been certified. Its frozen `(logtilt,logtau)` marks for `A=64,128,256,512` are respectively `(6.379403061027414,2.0138993946529746)`, `(6.351290070968021,1.9534147755932543)`, `(6.351575157103631,1.9543333762922757)`, and `(6.376799385723764,2.021968718064013)`.

A subsequent bounded global input-weight check for GRS[128,86] retained the same generic state-six packet cover and the finite regional conditioning factors from [state8_density_screen.py](theory/state8_density_screen.py). It did not use that module's state-eight-specific local improvements. The results expose a separate middle-density gap:

| Approximate input density | Exact input weight | Floating margin |
|---:|---:|---:|
| .03 | 46940 | 2089.094 |
| .08 | 125174 | -5881.714 |
| .15 | 234701 | -84388.104 |
| .20 | 312934 | -177436.881 |

The generic coarse cover is therefore not a route to a complete state-six certificate as it stands. This remains a deficiency of the evaluated upper bound, not a lower bound on bad-code probability. No further optimization of this cover was attempted.

### Audit of the proposed generic joint repair

The split-joint argument in [STATE8_DENSITY_CERTIFICATE.md](STATE8_DENSITY_CERTIFICATE.md) does extend to state six. It only requires superregular emission and injection, `CA=0`, nonzero eta, and fresh independent input/state diagonals. It does not require the special eight-state block identity. This audit establishes a valid local bound; the preceding coarse numerical table did not use it.

For general state size `s`, fix the incoming zero flat before exposing current diagonals. Let `H` be the remaining actual nonzero state support, with `h=|H|`, and let `I` have `w` live input coordinates. A feasible return requires `w>s-h`, since otherwise the superregular matrix `C[H^c,I]` forces every input value to zero. Its numerical return matrix `[E_H,C_I]` consequently has rank `s`.

Select `i` zero outputs on `I` and `z` outside `I`, and put `v=16-i-z`. Superregularity requires `z<h`; the output equations themselves have rank `i+z`. Under return, `R=CX` and `Y=(I+eta AC)X`. The latter map is invertible because `(AC)^2=0`. Hence the joint return/output kernel injects into the `v` allowed output coordinates, giving joint rank at least `h+w-v`. Together these facts give one uniform charge

\[
 \rho_s(w,i,z)=\max\{s,i+z,w-16+i+2z+1\}.
\]

The incoming zero flat has rank at most `min(w-1,s-1-z)`. Counting its canonical independent bases yields the complete local return cover

\[
 J_{w,s}=\sum_{\substack{0\le i\le w,\ 0\le z\le\min(s-1,16-w)\\
                        16-i-z\ge m_{w,s}}}
 \binom wi\binom{16-w}z y^{16-i-z}(1-y)^{i+z}
 \tau^{-\rho_s(w,i,z)}
 \sum_{k=0}^{\min(w-1,s-1-z)}\binom sk\tau^{-k},
 \qquad m_{w,s}=\max(1,17-s-w).
\]

The minimum with the corresponding pure-return cover is valid because both covers have the same actual-zero successor. The alive diagonal-pivot lemma applies once to the full joint matrix, permitting dependencies among the live source forms. A separate root factor must not be multiplied into this return entry. The generic middle screen and any eventual directed replay must implement and test this stronger bound separately from the coarse results above.

## Outer arithmetic, separated from runtime

GRS[128,87] can reuse the structured GRS[128,88] map and compute the missing systematic value at point 87 by interpolation from points 0 through 86. Its transpose takes the old 88th output `t` and adds `c_j t` to each of the first 87 outputs.

The interpolation coefficients are one seven times and each of `12,13,72,76,101,109,211,218,250,255` eight times, in the AES-field byte representation. Together with one, their GF2 span has basis `1,12,72,76,101,211`. Thus this extra shortening step needs five nonidentity fixed GF256 products per four-row SIMD tile, plus Boolean reconstruction and output updates. The test checks all 87 monomial interpolation identities exactly.

The outer-only counts from [the earlier structured cost comparison](OUTER_LENGTH_TRADEOFF.md) therefore extend as follows. That note's complete inner totals concerned an older four-state construction and must not be reused here.

| Outer | Coefficient-shared outer products / four rows | FFT-base plus shared-shortening products / four rows |
|---|---:|---:|
| GRS[64,44] | 101 | 132 |
| GRS[128,88] | 232 | 335 |
| GRS[128,87] | at most 237 | at most 340 |

These are arithmetic circuit upper counts, not instruction schedules or performance predictions. The odd dimension 87 changes record packing and adds shortening output updates. The existing dense inner uses 70 GFNI512 instructions per sixteen output records, including input diagonals; with that particular inner, the coefficient-shared total is `2477/348 = 7.11782` GFNI512 instructions per source record for GRS[128,87]. A separately developed exact FFT inner may change this comparison. No timing from another factorization is transferred to these estimates.

## Replay and remaining decision

Run from the repository root:

```text
python -B research/workstreams/char2_field_spin/theory/joint_outer_screen_tests.py
python -B research/workstreams/char2_field_spin/theory/joint_outer_screen.py --outer 128 87
python -B research/workstreams/char2_field_spin/theory/joint_outer_screen.py --outer 128 88 --active 160 192 224 240 288 320 384
python -B research/workstreams/char2_field_spin/theory/joint_outer_screen.py --outer 128 87 --shortening-cost
python -B research/workstreams/char2_field_spin/theory/superregular_return_screen.py --outer 128 87 --state 6 --width 16 --active 1 4 16 64 128 256 --iterations 30
python -B research/workstreams/char2_field_spin/theory/superregular_return_screen.py --outer 128 86 --state 6 --width 16 --active 64 128 256 512 --iterations 30
python -B research/workstreams/char2_field_spin/theory/joint_outer_screen.py --outer 128 86 --generic-six-density .03 .08 .15 .20
```

Seven tests pass. They compare the local operator to the directed density implementation, check the joint-versus-pure relation with exact rational marks, verify a two-packet coefficient product and finite normalization against exact fractions, replay both gap witnesses, reject invalid geometry, and check the odd-dimension interpolation.

The next proof step for state-eight GRS[128,87] is a complete sparse replay for active-row counts 1 through 1105, followed by a matching dense receipt and an exact class union. Its middle class is now directed-certified separately. State-six GRS[128,86] still needs a middle-density improvement. The smaller-outer screen itself does not select a faster implementation, establish a whole-code certificate, or improve the already-known half-rate 10% result.
