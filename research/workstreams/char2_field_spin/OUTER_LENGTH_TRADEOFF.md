# Larger outers at rate 11/16: arithmetic comparison

GRS[64,44] has a plausible structured outer with 101 fixed GF256 products per four rows. With the current sixteen-coordinate Cauchy inner, its complete product volume is `357/176 = 2.02841` GFNI512 instructions per source record. That is 5.62% above the current `169/88` estimate for GRS[32,22]. GRS[128,88] gives `93/44 = 2.11364`, or 10.06% above the current estimate.

These are arithmetic upper counts, not runtime results or distance certificates. The simple coefficient-sharing schedule has substantially more Boolean work at larger lengths. The separate sparse proof screens do not yet promote either larger outer or the width-eight alternative. No native implementation, benchmark, or dense-matrix fallback was used for this comparison.

## Cost convention and existing references

Each source record has 128 bits. Four outer rows fill the four record lanes of a ZMM register. One fixed GF256 multiplication site then costs one GFNI512 instruction. The forward and transpose linear circuits have the same number of multiplication sites. A Boolean count below treats a three-input XOR as one instruction.

The current Cauchy packet costs sixteen GFNI512 instructions for sixteen encoded records: five for each fixed map, one for the emission scale, one for the state diagonal, and four for the input diagonals. Its cost per source record is therefore `b/k` for an outer GRS[b,k]. All three compared outers have `b/k=16/11`.

The relevant existing references are:

- [The 34-product GRS[32,24] circuit](theory/grs32_24_parity_reference.py), with [its derivation](GRS32_24_COST.md).
- [The shortening circuits](theory/parity_state_cost_reference.py), including the optimized 41-product GRS[32,22] outer.
- [The existing additive-FFT generator](../../../spin/experiments/packet8_rs64_fused/generate_fft.py), which interpolates on 32 additive points and evaluates on another 32-point coset. Its associated native experiment is not the construction being costed here.
- [The current Cauchy packet cost](CAUCHY_DIAGONAL_PACKET_COST.md), which supplies the common inner count.

All new arithmetic below uses the same AES GF256 representation as these references. No general multiplication in GF128 is required.

## A larger systematic map from two structured maps

Let `n` be 16 or 32 and put `m=n/4`. Start with systematic GRS[4n,3n], whose source points are the byte labels `0,...,3n-1`. Its parity points are `3n,...,4n-1`. Restrict its polynomial degree to below `k=3n-m`, giving GRS[4n,k].

For the transpose, write `P` for the parity matrix of the base GRS[4n,3n] code. Given all `4n` input records `x`, first compute

\[
 v=P^T x_{3n:4n}.
\]

Let `L` interpolate at source points `0,...,k-1` and evaluate at the `m` missing source points `k,...,3n-1`. Form

\[
 t=x_{k:3n}+v_{k:3n},
 \qquad
 \operatorname{out}=x_{0:k}+v_{0:k}+L^Tt.
\]

This is the transpose of the desired systematic GRS map. It needs `m` binary XORs to form `t`, then `k` ternary XORs to combine the three output contributions. Both nontrivial maps below have additive-coset structure.

For the base map, the parity transpose has three symmetric `n`-by-`n` XOR-circulant blocks. The shortening map has eleven symmetric `m`-by-`m` XOR-circulant blocks. Translation by the corresponding additive subspace preserves both interpolation point sets, which gives these identities. The additive subspaces are byte-label vector spaces; they are not assumed to be subfields.

## Coefficient sharing in the circulant circuit

For a block of size `r=2^d`, use the existing superset zeta transform `M`. If `c` is the circulant's first row, put `gamma=Mc` and `z=Mx`. The output is `Mh`, where

\[
 h_S=\sum_{I\subseteq S}\gamma_I z_{S\setminus I}.
\]

For one transformed input `z_T`, the required coefficients across all output blocks are

\[
 \mathcal C_T=\{\gamma_I^{(j)}:I\cap T=\varnothing,\ j\text{ an output block}\}.
\]

Multiplication is GF2-linear. A basis of `span_GF2({1} union C_T)` containing one realizes all these coefficient products with

\[
 \dim_{\mathbb F_2}\operatorname{span}(\{1\}\cup\mathcal C_T)-1
\]

GF256 multiplication sites; the remaining combinations use XORs. This counts a concrete structured upper circuit. It is not a lower bound, since products of different transformed inputs may admit further sharing.

Read-only interpolation and GF2 elimination gave the following counts:

| Target outer | Base map | Base products | Shortening map | Shortening products | Complete products |
|---|---|---:|---|---:|---:|
| GRS[32,22], control | three blocks of size 8 | 34 | eleven blocks of size 2 | 8 generic, 7 optimized | 41 optimized |
| GRS[64,44] | three blocks of size 16 | 82 | eleven blocks of size 4 | 19 | 101 |
| GRS[128,88] | three blocks of size 32 | 186 | eleven blocks of size 8 | 46 | 232 |

The sixteen base counts for GRS[64,48], in transformed-input mask order, are

```text
7,7,7,6,7,7,5,2,7,7,7,2,7,2,2,0.
```

Its shortening counts are `7,6,5,1`. The thirty-two base counts for GRS[128,96] are

```text
7,7,7,7,7,7,7,6,7,7,7,6,7,5,7,2,
7,7,7,7,7,7,6,2,7,6,7,2,7,2,2,0.
```

Its shortening counts are `7,7,7,6,7,4,7,1`. The base constant coefficients are one in every block. The shortening constant coefficients are four copies each of `82,83`, then three copies of one, for GRS[64,44]; they are four copies each of `182,183`, then three copies of one, for GRS[128,88]. All numbers denote AES-field bytes.

The control reproduces the existing 34-product base count. Its generic shortening count is one above the existing seven-product optimized circuit. This difference is retained rather than extrapolating an unproved larger shortening optimization.

## Complete product volume and the Boolean tradeoff

For `P` complete outer products per four rows, the full product volume is

\[
 \frac{P}{4k}+\frac bk.
\]

| Outer | Outer products / four rows | Outer products / source record | Complete product volume | Increase over 169/88 |
|---|---:|---:|---:|---:|
| GRS[32,22] | 41 | 41/88 | 169/88 = 1.92045 | control |
| GRS[64,44] | 101 | 101/176 | 357/176 = 2.02841 | 5.62% |
| GRS[128,88] | 232 | 29/44 | 93/44 = 2.11364 | 10.06% |

To bound the Boolean work without designing a new optimized circuit, select the coefficient basis greedily in ascending byte order, with one already available. Reconstruct each required coefficient from this basis, share identical products, and use ternary XORs for sums. This gives the following conservative counts:

| Target outer | Base Boolean count | Shortening Boolean count | Final combinations | Complete Boolean upper count |
|---|---:|---:|---:|---:|
| GRS[64,44] | 434 | 165 | 4 + 44 | 647 |
| GRS[128,88] | 1427 | 569 | 8 + 88 | 2092 |

The optimized GRS[32,22] reference uses 170 Boolean instructions per four rows. The larger counts are not optimized, but they expose a real risk omitted by product volume alone. Shared coefficients also create more intermediate values and more fan-out. A register-bounded schedule may require cache-resident scratch or recomputation; neither cost is included in the multiplication-site count.

### Additive-FFT alternative

The full set of `4n` points is an additive subspace, whose vanishing polynomial has constant nonzero derivative. Lagrange duality therefore identifies the base parity transpose with interpolation on the `n` parity points followed by evaluation on the three source cosets. No extra diagonal scalings are needed.

Use one inverse additive FFT on `3n+{0,...,n-1}`, then three forward FFTs on the source cosets. Applying the existing normalized butterfly formulas gives:

| Base map | Common inverse products / XORs | Three forward products / XORs | Base total products / XORs |
|---|---:|---:|---:|
| GRS[64,48] | 32 / 64 | 81 / 177 | 113 / 241 |
| GRS[128,96] | 80 / 160 | 209 / 449 | 289 / 609 |

Combining this base circuit with the preceding coefficient-shared shortening gives:

| Target outer | Complete products | Complete Boolean count | Full product volume |
|---|---:|---:|---:|
| GRS[64,44] | 132 | 454 | 97/44 = 2.20455 |
| GRS[128,88] | 335 | 1274 | 77/32 = 2.40625 |

Thus the straightforward FFT trades additional GFNI products for fewer XORs and a more familiar staged schedule. The coefficient-shared circuit is the better product estimate, not an established faster implementation. The existing RS64 experiment already contains bounded sixteen-vector FFT staging; it is a scheduling reference, not a transferable timing result.

## Width eight is another bounded alternative

Keep state size four and restrict each packet to two additive GF4 cosets. The eight labels form an affine three-dimensional GF2 space. Every polynomial of degree at most six sums to zero over that space. Consequently `CA=0` still holds for the restricted cubic maps. The emission remains row-scaled Cauchy, with minimum nonzero state-output weight five. Every allowed-output subset is a subset of an old sixteen-label phase, so the existing eta-scaled all-subsets invertibility result still applies. A one-input return has output weight exactly eight.

The two maps separately cost five GFNI256 instructions each. They can instead share a full-width circuit: use the four 128-bit lanes for injection coset 0, injection coset 1, emission coset 0, and emission coset 1. Both maps use the same circulant coefficients. Only the two emission lanes receive the all-ones correction. The current input and refreshed state are both available before the state update, so this packing does not assume independence between successive packets.

The packed pair therefore needs five GFNI512 instructions. Add one for the eta-scaled emission state, one for the random state diagonal, and two for the eight input diagonals. The total is nine GFNI512 instructions per eight outputs. Packing, unpacking, and horizontal state reduction still need a separate schedule audit.

| Outer | Width-sixteen full volume | Width-eight full volume |
|---|---:|---:|
| GRS[32,22] | 169/88 = 1.92045 | 185/88 = 2.10227 |
| GRS[64,44] | 357/176 = 2.02841 | 389/176 = 2.21023 |
| GRS[128,88] | 93/44 = 2.11364 | 101/44 = 2.29545 |

The width-eight change adds `2/11` products per source record, independent of outer length. It is a 9.47% increase over the current GRS[32,22] estimate. Preserve all 64 labels by splitting each old phase into two eight-label packets. The number of state-diagonal coefficient bytes doubles, from `N/4` to `N/2`. The two-map packing shares only arithmetic inside one packet; it does not add a second route or parallel state stream.

## Traffic, validation, and recommendation

At a fixed padded source size, the leading record traffic depends on `b/k`, not outer length. The old minimum schedule therefore remains `944/11` bytes per source record; the regional temporary-buffer schedule remains `1456/11`. Larger outer blocks can add local scratch traffic. Padding to 64 rows also changes the exact source size: at K20 the padded dimensions are 1048960, 1050368, and 1053184 for the three lengths. At K16 the two larger outers require about 2.13% more padded records than GRS[32,22].

Read-only scalar checks compared every entry of the new base and shortening circulants against independent interpolation: 4720 new entries. The analogous controls add 236 entries. The additive-FFT map was checked on 56 scalar coordinate basis inputs across sizes 8, 16, and 32. A further 768 checks verified the width-eight `CA=0` identities, one-input return weights, and four-point coset maps. These checks establish the stated scalar identities and coefficient counts. No larger native circuit was generated or tested.

The parent's separate K20 regional-rank screens gave the following checkpoints:

| Alternative | Active rows | Floating upper-bound margin |
|---|---:|---:|
| GRS[64,44], width 16 | 8 | +80.56 bits |
| GRS[64,44], width 16 | 16 | -67.02 bits |
| GRS[128,88], width 16 | 16 | +64.38 bits |
| GRS[128,88], width 16 | 32 | -275.46 bits |
| GRS[32,22], width 8, old two-/three-input return bounds | 16 | -154.20 bits |
| GRS[32,22], width 16, same old return bounds | 16 | -163.22 bits |

These values are exploratory upper bounds, not directed certificates or bad-code lower bounds. Increasing outer length moves the unresolved sparse row-count range; the tested width-eight change gives only a small improvement. No alternative is promoted on these results.

The later [window obstruction](WINDOW_DIMENSION_SEARCH.md) parks the four-state GRS[32,22] profile. It does not depend on the packet width or coefficient distribution, provided the saved windows remain packet-aligned. If a larger outer closes the whole sparse class, GRS[64,44]'s product increase is small enough to revisit a register-bounded circuit design. Only a later complete-call benchmark can establish competitiveness.

A smaller geometry change is GRS[33,22], at rate 2/3. Appending one parity evaluation to the existing outer adds at most seven nontrivial GF256 products per four rows: its fixed coefficients span at most eight binary dimensions, including the identity product. The outer upper count is therefore 48, and the full product count is at most `48/88 + 33/22 = 45/22`, about 6.51% above `169/88`. This is a candidate cost bound, not an implemented circuit or a distance proof. See the [independent repair comparison](SPARSE_MULTIROW_REVIEW.md).
