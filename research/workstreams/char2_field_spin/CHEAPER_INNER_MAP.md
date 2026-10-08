# Cheaper inner maps: a basis obstruction and an exact FFT alternative

Changing the eight state coordinates to additive-polynomial coefficients does not preserve the current certificate. It preserves the emission code's MDS distance, but destroys two coordinate-rank properties used by the proof. The exact residue-coordinate map nevertheless admits a packed FFT implementation with fewer GFNI instructions. That implementation changes no construction parameter or randomness.

The scalar audit is [cheaper_inner_map.py](theory/cheaper_inner_map.py). The isolated native alternative is [Cauchy256InnerFft.h](../../../spin/experiments/char2_field_spin/Cauchy256InnerFft.h). Subsequent [native testing](FIELD256_NATIVE.md) passes correctness and sanitizer checks, but the FFT implementation is slower than the dense inner.

## Which ranks the certificate needs

Fix the current packet, its phase, and the common history before revealing its input and state diagonals. The certified map is

\[
R=D_sS,\qquad Y=X+2AR,\qquad S'=R+A^TX,
\]

where `A` is the `16 x 8` residue-coordinate matrix from [STATE_CAPACITY_REPAIR.md](STATE_CAPACITY_REPAIR.md). Its entries are scaled Cauchy coefficients over AES GF256. Every square coordinate minor is nonzero, and `A^TA=0`.

The [middle-class certificate](STATE8_DENSITY_CERTIFICATE.md) uses more than the `[16,8,9]` parameters of the image of `A`:

- For selected inactive output zeros `Z` and live state coordinates `H`, the matrix `A[Z,H]` must have rank `min(|Z|,|H|)`. This gives both the diagonal root moment and the implication `|Z|<|H|` on a surviving branch.
- For a return with active inputs `I`, `C[H-complement,I]` must have rank `min(8-|H|,|I|)`. A feasible return then has eight numerical diagonal pivots. This also limits the number of incoming zero coordinates to `|I|-1`.
- The MDS emission distance supplies a separate output-weight bound. The identity `CA=0` gives invertibility of the extended packet map and its causal stream. These two properties alone do not imply either preceding coordinate-rank statement.

Thus full superregularity is sufficient at exactly the places where the current formulas use root factors and eight-pivot returns. A replacement may use weaker, placement-dependent ranks, but it requires a different packet cover and a new whole-code replay.

## Exact factorization without changing state coordinates

Write

\[
W=(0,1,92,93,224,225,188,189),\qquad V=16+W.
\]

These are the ordered eight residue nodes after translation. Let `F_o` be normalized additive-FFT evaluation of degree-below-eight polynomials on `o+W`. Translating the polynomial variable by 16 gives

\[
A=\begin{pmatrix}F_{u_0+16}\\F_{u_8+16}\end{pmatrix}F_0^{-1},
\qquad
C=(F_0^{-1})^T
\begin{pmatrix}F_{u_0+16}^T&F_{u_8+16}^T\end{pmatrix}.
\]

Here `u_0` and `u_8` are the first labels in the packet's two `W` cosets. Addition in these formulas is byte XOR. The derivative of the vanishing polynomial of `W` is one, so these are the original scaled-Cauchy coefficients without an additional normalization diagonal.

The native alternative keeps four 128-bit records in each ZMM register. An eight-record transform uses two registers. Its four-pair stage crosses those registers without a shuffle. Its smaller stages use fixed 128-bit-lane broadcasts and masked XORs. All offsets, coefficients, masks, and stage orders are compile-time constants. There is no dynamic coefficient loop, heap allocation, or conversion of persistent state coordinates.

The field transpose reverses the butterfly stages and transposes each two-input map. It does not apply a binary transpose to the GFNI multiplication instruction. Eta is charged as two additional packed products. The backward state update remains `D_s*(h+2*A^T*v)`.

## Packed instruction accounting

The following counts include both fixed maps, the fresh input and state diagonals, eta, and the six recurrence XORs. They exclude loads, stores, coefficient materialization, and compiler register moves. A Boolean operation is one XOR, masked XOR, or ternary XOR instruction; it is not a claim of equal latency across instructions.

| Sixteen-record packet | GFNI512 | Boolean | 128-bit-lane shuffles | Status |
|---|---:|---:|---:|---|
| Existing dense residue maps | 70 | 38 | 24 | Same certified construction |
| Exact residue-coordinate FFT | 34 | 53 | 46 | Same certified construction |
| New coefficient basis, symmetric pair | 28 | 40 | 32 | Coordinate-rank bounds fail |
| True coefficient-state basis change | 28 | 45 | 40 | Coordinate-rank bounds fail |

The native FFT hooks exclude the four external input-diagonal products: each hook has 30 GFNI instructions. The two fixed maps alone have 26. These counts are identical in all four phases.

The exact FFT reduces the GFNI count by more than half, but it does not halve the total listed instruction count. Serial native comparisons show that extra lane shuffles and dependency chains outweigh the arithmetic saving: full encoding regresses about 8% at K16/K18 and 5% at K20. The dense inner remains selected.

## Concrete coefficient-basis failures

Let `M=F_0`. The faster emission matrix is `A'=AM`, which omits the interpolation of the eight state values. Consider fresh independent diagonal refresh in these new state coordinates. This changes the ensemble: an arbitrary basis change does not commute with a diagonal random matrix.

Two injection choices are inexpensive. The symmetric choice is `C'=(A')^T`; degree-fourteen power sums still give `C'A'=0`. The true state-coordinate change instead uses `C'=M^-1 C`. It retains the original product `AC`. Its extra coefficient-space map is

\[
G=M^{-1}M^{-T},\qquad G_{i,j}=\mathbf1_{i\mathbin{\&}j=0}.
\]

This is the third tensor power of `[[1,1],[1,0]]`. It needs only XORs and lane permutations. Neither injection choice repairs the following defects.

For every phase, two equal inputs in positions `(0,1)` produce an injection with exactly four zero coordinates. Under the symmetric choice these are `(0,2,4,6)`; the live coordinates are `(1,3,5,7)`. The return matrix on those four live state coordinates and the two active inputs has rank **five**, not eight. The true-coordinate-change injection also has a rank-five return witness, with the complementary coordinate pattern. Therefore the old uniform eight-pivot charge is false for both variants.

For phase zero, state coordinates `(0,4)` with values `(1,80)` produce zeros in output positions `(0,1,2,3)`. The other phases have the same four-zero pattern, with ratios `188`, `170`, and `55`. If those two incoming coordinates are proportional live forms, the four selected equations impose only one fresh diagonal ratio. Their survival probability is `1/255`, whereas the old rank-one/four-root formula bounds it by `1/255^4`. The new basis therefore invalidates the current diagonal root moment, not merely a convenient proof of that moment.

All entries of these emission matrices are nonzero. Moreover, the scalar check exhausts all `binom(16,8)=12870` maximal minors in each phase and finds full rank. The emission code still has distance nine. The failure concerns smaller coordinate minors.

For comparison, the systematic basis `A'=[I;T8]`, where `T8` is the half-coset interpolation matrix, is still cheaper. It has `C'A'=0` and the same MDS image. But an input in the first half injects into only one state coordinate, so a matching singleton return has rank one. This is not a certificate-compatible replacement either. These local failures do not by themselves prove the existence of a low-weight whole-code word.

## Validation and next gate

Run:

```sh
python -B research/workstreams/char2_field_spin/theory/cheaper_inner_map.py
```

The replay passes 51,576 geometry and original-map basis checks, plus 2,688 lane-level FFT, dense-map, and extended-adjoint checks. It also verifies the true basis change preserves `AC`, records the deficient small minors, and checks both variants' numerical return matrices. The lane interpreter implements the packed broadcasts and masks independently of the scalar shear circuit.

The native header exports `c256InnerTransposeFftPacked<Phase>` and `c256InnerForwardFftPacked<Phase>`, with the existing packed input, output, state, and diagonal arguments. Input and output may alias. Native whole-pipeline tests and ASAN/UBSAN pass; the measured performance gate fails.

A new coefficient-basis certificate is lower priority: its inexpensive arithmetic comes with explicit losses of three return pivots and three powers of diagonal root probability. A six-state alternative should first pass its own whole-code proof gate; restricting the original residue matrix to six columns preserves superregularity and `CA=0`, but changes the return exponent from eight to six. The [updated screens](CHEAPER_DESIGN_STATUS.md) still have a middle-class gap for that alternative.
