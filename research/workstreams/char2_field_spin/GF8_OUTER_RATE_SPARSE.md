# Variable-rate GF256 outer: sparse classes

Increasing the dimension of the length-128 outer reduces redundancy. This
replay checks which active-row classes retain a useful bound when the
width-16, state-8 inner and its setup distribution remain fixed.

At requested K18, dimension 88 now has a complete sparse-class bound below
`2^-72` for active-row counts `1..375`. The separate middle and dense
components are still required for a whole-code theorem. This component
does not establish a construction failure at dimensions where its bound
is too large.

## Fixed construction and changing geometry

The base field is GF256. Each outer row belongs to the fixed GRS code
`[128,k]`; its minimum distance is `129-k`. Each of 128 regions contains
one coordinate from every row. Setup independently samples a uniform row
permutation in each region and the nonzero input/state diagonals of the
previously certified state-eight inner.

For requested source length `2^ell`, let `L` be the least multiple of 64
with `kL >= 2^ell`. The padded source length is `kL`, the output length is
`N=128L`, and the forbidden output weight is `e=ceil(N/10)-1`.

Let `A` count the nonzero outer rows, and put `p=128-k`. For fixed rational
marks `0<y<=1` and `1<=tau<=255`, the replay evaluates

    F_h = [x^h] (sum_{w=0}^{16} binom(16,w) x^w T_w)^(L/16),
    R_A = sum_{h=0}^A binom(A,h)/binom(L,h) tau^h F_h,
    B_A = binom(L,A) tau^(-pA-1) y^(-e) e_0 R_A^128 1.

Here `T_w` is the established joint-return packet cover. The last vector
contains ones. This is the existing source-rank comparison, not an exact
word enumerator. Only the outer redundancy and finite geometry change.

For an interval `a..b`, every entry of `R_A` increases with `A`. The replay
uses `R_b`, maximizes the exact source factor over `a..b`, and multiplies
by the interval size. Adaptive subdivision proves a complete interval
union; sampled points are not used as coverage evidence. The final bound
is the exact sum of outward-rounded dyadic upper endpoints.

## K18 dimension 88: completed sparse component

The finite geometry is `L=3008`, `K_pad=264704`, `N=385024`, and `e=38502`.
The replay covers `A=1..375` in 67 disjoint adjacent intervals. Its exact
dyadic sum is below `2^-72`. The partition hash is

`0f633ee9ea876ec634d2fcf0916117c6610df7bd75c9f9ddb504a10f5f37ea38`.

Outside this sparse range, the outer minimum distance gives
`W>=41*376=15416`. Therefore a middle component beginning at `W=15401`
overlaps it. The whole-code aggregator must check this inequality, the
middle/dense join, and the exact sum of all three component bounds.

```text
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_sparse.py --dimension 88 --k-log 18 --cover 375 --precision 4096 --cell-bits 64
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_sparse_tests.py
```

The seven tests include the complete replay, agreement with the frozen
dimension-87 bound, unchanged finite geometry, polynomial and matrix-power
controls, the exact continuation cone, and receipt scope.

## K20 obstruction: intermediate active-row counts

The first sparse row is not the rate bottleneck. For dimension 88 at
K20, its directed bound is below `2^-200`. The difficult classes have
roughly two hundred active rows. Optimizing the established two-state
comparison gives these floating point margins:

| Active rows | 128 | 192 | 256 | 320 |
|---:|---:|---:|---:|---:|
| Margin in bits | 154.17 | 26.55 | 4.59 | 72.34 |

These are valid parameter searches, not complete certificates or lower
bounds on the construction's failure probability. They show that merely
retargeting the two existing marks does not reach 40 bits near the trough.

A wider directed screen uses the frozen dimension-87 marks. At K20,
dimension 89 gives negative margins at `A=128,256,512`; higher sampled
dimensions are worse. At K18, dimension 89 passes the sampled sparse
classes, but it has the same padded output length as dimension 88.
Dimension 90 at K18 already has a negative frozen bound at `A=128`.

## Proof-only rank-history control

The optional nine-state transfer retains the incoming rank cap. State
zero means the actual zero witness state; states `1..8` mean a nonzero
witness state with that rank cap. This changes the proof, not the encoder.

Its domination contract applies to nonnegative continuation functions
that increase on states `1..8`; the value at state zero is independent.
Each local packet operator preserves this cone. Every coefficient `F_h`
is a positive sum of chronological products of such operators. The
regional weighting also uses nonnegative coefficients, so regional
composition preserves the contract. The terminal all-ones vector belongs
to the cone. No entrywise domination of an exact rank-transition law is
assumed.

The generic replay uses nine-by-nine matrices of Arb polynomials. It
truncates only degrees above the active-row endpoint. This leaves every
coefficient used by the regional comparison unchanged. Independent reviews
checked the composition order, cone use, and interval monotonicity.

At K20, dimension 88, and the frozen marks, nine-state directed bounds
give 44 bits at `A=192` and 23 bits at `A=256`. The corresponding two-state
bound at `A=256` has floor -1 bit. The gain is about 24 bits, but the
remaining gap prevents a 40-bit certificate. No full rank-history sparse
cover is claimed.

One final control reuses the optimized two-state marks
`(logtilt,logtau)=(6.307819162762217,1.8517167647177386)` at `A=256`.
The directed nine-state margin rises to 30 bits. This is a bounded proof
improvement, but it still misses the target before any union over classes.
The rational marks in this replay are
`y=1923721458/2108508031` and `tau=617898613/96989975`.

```text
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_sparse.py --dimension 88 --k-log 20 --active 192 256 --rank-history --precision 4096
```

The next useful task is the K18 whole-code aggregate. Further K20 work
needs a sharper bound near the sparse trough, not another check of the
single-row class. No encoder implementation or timing changed here.
