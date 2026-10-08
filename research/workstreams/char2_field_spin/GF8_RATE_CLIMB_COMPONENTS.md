# Middle and dense certificates near the GF(256) rate limit

The first punctured checkpoint used a `[128,84]` outer and retained seven of eight packet outputs, giving nominal rate 3/4.
This round raises the rate without changing the full-Cauchy inner or adding another permutation.
It changes the outer dimension and the deterministic output-puncturing pattern.

All guarantees below concern GF(256) symbol distance.
They do not assert 10% binary bit distance.
Sixteen parallel byte lanes may share one sampled scalar code; the corresponding 128-bit tuple code inherits its nonzero-coordinate distance.

## Construction and scope

Each outer row uses a fixed GF(256) MDS code `[128,d]`.
One independent uniform permutation of the rows acts within each of the 128 outer-coordinate regions.
The inner processes eight symbols at a time, with eight persistent state symbols.
Its full maps and setup distribution are those in [the punctured checkpoint](GF8_PUNCTURED_COMPONENTS.md).
The matrix digest remains
`6e3cf3e5d6b5876838afc5b2d7230f931974f363f9f67b402843f677521450a1`.

A puncturing pattern `(m_0,...,m_(P-1))` retains the first `m_j` outputs of packet phase `j`.
The pattern repeats in chronological order from packet zero.
The state is not reset at phase or region boundaries.
Input and state diagonals are independent uniform nonzero field scalars, sampled once and fixed for every message.
The initial state is zero and the final state is discarded.

For `L` outer rows, the raw and retained lengths are

\[
N_{\rm in}=128L,\qquad
N_{\rm out}=\frac{N_{\rm in}}{8P}\sum_{j=0}^{P-1}m_j.
\]

The nominal rate is the padded dimension `dL` divided by `N_out`.
The requested-dimension rate can be smaller because natural row counts require padding.
The certificates cover the entire padded code, so shortening to the requested dimension preserves the distance guarantee.

## Closed component intervals at requested K = 2^20

Each component below has been replayed with directed arithmetic.
The final column lists sparse, middle, and dense component margins, not a single aggregate margin.
The aggregate must add their failure bounds.

| Outer dimension | Kept-output pattern | Nominal rate | Padded dimension | Retained length | Component margins, bits |
|---:|---|---:|---:|---:|---|
| 80 | `(6)` | `5/6` | 1,049,600 | 1,259,520 | 264 / 17,936 / 5,480 |
| 67 | `(5)` | `67/80` | 1,050,560 | 1,254,400 | 368 / 51,673 / 1,395 |
| 84 | `(6,6,6,7)` | `21/25` | 1,053,696 | 1,254,400 | 115 / 8,825 / 1,734 |
| 74 | `(5,6)` | `37/44` | 1,051,392 | 1,250,304 | 312 / 32,550 / 184 |

The last profile has nominal rate about 84.091%.
Its requested-dimension rate is `1024/1221`, about 83.866%.
The GF(256) random-code reference `1-H_256(0.1)` is about 84.145%; this comparison uses nominal rate.
The result is a concrete construction and failure bound near that reference, not a claim of optimality.

The exact coverage interfaces are:

| Profile | Sparse active rows | Middle raw weights | Dense raw weights |
|---|---|---|---|
| `80/(6)` | 1–1,028 | 50,381–1,343,487 | 1,343,488–1,679,360 |
| `67/(5)` | 1–971 | 60,212–1,806,335 | 1,806,336–2,007,040 |
| `84/(6,6,6,7)` | 1–1,070 | 48,169–1,364,787 | 1,364,788–1,605,632 |
| `74/(5,6)` | 1–991 | 54,559–1,636,761 | 1,636,762–1,818,624 |

Every active outer row contributes at least `129-d` raw nonzero symbols.
For each row of the table, `(129-d)(A_max+1)` is at least the first middle weight.
Thus the sparse receipt covers every remaining nonzero source below the middle interval.
The middle and dense intervals meet at consecutive integers.

Punctured kernel messages are included.
Their retained output is zero, which belongs to the forbidden output ball.
The dense proof includes all full-output completions of that zero retained output.
The all-zero full completion has only zero raw input because the full causal inner is invertible.
It therefore cannot represent a nonzero source.

## Directed middle replay

Let `M_m(x,y,tau)` denote the positive two-state transfer for a packet retaining `m` outputs.
The period transfer is

\[
M_{m_0}(x,y,\tau)\cdots M_{m_{P-1}}(x,y,\tau).
\]

The replay forms this product over exact rational numbers, then raises it to `N_in/(8P)` with outward rounding.
The matrices generally do not commute.
An averaged matrix, or a rearranged pattern, would describe a different bound.
The output mark uses the retained length and the strict forbidden-weight threshold `ceil(N_out/10)-1`.

The source polynomial uses the actual outer parity `128-d` and minimum support `129-d`.
Its marked support sequence is discretely convex when the combined support mark is at least two.
The existing relaxed polynomial is used otherwise.
Each density cell maximizes its positive source terms separately, including their possible interior maxima.
Consequently sampled optimizer success is not an acceptance condition for the directed proof.

The helper `gf8_rate_climb_middle_certificate.py` supports fixed patterns and exact chronological products.
The wrapper `gf8_periodic_climb_middle_certificate.py` supplies frozen rational witnesses for both periodic profiles.
It covers 820 cells for `84/(6,6,6,7)` and 870 cells for `74/(5,6)`.
Their respective failure upper bounds are

\[
2.274753844331692\cdot10^{-2657}
\quad\text{and}\quad
1.995482685218021\cdot10^{-9799}.
\]

These displayed decimals are rounded upward from the full directed receipts.
The corresponding partition digests are
`4168365351436376ec9c73ff4775eb3eb56a7b2fe60c80c71bb14e6fdfc44526`
and
`aa679e45f42e36af5a41ba4d386dd207c1a7737e7005da429c08bdc7db8cbd44`.

## Dense replay and the final refinement

For a fixed retained output, the dense proof includes `256^(N_in-N_out)` full-output completions.
It conditions on the state diagonals, inverts the full causal inner, and leaves the input diagonals and regional route unexposed.
Each completion has one raw inverse-input weight under this conditioning.
The proof therefore takes the largest verified weight-cell bound, then applies one output-ball and completion union.
It does not add a second union over weight cells.

`gf8_rate_climb_dense.py` and `gf8_periodic_climb_dense.py` close the first three profiles with the existing relaxed membership theorem.
The last profile needs the exact MDS support spectrum to remove the accumulated nonzero-diagonal relaxation.
The separate helper `gf8_exact_spectrum_dense.py` supplies that refinement.
It bounds the positive excess above the constant random-subspace membership term with Maclaurin's inequality on forced-zero coordinates.
This bound applies to unequal regional support probabilities; it does not assume the actual routed symbols are independent.

At `74/(5,6)`, the exact-spectrum dense receipt starts at raw density `.9` and covers all 181,863 remaining weights.
Its 21 cells give 184 bits after accounting for the retained-output ball and all omitted-output completions.
The proof and its exact small-field/profile tests received independent audits.

## What changed in the search

Increasing only the outer dimension while retaining seven outputs reached a low-density obstruction in the available bound.
For outer dimension 88, density `.02` had approximately -78 bits, although density `.03` had about 2,692 bits.
At dimension 90 those values were approximately -3,067 and -1,840 bits.
The sparse dip also worsened between dimensions 88 and 89.
These negative first-moment bounds do not establish low actual distance.

Stronger outers followed by deterministic puncturing gave substantially better proved rates.
Periodic puncturing supplied the intermediate rates between fixed output counts without extra routing or randomness.
Near the final rate, exact outer membership—not another inner change—removed the remaining dense-bound loss.

The frozen 3/4 checkpoint and its scripts remain unchanged.
No timing claims follow from these certificate replays.
The next implementation step is a correctness-checked periodic output schedule, followed by a serialized end-to-end benchmark against the fixed-output checkpoints.

```text
python research/workstreams/char2_field_spin/theory/gf8_rate_climb_middle_dense_tests.py
python research/workstreams/char2_field_spin/theory/gf8_periodic_climb_middle_dense_tests.py
python research/workstreams/char2_field_spin/theory/gf8_rate_climb_certificate.py --profile periodic74
```
