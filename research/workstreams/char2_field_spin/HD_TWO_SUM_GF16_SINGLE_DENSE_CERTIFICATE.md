# Directed r1 dense certificate with full state matrices

The small-cohort certificate excludes one class: a single active source group containing at least three nonzero constituent words.
For this class, retaining the two-state matrices through all outer columns avoids the loss from a common-potential bound.

The construction is `two-sum-gf16-coset-g16-rs16x8-h4-swapped-pairs-k20-rate-half`.
It uses the g16 two-sum mixer with labels $06+\operatorname{GF}(16)$ inside AES $\operatorname{GF}(256)$, sixteen RS$[16,8]$ constituents, the alternating independently swapped pair route, and the existing h4 DHD inner map.
The probability space is its ideal random setup distribution.
The selected scalar message class has exactly one active source group and a fixed nonzero constituent row set $A\subseteq[16]$ of arity 3 through 16.
The source group and row set are fixed before setup randomness is averaged.

The calculation covers all 8,192 source groups and all 65,399 dense row sets.
It does not cover two active groups or change the scope of any earlier receipt.
It is not a whole-code certificate.

## Ordered regional matrices

In each region, this class has either zero or one active position.
There are no full pairs under either matching.
Let $R_i(0)$ and $R_i(1)$ be the directed two-state regional matrices for those occupancies.
They average the independent regional pair permutation, fair swap, nonzero output amplitude, and inner setup as in the established paired-route operator.
The source group is not resampled between regions or outer columns.

For $w=0,\ldots,16$, define

$$
Q_w=\sum_{\substack{T\subseteq[16]\\|T|=w}}
R_0(1_T(0))R_1(1_T(1))\cdots R_{15}(1_T(15)).
$$

The coefficient recurrence multiplies matrices in this chronological order.
It does not commute them.
Each labelled mask occurs once, so there is no additional factor $\binom{16}{w}$.

## Fixed-arity comparison and column floor

For a fixed row set $A$, let $N_A(T)$ be the exact uniform-input output-support count of the two-sum mixer.
The comparison allows zero coordinates within $A$.
The closed support kernel gives

$$
S_a(w)=\sum_{|A|=a}N_A(T)^{16},\qquad |T|=w,
\qquad
W_a(w)=\left(\frac{S_a(w)}{255^{8a}}\right)^{1/16}.
$$

Hölder's inequality bounds the sum over fixed row sets by a product of the $W_a$ weights, while keeping each row set fixed across all 16 columns.
The denominator includes the RS prescribed-support comparison.
No additional row-set count or density factor is inserted.

For arity $a$, form two nonnegative matrix coefficients:

$$
B_{a,0}=W_a(0)Q_0,\qquad
B_{a,1}=\sum_{w=1}^{16}W_a(w)Q_w.
$$

The dense group must be nonempty in at least nine outer columns.
Indeed, it contains a nonzero RS$[16,8]$ word, and the column mixer and diagonals are invertible.
Let $e_0=(1,0)$ and $\mathbf1=(1,1)^T$.
The retained marked transfer is

$$
T_a(z)=e_0\left(\sum_{j=9}^{16}[Y^j](B_{a,0}+YB_{a,1})^{16}\right)\mathbf1.
$$

The implementation expands this positive matrix polynomial by ordered multiplication.
It does not replace matrix words by a scalar binomial formula.
The initial state is zero, and both terminal states are counted directly.
There is no common state potential and no endpoint-potential factor.

For each fixed arity, the failure bound for output weight at most 734,003 is

$$
\frac{8192}{255}\,z^{-734003}T_a(z).
$$

The source group count occurs once.
Projectivization contributes $255^{-1}$.
Summing these bounds over the 14 arities covers the selected class.

## Frozen marks and directed replay

The selected proposal uses this deterministic schedule:

| Arity | Mark $z$ |
| --- | --- |
| 3–5 | $1999/2000$ |
| 6–10 | $1199/1200$ |
| 11–16 | $499451537/500000000$ |

The mark is fixed for the whole arity class, before the setup is sampled.
No extra mark-bank factor is needed in the final union.
The proposal `hd_two_sum_gf16_single_dense_selected_float.json` records the selected rows and pins both source artifacts.
Earlier failed potential and mark controls remain separate evidence; they are not treated as certificates.

`theory/hd_two_sum_gf16_single_dense_certificate.py` pins this proposal and its reviewed dependencies.
It rebuilds the three regional banks from exact rational marks, computes the positive ordered matrix coefficients, and evaluates every selected arity.
Each dyadic upper endpoint $u$ for a root is checked exactly by
$u^{16}255^{8a}\geq S_a(w)$.
All precision-sensitive operations execute in one Arb context.
Saved FLOAT bounds are not consumed by the replay.

Eight focused tests pass. They include exact labelled-mask enumeration, independent noncommuting matrix words, an unreachable-live-state control, once-only normalization, strict scope and schedule checks, and an independent rounded-union check.
The complete lightweight `hd_two_sum*_tests.py` regression also passes: 142 tests in 6.536 seconds.
Two independent source reviews cleared the frozen adapter before serial replay.

| Replay | Arithmetic time | Log₂ upper bound |
| --- | ---: | ---: |
| 256-bit fresh replay | 0.6688 s | $-26611526601709/(5\cdot10^{11})$ |
| 384-bit fresh replay | 0.6617 s | $-26611526601709/(5\cdot10^{11})$ |
| 512-bit sum of rounded arity ceilings | — | $-26611526601709/(5\cdot10^{11})$ |

Every arity ceiling and the raw aggregate agree between the two fresh replays.
The independent sum consumes the 14 rounded arity bounds, not either saved aggregate.
Thus this selected class has failure probability at most $2^{-53.223053203418}$ for output weight at most 734,003.
The probability is over the ideal construction's setup; this is not an empirical distance measurement.

The artifacts are in `theory/`:

| Artifact | SHA-256 |
| --- | --- |
| `hd_two_sum_gf16_single_dense_certificate.py` | `e08f857a8c7a888b821a20a962600e718e0768290e21c6502abedd81629fac10` |
| `hd_two_sum_gf16_single_dense_certificate_tests.py` | `01c552e2a4fd8a726d00f2785eed7cd4d67ca0b3c23363c53af7f9e58ad4a765` |
| `hd_two_sum_gf16_single_dense_receipt256.json` | `028192f0c09a303d8985fe3d5a85a7bb5f996a78998170b8ad26b02d4abb842f` |
| `hd_two_sum_gf16_single_dense_receipt384.json` | `c0ab40feeee3d4efd11098371f0bbf0df513db3728c04725786ec3286a0d4da4` |
| `hd_two_sum_gf16_single_dense_dual_verified.json` | `3cd4b536f4c2d0bacd417cd9da4b6bb8603e311f8da95fb925f40657a63e2732` |

Hashes record provenance, not an independent proof of arithmetic execution.
An independent 1024-bit calculation sums the 14 rounded arity ceilings and adds the disjoint 322-cell rounded union.
Both the $r=1$ dense sum alone and this combined sum have ceiling $-26611526601709/(5\cdot10^{11})$.
The combination covers all zero/one-dense cohorts for $1\leq r\leq17$, but still excludes cohorts with two or more dense groups and active-group counts above 17.
The next proof step is to retain ordered matrices for broader cohorts before considering another construction change.
