# Sparse certificates during the GF256 rate climb

This record concerns setup failure for the fixed full-Cauchy width-eight,
state-eight construction at requested dimension `K=2^20`.
The target remains 10% symbol distance and at most `2^-40` setup failure.
Each 128-bit record carries sixteen parallel GF256 lanes of the same sampled
scalar code. The argument does not add sixteen setup-failure probabilities.

The construction samples one independent uniform regional permutation per
outer coordinate. Input and state diagonals contain independent uniform
nonzero GF256 entries. Each packet retains a fixed prefix of its eight outputs.
A periodic variant varies the prefix length in a fixed chronological pattern;
it adds no permutation or randomness.

The directed sparse calculations below are complete for their stated active-row
ranges. They are not whole-code certificates by themselves. The generic
aggregate must also verify the middle and dense components, including the
omitted-output completion factor.

## Directed results

The outer constituent is GRS `[128,d]`. Write `A` for its number of active rows.
All results use the same concrete packet matrices, whose SHA-256 digest is
`6e3cf3e5d6b5876838afc5b2d7230f931974f363f9f67b402843f677521450a1`.

| Outer dimension `d` | Retained outputs per packet | Nominal rate | Complete `A` range | Sparse margin | Directed intervals |
|---:|---|---:|---:|---:|---:|
| 88 | 7 | `11/14` | 1–1120 | 79 bits | 317 |
| 80 | 6 | `5/6` | 1–1028 | 264 bits | 44 |
| 67 | 5 | `67/80` | 1–971 | 368 bits | 40 |
| 84 | repeat `(6,6,6,7)` | `21/25` | 1–1070 | 115 bits | 57 |
| 74 | repeat `(5,6)` | `37/44` | 1–991 | 312 bits | 40 |

The margins bound the sum over every active-row count in the stated range.
The interval partition digests, in table order, are:

- `ce15f097cae15183872a14ae8fa3752a6731fa52eebab02c1092c3f8ca33f64b`
- `f76f8cbb7482c1fc562ee17b0c35e0f8d611347751e81a5edab9e8eadb81f39b`
- `fd4a6612693abf0436ab6ad3a481ba62e9852b0b8ccb41c2af6b021a427f2eb4`
- `e6066f398a2433f20f9b761361515b3c89fe4228d05b2b90c5ee3b8928618d83`
- `bff62cfe77bc98c86f283dd6d5e3406e5afb0d12c023c659bc508e1d3b01b9ee`

Padding changes the requested-dimension rate. For example, the 84% periodic
construction has padded dimension 1,053,696 and output length 1,254,400.
Its requested-dimension rate is `1024/1225`, while the padded code has rate `21/25`.
The requested-dimension subcode inherits the padded code's distance bound.

For the period `(5,6)` construction, the nominal rate is `37/44`, approximately
84.0909%. The GF256 random-code/Gilbert–Varshamov benchmark at relative distance
0.1 is `1-H_256(0.1)`, approximately 84.1446%. The nominal gap is about 0.0537
percentage points. This benchmark is not a universal upper bound on code rate.
At the exact requested K20 size, 2,816 padded input symbols reduce the rate to
`1024/1221`, approximately 83.8657%. Its gap to the same benchmark is about
0.2789 percentage points. These finite-size and nominal comparisons should not
be interchanged.

## Sparse-to-middle coverage

Each active GRS row has at least `129-d` nonzero symbols. Therefore, a sparse
cover through `A_max` leaves only words of weight at least
`(129-d)(A_max+1)` for the middle or dense components.

| Profile | Raw input length | Proposed first middle weight | First weight not covered by the sparse row range |
|---|---:|---:|---:|
| `d=88`, keep 7 | 1,531,904 | 45,958 | 45,961 |
| `d=80`, keep 6 | 1,679,360 | 50,381 | 50,421 |
| `d=67`, keep 5 | 2,007,040 | 60,212 | 60,264 |
| `d=84`, period `(6,6,6,7)` | 1,605,632 | 48,169 | 48,195 |
| `d=74`, period `(5,6)` | 1,818,624 | 54,559 | 54,560 |

These joins have overlap rather than a gap. Every proposed first middle weight
is the ceiling of 3% of the raw input length.

## What the periodic replay changes

Fix rational marks `0<y<=1` and `1<=tau<=255`. Let `T_{h,w}` be the audited
two-state local transfer matrix for a packet with input weight `w` and `h`
retained outputs. Its states distinguish a zero and nonzero recursive state.
Define its matrix polynomial by

\[
P_h(z)=\sum_{w=0}^{8}\binom8w T_{h,w}z^w.
\]

For a public pattern `(h_0,...,h_{p-1})`, the replay multiplies these polynomials
in chronological order. A region with `L` positions uses

\[
P(z)=\bigl(P_{h_0}(z)\cdots P_{h_{p-1}}(z)\bigr)^{L/(8p)}.
\]

The metadata requires `8p` to divide `L`. Thus every region starts at the same
public phase, without resetting the recursive state. The proof neither
averages the phases nor interchanges their matrix factors.

For `A` active outer rows, the regional matrix is

\[
R_A=\sum_{h=0}^{A}
\frac{\binom Ah}{\binom Lh}\tau^h[z^h]P(z).
\]

The existing source and terminal-state argument bounds the class failure by

\[
\binom LA\tau^{-(128-d)A-1}y^{-D}
e_0 R_A^{128}\mathbf1,
\]

where `D=ceil(N_out/10)-1`. The bound therefore excludes every codeword whose
weight is strictly below 10% of the retained output length.

All polynomial coefficients are nonnegative. Truncating degrees above the
largest tested `A` retains every coefficient used by the regional sum.
For an interval of active-row counts, the regional factor is evaluated at its
upper endpoint, while the source factor uses its own exact interval maximum.
Rational local matrices and outward Arb arithmetic determine acceptance.
Floating-point optimization only proposes marks.

## The keep-seven stopping point

With seven outputs retained in every packet, increasing the outer dimension
from 88 to 89 makes the current sparse bound fail. Reoptimization gave:

| Outer dimension | `A=128` | `A=256` | `A=512` |
|---:|---:|---:|---:|
| 88 | +227.16 bits | +156.02 bits | +1078.43 bits |
| 89 | -112.50 bits | -510.03 bits | -266.84 bits |
| 90 | -446.80 bits | -1165.24 bits | not rerun |

These are floating-point screens of one upper bound. A negative result is not
a low-distance codeword or a proof that the construction fails. Retaining fewer
outputs permits a stronger outer constituent at the same final rate and removes
this particular sparse obstruction in the tested profiles.

## Reproduction

Run these commands from the adjacent `theory` directory:

```text
python -B gf8_rate_climb_sparse.py --dimension 88 --outputs 7 --last 1120
python -B gf8_rate_climb_sparse.py --dimension 80 --outputs 6 --last 1028
python -B gf8_rate_climb_sparse.py --dimension 67 --outputs 5 --last 971
python -B gf8_periodic_climb_sparse.py
python -B gf8_periodic_climb_sparse.py --dimension 74 --pattern 5 6 --last 991
python -B gf8_rate_climb_sparse_tests.py
python -B gf8_periodic_climb_sparse_tests.py
python -B gf8_rate_climb_certificate_audit_tests.py
```

The five fixed-profile tests, four periodic tests, and seven aggregate audit
tests pass. The periodic tests compare the directed polynomial coefficients
against an independent exact rational product. Aggregate tests reject actual
coverage gaps, omitted completion factors, phase mismatches, and unsupported
failure margins. Their synthetic fixtures test the combiner, not code distance.

The frozen outer-84/keep-seven helpers and certificates remain unchanged.
No native implementation, timing, or performance claim is part of this record.
