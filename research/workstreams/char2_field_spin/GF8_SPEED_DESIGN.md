# Smaller outer and packet for the high-rate GF256 code

The speed candidate halves both the outer length and the inner packet size.
It uses GRS[64,40], a four-input/four-state Cauchy inner, and three retained
outputs per packet. Its nominal rate is 5/6. This is a new research profile;
the earlier [rate checkpoints](GF8_RATE_FRONTIER.md) remain unchanged.

At requested dimension `K=2^20`, the complete replay certifies 10% symbol
distance with 57 bits of whole-code setup-failure margin. The sparse component
determines that margin; middle and dense components have much larger margins.
The certificate concerns GF256-symbol distance, not binary Hamming distance.
Native correctness and timing are separate from these proof receipts.

## Fixed construction

Use GF256 with polynomial `0x11b`. Pad the source to 26,240 rows of 40
symbols, then encode each row with the same fixed GRS[64,40] constituent.
Collect each of the 64 outer coordinates into a region. Independently sample
a uniform permutation of the row labels in each region. Concatenate the
regions, then process four symbols per inner packet.

All randomness is sampled once and reused for every message. Each packet
has independent diagonal matrices `D_x` and `D_s`, whose entries are uniform
over the 255 nonzero field elements. These diagonals and the regional
permutations are mutually independent. The state starts at zero and continues
across region boundaries.

For raw packet input `x` and incoming state `q`, both in GF256^4, define

\[
 X=D_xx,\qquad R=D_sq,\qquad
 Y=FX+A_eR,\qquad q'=R+CX.
\]

Retain coordinates 0, 1, and 2 of `Y`. Discard coordinate 3 and the terminal
state. Each matrix below is a four-point XOR-circulant: entry `(i,j)` is
entry `i XOR j` of its first row.

| Matrix | First row, in decimal field-byte representation |
|---|---|
| `F` | (86, 90, 226, 133) |
| `A_e` | (30, 235, 104, 156) |
| `C` | (27, 28, 18, 20) |

The implementation may use the exact identity

\[
 F=107A_eC,\qquad Y=A_e(R+107CX).
\]

Here 107 denotes a field element, not an integer multiplication. Both state
maps have unit row sums. Their four-point circulant circuits replace the
larger eight-point maps without changing the stated packet map.

The maps come from the full 8-by-8 Cauchy matrix
`T[i,j] = 1 / (i XOR (128+j))`, for `0 <= i,j < 8`.
Partition `T` into four 4-by-4 blocks. Normalize the bottom rows by the
inverse of the bottom-right block, then scale the input and output blocks.
The resulting matrix is `[[F,A_e],[C,I]]`.
The bottom-row normalization preserves rank when all four state rows are
selected; it does not assert that every minor of the normalized matrix is
nonzero.

The receipt stores `A=A_e/2` and `eta=2`. Its map digest covers row-major
`F,A,C`, followed by `eta`:

`ecd9399dee2d87cde4c12a66704a7884d7a590d38603a3de18cfeb07292fabd3`.

## Geometry and distance event

| Quantity | K20 profile |
|---|---:|
| Requested dimension | 1,048,576 |
| Padded dimension | 1,049,600 |
| Outer rows | 26,240 |
| Raw outer length | 1,679,360 |
| Retained output length | 1,259,520 |
| Nominal rate | 83.3333% |
| Requested-dimension rate | 83.2520% |
| Forbidden output weights | 0 through 125,951 |
| Certified minimum distance | 125,952, exactly 10% |

The bad setup event is the existence of a nonzero padded message whose
retained output has forbidden weight. In particular, this event includes
nonzero messages in the punctured encoder's kernel. No injectivity assumption
is made after puncturing. The requested-dimension code is a shortened subcode,
so the padded-code guarantee also applies to it.

For the once-sampled distribution above, the directed certificate proves
that this bad event has probability at most `2^-57`. Equivalently, with
probability at least `1-2^-57`, every nonzero padded message produces at least
125,952 nonzero retained symbols.

A 128-bit record contains sixteen lanes of the same sampled scalar code.
Every nonzero record message has a nonzero scalar lane. Its encoded lane
support is contained in the encoded record support. Thus the same setup
event controls 10% nonzero-record distance, without a sixteen-fold union.

## Complete proof coverage

Let `A` count active outer rows and let `W` be the raw outer word's weight.
The GRS constituent has minimum distance 25, so `W >= 25A`.

| Component | Complete covered range | Margin | Directed cells |
|---|---|---:|---:|
| Sparse | `1 <= A <= 2015` | 57 bits | 255 |
| Middle | `50,381 <= W <= 1,343,487` | 10,101 bits | 770 |
| Dense | `1,343,488 <= W <= 1,679,360` | 8,052 bits | 40 |

Every source outside the sparse range has `W >= 25*2016 = 50,400`.
The middle range starts at 50,381, and its upper endpoint immediately
precedes the dense range. These ranges exhaust all nonzero padded sources.

The sparse and middle proofs use the positive two-state transfer for the
Cauchy packet map. For any selected `z` output rows and any input/state
column support of size `w+h`, the selected-output rank is `min(z,w+h)`.
After adding all four return-to-zero equations, the rank is
`min(z+4,w+h)`. Exact enumeration verifies all 8,192 such rank cases.
The certificate also checks every square minor of the birth, return,
emission, and injection maps.

The sparse replay uses exact rational local coefficients and outward Arb
arithmetic at 4,096-bit precision. It retains the active-row count and uses
64 regional factors. Separate fixed witnesses around the sparse dip avoid
unnecessary interval slack; they do not change the proof or the construction.
The middle replay uses 80-digit directed decimal arithmetic, 64 central-mode
conditioning factors, and the length-64 source polynomial.

For the dense proof, first fix the state diagonals and a completed full
output. The full causal inner is invertible because `F` is invertible.
Its inverse determines a word independently of the input diagonals and
regional permutations. The exact-MDS membership bound applies to that word.
Puncturing omits 419,840 symbols, so the union over completions contributes
`256^419840 = 2^3358720`. The retained output ball contributes at most
`2^(1.27*1259520)`. The dense bound is uniform over its weight range; the
cell partition does not introduce another union over possible weights.

The aggregator adds the three directed failure bounds and verifies a total
of at most `2^-57`. It checks the matrix
identity, all common metadata, both coverage joins, the completion factor,
and the strict output threshold. Its replay is the authority for the
whole-code margin, rather than the minimum of three rounded component margins.

Component partition hashes are:

| Component | SHA-256 |
|---|---|
| Sparse | `ffa7d038c7e68bd2400c3847bd7cb93459cfcf9da0cade842cb6122584cf2c01` |
| Middle | `2331ad47429d48eebbc27894fab66ded5014bbff41d0894b52b3fb808b7ea20b` |
| Dense | `4fa1aff001b87e1a0175f8a852cb91fe44a1fd44b1478d53d5907986046771a4` |

## Why this profile

Floating screens compared the outer and inner changes separately. Positive
entries below are feasible diagnostic witnesses, not complete certificates.
Negative entries indicate failure of the tested bound, not low-distance
codewords. All rows use requested `K=2^20` and 10% retained-output distance.

| Outer | Inner width/state; outputs retained | Nominal rate | Sparse `A=128` | Middle `W/N=0.03` |
|---|---|---:|---:|---:|
| [128,80] | 4/4; 3 | 83.3333% | 1,666 bits | 16,402 bits |
| [64,40] | 8/8; 6 | 83.3333% | 870 bits | 17,086 bits |
| [64,40] | 4/4; 3 | 83.3333% | 122 bits | 12,730 bits |
| [64,39] | 4/4; 3 | 81.2500% | not needed | not replayed |
| [64,45] | 8/8; 7 | 80.3571% | -898 bits | -5,732 bits |

The selected [64,40] profile combines the smaller outer with smaller state
maps. Its tighter sparse margin remains above the 40-bit target.
The [64,39] fallback had over 200 bits at all sampled active counts from
32 through 512, but was not needed. The failed [64,45] control also shows
why rate alone does not predict this proof's sparse behavior.

## Replay and scope

From the repository root, run:

```text
python research/workstreams/char2_field_spin/theory/gf8_speed_certificate.py
python research/workstreams/char2_field_spin/theory/gf8_speed_certificate_tests.py
python research/workstreams/char2_field_spin/theory/gf8_speed_sparse_certificate_tests.py
python research/workstreams/char2_field_spin/theory/gf8_speed_geometry_tests.py
```

The canonical profile and exact maps are in
[`gf8_speed_profile.py`](theory/gf8_speed_profile.py) and
[`gf8_speed_maps.py`](theory/gf8_speed_maps.py).
The complete replay is
[`gf8_speed_certificate.py`](theory/gf8_speed_certificate.py).
Its helpers do not modify the earlier frozen certificates.

The complete aggregate replay passed. All twelve speed geometry, certificate,
and sparse tests also passed. A separate agent audited the component bounds,
length-64 conditioning factors, completion union, and coverage joins.

The scripts accept some other dimensions and sizes for experiments. This
note certifies only the stated K20 [64,40] profile.
A deterministic benchmark seed is a reproducibility fixture, not a proof
that its generated setup has the uniform distribution used by the theorem.
The [native follow-up](GF8_SPEED_NATIVE.md) validates the implementation and
selects cached regional stores with eight-packet batching at about4.1ms.
Its matched control is the earlier84% profile; no timing claim is part of
the mathematical certificate.
