# Dense classes for the GF256 outer-rate sweep

The dense bound permits outer dimensions 88, 89, and 90 at the existing
24% inverse-input boundary. Dimension 92 passes from 25%, and dimension 96
passes from 28%. These are component certificates for 10% output distance;
the sparse and middle classes must be joined separately.

The [new replay](theory/gf8_outer_rate_dense.py) keeps length 128, width 16,
state size eight, and GF256. The outer is GRS[128,d]. For requested
`K=2^k_log`, put

\[
L=64\left\lceil\frac{K}{64d}\right\rceil,\qquad
N=128L,\qquad e=\lceil N/10\rceil-1.
\]

The route samples an independent permutation in every coordinate region.
Input diagonals are independent nonzero GF256 scalars. The inner is
independent of this randomness and invertible. Its intended realization
is the unchanged symmetric Cauchy recurrence. The dense argument requires
only independence and invertibility.

At K20, the complete directed replays give:

| Dimension d | Dense density cutoff | Integer weights covered | Cells | Margin, bits |
|---:|---:|---:|---:|---:|
| 88 | .24 | 367657 through 1531904 | 152 | 63954 |
| 89 | .24 | 363725 through 1515520 | 152 | 35507 |
| 90 | .24 | 359793 through 1499136 | 152 | 7389 |
| 92 | .25 | 366592 through 1466368 | 150 | 3714 |
| 96 | .28 | 392233 through 1400832 | 144 | 17894 |

For dimension 88, K16 covers weights 23593 through 98304 with 2988 bits.
K18 covers weights 92406 through 385024 with 15031 bits. Both use 152 cells.
The finite requested rates differ from the design rate because of padding.

## Endpoint bound and acceptance

Condition on the inner map B and fix a forbidden output y. Write
`x=B^-1 y` and `alpha=wt(x)/N`. The existing
[regional membership theorem](FIELD128_DENSE_CERTIFICATE.md) bounds the
probability that x belongs to the routed, diagonally scaled outer code.
For a cell `[a,b]` and a rational tilt `lambda>=1`, define

\[
Q_d(a,b,\lambda)
=(1-a)^{128}
+256^{-(128-d)}
\left(1+\left(\frac{256}{255}\lambda-1\right)b\right)^{128}.
\]

The first term decreases with density; the second increases. Thus
`Q_d(a,b,lambda)` bounds the membership bracket throughout the cell.
The extraction factor is at most `lambda^(-Na)`.

The output ball contains fewer than `2^(127N/100)` vectors. The replay
checks the existing exact rational comparison

\[
(10/9)^{100}2295^{10}<2^{127}.
\]

For an integer cell `[u,v]`, set `a=u/N` and `b=v/N`. A lower bound on
the resulting failure margin is

\[
u\log_2\lambda-L\log_2 Q_d(a,b,\lambda)
-128\log_2(L+1)-\frac{127N}{100}.
\]

The replay forms Q with exact fractions and evaluates logarithms with
192-bit outward Arb intervals. It rounds the lower margin downward.
Floating optimization only proposes the rational tilt; acceptance uses
the directed bound. A failing cell is bisected until it passes or a
single failing weight is reached.

The cells partition every integer weight from the cutoff through N.
Their minimum margin gives one uniform membership bound for each fixed
forbidden output. The output union therefore adds no factor for the
number of density cells.

An approximate point search places the K20 40-bit cutoffs near .227641,
.232140, .236914, .247419, and .273590 for dimensions 88, 89, 90, 92,
and 96. These estimates concern this bound only. They are neither complete
coverage certificates nor lower bounds on the code's failure probability.

## Aggregate interface and replay

Each receipt identifies the outer, field, packet width, state size,
construction, requested and padded source sizes, regional length,
output length, and forbidden threshold. It returns `weights=[first,N]`,
`complete=true`, `rigorous=true`, and an integer
`failure_upper_dyadic_exponent=m`, meaning failure probability at most
`2^-m` for this dense class.

An aggregate must match every geometry field. It must check
`middle_last+1=dense_first` and
`(129-d)*(sparse_last+1)>=middle_first`, then sum the three actual
probability bounds. Component bit margins must not be added.

```text
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_dense_tests.py
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_dense.py --dimension 88 89 90
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_dense.py --dimension 88 --k-log 16 18
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_dense.py --dimension 92 --first 1/4
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_dense.py --dimension 96 --first 7/25
```

Five tests pass. They check geometry, frozen-bracket agreement, endpoint
monotonicity, an independent exact integer-power inequality, full dense
coverage, and rejection of a failing or invalid class. Existing certificate
files remain unchanged. No native work or benchmark accompanies this replay.

The [dimension-88 K18 aggregate](theory/gf8_outer_rate_certificate.py) was
independently replayed after these component checks. Its sparse, middle,
and dense margins are 72, 129, and 15031 bits; their exact probability sum
retains 72 bits. The sparse class ends at 375 active rows. Every remaining
source has weight at least `41*376=15416`, above the middle lower endpoint
15401. The middle ends at 92405 and the dense class begins at 92406.
This closes the ideal K18 ensemble with requested rate `32/47` and minimum
distance 38503. The next rate decision is the unresolved K20 sparse trough;
the dense component is already sufficient there.
