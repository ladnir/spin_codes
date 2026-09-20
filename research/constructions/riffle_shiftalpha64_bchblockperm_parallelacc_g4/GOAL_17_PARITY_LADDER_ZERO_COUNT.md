# Goal 17: begin the parity ladder with zero parity

## Question

Goal 16 reduced the accumulator input to one packet histogram

\[
\mathbf h=(h_0,h_1,h_2,h_3,h_4).
\]

The remaining outer quantity is the expected number \(C(\mathbf h)\) of
nonzero messages that produce each histogram. This goal computes that
quantity for the zero-parity baseline. Later goals will add one parity block
and then two parity blocks.

This order separates two effects. The accumulator kernel remains fixed. Only
the outer constraints and the number of packet positions change.

## Zero-parity baseline

Let \(B=16384\). The baseline encodes each of the \(B\) field symbols with
the extended binary BCH code \([128,64,22]\). It independently permutes the
128 bits in each physical BCH block. It adds no field-valued parity block.

Each BCH block contains 32 four-bit packets. Hence the zero-parity baseline
has

\[
N_0=32B=524288
\]

packet positions.

Let \(A_r\) denote the BCH weight enumerator. For a local packet histogram
\(\mathbf a=(a_0,\ldots,a_4)\), define

\[
|\mathbf a|:=\sum_{k=0}^4 a_k,
\qquad
\|\mathbf a\|:=\sum_{k=0}^4 k a_k.
\]

When \(|\mathbf a|=32\) and \(\|\mathbf a\|=r\), define

\[
P_r(\mathbf a)
=
\frac{32!}{\prod_{k=0}^4a_k!}
\prod_{k=0}^4\binom4k^{a_k}.
\tag{1}
\]

Equation (1) counts the weight-\(r\) binary strings with histogram
\(\mathbf a\). Therefore

\[
\sum_{\substack{|\mathbf a|=32\\\|\mathbf a\|=r}}
P_r(\mathbf a)
=
\binom{128}{r}.
\tag{2}
\]

## One-block polynomial

Fix one physical BCH block. Sum over its \(2^{64}\) possible field inputs,
and average over its uniform bit permutation. The resulting histogram
polynomial is

\[
F(\mathbf t)
:=
\sum_r \frac{A_r}{\binom{128}{r}}
\sum_{\substack{|\mathbf a|=32\\\|\mathbf a\|=r}}
P_r(\mathbf a)\mathbf t^{\mathbf a},
\tag{3}
\]

where \(\mathbf t^{\mathbf a}:=\prod_{k=0}^4t_k^{a_k}\).
Equation (2) gives \(F(\mathbf 1)=2^{64}\).

The bit permutations of distinct physical blocks are independent. Thus the
complete expected histogram polynomial is

\[
F(\mathbf t)^B.
\tag{4}
\]

The all-zero message contributes \(t_0^{N_0}\). The exact zero-parity count
for nonzero messages is therefore

\[
C_0(\mathbf h)
=
[\mathbf t^{\mathbf h}]
\left(F(\mathbf t)^B-t_0^{N_0}\right).
\tag{5}
\]

The coefficient in (5) is an ensemble expectation. A fixed setup can have a
different integer count. Linearity of expectation makes (5) the quantity
needed by the first-moment distance proof.

The global packet permutation does not occur in (3)--(5). It preserves the
global packet histogram.

## Minimum shell

The zero-parity outer code has field-symbol distance one. Its minimum binary
shell contains one active BCH block of weight 22. The shell has

\[
B A_{22}
=16384\mathbin{\cdot}243840
=3{,}995{,}074{,}560
\tag{6}
\]

labeled outer words.

Only 136 local packet histograms satisfy

\[
|\mathbf a|=32,
\qquad
\|\mathbf a\|=22.
\]

For each such histogram, the expected global count is

\[
C_{0,22}(N_0-32+a_0,a_1,a_2,a_3,a_4)
=
B A_{22}\frac{P_{22}(\mathbf a)}{\binom{128}{22}}.
\tag{7}
\]

The active packet support ranges from 6 through 22. The most likely local
histogram is

\[
(a_0,a_1,a_2,a_3,a_4)=(14,14,4,0,0).
\]

It has 18 active packets and probability approximately \(0.16768187\).
The receipt records every coefficient in (7) as an exact rational number.

## Accumulator interface

Let \(A_{N,4}(\mathbf h,w)\) and \(T_{N,4}(\mathbf h)\) be the quantities
from Goal 16. The zero-parity expected output enumerator is

\[
\overline A^{(0)}_w
=
\sum_{\mathbf h} C_0(\mathbf h)
\frac{A_{N_0,4}(\mathbf h,w)}{T_{N_0,4}(\mathbf h)}.
\tag{8}
\]

Equation (8) is exact. This goal computes the first factor. It does not yet
expand (8) for every zero-parity message. The minimum shell can nevertheless
be composed exactly, as described next.

## Exact minimum-shell composition

Fix an active packet sequence. Let \(m_b\) count its prefix accumulator
states of weight \(b\). The inactive packets occupy the gaps before, between,
and after the active packets. A state of weight zero makes its following gap
free. A state of weight \(b>0\) charges \(b\) output bits per packet in its
following gap.

The exact evaluator groups active paths by

\[
(H,m_0,m_1,m_2,m_3,m_4).
\]

It then counts every compatible gap placement. The gap generating function
has poles only at roots of unity of orders 1 through 4. Its coefficients are
therefore quasipolynomials of period 12. Finite differences evaluate those
quasipolynomials at the global length without iterating over all 524288
packet positions.

For the complete weight-22 occupation-one shell, the exact results are:

| Bad output weight \(D\) | Relative weight | Per-word log probability | Expected-count log | Margin |
|---:|---:|---:|---:|---:|
| 76000 | 0.0362396 | -32.0397 | -0.1442 | 0.1442 bits |
| 77000 | 0.0367165 | -31.8612 | 0.0344 | -0.0344 bits |
| 188743 | 0.0900000 | -19.4650 | 12.4306 | -12.4306 bits |

Thus the minimum shell crosses expected count one between 76000 and 77000.
At the 9% target, this shell alone has about \(2^{12.43}\) expected bad
words. Consequently, the zero-parity first-moment proof cannot reach 9%.

The positive margin at 76000 concerns only the minimum shell. Higher BCH
weights and messages with several active blocks still contribute to (8).
Therefore this goal does not prove distance 76001 for the complete
zero-parity code.

## What changes when parity is added

The next stage adds the parity symbol \(p_0=\sum_i x_i\). Its
minimum-occupation words use two field coordinates. Both coordinates contain
the same nonzero field value. Consequently, the one-parity minimum shell
will use products of two copies of the conditional polynomial in (3).

Double parity raises the minimum occupation to three. The three field values
then depend on the selected outer support. That dependence is the first
stage that requires the shifted coefficient schedule.

Thus the ladder isolates the source of the previous complexity:

\[
\text{one local factor}
\longrightarrow
\text{one repeated-value pair}
\longrightarrow
\text{one support-dependent triple}.
\]

## Verification

Run

```powershell
python scripts/analyze_riffle_parity_ladder_zero.py
```

The script checks the complete BCH spectrum, complement symmetry, every
weight-22 packetization coefficient, and all probability masses. It writes
`receipts/goal17_zero_parity_histogram_count.json`.

Run the exact accumulator composition with

```powershell
python scripts/analyze_riffle_parity_ladder_zero_accumulator.py
```

This second audit performs 15,257 direct gap-coefficient comparisons and
2,008 exhaustive all-gap comparisons. It writes
`receipts/goal17_zero_parity_accumulator_exact.json`.

The next goal should compute and compose the one-parity minimum shell with
the same exact interface. The complete-code sums remain separate proof
obligations at every parity level.
