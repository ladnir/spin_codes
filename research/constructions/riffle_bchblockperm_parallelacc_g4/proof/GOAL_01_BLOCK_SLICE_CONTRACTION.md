# Goal 01 report: one-block slice contraction

## Result

The within-block permutation removes most of the worst-case grouping loss for
one minimum-weight BCH word. At global binary length \(2{,}097{,}408\) and
distance threshold 40, the bounds are:

| Analysis | Log probability bound |
|:---|---:|
| Old one-lane packet bound | \(-28.1478\) |
| Block-averaged joint-lane moment | \(-95.0755\) |
| Exact-zero-gap bound | \(-93.5939\) |
| Full global bit permutation, exact | \(-155.2171\) |

The joint-lane moment recovers 66.93 of the 127.07 bits between the old bound
and the full-bit benchmark. This is a one-block statement. It does not prove a
distance bound for the complete code.

## Exact packet-support law

Fix a nonzero BCH word \(c\in\mathbb F_2^{128}\) of weight \(h\). If
\(\sigma\gets S_{128}\), then \(\sigma(c)\) is uniform over the binary
weight-\(h\) slice. Partition \(\sigma(c)\) into 32 ordered four-bit packets.

Let \(H\) denote the number of nonzero packets. Inclusion-exclusion gives

\[
\Pr[H=s]
=
\frac{\binom{32}{s}}{\binom{128}{h}}
\sum_{j=0}^{s}
(-1)^{s-j}\binom{s}{j}\binom{4j}{h}.
\tag{1}
\]

For \(h=22\), equation (1) gives

\[
\mathbb E[H]=17.0987926509.
\]

The mode is \(H=17\), with probability 0.27081. The minimum possible support
is six packets, but its probability is only

\[
8.3565445613\times10^{-17}<2^{-53.41}.
\]

The old fixed-word proof effectively assigned probability one to such maximal
packing.

## Block-averaged gap moment

Fix \(z,\tau\in(0,1)\). For the ordered nonzero packet values
\(v_1,\ldots,v_H\), define prefix states and weights

\[
q_j:=\sum_{i=1}^jv_i,
\qquad
a_j:=\operatorname{wt}(q_j).
\]

The global packet permutation maps their support to a uniform \(H\)-subset of
the \(N=524{,}352\) packet positions. The Goal 02 gap calculation gives

\[
\Pr[W\le D\mid v_1,\ldots,v_H]
\le
z^{-D}
\frac{\tau^{-(N-H)}}{\binom NH}
\prod_{j=1}^{H}\frac{z^{a_j}}{1-\tau z^{a_j}}.
\tag{2}
\]

The weight-\(h\) slice is invariant under permutations of the 32 packet slots.
Therefore, averaging the right side of (2) over the local block permutation
also averages the active packet order induced by the global permutation.

The replay script evaluates that average by dynamic programming. Its state is

\[
(t,r,H,q),
\]

where \(t\) is the processed packet count, \(r\) is the used binary weight,
and \(q\in\mathbb F_2^4\) is the current prefix state. The computation retains
all 16 prefix states. It does not select one lane or replace \(h\) by \(h/4\).

At \(h=22\) and \(D=40\), optimization of \(z\) and \(\tau\) gives

\[
\Pr[W\le40]\le2^{-95.0754829299}.
\tag{3}
\]

The effective twenty-second root in (3) is 0.05001093. The old one-lane root
is 0.41195421.

## Intended distance scale

At the previous target \(D=188{,}766\), approximately \(0.09L\), the averaged
moment in (2) becomes trivial. Its optimized upper bound is one. For
comparison, a full global bit permutation has exact one-block probability

\[
2^{-20.1472349681}
\]

at the same \(h=22\) and \(D\).

The trivial moment does not give a counterexample. Exact slice enumeration
shows that zero returns are common:

\[
\Pr[\text{at least one zero prefix state}]=0.6814216425,
\]

\[
\mathbb E[\text{zero prefix states}]=1.1301729266,
\qquad
\Pr[q_H=0]=0.1249999998.
\]

The geometric factor in (2) assigns a singular cost to each zero-state gap as
\(\tau\) approaches one. The next proof should sum all zero-cost gaps exactly
and apply a moment only to positive-state gaps. A split that treats zero
returns as rare cannot close the target.

This exact zero-gap summation is possible without the full state-weight
histogram. Fix a path with \(H\) active packets and \(m_0\) zero prefix states.
There are

\[
K:=H-m_0
\]

positive-state gaps. If \(W\le D\), the total extra length assigned to those
gaps is at most \(D-K\). The number of positive-gap allocations is therefore
at most \(\binom DK\). The leading gap and the \(m_0\) zero-state gaps absorb
the remaining positions in at most

\[
\binom{N-H+m_0}{m_0}
\]

ways. Conditional on \((H,m_0)\), this gives

\[
\Pr[W\le D\mid H,m_0]
\le
\min\left\{
1,
\frac{\binom D{H-m_0}\binom{N-H+m_0}{m_0}}
{\binom NH}
\right\}.
\tag{4}
\]

The exact dynamic program also computes the joint slice distribution of
\((H,m_0)\). Averaging (4) at \(D=188{,}766\) gives

\[
\Pr[W\le188{,}766]\le2^{-11.7208369330}.
\tag{5}
\]

Thus exact treatment of zero gaps restores a nontrivial target-scale bound.
The dominant contributions have 15--18 active packets and 4--6 zero prefix
states. Equation (4) ignores that positive states often have weight greater
than one, so (5) retains identifiable slack.

## Scope

Goal 01 proves the slice law in (1) and a valid one-block moment bound. It does
not compose the field double-parity outer, interleavings of several active BCH
blocks, or the complete BCH spectrum. Goal 02 begins with the three-block
minimum forced by the field outer. Goal 01 makes no performance claim for the
128-bit permutations.
