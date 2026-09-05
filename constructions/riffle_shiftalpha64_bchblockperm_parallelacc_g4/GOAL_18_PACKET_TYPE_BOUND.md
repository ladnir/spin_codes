# Goal 18: replace the packet-type table by coefficient tilting

## Question

Goal 16 indexes the accumulator law by four packet counts. An eager table at
the target length has

\[
\binom{N+4}{4}
\]

entries. For \(N=524288\), this table is too large. This goal asks whether a
bound can retain the packet information without expanding that table.

## Coefficient bound

For packet weight \(k\), define the matrix

\[
K_k(z)_{a,b}:=c_4(a,b,k)z^b.
\]

For positive \(\mathbf x=(x_0,\ldots,x_4)\), define

\[
M(\mathbf x,z):=\sum_{k=0}^4x_kK_k(z).
\]

Let \(A_{N,4}(\mathbf h,w)\) and \(T_{N,4}(\mathbf h)\) be the quantities
from Goal 16. Fix \(0<z<1\). Nonnegativity of every coefficient gives

\[
\sum_{w=0}^D A_{N,4}(\mathbf h,w)
\le
z^{-D}\mathbf x^{-\mathbf h}
e_0^{\mathsf T}M(\mathbf x,z)^N\mathbf 1.
\tag{1}
\]

Consequently,

\[
\Pr[\operatorname{wt}(y)\le D\mid\mathbf h]
\le
\frac{z^{-D}\mathbf x^{-\mathbf h}
e_0^{\mathsf T}M(\mathbf x,z)^N\mathbf 1}
{T_{N,4}(\mathbf h)}.
\tag{2}
\]

Equation (2) replaces multivariate coefficient extraction by a five-state
matrix power and five positive scalar parameters. One scale in \(\mathbf x\)
is redundant, so the numerical optimization has five effective variables:
four packet tilts and the output tilt \(z\).

The bound retains the random order created by the global packet permutation.
This point is essential. A product bound applied separately to each packet
weight has \(\tau_0=1\), because a zero packet can occur while the accumulator
state is zero. Such a bound permits all zero packets to be harmless and loses
the random-placement effect.

## Minimum-shell probe

The numerical probe applies (2) to the complete zero-parity, weight-22,
occupation-one shell. It averages all 136 local packet histograms before
optimizing the common tilts. The table compares the resulting upper bound
with the exact calculation from Goal 17.

| \(D\) | Exact expected-count log | Bound log | Loss |
|---:|---:|---:|---:|
| 76000 | -0.1442 | 12.6403 | 12.7844 bits |
| 77000 | 0.0344 | 12.8144 | 12.7800 bits |
| 188743 | 12.4306 | 24.8610 | 12.4304 bits |

The optimizer used both deterministic starts and a sequential global search.
Both searches returned the same values. These values are numerical evidence,
not proof certificates.

The bound removes the \(N^4\) expansion, but it cannot replace the exact
sparse-shell calculation near its threshold. At \(D=76000\), the exact shell
has only 0.144 bits of margin, whereas coefficient tilting loses 12.784 bits.

## Summing the outer types

Suppose \(C(\mathbf h)\) is the expected number of outer words with packet
histogram \(\mathbf h\). Equation (2) bounds each summand in

\[
\sum_{\mathbf h} C(\mathbf h)
\Pr[\operatorname{wt}(y)\le D\mid\mathbf h].
\tag{3}
\]

A crude replacement of (3) by the largest summand times the number of types
would cost

\[
\log_2\binom{524292}{4}=71.4151
\]

bits. This finite-length loss is not acceptable for the current margins.
The fact that this cost is only \(O(\log N)\) asymptotically does not make it
small at the target length.

The viable use of (2) is therefore hybrid.

1. Treat low outer occupation exactly. Goal 17 already supplies this method
   for occupation one.
2. Use coefficient tilting for a bulk region where individual terms have a
   substantial negative exponent.
3. Sum the bulk with adaptive boxes or a certified log-sum-exp calculation.
   Do not charge every formal packet type at the worst-case value.
4. Use a smaller exact construction to locate the transition between the
   sparse and bulk regions.

This division removes the need for a complete \(N^4\) table. It does not
reduce the exact sufficient statistic below four packet counts.

## Reproduction

Run

```powershell
python scripts/explore_riffle_parallelacc_type_bound.py --global-search
```

The script writes
`receipts/goal18_packet_type_bound_probe.json`. It checks the local type mass
and compares scaled matrix exponentiation against direct iteration through
length 16.

The next experiment should compute a complete small-parameter spectrum. The
existing extended BCH \([8,4,4]\) instance over \(GF(16)\) supports exact
zero-, one-, and two-parity calculations. That spectrum can measure the
sparse-to-bulk transition and the loss from (2) across all occupations.
