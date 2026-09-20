# Riffle DP-2Lap g=4: 24-node certificate for \(C_{18}\oplus C_{20}\)

## Claim

Let \(T\) be the 64-bit zero-input transpose recurrence. Fix a nonzero packet
value \(v\in\mathbb F_2^4\). For a character
\(\chi\in C_{18}\oplus C_{20}\), define the 384-bit observation

\[
c^{(24)}_{\chi,v}(s,t)
=\langle\chi,T^{t+1}(v\mathbin{\ll}4s)\rangle,
\qquad 0\le s<16,\quad 0\le t<24.
\]

The exhaustive certificate proves the following statement for every nonzero
\(\chi\in C_{18}\oplus C_{20}\):

\[
\min\{\operatorname{wt}(c^{(24)}_{\chi,15}),
384-\operatorname{wt}(c^{(24)}_{\chi,15})\}\ge73,
\]

and, for every \(v\ne15\),

\[
\min\{\operatorname{wt}(c^{(24)}_{\chi,v}),
384-\operatorname{wt}(c^{(24)}_{\chi,v})\}\ge97.
\]

## Finite certificate

Each 24-node code has length 384 and dimension 38. For every packet value,
the certificate constructs ten disjoint 38-coordinate information sets.
They cover 380 coordinates.

Suppose a word or its complement has weight below the required distance
\(d_v\). Its restrictions to the ten information sets have total weight at
most \(d_v-1\). At least one restriction therefore has weight at most

\[
r_v=\left\lfloor\frac{d_v-1}{10}\right\rfloor.
\]

Thus, \(r_{15}=7\) and \(r_v=9\) for \(v\ne15\). Each restriction determines
one codeword. It suffices to enumerate the Hamming balls of radius \(r_v\)
around both zero and one for every information set.

The primary implementation uses fixed-width 384-bit words and split-19
lookup tables. It checks 32,062,999,280 information vectors across all packet
values. The independent implementation uses distinct information partitions
and direct recursive combination enumeration. It checks the same number of
vectors. Both implementations pass every packet value.

A Python audit reconstructs the intended recurrence and component kernel.
It verifies the rank and disjointness of all information sets in both
certificate families. The audit also checks every receipt, candidate count,
threshold, and artifact hash.

The smallest side encountered by the primary search is 120 for
\(v\ne15\) and 119 for \(v=15\). These values are observed upper bounds on
the exact distances. The certificate proves only the stated lower bounds 97
and 73; it does not determine either exact minimum.

## Full-orbit consequence within the subcode

The full orbit contains 32,772 nodes. Partition its first 32,760 nodes into
1,365 consecutive 24-node blocks. The last 12 nodes remain unpaired.

The transpose recurrence is invertible and preserves
\(C_{18}\oplus C_{20}\). Hence every complete block starts from a nonzero
character in the same subcode. For concatenated binary blocks, the two-sided
weight of the concatenation is at least the sum of the blockwise two-sided
weights.

For \(v\ne15\), the complete blocks contribute at least

\[
1{,}365\cdot97=132{,}405>131{,}088=524{,}352/4.
\]

The unmatched block can contribute zero. The resulting full-orbit character
bias is at most

\[
\frac{524{,}352-2\cdot132{,}405}{524{,}352}
=\frac{259{,}542}{524{,}352}
\approx0.494977<\frac12.
\]

For \(v=15\), the complete blocks contribute at least

\[
1{,}365\cdot73=99{,}645>98{,}316=3\cdot524{,}352/16.
\]

The corresponding bias is at most

\[
\frac{325{,}062}{524{,}352}
\approx0.619931<\frac58.
\]

Therefore, the desired global character caps hold for every nonzero
character in \(C_{18}\oplus C_{20}\).

## Scope

The certificate resolves the degree-18 plus degree-20 mixed-component
obstruction. It does not prove the 24-node bounds for the other uncovered
component supports or for the full 64-dimensional character code. Therefore,
it does not yet close the support-33 terminal-zero row or the complete
Riffle DP-2Lap g=4 proof.

## Artifacts

- Primary source: `scripts/certify_riffle_dp_2lap_g4_c18_c20_24node.cpp`
- Independent source: `scripts/verify_riffle_dp_2lap_g4_c18_c20_24node.cpp`
- Aggregate audit source: `scripts/audit_riffle_dp_2lap_g4_c18_c20_24node.py`
- Aggregate audit receipt:
  `explorations/riffle_dp_2lap_g4_c18_c20_24node_audit.json`
- Primary receipts:
  `explorations/riffle_dp_2lap_g4_c18_c20_24node_v01.json` through
  `explorations/riffle_dp_2lap_g4_c18_c20_24node_v15.json`
- Independent receipts:
  `explorations/riffle_dp_2lap_g4_c18_c20_24node_v01_independent.json`
  through
  `explorations/riffle_dp_2lap_g4_c18_c20_24node_v15_independent.json`
