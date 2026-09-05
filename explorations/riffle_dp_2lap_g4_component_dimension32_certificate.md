# Riffle DP-2Lap g=4: endpoint-aware certificate through dimension 32

## Result

The dimension-31 certificate covers every nonempty irreducible-component
support of total dimension at most 31. Exactly four new supports occur at
dimension 32:

| mask | component degrees |
|---:|:---|
| `0x2d` | 1, 4, 9, 18 |
| `0x34` | 4, 10, 18 |
| `0x4b` | 1, 2, 9, 20 |
| `0x52` | 2, 10, 20 |

For each listed support (S), exact exhaustive certificates prove

\[
\begin{aligned}
d_{24}(S,v)&\ge97 &&\text{for }v\in\{1,\ldots,14\},\\
d_{24}(S,15)&\ge72,\\
d_{12}(S,15)&\ge36.
\end{aligned}
\]

The quantity (d_w(S,v)) is the minimum two-sided weight in the
(16w)-bit observation code. It minimizes weight and complement weight over
all nonzero characters supported on (S).

The result extends the complete endpoint-aware tranche to all 65 supports of
total dimension at most 32. The separate (C_{18}\oplus C_{20}) certificate
adds mask `0x60`. Exact merged coverage is therefore 66 of 127 nonempty
supports, leaving 61 supports open.

## Exhaustive argument

Each 24-node code has length 384 and dimension 32. Its 12 disjoint
information sets partition all coordinates. A word with two-sided weight
below 97 restricts to weight at most eight on one information set. A
value-15 word with two-sided weight below 72 restricts to weight at most five.

Each 12-node endpoint code has length 192. Its six disjoint information sets
also partition all coordinates. A value-15 word with two-sided weight below
36 restricts to weight at most five on one information set.

Each restriction determines one codeword. Exhausting the specified Hamming
balls around the zero and all-one restrictions excludes every violating
character. Each implementation checks

\[
10{,}119{,}775{,}656
\]

information vectors across 64 cases.

The primary implementation uses fixed-dimension split-half tables. The
independent verifier uses different information-set partitions and direct
recursive fixed-weight enumeration. Both implementations pass every case.

A Python audit reconstructs the BCH state map, transpose recurrence,
component kernels, and observation columns. The audit verifies 128 full
receipts and 16 preflight receipts. It also checks each information-set rank,
partition, receipt hash, executable hash, and exhaustive vector count.

## Global implication and scope

The local bounds imply the required full-orbit character caps. Values 1
through 14 receive weight at least

\[
1{,}365\cdot97=132{,}405>131{,}088.
\]

For value 15, the complete blocks and endpoint give

\[
1{,}365\cdot72+36=98{,}316.
\]

The same bounds apply to complement weight. Consequently, every character in
the 65-support tranche satisfies the required bias cap.

This certificate does not cover the other 61 supports. It also closes only
the full-character terminal-zero step in the larger proof ledger.

Evidence:

- `explorations/riffle_dp_2lap_g4_component_dimension32_primary_probe.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension32_primary_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension32_independent_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension32_audit.json`;
- `scripts/run_riffle_dp_2lap_g4_component_dimension32.py`;
- `scripts/audit_riffle_dp_2lap_g4_component_dimension32.py`.

## Next proof decision

The next finite frontier is dimension 33. It contains masks `0x2e`, `0x35`,
`0x4c`, and `0x53`. A complete sweep would check 12,216,115,184 information
vectors per implementation.

Dimension 32 is also the planned checkpoint for structural analysis. The
non-15 certificates have substantial observed margin: the smallest side
among the enumerated representatives is 120, compared with the required 97.
The value-15 equality cases arise from the known degree-one component. Before
continuing through dimension 33, the next goal should test whether these two
patterns yield a uniform anti-cancellation lemma for the remaining supports.
