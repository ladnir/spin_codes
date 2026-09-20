# Riffle DP-2Lap g=4: endpoint-aware certificate through dimension 31

## Result

The previous certificate covers every nonempty irreducible-component support
of total dimension at most 30. Exactly four new supports occur at dimension
31:

| mask | component degrees |
|---:|:---|
| `0x2c` | 4, 9, 18 |
| `0x33` | 1, 2, 10, 18 |
| `0x4a` | 2, 9, 20 |
| `0x51` | 1, 10, 20 |

Exact exhaustive certificates prove the corrected local bounds on all four
supports. For each support (S),

\[
\begin{aligned}
d_{24}(S,v)&\ge97 &&\text{for }v\in\{1,\ldots,14\},\\
d_{24}(S,15)&\ge72,\\
d_{12}(S,15)&\ge36.
\end{aligned}
\]

The symbol (d_w(S,v)) denotes the minimum two-sided weight in the
(16w)-bit observation code. Thus, it minimizes the weight and complement
weight over all nonzero characters supported on (S).

The result extends the complete endpoint-aware tranche from 57 supports
through dimension 30 to 61 supports through dimension 31. The separate
(C_{18}\oplus C_{20}) certificate adds mask `0x60`. Exact merged coverage is
therefore 62 of the 127 nonempty supports, leaving 65 supports open.

## Exhaustive argument

Each dimension-31 code has length 384 at 24 nodes. It contains 12 disjoint
information sets. A word with two-sided weight below 97 restricts to weight
at most eight on one set. A value-15 word with two-sided weight below 72
restricts to weight at most five.

The 12-node endpoint code has length 192 and contains six disjoint
information sets. A value-15 word with two-sided weight below 36 restricts to
weight at most five on one set.

Each restriction determines one codeword. Exhausting the corresponding balls
around the zero and all-one restrictions excludes every violating character.
The complete workload per implementation is

\[
7{,}716{,}616{,}224
\]

information vectors across 64 cases.

The primary implementation uses compile-time dimension specialization and
split-half lookup tables. The independent verifier changes the information
sets and uses direct recursive fixed-weight enumeration. Both implementations
pass all 64 cases.

A Python audit independently reconstructs the BCH state map, transpose
recurrence, component kernels, and observation columns. It verifies all 128
full receipts and 16 preflight receipts. It also checks each information-set
rank, receipt hash, exhaustive candidate count, and primary-independent
partition difference.

## Global implication and scope

The dimension-31 result uses the same global implication as the dimension-30
certificate. The 1,365 complete 24-node blocks give weight at least 132,405
for values 1 through 14. For value 15, the complete blocks and endpoint give

\[
1{,}365\cdot72+36=98{,}316.
\]

The inequalities also apply to complement weight. Consequently, every
character in the 61-support tranche satisfies the required full-orbit bias
cap.

This certificate does not cover the other 65 supports. It also closes only
the full-character terminal-zero step in the larger proof ledger.

Evidence:

- `explorations/riffle_dp_2lap_g4_component_dimension31_primary_probe.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension31_primary_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension31_independent_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension31_audit.json`;
- `scripts/run_riffle_dp_2lap_g4_component_dimension31.py`;
- `scripts/audit_riffle_dp_2lap_g4_component_dimension31.py`.

The next finite frontier is dimension 32. It consists of masks `0x2d`,
`0x34`, `0x4b`, and `0x52`. The same certificate method predicts
10,119,775,656 information vectors per implementation for that tranche.
