# Riffle DP-2Lap g=4: endpoint-aware certificate through dimension 33

## Result

The dimension-32 certificate covers every nonempty irreducible-component
support of total dimension at most 32. Exactly four new supports occur at
dimension 33:

| mask | component degrees |
|---:|:---|
| `0x2e` | 2, 4, 9, 18 |
| `0x35` | 1, 4, 10, 18 |
| `0x4c` | 4, 9, 20 |
| `0x53` | 1, 2, 10, 20 |

For each listed support \(S\), exact exhaustive certificates prove

\[
\begin{aligned}
d_{24}(S,v)&\ge97 &&\text{for }v\in\{1,\ldots,14\},\\
d_{24}(S,15)&\ge72,\\
d_{12}(S,15)&\ge36.
\end{aligned}
\]

The quantity \(d_w(S,v)\) is the minimum two-sided weight in the
\(16w\)-bit observation code. It minimizes weight and complement weight over
all nonzero characters supported on \(S\).

The result extends the endpoint-aware tranche to all 69 supports of total
dimension at most 33. The separate \(C_{18}\oplus C_{20}\) certificate adds
mask `0x60`. Exact merged coverage is therefore 70 of 127 nonempty supports,
leaving 57 supports open.

## Exhaustive argument

Each 24-node code has length 384 and dimension 33. The certificate constructs
11 disjoint information sets. A word with two-sided weight below 97 restricts
to weight at most eight on one information set. A value-15 word with
two-sided weight below 72 restricts to weight at most six.

Each 12-node endpoint code has length 192 and dimension 33. The certificate
constructs five disjoint information sets. A value-15 word with two-sided
weight below 36 restricts to weight at most seven on one information set.

Each restriction determines one codeword. Exhausting the specified Hamming
balls around the zero and all-one restrictions excludes every violating
character. Each full implementation checks

\[
12{,}216{,}115{,}184
\]

information vectors across 64 cases.

The primary implementation uses fixed-dimension split-half tables. The
independent verifier uses different information-set partitions and recursive
fixed-weight enumeration. Both implementations pass every case.

A Python audit reconstructs the BCH state map, transpose recurrence,
component kernels, and observation columns. The audit verifies 128 full
receipts and 16 preflight receipts. It also checks every information-set
rank, partition, receipt hash, executable hash, and exhaustive vector count.

## Global implication and scope

For values 1 through 14, the complete 24-node blocks contribute weight at
least

\[
1{,}365\cdot97=132{,}405>131{,}088.
\]

For value 15, the complete blocks and endpoint contribute weight at least

\[
1{,}365\cdot72+36=98{,}316.
\]

The same bounds apply to complement weight. Consequently, every character in
the 69-support tranche satisfies the required full-orbit bias cap.

This certificate does not cover the other 57 supports. It also closes only
the full-character terminal-zero step in the larger proof ledger.

Evidence:

- `explorations/riffle_dp_2lap_g4_component_dimension33_primary_probe.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension33_primary_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension33_independent_run.json`;
- `explorations/riffle_dp_2lap_g4_component_dimension33_audit.json`;
- `scripts/run_riffle_dp_2lap_g4_component_dimension33.py`;
- `scripts/audit_riffle_dp_2lap_g4_component_dimension33.py`.

## Next proof decision

The dimension-34 frontier contains four supports:

- `0x2f=(1,2,4,9,18)`;
- `0x36=(2,4,10,18)`;
- `0x4d=(1,4,9,20)`;
- `0x54=(4,10,20)`.

Its predicted workload is 15,745,416,320 information vectors per
implementation. Dimension 35 introduces an information-set-radius increase
and is substantially more expensive. Dimension 34 is therefore the next
reasonable exact frontier and the natural checkpoint before that cost jump.
