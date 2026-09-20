# Riffle DP-2Lap g=4: dimension-30 tranche counterexample

## Target

The ambitious tranche contained every nonempty irreducible-component support
of total dimension at most 30. The component degrees are

\[
1,2,4,9,10,18,20.
\]

There are 57 target supports. Eight maximal supports cover them under
subcode inclusion:

\[
\begin{gathered}
(10,20),\ (1,9,20),\ (1,10,18),\ (2,10,18),\\
(1,2,4,18),\ (1,2,4,20),\ (1,2,9,18),
(1,2,4,9,10).
\end{gathered}
\]

For each maximal subcode, the proposed certificate required 24-node
two-sided distance 73 for packet value 15 and distance 97 for every other
nonzero value.

## Exact counterexample

The value-15 claim is false. In the subcode with component degrees
\((1,9,20)\), the character

\[
\chi=\mathtt{0xa685aac60e5acfb3}
\]

has 24-node weight 72. Hence

\[
\min\{72,384-72\}=72<73.
\]

The character is the known pure degree-one character. Applying the degree-one
factor to \(\chi\) gives zero. Applying either the degree-nine or degree-20
factor gives a nonzero result. Thus, the failure does not arise from new
mixed-component cancellation.

The same witness refutes the distance-73 claim for every target subcode that
contains the degree-one component. This accounts for 28 of the 57 target
supports.

## Full-orbit interpretation

An independent implementation reconstructs the recurrence and replays all
32,772 nodes. Every 12-node block has weight 36. Every complete 24-node block
has weight 72. The final 12-node block also has weight 36.

The complete pairs and the remainder therefore contribute

\[
1{,}365\cdot72+36=98{,}316.
\]

The full character sum is

\[
524{,}352-2\cdot98{,}316=327{,}720,
\]

so the absolute bias is exactly

\[
\frac{327{,}720}{524{,}352}=\frac58.
\]

The witness refutes the stronger distance-73 device. It does not refute the
desired value-15 global cap. Instead, it attains that cap and shows why the
last 12-node block cannot be discarded for value 15.

## Audited partial results

The primary run checked 2,917,959,325 information vectors before it found the
counterexample. The audit independently reconstructs the intended codes and
validates all 30 generated receipts and their information sets.

The maximal subcode \((10,20)\) passed every packet value. It certifies the
strong thresholds for the three nonempty supports \((10)\), \((20)\), and
\((10,20)\).

For \((1,9,20)\), packet values 1 through 14 passed distance 97 before value
15 refuted distance 73. Subcode containment gives non-15 distance-97
certificates for nine target supports after combining this result with the
\((10,20)\) certificate.

The exhaustive verifier planned for a universal passing theorem was not run.
Once the primary implementation produced an exact witness, direct independent
replay was both sufficient and stronger evidence for the refutation.

## Revised proof target

The value-15 accounting must retain the endpoint:

- prove 24-node two-sided distance at least 72;
- prove two-sided distance at least 36 for the final 12-node block.

These bounds give exactly 98,316 over the full orbit. The degree-one witness
shows that both bounds are tight.

For packet values 1 through 14, retain the 24-node distance-97 target. That
target exceeds the global requirement without using the final 12 nodes.

The next certificate should apply this asymmetric rule to the eight maximal
supports. It should stop on any counterexample and otherwise cover all 57
target supports by inclusion.

## Artifacts

- Generalized primary source:
  `scripts/certify_riffle_dp_2lap_g4_component_24node.cpp`
- Sequential runner:
  `scripts/run_riffle_dp_2lap_g4_component_24node.py`
- Primary run receipt:
  `explorations/riffle_dp_2lap_g4_component24_primary_run.json`
- Counterexample receipt:
  `explorations/riffle_dp_2lap_g4_component24_primary_s49_v15.json`
- Independent verifier:
  `scripts/verify_riffle_dp_2lap_g4_component24_s49_v15_counterexample.py`
- Independent verification receipt:
  `explorations/riffle_dp_2lap_g4_component24_s49_v15_counterexample_verification.json`
- Audit source:
  `scripts/audit_riffle_dp_2lap_g4_component24_dimension30_counterexample.py`
- Audited coverage ledger:
  `explorations/riffle_dp_2lap_g4_component24_dimension30_counterexample_ledger.json`
