# Riffle DP-2Lap g=4: endpoint-aware component certificate through dimension 30

## Result

Let (S) be a nonempty set of irreducible character components whose total
dimension is at most 30. For packet value (v), form the binary observation
code obtained by evaluating the characters in (S) on all 16 packet slots
through consecutive states of the autonomous transpose recurrence.

Exact exhaustive certificates prove the following two-sided distances:

\[
\begin{aligned}
d_{24}(S,v)&\ge 97 &&\text{for }v\in\{1,\ldots,14\},\\
d_{24}(S,15)&\ge 72,\\
d_{12}(S,15)&\ge 36.
\end{aligned}
\]

Here (d_w) is the minimum of the weight and complement weight over all
nonzero characters in the (16w)-bit observation code. The value-15 bounds
are tight. The pure degree-one character
`0xa685aac60e5acfb3` has weight 72 in the 24-node code and weight 36 in the
12-node code.

This result proves the desired full-orbit character caps for all 57 component
supports of total dimension at most 30. Together with the earlier exact
(C_{18}\oplus C_{20}) certificate, it covers 58 of the 127 nonempty
component supports. The other 69 supports remain open.

## Why the local bounds imply the orbit bounds

The full orbit contains 32,772 nodes and 16 bits per node. Partition it into
1,365 complete 24-node blocks and one final 12-node block. The transpose
recurrence is invertible, so a nonzero character remains nonzero at every
block boundary. The same local distance bound therefore applies after every
transition.

For values 1 through 14, the complete blocks alone give

\[
1{,}365\cdot97=132{,}405>131{,}088.
\]

For value 15, exact endpoint accounting gives

\[
1{,}365\cdot72+36=98{,}316.
\]

The inequalities apply to both weight and complement weight. Thus every
covered nonzero character has full-orbit weight between (3M/16) and
(13M/16), where (M=524{,}352). Equivalently,

\[
|b_{15}(\chi)|\le\frac58,
\qquad
|b_v(\chi)|\le\frac12\quad(v\ne15).
\]

The degree-one witness attains the value-15 endpoint equality and has
full-orbit bias exactly (5/8). The endpoint term is therefore necessary.

## Finite certificate

The 57 target supports are contained in eight maximal supports:

| mask | component degrees | dimension |
|---:|:---|---:|
| `0x50` | 10, 20 | 30 |
| `0x49` | 1, 9, 20 | 30 |
| `0x31` | 1, 10, 18 | 29 |
| `0x32` | 2, 10, 18 | 30 |
| `0x27` | 1, 2, 4, 18 | 25 |
| `0x47` | 1, 2, 4, 20 | 27 |
| `0x2b` | 1, 2, 9, 18 | 30 |
| `0x1f` | 1, 2, 4, 9, 10 | 26 |

A distance lower bound for a maximal subcode also holds for each contained
subcode. It is therefore enough to certify 16 cases for each row: all 15
packet values in the 24-node code and packet value 15 in the 12-node code.

For a code of length (n) and dimension (k), the implementation constructs
(q=\lfloor n/k\rfloor) pairwise disjoint information sets. Suppose a
codeword has weight below (d). One information set then contains at most

\[
r=\left\lfloor\frac{d-1}{q}\right\rfloor
\]

ones. Its restriction determines the codeword uniquely. Exhausting the
radius-(r) Hamming ball on every information set therefore excludes every
such word. Applying the same argument around the all-one restriction excludes
complement weight below (d).

The primary implementation uses fixed-dimension split-half lookup tables.
It passed all 128 cases after checking 6,425,612,528 information vectors. A
second implementation uses independently seeded reversed partitions and
direct recursive fixed-weight enumeration. It passed the same 128 cases and
checked the same number of vectors.

A Python audit independently reconstructs the systematic BCH map, autonomous
transpose recurrence, component kernels, and observation columns. It checks
all 256 receipts, the rank and disjointness of every listed information set,
all candidate counts, all receipt hashes, the 57-support inclusion cover, and
the tight 72/36 witness. The primary and independent partitions differ in
every certified case.

## Evidence and scope

- Primary aggregate receipt:
  `explorations/riffle_dp_2lap_g4_component_endpoint_primary_run.json`
- Independent aggregate receipt:
  `explorations/riffle_dp_2lap_g4_component_endpoint_independent_run.json`
- Independent audit and support ledger:
  `explorations/riffle_dp_2lap_g4_component_endpoint_dimension30_audit.json`
- Certificate and verifier source:
  `scripts/certify_riffle_dp_2lap_g4_component_24node.cpp`
- Sequential runner:
  `scripts/run_riffle_dp_2lap_g4_component_24node.py`
- Audit source:
  `scripts/audit_riffle_dp_2lap_g4_component_endpoint_dimension30.py`

This is an exact finite certificate for the listed component subcodes. It is
not yet a proof of the character caps for arbitrary 64-dimensional
characters. The earlier (C_{18}\oplus C_{20}) result adds mask `0x60`; the
remaining 69 supports require additional certificates or a structural
argument that avoids support-by-support enumeration.

The next bounded tranche is dimension 31. It adds four frontier supports:
`0x2c = (4,9,18)`, `0x33 = (1,2,10,18)`, `0x4a = (2,9,20)`, and
`0x51 = (1,10,20)`. The same information-set calculation predicts
7,716,616,224 vectors per implementation for the complete four-support
endpoint-aware sweep.
