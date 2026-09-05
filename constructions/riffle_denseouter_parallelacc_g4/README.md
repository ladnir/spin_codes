# Riffle DenseOuter-ParallelAcc g=4

**Status:** `PROVED_ASYMPTOTIC_WARMUP`

This candidate packetizes the proved dense-outer-plus-accumulator warmup. It
uses a terminated random sliding dense outer, a uniform permutation of
four-bit packets, and four lane-parallel accumulators.

Goal 01 proves a complete asymptotic linear-distance theorem:

- construction: `CONSTRUCTION.md`;
- goal: `GOAL_01_PACKETIZED_WARMUP.md`;
- proof: `proof/GOAL_01_PACKETIZED_WARMUP.md`.

Goal 02 strengthens the fixed-length inner contraction by using all four state
weights:

- goal: `GOAL_02_JOINT_LANE_CONTRACTION.md`;
- proof: `proof/GOAL_02_JOINT_LANE_CONTRACTION.md`;
- receipt: `receipts/goal02_joint_lane_summary.json`;
- replay script: `../../scripts/analyze_riffle_denseouter_parallelacc_g4_goal02.py`.

Goal 03 locates and improves the finite distance frontier without changing the
construction:

- goal: `GOAL_03_DISTANCE_FRONTIER.md`;
- proof: `proof/GOAL_03_DISTANCE_FRONTIER.md`;
- receipt: `receipts/goal03_distance204_summary.json`.

At binary length 2,097,408, exact counting for one- and two-packet supports
followed by the general moment bound certifies distance greater than 204. The
required outer memories are 97 for 20 failure bits and 137 for 40 failure
bits. The current proof relaxation stops at 205 on a four-packet return; this
is not an upper bound on the construction's true distance.

Goal 01 uses one scalar lane of the parallel accumulator. Goals 02 and 03 use
all four state weights, while retaining a uniform state-path relaxation for
supports of size at least three.

The inherited `outerDense.tex` has two inconsistencies. It claims the benefit
of zero termination but displays only the first \(k\) parity outputs. It also
claims that a span of length \(\ell\) activates \(\ell+M\) parity windows,
which fails when the span contains a zero gap longer than \(M\). This candidate
emits all \(k+M\) parity symbols and replaces the span count by a support-cluster
argument. Its theorem therefore applies to the construction in
`CONSTRUCTION.md`, not literally to the proof in `outerDense.tex`.
