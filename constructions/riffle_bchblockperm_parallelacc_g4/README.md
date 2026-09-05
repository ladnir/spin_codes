# Riffle BCHBlockPerm-ParallelAcc g=4

**Status:** `PROOF_GYM_EXPLORATION`

This candidate adds an independent 128-bit permutation inside every inherited
extended-BCH block. It then applies the inherited global four-bit packet
permutation and a four-lane parallel accumulator.

Goal 01 proves the exact one-block slice law. Goal 02 transfers the first
genuine double-parity profile to a uniform 384-bit superblock:

- construction: `CONSTRUCTION.md`;
- goal: `GOAL_01_BLOCK_SLICE_CONTRACTION.md`;
- report: `proof/GOAL_01_BLOCK_SLICE_CONTRACTION.md`;
- receipt: `receipts/goal01_block_slice_summary.json`;
- replay: `../../scripts/analyze_riffle_bchblockperm_parallelacc_g4_goal01.py`.
- three-block goal: `GOAL_02_THREE_BLOCK_ZERO_GAP_TRANSFER.md`;
- three-block proof: `proof/GOAL_02_THREE_BLOCK_ZERO_GAP_TRANSFER.md`;
- three-block receipt: `receipts/goal02_three_block_summary.json`;
- three-block replay:
  `../../scripts/analyze_riffle_bchblockperm_parallelacc_g4_goal02.py`.

At BCH weight 22 and output threshold 40, the local block permutation improves
the old one-lane bound by 66.93 probability bits. At relative distance 0.09,
the generic moment is trivial, but exact zero-gap summation gives a nontrivial
\(2^{-11.72}\) bound.

For the three-block profile \((22,22,22)\), Goal 02 gives a rigorous
\(2^{-24.84}\) zero-gap bound. Retaining positive accumulator-state weights
improves it to \(2^{-36.45}\). The next proof step must cover all attainable
BCH weight triples before summing the field-outer spectrum.

Goal 03 measures the deterministic alpha schedule. Uniform random coefficient
indices behave independently at ordinary BCH weights, but the initial lags
preserve many low-weight words. The authenticated support-at-most-13 tail
already contains 3,622 one-data-symbol inputs with BCH profile
\((22,22,22)\). See:

- `GOAL_03_ALPHA_JOINT_SPECTRUM.md`;
- `proof/GOAL_03_ALPHA_JOINT_SPECTRUM.md`;
- `receipts/goal03_alpha_joint_spectrum_summary.json`;
- `../../scripts/analyze_riffle_bch_alpha_joint_spectrum.py`.

A separately registered shifted geometric candidate removes every weight-22
pair in the complete minimum shell at negligible arithmetic cost. The current
candidate remains unchanged.

Goal 04 tests the construction-independent Hölder route. The inequality does
remove every deterministic nonzero alpha coefficient, but the current
\((22,22,22)\) certificate leaves a 15.448-bit shortfall for a 20-bit target.
The route remains viable after a stronger local bound, or as a high-weight
remainder bound combined with exact shifted low-shell exclusions. See:

- `GOAL_04_HOLDER_PRODUCT_ENVELOPE.md`;
- `receipts/goal04_holder_minimum_shell_gate.json`;
- `../../scripts/analyze_riffle_bch_holder_gate.py`.
