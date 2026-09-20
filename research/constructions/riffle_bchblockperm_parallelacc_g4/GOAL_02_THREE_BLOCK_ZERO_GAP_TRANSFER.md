# Goal 02: three-block zero-gap transfer

## Question

Does the field double-parity minimum of three active BCH blocks give enough
accumulator contraction at relative distance 0.09?

## Target

1. Fix three BCH block weights \((h_1,h_2,h_3)\), beginning with
   \((22,22,22)\).
2. Average over the three independent 128-bit permutations and the global
   packet interleaving.
3. Bound the joint distribution of packet support and zero prefix states.
4. Sum all zero-cost gaps exactly.
5. Use the positive prefix-state weights to sharpen the remaining gap count.

The proof calculation and a refutation search must run concurrently. The
refutation search should seek interleavings with more zero returns than the
independent-state heuristic predicts.

## Pass condition

Goal 02 passes if it gives a nontrivial, reproducible upper bound for the
three-block low-output event at threshold 188,766. It need not yet sum the
complete field-outer spectrum.

## Status

**PASS for the minimum profile.**  The exact superblock transfer gives

\[
\Pr[W\le 188{,}766\mid (22,22,22)]\le 2^{-24.8379}.
\]

The positive-state-weight refinement improves this bound to

\[
\Pr[W\le 188{,}766\mid (22,22,22)]\le 2^{-36.4472}.
\]

The result is not a full distance proof.  The remaining proof must cover all
attainable BCH weight triples and combine them with the field-outer spectrum.

See `proof/GOAL_02_THREE_BLOCK_ZERO_GAP_TRANSFER.md` for the argument and
`receipts/goal02_three_block_summary.json` for the numerical receipt.
