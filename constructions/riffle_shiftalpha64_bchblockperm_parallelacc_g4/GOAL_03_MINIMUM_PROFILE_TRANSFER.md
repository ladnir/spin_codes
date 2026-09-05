# Goal 03: transfer the actual minimum profile

## Question

Does the exact shifted minimum profile \((24,24,24)\) contract sufficiently
at output threshold \(D=188{,}766\), approximately \(0.09N\)?

## Result

**PASS for the complete one-data minimum shell.**

Three independently permuted weight-24 BCH blocks equal a uniform weight-72
slice in 384 coordinates, conditioned on each labeled 128-bit part having
weight 24. The conditioning event has probability

\[
2^{-6.1583525248}.
\]

The positive-state-weight dynamic program treats zero accumulator states
exactly and uses the actual Hamming weight of every positive state. At scaled
cost 33.7, it gives

\[
\Pr[W\le188{,}766\mid(24,24,24)]
\le 2^{-40.2949878173}.
\]

The shifted spectrum contains exactly 1,175 such profiles. Their complete
union-bound contribution is therefore

\[
1{,}175\cdot2^{-40.2949878173}
=2^{-30.0965427758}.
\]

This contribution passes a 20-bit failure target by 10.0965 bits.

## Comparison with the unshifted shell

The unshifted schedule contains 1,364,834 minimum profiles
\((22,22,22)\). The existing profile certificate gives their contribution as
\(2^{-16.0669222272}\). The shifted minimum-shell accounting gains 14.0296
probability bits and changes this shell from a blocker to a passing term.

## Scope

This result closes only the complete one-data minimum shell. A full distance
proof must sum larger one-data profiles and inputs with several nonzero data
symbols. The fixed scaled cost need not be optimal for the bound to hold.

## Reproduction

- `receipts/goal03_minimum_profile_transfer.json`;
- `../../scripts/analyze_riffle_bchblockperm_parallelacc_g4_goal02.py` with
  `--part-weight 24`.
