# Goal 01: wrap-state uniformity

**Status:** complete. Conditional uniformity and the exact list budget are
proved. The list-size bound itself remains open. See
proof/GOAL_01_WRAP_UNIFORMITY_PROOF.md.

## Objective

Prove the exact distributional effect of the one-time wrap multiplier and
derive the list-size target for the authenticated support-33 population.

## Completion criteria

Goal 01 is complete when:

1. the field modulus and setup order are authenticated;
2. the wrapped state is proved uniform after fixing all inherited randomness
   and any nonzero first-lap terminal state;
3. the zero-terminal and autonomous-prefix transfers are classified;
4. the exact maximum admissible support-33 list size is recorded; and
5. an independent audit reproduces the arithmetic.

## Excluded work

Goal 01 does not prove the required list-size cap. It does not prove the full
construction and does not benchmark the integrated candidate.
