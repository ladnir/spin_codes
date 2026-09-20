# Workstream brief: linear-time audit

## Objective

Audit the asymptotic work of Structured SPIN and determine whether replacing
only the outer constituent can yield a linear-time family.

## Inputs

- frozen implementation and operation receipt;
- current structured outer and RM2Sub inner descriptions;
- BA paper and existing BA spectrum tools;
- SPIN roadmap and theory assumptions.

## Owned scope

Write only within `workstreams/linear_time_audit/`.

This first wave is read-only with respect to implementations and constructions.
Do not benchmark.

## Tasks

1. Define the complexity model, including setup, ordinary encoding, transposed
   encoding, word operations, and bit operations.
2. Derive

   `n/B * C_out(B) + n/t * C_in(t,s) + O(n)`

   for the complete family.
3. Audit how `B`, `t`, and `s` must scale in the available proof arguments.
4. Determine whether constant `t,s` remain plausible as `B=Theta(log n)`.
5. Classify the frozen outer and inner as linear, quasilinear, or superlinear
   local circuits under each cost model.
6. List candidate linear-size outer and inner constituents, but do not design a
   new code until the theory requirements are known.
7. State the exact decision that would show an outer-only replacement is
   sufficient.

## Deliverables

- `COST_MODEL.md`;
- `CURRENT_SCALING_AUDIT.md`;
- `CONSTITUENT_MENU.md`;
- `DECISION_TREE.md`;
- `MERGE_SUMMARY.md`.
