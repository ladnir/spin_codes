# Workstream brief: finite and asymptotic theory

## Objective

Formulate theorem targets that cover arbitrary finite lengths and asymptotic
SPIN families. Make assumptions and probability spaces explicit.

## Inputs

- `framework.tex` and its flexible corollaries;
- accumulator and random-convolution proof material;
- Structured SPIN proof receipts and proof-status files;
- `RIFFLE_NEXT_WORK_ROADMAP.md`;
- the Controlled Writing for Cryptography skill.

## Owned scope

Write only within `workstreams/finite_asymptotic_theory/`.

## Tasks

1. State the exact finite-length first-moment theorem for general outer
   spectrum, interleaver, and inner transfer law.
2. Define admissible lengths and formulate one wrapper for all other finite
   lengths using padding, shortening, puncturing, or adjacent block sizes.
3. Separate construction existence from a requested distance or margin.
4. Formulate asymptotic parameter schedules for Accumulator SPIN, Random SPIN,
   and Structured SPIN.
5. Determine the expected global margin scale for one-block and many-block
   message classes.
6. State a conditional Structured SPIN theorem using the weakest sufficient
   outer-spectrum and inner-transfer hypotheses.
7. Identify exactly what must be proved for the actual constituents.

## Deliverables

- `FINITE_LENGTH_FRAMEWORK.md`;
- `ARBITRARY_LENGTH_WRAPPER.md`;
- `ASYMPTOTIC_SCALING.md`;
- `STRUCTURED_SPIN_THEOREM_TARGET.md`;
- `MERGE_SUMMARY.md`.

Do not replace a missing proof with empirical evidence.
