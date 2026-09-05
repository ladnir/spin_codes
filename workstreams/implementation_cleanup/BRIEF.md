# Workstream brief: implementation cleanup

## Objective

Package and clean the exact Structured SPIN implementation while preserving its
semantics and optimized hot-path structure.

## Inputs

- frozen source manifest and source archive;
- end-to-end performance receipt;
- frozen correctness checksum;
- legacy generators and Peach patch scripts;
- repository `AGENTS.md` performance and benchmark rules.

## Owned scope

Write only within:

- `workstreams/implementation_cleanup/`;
- a new, non-frozen `implementation/` subdirectory beneath the frozen
  construction, if a buildable source layout is required.

Never edit `frozen_source/` or existing receipts.

## Tasks

1. Inventory the exact source and its dependencies in the local libOTe tree.
2. Define the smallest canonical source layout and build entry point.
3. Establish a frozen-vs-clean correctness oracle using checksum
   `0x95c9d722a9539fef` and independent staged checks.
4. Separate generated schedules, setup data, hot-path code, tests, and
   benchmarks without adding runtime abstraction overhead.
5. Preserve fixed-width operations, explicit batching, paired fanout ILP, and
   predictable memory access.
6. Produce a reproducible build/test/benchmark guide and an operation ledger.
7. Run benchmarks serially on Peach only after checking that no benchmark is
   active.

## Acceptance criteria

- construction and randomness distribution unchanged;
- independent correctness checks pass;
- checksum unchanged;
- every hot-path change explained at the code-generation level;
- no unexplained material regression from the 10.823872 ms frozen median;
- exact commands and hashes recorded;
- `MERGE_SUMMARY.md` completed.
