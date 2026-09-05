# Integration contract

## Shared baseline

The immutable baseline is **Structured SPIN (B=256, t=128, s=19)**. Its frozen
sources and hashes are recorded under
`constructions/riffle_parityfanout31x33_bchperm_transpose_bitshuffle_splitstate_preaddmul_rm2sub_t128_s19/`.

Do not modify:

- `frozen_source/`;
- `MAIN_CODE_FREEZE.md`;
- existing benchmark and proof receipts;
- `SPIN_NAMING.md`;
- another workstream's owned directory.

A mathematical or algorithmic change creates a separately named SPIN variant.

## Workstream behavior

1. Read the assigned brief and cited source material before acting.
2. Preserve unrelated user changes.
3. Write only within the assigned scope.
4. Distinguish proved claims, numerical diagnostics, conjectures, and proposed
   proof routes.
5. Do not silently promote the modeled outer spectrum or binary64 ledger to a
   theorem.
6. Do not change SPIN terminology. Propose naming changes in the merge summary.
7. Do not merge, rebase, reset, or discard work from another task.

## Benchmark serialization

Only `implementation_cleanup` may run performance benchmarks in the first
wave. It must check for an existing benchmark process before every run. No two
benchmarks may run simultaneously.

## Required handoff

Every workstream must produce `MERGE_SUMMARY.md` containing:

- objective and scope;
- files created or changed;
- decisions recommended to the integration owner;
- claims established and their evidence;
- assumptions and open issues;
- conflicts or dependencies affecting another workstream;
- the next smallest useful task.

The integration owner decides which recommendations enter shared paper,
implementation, and theorem files.
