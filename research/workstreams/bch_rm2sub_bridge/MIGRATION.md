# BCH migration onto the GitHub history

The working branch is `codex/bch-github-migration`, based directly on
`cf78d822`. Its workspace is `C:/Users/peter/repo/permute_conv-github-bch`.
The old `ba80` worktree remains intact as an archive and recovery source.
No unrelated Git history or earlier estimator snapshot was merged.

## What moved

The migration preserves 345 source and documentation files, about 1.69 MB,
from the BCH bridge and the underlying BCH spectrum work. Three of these are
frozen Python dependencies under a historical `generated/` source directory.
Their locations are preserved because the audited modules import them there.
Existing BCH prose on main matched the old notes after line-ending normalization.
No proof statement or frozen arithmetic kernel was edited.

`MIGRATION_MANIFEST.json` records 804 files by path, byte size, and SHA-256.
The tracked set includes code, notes, the manifest, selected map snapshots,
the final closure audit and union ledger, and the migration validation record.
Approximately 638 MB of runtime inputs were copied locally. Most of those
bytes are repeated LP exports; the bulk inputs remain ignored and uncommitted.
Raw-byte Git attributes preserve hashes across checkouts.

The source commits do not fully identify the migrated content: the bridge was
previously untracked. The per-file hashes identify that content. Original
source-worktree hashes in map manifests remain provenance records; the
snapshot hashes authenticate the copied inputs actually consumed by the proof.

## Validation and its scope

The frozen final audit was replayed in this workspace. It reaggregated all
8,192 occupancies, reconstructed the 1,163-row outer LP, replayed 46 exact
shell caps and the OA29 caps, and passed nine tests in six modules.
The resulting margin remains 50.43906854330953 bits. The exact rational
setup-failure bound remains below `2^-50` under the model stated in
`T64_S20_FULL_CLOSURE.md`.

`MIGRATION_VALIDATION.json` records that replay. The wrapper compares every
deterministic field with the frozen receipt. Only recorded subprocess test
timings may differ; the original audit still requires each test subprocess
to succeed. This was not a fresh 512-bit recomputation of all inner bounds.
It does not extend the theorem to other parameters or the full SPIN protocol.

All 804 migrated files passed byte-for-byte checks. Two migration-tool tests
passed, including refusal to overwrite changed files and rejection of path
escapes. Of 24 selected estimator tests, 23 passed. The remaining test,
`test_producer_authenticates_each_refined_occupation` in
`landscape_db/test_composition_boxes.py`, requires omitted sweep fixtures;
the first missing file is `activation_pilot_v1/manifest.json`. The estimator's
full database and data-dependent suite have not been validated here.
No benchmarks or parameter sweeps were run during migration.

## Restore local proof inputs and check the migration

From the GitHub workspace root, on the existing Windows installation:

```text
python -B workstreams/bch_rm2sub_bridge/migrate_legacy_workspace.py --restore --source-root C:/Users/peter/.codex/worktrees/ba80/permute_conv
python -B workstreams/bch_rm2sub_bridge/migrate_legacy_workspace.py --verify
python -B workstreams/bch_rm2sub_bridge/verify_migrated_closure.py
```

The restore helper verifies hashes and refuses to overwrite different content.
It needs the retained source workspace, or an equivalent unpacked directory
with the manifest's relative paths. No public bulk-artifact download is
provided by this commit. A fresh clone alone cannot replay the complete audit.
The existing frozen receipts also retain Windows-style relative paths;
cross-platform execution was not established by this migration.

Python dependencies include NumPy, SciPy, SymPy, and python-flint. Running new
LP solves additionally requires the legacy solver setup; the replay uses
existing solutions and does not invoke an LP solver.

## Next work

Integrate this branch through a normal GitHub merge. Future BCH work should
continue here, not on the disconnected old branch. Coordinate an estimator
fixture export or regeneration procedure with its owner, then add the BCH-256
certified and heuristic comparisons across its message-length grid.
Keep rigorous bound improvements separate from growth-curve evidence.
