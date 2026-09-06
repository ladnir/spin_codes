# GitHub-first coordination handoff

## User decision

GitHub is authoritative for both the manuscript and the research workspace.
Overleaf may remain a paper-only publishing copy. Its limits must not drive
research history rewrites, disconnected snapshots, or repeated file imports.

Canonical repository: https://github.com/ladnir/permute_conv

Shared integration branch: `origin/main`. Fetch it before choosing a base;
the commit listed below is a checkpoint, not a permanently frozen main branch.

## Instructions to the other agents

1. Preserve your current changes and continue owning only your workstream.
2. Use a named `codex/` branch descended from the curated GitHub history.
   If your checkout has a detached HEAD, create a branch at that HEAD first.
3. Commit focused changes and share the branch name, commit ID, tests, and
   overlapping paths. Use reviewed GitHub merges or pull requests to integrate.
4. Do not switch, reset, rebase, or clean another agent's active worktree.
   Existing scope ownership and frozen-proof restrictions still apply.
5. Do not merge Overleaf history, force-push shared history, or repeat orphan
   snapshots. Older unrelated branches need a one-time scoped migration.
6. Keep code, proof notes, compact spectra, and essential audit evidence in Git.
   Keep bulk experiments and build outputs outside tracked files. Run
   `python scripts/check_repository_hygiene.py` before pushing.
7. Preserve byte hashes in existing receipts. Correct checkout line endings
   where necessary; do not regenerate expected hashes merely to hide a mismatch.
8. Never run two benchmarks concurrently. Coordinate before starting one.

The manuscript under `paper/` is authoritative after integration into main.
Paper contributions from other branches still require review and integration;
this decision does not automatically choose between conflicting drafts.
Preserve any unique Overleaf edits for a one-time content review. Subsequent
exports should include only TeX, bibliography, and required figures/assets and
record the source GitHub commit. The export mechanism is not implemented yet.

## Current checkpoints and ownership (updated 2026-09-06)

- The shared GitHub checkpoint is now `cf78d822`: finite analysis code without
  the experiment history. Work continuing in `0add` remains owned by the
  finite-theory agent. Its omitted sweep fixtures are not automatically available
  in a fresh clone.
- `codex/bch-github-migration` starts directly at `cf78d822` in
  `C:/Users/peter/repo/permute_conv-github-bch`.
  It carries the compact BCH sources, map snapshots, closure summaries, and
  a checksum manifest. The required bulk audit inputs are local and ignored.
  See `workstreams/bch_rm2sub_bridge/MIGRATION.md` for replay and restore commands.
- Our old BCH worktree is
  `C:/Users/peter/.codex/worktrees/ba80/permute_conv`, on
  `codex/bch-m22-proof-notes`. Its `c3d8225` estimator import belongs to the
  old history and should not be replayed onto the curated branch.
- The BCH bridge has now been migrated and its frozen final audit replayed.
  The exact union, outer LP, 46 shell caps, and closure tests agree with the
  original result. The old worktree is still the recovery source for bulk inputs;
  do not delete it until those inputs have another verified archive.

Next integration task: merge the BCH migration branch, coordinate the estimator
fixture export with its owner, and add BCH-256 evidence across the parameter grid.
Continue research on branches sharing the curated GitHub history. The earlier
`codex/bch-github-integration` branch need not be merged: its policy changes are
carried here without its old experiment tranche.
