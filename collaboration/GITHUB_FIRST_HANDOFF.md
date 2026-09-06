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

## Current checkpoints and ownership

- At preparation, fetched `origin/main` was `705f4a3`.
- The finite-theory activation tranche is `88aa11a`, descended from that main.
  Work continuing in `0add` remains owned by the finite-theory agent.
- `codex/bch-github-integration` starts at `88aa11a` in
  `C:/Users/peter/repo/permute_conv-github-bch`.
  It carries this handoff, the publication-policy update, and the two-line
  CRLF correction for the activation-pilot JSON receipt dependencies.
- Our old BCH worktree is
  `C:/Users/peter/.codex/worktrees/ba80/permute_conv`, on
  `codex/bch-m22-proof-notes`. Its `c3d8225` estimator import belongs to the
  old history and should not be replayed onto the curated branch.
- The current BCH bridge under that old worktree's
  `workstreams/bch_rm2sub_bridge/` is still untracked and has NOT been migrated.
  The BCH agent owns its checkpoint, compact artifact selection, and validation.
  Do not delete or retire the old worktree before that process completes.

Next integration task: migrate the BCH bridge once, retain the evidence needed
to reproduce or audit its claims, and test it beside the activation estimator.
Then continue research on branches sharing the curated GitHub history.
