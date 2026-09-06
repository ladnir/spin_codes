# GitHub publication policy

GitHub is the authoritative repository for the manuscript, source code, proof
prose, compact manifests, spectra, and the small machine-readable artifacts
needed to audit a claimed result. The shared integration branch is `origin/main`
at https://github.com/ladnir/permute_conv.

Overleaf is an optional paper-only publishing mirror, not a competing source
of truth. Export the manuscript and required assets from a recorded GitHub
commit. Do not import Overleaf history into research branches. If edits are
made in Overleaf, review and commit their content in GitHub before the next
export. Existing Overleaf-only edits must be preserved for that review.

All active research branches should descend from the curated GitHub history.
Use separate worktrees and ordinary scoped commits, then integrate through
reviewed merges or pull requests. Do not rewrite research history to satisfy
Overleaf limits. See `collaboration/GITHUB_FIRST_HANDOFF.md` for the migration
handoff and remaining work.

Do not commit raw experiment dumps, generated receipt trees, build products,
profiling output, or large binary tables. A compact JSON, CSV, or text result
may be committed when it has a stable descriptive name, is referenced by the
proof documentation, and is needed to reproduce or audit a stated claim.

The repository hygiene check rejects tracked files larger than 5 MiB and known
bulk-output formats or directories. If a future proof artifact genuinely needs
more space, store its generator and a checksum in Git, then publish the artifact
through an external release or artifact store after review.

Before pushing, run:

```text
python scripts/check_repository_hygiene.py
```
