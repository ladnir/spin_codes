# GitHub publication policy

GitHub is the canonical repository for source code, proof prose, compact
manifests, spectra, and the small machine-readable artifacts needed to audit a
claimed result. Overleaf may remain a thin paper-editing mirror.

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
