# Storage limits for this workstream

The local disk filled during the proof search. Storage is a hard resource
constraint, not a reason to keep lowering a job's free-space guard.

- Before starting a local numerical run or build, check free space on the
  actual output volume. Keep at least **20 GiB free**, including the estimated
  peak output and temporary-file footprint. If that cannot be met, use the
  configured remote research host or stop; do not launch the local job.
- Read-only inspections and bounded maintenance of existing result indexes
  may run below that floor. Such maintenance must not generate coefficients
  or launch a search, and must preflight at most **5 MiB of total writes**.
- Use bounded remote scratch directories for bulk coefficient searches and
  intermediate trees. Set time, node, memory, and output-size limits before
  launching. Remote storage is not unlimited either.
- Bulk generated proof data belongs on the remote research host by default;
  local storage is for source and compact selected results, not accumulating
  search histories. Include transport archives and extracted duplicates in
  each batch's storage budget. If cleanup is blocked, stop accumulating local
  output rather than treating the unused import allowance as a target.
- After the 2026-10-10 disk-cleanup request, do not add new local bulk proof
  data, randomness tapes, transport archives, or duplicate build trees without
  explicit user approval. Stream preservation transfers directly to durable
  remote storage; do not stage a second local archive. A new turn is not a
  reset of the storage budget. Keep the cumulative local generated-data
  footprint flat or decreasing; source edits and compact result indexes are
  the bounded exception, not permission to keep importing raw results.
- Keep bulk proof trees and audit dependencies remote. Verify them there,
  preserve a durable archive outside `/tmp`, and retain only compact results,
  indexes, hashes, and archive/restore references locally. Do not download
  entire search trees merely to verify or index them.
- Budget at most **16 MiB of new local generated data per work batch**,
  including temporary transport copies. This is a ceiling, not a recurring
  allowance: report cumulative local growth and reclaim disposable outputs
  before another batch. Larger imports require an explicit user request.
  Do not repeatedly fetch full source/build/checkpoint archives when only
  a result changed.
- Preserve source, certified checkpoints, and all dependencies of the selected
  proof index. Ignored files are not automatically disposable. Before evicting
  unique evidence, verify a durable archive and record its location, contents,
  hashes, and restore procedure. Remote /tmp is not a durable archive.
- Remove only verified rebuildable caches or verified redundant copies, using
  exact validated paths. Never remove a whole checkout, rewrite Git history,
  or bypass an execution-policy rejection to make space.
- Do not edit frozen, hash-pinned proof helpers just to change their historical
  storage guards. Apply these stricter checks in new orchestration and before
  invocation; their legacy guards alone are insufficient.
- End each substantial run by reporting the output footprint and reclaiming
  disposable scratch. Never run two benchmarks simultaneously.

On 2026-10-10, recursive cache deletion was rejected by the execution policy.
Lossless NTFS compression was used instead. No certificate data was deleted.
The compiled SPIN proof objects and isolated replay objects remain present;
compression changes storage allocation, not file contents or proof status.
