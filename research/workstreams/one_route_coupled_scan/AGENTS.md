# Research constraints

This workstream studies a new code family. Preserve the existing SPIN library,
paper, certified GF256 reference, and cheap-lane proof work.

- Read README.md, CONSTRUCTION.md, PROOF.md, and GOAL.md before continuing.
- Use Controlled Writing for Cryptography for mathematical prose.
- Distinguish proved local statements, finite toy checks, whole-code
  certificates, heuristic predictions, and measured performance.
- Setup is sampled once and fixed for all messages. A fixed-message
  probability is not a minimum-distance certificate.
- Keep one structured global routing stage. Register shuffles and arithmetic
  inside an inner step are not additional global routing passes.
- Never run two benchmarks at once. The root agent owns all numerical runs;
  subagents may derive, review, or edit disjoint source files without running.
- Before numerical runs or builds, keep at least 20 GiB free on the output
  volume. Set explicit time, memory, and output bounds before large jobs.
- Keep only source and compact selected results locally. Do not create raw
  coefficient dumps, randomness tapes, search trees, archives, duplicate build
  trees, or generated reports without explicit approval. Use Python -B.
- Preserve dirty user work and frozen proof helpers. No promotion, paper
  changes, commits, or pushes without authorization.
- Update the progress ledger with what changed, what was checked, the exact
  unresolved obstruction, and the next discriminating test. Do not replace
  conceptual analysis with unbounded parameter sweeps.
