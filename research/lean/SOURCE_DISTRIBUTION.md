# Source distribution

The Git snapshot contains the complete Lean source, including the generated
numerical certificates, pinned dependencies, generators, documentation, and
compact verification records. `SpinCodes.Native` is the public import;
`SpinCodes.NativePin` checks the four main statements through that import.

Build objects, transfer archives, compiler logs, temporary replay directories,
and local caches are excluded. Obtain the pinned external dependencies and
rebuild the project sources as described in `FINAL_REPRODUCTION.md`. Existing
source/object reports document the original run; they do not provide the
omitted compiled objects or make a fresh checkout an audited cached assembly.
The cached-assembly collectors also require local evidence outside this source
distribution. Use a source build or project-source replay for a fresh checkout.

The detailed `encoder_native_evidence_coverage.json` exceeds the repository's
5 MiB file limit and remains local. A compact, explicitly derived summary is
included as `scripts/map_data/encoder_native_evidence_coverage_summary.json`,
with the original report's hash. It is not a replacement compiler record.
The human-readable evidence review retains the provenance distinctions.

The paper revision record refers to the complete working manuscript built
during the revision. The scoped Git paper commit contains the approved Lean
presentation edits applied to the committed manuscript; unrelated working-tree
revisions are left for their owners. These are distinct manuscript snapshots.

The `.gitattributes` file preserves exact bytes under this directory so that
Git line-ending conversion does not alter recorded source hashes. Do not
interpret a historical report as certifying later edited sources or a different
paper snapshot.
