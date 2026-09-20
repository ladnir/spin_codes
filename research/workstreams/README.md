# Research workstreams

For the current paper, start with the [artifact guide](../artifact/README.md),
not the chronological research notes. These are the primary entry points:

- [Selected BCH-256 encoder](bare_bch_rm2sub/README.md): implementation,
  correctness tests, and performance methodology.
- [IMT inner](inner_design/IMT.md): the new inner's interface, certified
  instances, and source-name mapping; separate from the supported RM2Sub baseline.
- [Finite proof handoff](bch_rm2sub_bridge/PAPER_HANDOFF.md): selected maps,
  exact ledgers, and size-specific replay notes.
- [Asymptotic certificate snapshot](paper_architecture/certificates/single_sampled_ba_rm2sub/README.md):
  growing-outer Structured SPIN evidence.
- [Small-code engineering results](finite_asymptotic_theory/landscape_db/CURRENT_ENGINEERING_RESULTS.md):
  matched/nested parameter sweeps and their interpretation.

Other directories retain supporting studies and earlier designs. Their
conclusions are scoped to their own maps and parameters. Historical `Next
step` paragraphs are not the current manuscript plan. Frozen sources remain
at their original paths to preserve imports and provenance.
