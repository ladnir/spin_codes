# SPIN blog draft

Edit [post.md](post.md). [preview.html](preview.html) is a local reading preview, generated from the Markdown. Nothing has been published.

The five figures are editable SVGs in [figures](figures/), with high-resolution PNG exports. They use the talk's orange highlights, binary examples, and construction diagrams, adapted for a light blog page. Alt text and captions are included in the post.

## Figures and evidence

| Figure | Content | Evidence |
|---|---|---|
| 1 | Witness, encoding, commitment, and proof claims | Conceptual diagram, based on the PCS interface in `research/paper/pcs.tex` |
| 2 | Distance and sampled disagreements | Two illustrative word pairs; not code-family distance claims |
| 3 | Sparse input through an outer, permutation, and recursive inner | Actual small linear encoder execution. The repeated toy outer has parameters [8,4,4]; the inner has memory 6. Full inputs, permutation, feedback masks, and output are in `figures/sparse-example.json`. This is not a production SPIN parameter set. |
| 4 | Structured routing | Transpose and regional shuffle from `research/paper/structured_spin.tex`; separate toy dimensions from Figure 3 |
| 5 | FLOCK time breakdown at 16,384 compressions | Loaded from `research/paper/data/application_results.json`, using the September 29 integrated implementation measurements. The builder exports `performance.json` alongside the post. Remaining time is total minus commitment minus opening. Unrounded values are in `figures/flock-breakdown.json`. |

SPIN construction and theorem claims follow the current local manuscript, especially `abstract.tex`, `spin_intuition.tex`, and `structured_spin.tex`. Finite guarantees refer to the listed message lengths and probability over code setup, not full SNARK soundness. The measured PCS integration has conditional fixed-matrix bounds; the manuscript separately discusses extraction and Fiat–Shamir composition.

The larger-workload discussion uses the component data, not an extrapolation from the 16k case. It distinguishes SPIN's improved cost per unit of data from Ligerito's favorable opening scaling. The code-switching integration is future work, with no claimed measured compression ratio or retained speedup.

The updated performance data comes from `hypercat-spin-integration/results/flock-spin/commit-integration-20260929/summary.json`. That directory retains all raw samples, per-process receipts, validation logs, source hashes, and build profiles. The baseline source revisions and precise warmup/aggregation method are also recorded in `performance.json`. The paper and blog now use these same measurements. Earlier paper FLOCK results remain in the application data's `historical_flock` entry.

## Before publication

- Confirm the final paper URL and replace the relative manuscript link with a stable public link and publish `performance.json` alongside the post. The repository link is already public.
- Reconcile final manuscript revisions with the draft's benchmark and distance statements.
- Decide whether to keep the structured-routing detail in the main post or move it to an expandable section.

## Rebuild

`build.mjs` recreates the figures and HTML from the Markdown and measurement data. It uses the bundled Node runtime's `sharp` and `marked` packages. It does not run benchmarks, modify the slides, or modify the paper.

```powershell
& 'C:\Users\peter\.cache\codex-runtimes\codex-primary-runtime\dependencies\node\bin\node.exe' .\build.mjs
```

The preview reflects the Markdown at the last rebuild. SVGs remain editable, but direct edits to generated figures should also be applied to `build.mjs` to preserve them on rebuild.
