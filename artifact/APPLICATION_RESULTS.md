# Ordinary encoding, SPIN–Brakedown, and Flock

The paper now presents the application chain in three stages: ordinary
encoding, a standalone SPIN–Brakedown PCS, and integration into Flock.
The application prose focuses on the changes, their purpose, and the result.
Kernel tuning and experiment history remain in the implementation records.

## Tables and reproduction

`paper/data/application_results.json` retains compact measurements extracted
from Hypercat revision `015a4f8`. It includes the SHA-256 and relative path of
each source summary. The Flock implementation change is at revision `2d667ca`;
the measured executable/source receipts remain authoritative for the runs.

From the paper repository root:

```
python -B paper/build_application_tables.py --check
```

This checks generated table contents and both conditional parameter budgets
with integer arithmetic. It does not run benchmarks or replay distance proofs.
Without `--check`, it regenerates the three TeX table bodies. With
`--source PATH`, it imports the selected summaries from a Hypercat checkout;
review the resulting manifest before accepting new measurements.

| Paper table | Source summary in Hypercat | Aggregation |
| --- | --- | --- |
| Ordinary encoding (`tab:ordinary-encoding`) | `results/spin-brakedown/peach-paired/summary.json` | Median of 31 trials for each `fused-k16/18/20` run. |
| Standalone PCS (`tab:spin-pcs-standalone`) | `results/spin-brakedown/peach-security/summary.json` | `pooled_optimized`, ten trials from two processes per shape. |
| Flock (`tab:spin-flock`) | `results/flock-spin/native-current-comparison/summary.json` | Mean of four process medians; four measured trials after five warmups per process. |

The Bolt commitment comparison comes from `results/bolt-one-thread/summary.json`.
Its SPIN commitment value is 413.46894 ms from the earlier paired comparison.
It is not obtained by subtracting independently aggregated phase medians from
the newer 519.602439 ms complete prover result. Bolt timings are medians of five
trials after one warmup. The fastest Bolt-max run uses SHA-256; the alternative
BLAKE3 run takes 1928.190524 ms. Both use the same one-core hardware conditions.
The builds have different compiler revisions and different layouts; no claim of
a controlled encoder-only replacement is made for this comparison.

Ordinary encoding processes 128 parallel binary instances. The call maps K
128-bit blocks to 2K such blocks. It excludes commitment, opening, setup,
allocation, and correctness checks. The paired experiment alternates direction
order and uses independent buffers. The largest initial revised run reports
11.085 ms forward; the two later revised process medians are 11.489 and 11.317 ms.
The paper reports the complete three-size experiment and the confirmation range,
without selecting the fastest run as a stable advantage over the transpose.

The existing transposed tables remain separate experiments. In particular,
11.259, 11.104, and 11.352 ms are measurements from different experiments, not
three summaries of the same samples.

## Parameter calculation behind the concise paper paragraph

Condition once on the selected SPIN code having the stated distance. Reuse
that code across rows and proofs. Extending its binary generator to GF(2^128)
preserves distance: each nonzero extension-field codeword has a nonzero binary
basis component whose support is contained in its support. Binary codewords
also embed in the extension, giving equality of minimum distances.

For a complete matrix fixed before the challenge, let d be the code distance,
N its block length, e=floor((d-1)/3), and a=e+1. One uniform random testing
combination followed by t independent column samples gives the sufficient bound

    a / 2^128 + (1 - a/N)^t.

The event is acceptance of a matrix farther than e columns from the interleaved
code, or an incorrect evaluation of its unique nearby codeword matrix. All
folded messages are fixed before column sampling. The values K=65536 and
262144 have N=2K and assumed distances 13108 and 52429. At both shapes, t=2045
makes this expression smaller than 2^-100. Two Flock openings use K=65536 and
t=2110 each; the sum of their fixed-matrix bounds is below 2^-102.

The detailed derivation is in Hypercat's `docs/spin-folding-security.md` and
`docs/spin-fold-budget.md`. These bounds condition on distance and do not include
commitment extraction, Fiat–Shamir losses, or composed Flock knowledge error.
The code's finite certificate covers its ideal setup ensemble. The fixed
benchmark setup is used under the standing distance assumption; the experiments
do not certify that particular seed's distance.

## Integration scope and remaining work

The Flock adapter supports its two weighted witness functionals. Its verifier
currently receives the circuit and proof, not separate public hash inputs and
outputs. The paper therefore describes compression-constraint proofs. A complete
public hash application needs the public-I/O wrapper and rejection tests.
Composed extraction and transcript accounting also remain to be completed.

The final comparison uses the native Fast100 configuration, with the direct
linear-check evaluator and applicable shared arithmetic optimizations. C-column
bank reuse and the dense-tail dispatch bring that native baseline up to date.
The prior direct RowMajor comparison in `docs/paper/spin-pcs.tex` is historical
and is not used in the new paper tables.

Primary references added in this pass were checked at:

- https://eprint.iacr.org/2026/310 (Bolt).
- https://eprint.iacr.org/2026/1329 (Flock).
- https://eprint.iacr.org/2025/1187 (Ligerito).
- https://cic.iacr.org/p/1/1/2/pdf (generic-code proximity).

No new performance measurement was taken during this writing pass. A future
ordinary-encoder comparison with other code families must benchmark their forward
implementations; the current transposed comparison cannot supply those numbers.

## Writing-pass validation

The application-table check and existing `paper/check_finite_integration.py`
pass. The latter preserves the five finite margins, selected-map checks, and
previous transposed timing entries; it is not a full numerical proof replay.
The author draft builds with TeX Live 2026 and has no unresolved references or
overfull boxes. The changed front matter and application pages were rendered
and visually inspected. Existing underfull-box and class/package warnings remain.
The new encoding subsection begins on page 32; the PCS and Flock sections begin
on pages 35 and 36 of this 61-page author draft. The PDF is a build product under
`output/pdf/spin_codes_draft.pdf` and is not committed.
