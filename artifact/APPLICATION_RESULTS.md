# Ordinary encoding, SPIN–Brakedown, and Flock

The paper now presents the application chain in three stages: ordinary
encoding, a standalone SPIN–Brakedown PCS, and integration into Flock.
The application prose focuses on the changes, their purpose, and the result.
Kernel tuning and experiment history remain in the implementation records.
PCS and Flock tables display whole milliseconds. Their proof-size columns
use one decimal place consistently because some entries are below 1 MiB.
The pinned measurements and calculations retain full precision.

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
Without `--check`, it regenerates the four TeX table bodies. With
`--source PATH`, it imports the selected summaries from a Hypercat checkout;
review the resulting manifest before accepting new measurements.
`--bolt-source PATH` separately imports the pinned Bolt opening projection
from the Bolt worktree; importing Hypercat data leaves that record unchanged.
`--ligerito-source SUMMARY_JSON` imports the standalone Ligerito run summary.
It is pinned separately in `paper/data/ligerito_standalone.json`, with its
source-summary hash and raw-run receipts.

| Paper table | Source summary in Hypercat | Aggregation |
| --- | --- | --- |
| Ordinary encoding (`tab:ordinary-encoding`) | `results/spin-brakedown/peach-paired/summary.json` | Median of 31 trials for each `fused-k16/18/20` run. |
| Standalone PCS (`tab:spin-pcs-standalone`) | SPIN `results/spin-brakedown/peach-security/summary.json`; standalone Ligerito below | Ten pooled trials from two processes per configuration. |
| Opening only (`tab:pcs-opening`) | Same SPIN and Ligerito summaries; Bolt source below | Measured `open_ms` medians; Bolt sum of component medians and work proxies. |
| Flock (`tab:spin-flock`) | `results/flock-spin/native-current-comparison/summary.json` | Mean of four process medians; four measured trials after five warmups per process. |

The Bolt commitment comparison comes from `results/bolt-one-thread/summary.json`.
Its SPIN commitment value is 413.46894 ms from the earlier paired comparison.
It is not obtained by subtracting independently aggregated phase medians from
the newer 519.602439 ms complete prover result. Bolt timings are medians of five
trials after one warmup. The fastest Bolt-max run uses SHA-256; the alternative
BLAKE3 run takes 1928.190524 ms. Both use the same one-core hardware conditions.
The builds have different compiler revisions and different layouts; no claim of
a controlled encoder-only replacement is made for this comparison.

### Standalone Ligerito rows

The Ligerito rows are newly measured standalone binary-MLE commitments and
openings, not extracted Flock phase times. The implementation base is Flock
`2d667ca`, the same optimized source used by the integrated comparison.
Revision `cd4189e` adds `crates/flock-prover/examples/pcs_standalone.rs`.
The runner, raw CSV files, machine/compiler/configuration hashes, build log,
summary validator, and timing contract are in `benchmarks/pcs-standalone/`.
The source summary is `results-paper/summary.json` within that directory.
The run records were committed in Flock revision `f73a91f`.
The recorded source hashes are for the Windows-exported archive (CRLF);
the harness, lockfile, and both m32 configurations were checked against
that archive and against the pinned Git blobs after newline normalization.

Both rows use 512 MiB of packed bits and one ordinary 32-coordinate binary
multilinear evaluation. Ring switching receives equality weights on the
first six coordinates and the remaining coordinates through the existing
API. No precomputed partial evaluation is passed to the opener. Recursive
Ligerito uses F256 challenges with the shipped m32 Fast100/Slim100 profiles;
the initial rates are 1/2 and 1/4. The profile names identify parameter
targets, not a new complete composed-security claim.

| Profile | Commit (ms) | Open (ms) | Commit + open (ms) | Verify (ms) | Opening bytes |
| --- | ---: | ---: | ---: | ---: | ---: |
| Fast100 | 942.91 | 3388.98 | 4331.15 | 2.08 | 508496 |
| Slim100 | 1899.70 | 3667.51 | 5567.25 | 1.47 | 271328 |

Timings pool ten trials per profile from two processes, one warmup each,
on Peach CPU15, one worker, requested 4.5 GHz, boost disabled. Rust 1.94.1
uses native CPU features and the bench profile (thin LTO, one codegen unit).
All runs hold the common serial benchmark lock. Internal allocations remain
timed; input generation/copy for the consuming API, setup, independent claim
evaluation, statement serialization, outer serialization, and decoding are
excluded. Each total is the sum of commit and open for one trial before
taking the median. Proof size excludes the serialized commitment (4137 and
2089 bytes, respectively), which is retained separately in the summary.

All 26 proofs verified (24 large proofs including warmups and two smoke
proofs). Each of the six processes also rejected an altered claim. The
m22 smoke runs checked the evaluation against a scalar Boolean-MLE oracle;
large repetitions checked stable proof hashes. These standalone boundaries
differ from Flock, which can supply partial evaluations computed earlier.

### Opening-only comparison

The authors' public Bolt implementation does not include the complete
opening prover. We measured dominant opening computations to obtain
preliminary estimates; the indicated gap did not motivate completing a
Bolt opening implementation for this comparison. This is why the paper
reports calibrated projections alongside measured SPIN openings.

The measured SPIN opening medians are 106.815013 ms (square) and
107.0409855 ms (longer rows). Both come directly from `open_ms` in the
conditional 100-bit runs, not from subtraction of unrelated measurements.
Opening includes folds, transcript work, and authenticated query preparation;
outer serialization and independent evaluation for validation are excluded.
Independently aggregated phase medians need not add to the total median.

`paper/data/bolt_opening_projection.json` retains the 512 MiB standalone
case from Bolt worktree revision `9b45089`, including the SHA-256 of
`tools/standalone/opening-results/projection.json`. The worktree uses upstream
`bcc-research/bolt-rs` revision
`3832e47b24e7b3e10525c9c5bcfc1cfe66d525f2`. Detailed scope, raw CSV files,
source/build receipts, and the component validator are under
`tools/standalone/`; see `CALIBRATED_PROJECTION.md` there.

The amortized-limit model adds row evaluation (636.655120 ms), independent
random column folding (372.763259 ms), two sumchecks (36.095448 ms), and
inner proximity proxies (227.253379 ms), totaling 1272.767206 ms. The
non-amortized model adds 455.895924 ms of leading Mulperm work calibrated
with streaming product scans, totaling 1728.663130 ms. The paper rounds
these projections to integer milliseconds and does not report a measured
opening speedup factor. The generator checks the sums before rendering.

Each component is a median of five trials after one warmup on Peach CPU 15,
one worker, requested 4.5 GHz and boost disabled, using Rust 1.94.1 with native
CPU features. Custom folds were cross-checked, sumcheck endpoints checked,
and all 36 inner-proxy proofs verified outside the prover timer. These are
component measurements, not a complete opening run. The kernels follow
upstream arithmetic and are not claimed to be optimal.

The inner proxies are separate complete Ligerito proximity proofs, not
Bolt's constrained, virtual-syndrome, or batched proofs. Some work is
omitted and other work is substituted, so neither projection is a runtime
bound. A finite amortization batch must pay its share of the matrix proof
and batching work. These runs also do not establish equal composed security
between implementations. The paper retains the equal-input-volume scope
and reports the measured Flock integration separately. The additional
Bolt/Flock table entries use the explicitly hypothetical assumptions below.

### Optimistic Bolt--Flock projection

Table `tab:spin-flock` includes measured SPIN and Ligerito rows alongside
explicitly labeled Bolt projection rows. The latter are a cost
substitution, not an implemented integration or a runtime bound. The pinned
Bolt JSON now also retains the 32 and 128 MiB calibration cases.

We give Bolt zero cost for converting Flock's weighted claims, assume the
proposed joint opening applies, and use the matrix-proof amortized limit.
The joint model retains both row evaluations but shares one independent
random fold, the pair of sumchecks, and the two inner proximity proxies.
It assigns no cost to additional batching preparation or claim adaptation.
It does not assume that Flock's cached evaluations eliminate the row
evaluations at any new challenge point. These choices follow the component
accounting in `tools/standalone/SHARED_OPENINGS.md` and the interface audit
in `CACHE_REUSE.md`, with free conversion as an explicit favorable assumption.

For each workload let O be the SPIN/Flock total minus its commitment and
opening phase summaries. We retain O and add measured Bolt commitment C
and projected shared opening J:

    J = 2 row evaluations + random fold + 2 sumchecks + inner proxies
    projected total = O + C + J

| Witness | Retained O (ms) | Bolt C (ms) | Shared J (ms) | Projected total (ms) |
| --- | ---: | ---: | ---: | ---: |
| 32 MiB | 60.1349785 | 79.717826 | 104.021498 | 243.8743025 |
| 128 MiB | 250.69794925 | 392.690179 | 417.035058 | 1060.42318625 |

O is an accounting residual of aggregated measurements, not a newly timed
phase. The surrounding computation is held fixed by assumption; another
integration could change it. The generator reconstructs J from component
medians and checks both it and the total against the pinned source record.
The Bolt table entries round milliseconds to integers. It does not project verification
time or proof size, and it does not claim a complete matched-security
Bolt/Flock integration. The existing standalone-proxy omissions still apply.

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

### Blaze compatibility

The PCS section now distinguishes the current explicit-fold Brakedown
implementation from the theoretical SPIN composition with Blaze's code-switching
framework. Blaze Sections 1.2.1 and 6 give the generic interleaving interface;
Section 8.2 explains that its inner proof avoids transmitting the full linear
combination. The outer code must be paired with a multilinear-evaluation
IOPP, which checks proximity and the decoded message's evaluation jointly.
Linearity and distance alone do not provide a fast concrete prover for that
interface. The SPIN block maps, routing, and recurrence would need such a
proof; substituting them into the RAA-specific proof is not an implemented
change. Generic proof machinery establishes the theoretical compatibility
and asymptotic proof-size reduction. The paper states that observation;
it does not claim a concrete SPIN--Blaze instantiation, parameter selection,
finite-size reduction factor, or measured performance.

At the current longer-row 512 MiB shape, the two explicit folded messages
contribute 8 MiB and the sampled columns about 4 MiB before authentication.
Replacing the former with an inner proof still leaves the latter. The
additional opportunity is to use longer rows and fewer rows: sampled
columns become shorter while the inner proof handles the longer folded
instance. Total size must include row-evaluation messages, authenticated
columns, Merkle paths, and the inner proof. No SPIN--Blaze proof-size or
runtime projection is included in the benchmark tables.

Source: https://eprint.iacr.org/2024/1609 (Sections 1.2.1, 6, and 8.2).

### Flock integration

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
After adding the opening-only comparison, the draft has 62 pages; the
opening table is on page 36 and Flock begins on page 37. Pages 36--38 were
rendered and inspected. The table arithmetic, finite-integration, repository
hygiene, and whitespace checks pass; the build has no unresolved references
or overfull boxes. No new benchmarks were run.
The PDF is a build product under
`output/pdf/spin_codes_draft.pdf` and is not committed.

The Bolt/Flock projection now appears directly in the Flock table; the
bar chart and its generator were removed. Verification time and proof size
are left unestimated for Bolt. The table arithmetic and PDF layout were
checked after this change. No benchmarks were run.

The standalone Ligerito pass adds measured Fast100 and Slim100 rows to
both PCS tables. The run-summary validation, application-table check,
finite-integration check, and repository hygiene check pass. The changed
PCS pages were rendered and inspected after rebuilding the draft.

### Original BAA comparison in the abstract

The abstract's approximately 3x encoding speedup compares SPIN's 11.104 ms
with approximately 32 ms for original rate-1/2 BAA, not with chosen-block
BAA. Peter confirmed on 2026-09-16 that separate measurements find comparable
original-BAA latency on Ryzen. Those separate run logs are not archived in
this worktree. The companion manuscript's original-BAA rerun records
32.769 ms on Intel; that exact value is not relabeled as a Ryzen measurement.
The archived same-host comparison remains unchanged: chosen Golay and RM
BAA take 23.997 and 27.660 ms, or 2.16x and 2.49x the SPIN latency.
