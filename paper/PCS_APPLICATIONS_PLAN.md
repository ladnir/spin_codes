# Adding the SPIN PCS and Flock applications

Review base: `3004281944ee75b918c17a62fc1aa164c5a4839a` on `main`.
Working branch: `codex/spin-paper-pcs-applications`.
Date: 2026-09-15.

Drafting update: the subsequent pass adds `pcs.tex`, `flock.tex`, and a separate
ordinary-encoding subsection, with concise explanations of the changes,
motivation, and results. `artifact/APPLICATION_RESULTS.md` records the evidence
and remaining scope. The review and proposal below describe the starting point.

This is a review and proposed revision, not a manuscript rewrite. The review
covers the manuscript's organization, main statements, implementation section,
and the retained application evidence. It is not a new audit of the numerical
distance certificates or every appendix proof. No new benchmarks were run.

## Recommended contribution narrative

The application chain is **SPIN -> SPIN–Brakedown PCS -> Flock**.

The code construction and distance analysis remain the central contribution.
The standalone PCS is the first application: it translates SPIN's forward
encoding performance and finite distance guarantee into commitment and opening
costs. The Flock integration then establishes the benefit inside a complete
prover pipeline. Implementation improvements made during integration are an
additional contribution; the final Flock comparison gives both PCS backends
the applicable improvements.

Use two distinct application contributions in the introduction:

1. A SPIN–Brakedown PCS, with its binary input interface, concrete parameters,
   forward implementation, and standalone comparison, including Bolt.
2. Integration into Flock, including shared implementation improvements and
   the final comparison of SPIN–Brakedown and Ligerito within Flock.

The intended performance claim is fastest proving among comparable PCS
implementations. The retained Bolt experiment currently establishes a faster
commitment path, not a measured complete Bolt opening. Preserve that distinction
until the comparison evidence supports the broader wording.

## Review of the current draft

| Location at the review base | Finding | Proposed revision |
| --- | --- | --- |
| `introduction.tex:205` | The draft explicitly says that SPIN does not instantiate a proof system. This conflicts with the proposed application contribution. | Replace the sentence with the PCS-to-Flock chain and expand the surrounding related-work paragraph. |
| `abstract.tex:33–34`, `introduction.tex:118–158` | The application story ends at transposed encoding; the five contribution bullets omit the PCS and Flock. | Add two application contributions and finish the abstract with standalone and integrated results. |
| `preliminaries.tex:7` | “All vector spaces are over F_2” would conflict with the extension-field challenges and folded messages. | Scope the binary convention to the coding sections; introduce the extension field locally in the PCS section. Preserve the row-vector convention. |
| `implementation.tex:4–8,79–86` | The measured workload is 128 parallel transposed encodings, not ordinary encoding or PCS proving. | Keep those results and their workload definitions. Add a separate ordinary-encoder subsection with its own paired experiment. |
| `finite_certificates.tex:115–163` | The finite theorem gives a setup event for five lengths, with one code shared across messages. It does not itself assign application security. | Reference this event once in the PCS instantiation, then state the conditional protocol bound. Do not charge setup failure once per row. |
| `main.tex:42–50` | The body ends after the encoder implementation section. | Append dedicated PCS and Flock sections before the bibliography. |
| `references.bib`, `artifact/PAPER_MAP.md` | The application references and benchmark evidence are absent. | Add primary references and a compact map from each new claim to its theorem, source revision, and measurement receipt. |

The existing separation between asymptotic codes, finite certificates, and
measured implementation is useful and should remain. The application uses the
finite BCH-256/RM2Sub construction. Do not present its measurements as timings
of the growing Golay–BA constituent used in the asymptotic theorem.

The new sections must distinguish code message length K, codeword length N=2K,
PCS matrix height R, total input bits RK, and Flock compression count. In
particular, 65,536 Flock compressions is not a claim about a 65,536-bit witness.

## Proposed manuscript structure

Keep Sections 2–8 and the coding proof appendices structurally unchanged.
Keep the title. Add a short application preview early in the introduction so
that readers know why both ordinary encoding and transposed encoding matter.

### Section 9: Encoder implementation and performance

Rename the present heading to clarify its scope. Retain the current transposed
tables and comparisons. Add an ordinary-encoder subsection before the comparison
with other transposed encoders:

- State that the ordinary implementation evaluates the same selected map.
- Explain the forward recurrence optimization using the manuscript's notation:
  from Y=X+A(Q) and CA=0, obtain C(Y)=C(X). Gathering and emission can therefore
  share a pass before computing the state update.
- Report the controlled paired forward/transpose experiment. At K=2^20, the
  fresh-process confirmation supports parity within about 1%; do not substitute
  the older unbalanced laptop measurements.
- Keep the 128-parallel-instance workload explicit. It is not a latency for one
  binary word of K bits.

Suggested addition: a compact paired table, with detailed kernel and build
information in the artifact. Do not pool this experiment with the existing
11.259 ms or 11.104 ms transpose experiments.

### Section 10: Polynomial commitments from SPIN

New file: `paper/pcs.tex`, label `sec:spin-pcs`.

1. **Interface and construction.** Define the committed Boolean evaluation
   table, public matrix/index mapping, extension-field evaluation points,
   commitment, and opening. Explain the Brakedown matrix protocol and the
   specific SPIN instantiation. Credit the inherited protocol explicitly.
2. **Distance and parameters.** Show that extending a binary generator to
   F_(2^128) preserves its minimum distance. Connect the finite theorem to the
   selected K=2^16 and K=2^18 shapes. Define the fixed-matrix bad event before
   stating its conditional bound. Specify message/challenge ordering and the
   one random testing fold plus evaluation fold.
3. **Implementation.** Explain bounded row batches, packed storage, the natural
   input mapping, batched column processing, and fast row combinations. Include
   only mechanisms that explain the measured cost.
4. **Standalone evaluation.** Report commit, open, verify, proof bytes, matrix
   dimensions, and memory. Present the fast and smaller-proof shapes. Follow
   this with the comparison to prior PCS implementations, including Bolt.

Use separate tables for complete SPIN opening results and the measured Bolt
commitment comparison. A missing Bolt opening measurement is marked unavailable,
not zero. Normalize the existing Bolt comparison by input bytes and disclose
the different coefficient fields and layouts.

The conditional calculation currently has the form

    a / 2^128 + (1 - a/N)^t,  where a = floor((d-1)/3) + 1.

At the two standalone shapes, t=2045 makes this fixed-matrix expression below
2^-100. A PCS security statement also needs the applicable extraction and
transcript argument. Import and specialize established results where possible;
do not treat the fixed-matrix calculation alone as the complete theorem.
Put this scope in one parameter paragraph and in table terminology, without
repeating it throughout the evaluation.

### Section 11: Application to Flock

New file: `paper/flock.tex`, label `sec:spin-flock`.

1. **PCS interface required by Flock.** Explain the two witness functionals,
   including the low-coordinate interpolation weights. State how the
   SPIN–Brakedown opening supports them and how the layout/basis conversion
   preserves the committed witness.
2. **Integration and implementation improvements.** Separate the PCS adapter
   from the direct linear-check evaluator, zerocheck arithmetic, buffer reuse,
   and other improvements that also apply to the Ligerito backend. Explain
   native C-column-bank reuse briefly; place detailed kernel accounting in
   an implementation appendix or artifact note.
3. **Matched final comparison.** Lead with the current Fast100 versus SPIN
   measurements in one executable. State shared input, hardware, worker count,
   applicable kernels, compiler settings, warmup, and aggregation. Retain each
   backend's native layout and opening protocol. Compare complete integrations.
4. **Performance and communication.** Report total prover time, verification,
   and serialized proof size at both workloads. Give the large-workload phase
   breakdown, identifying commitment/opening as the main difference.

Use one principal timing/proof-size table and one compact phase figure or table.
Put Fast's stricter profile in an appendix or supplementary table. The original
unoptimized Flock can appear as a clearly separate implementation ablation if
its conditions are made comparable; it must not supply the PCS-replacement
speedup denominator.

Describe the current benchmark as proofs of BLAKE3-compression constraints.
The adapter's verifier currently takes the circuit and proof, not public hash
inputs/outputs. A claim about a complete public hash application requires that
binding and corresponding rejection tests. The two functional openings also
need the composed extraction argument before a complete security claim.

### Front matter and supporting material

- Add the two application bullets after the finite-code contribution and update
  the contribution count from five to seven.
- Preserve the principal code theorem in the abstract; shorten the construction
  detail there to make room for two application sentences.
- Extend the related-work discussion to Bolt, Flock, Ligerito, and the generic
  proximity result used in parameter selection. Verify primary bibliographic
  records when drafting; implementation revision citations are supplementary.
- Update the organization paragraph, `main.tex`, `paper/README.md`, and
  `artifact/PAPER_MAP.md` together.
- Add a PCS/protocol appendix for the full parameter and composition argument
  if it would interrupt the application narrative.
- Budget roughly 4–6 additional main-text pages as a first draft, then shorten
  against the submission constraints. Do not assume the current draft has room.

## Results available for the new sections

These figures come from retained local records, not new runs in this worktree.

| Experiment | Result | Scope |
| --- | --- | --- |
| SPIN–Brakedown, 512 MiB, square shape | Commit + open 519.60 ms; verify 15.31 ms; proof 18.198 MiB | One opening; conditional 100-bit parameter target. |
| SPIN–Brakedown, same input, longer rows | Commit + open 548.43 ms; verify 13.42 ms; proof 12.432 MiB | Different shape, same input volume. |
| Bolt-max commitment, 512 MiB | 1719.90 ms versus SPIN 413.47 ms; 4.16x ratio | Fastest measured Bolt variant; commitment only. Fields, layouts, and implementations differ. |
| Flock, 16,384 compressions | SPIN 111.994 ms; Ligerito Fast100 216.479 ms; 1.933x throughput | Current complete prover comparison. |
| Flock, 65,536 compressions | SPIN 412.836 ms; Ligerito Fast100 567.928 ms; 1.376x throughput | Current complete prover comparison. |
| Large Flock proofs | SPIN 13,457,008 bytes; Ligerito Fast100 415,657 bytes | Measured serialization, not estimates. |
| Large Flock verification | SPIN 45.90 ms; Ligerito Fast100 27.73 ms | Retain alongside proving time. |

For the large Flock workload, commitment plus opening accounts for about
139 ms of the 155 ms total difference. This locates the implementation benefit;
it does not attribute every millisecond to the encoder in isolation.

Choose one proof-size unit for the manuscript. The standalone table above uses
MiB; the Flock byte counts are exact and can be converted without ambiguity.
Do not combine them with decimal MB labels.

## Evidence and reusable material

Application source checkout:
`C:/Users/peter/.codex/worktrees/flock-integration/hypercat`, revision `015a4f8`.
Flock implementation checkout: sibling `flock`, revision `2d667ca`.
Raw benchmark records identify their measured source bytes and executables;
the later archival commits do not replace those receipts.

Paths below are relative to that Hypercat checkout:

| Topic | Evidence |
| --- | --- |
| Forward encoder | `docs/spin-forward-peach.md`; `results/spin-brakedown/peach-paired/` |
| Standalone PCS | `docs/spin-100-bit-profile.md`; `results/spin-brakedown/peach-security/` |
| Conditional bound | `docs/spin-folding-security.md`; `docs/spin-fold-budget.md` |
| Bolt timings | `docs/bolt-one-thread-comparison.md`; `results/bolt-one-thread/` |
| Bolt communication accounting | `docs/bolt-proof-size-accounting.md`; `results/bolt-proof-size/` |
| Current Flock comparison | `results/flock-spin/native-current-comparison/README.md` and `summary.json`; `results/flock-spin/native-current-final/` |
| Flock interface | `hypercat-flock/src/lib.rs`; `docs/flock-spin-integration.md` (historical timings) |
| Earlier application prose | `docs/paper/spin-pcs.tex`; `docs/paper/spin-pcs-references.bib` |

The earlier application fragment is a source of protocol prose, not a section
to paste unchanged. Its Flock table, direct RowMajor baseline, implementation
description, and measurement method predate the current comparison. Its
standalone method says five trials, while the current standalone note pools
ten trials from two processes. Resolve these details from the final records.

The Bolt full-proof sizes in the accounting note are paper-anchored estimates;
they are not generated complete proofs. Keep them separate from measured sizes.

Follow the paper repository's publication policy: import compact result tables,
source manifests, and generators, not raw trace directories or binary archives.
Pin the external implementation dependencies rather than duplicating their
entire source trees into this manuscript worktree.

## Concrete next revision

1. Draft `pcs.tex` first, including the interface, theorem dependencies, selected
   shapes, and standalone/Bolt tables. Establish exactly what each comparison
   measures before selecting a headline about prior PCS performance.
2. Draft `flock.tex` from the current receipts, rewriting the earlier fragment
   around the optimized native comparison. Keep the two application stages
   visibly distinct.
3. Add ordinary encoding to Section 9 and connect the binary distance theorem
   to the PCS parameter paragraph. Complete the protocol obligations in a
   focused appendix or explicitly delimit the implementation claim.
4. Revise abstract, contributions, and related work after these sections settle.
5. Generate compact tables from pinned results, compile the paper, and inspect
   the rendered pages. No benchmark rerun is needed just to draft these sections.

The requested worktree and this proposal are ready for that drafting pass.
