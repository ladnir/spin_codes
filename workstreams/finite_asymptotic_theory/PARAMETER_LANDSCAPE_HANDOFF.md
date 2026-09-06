# Handoff: finite Structured SPIN landscape and BCH--RM2Sub closure

Date: 2026-09-06.

Active goal: update the BCH-64/128 engineering grid with higher-occupation
evidence. Read `landscape_db/BCH_DOMINANCE_ANALYSIS.md`. The Q1 curves do
not establish full-margin scaling or implementation state-size choices
until the remaining occupation sum is controlled. The new multi-bit
four-state evaluator retains zero syndromes explicitly. A wider sparse
counting search and complete dense type covers are being evaluated;
partial coverage must remain labeled. BCH-256 is excluded from this goal.

The active larger-t calculation is BCH-128, K=2^20, t128/s20. The
`close_bch_dense_v10.py` producer retains a complete Q257..16384 cover
in `landscape_db/bch_dense_v10_b128_t128_s20_e20_q257/`; its bound is
still being tightened. It separates the all-one outer word, optimizes
counting probabilities jointly, uses the new character-based activation
density bound, and isolates zero-count faces when subdividing types.
Do not restart while its process is live. Earlier v5..v9 refinement jobs
were explicitly stopped after retaining their checkpoints.

Read the activation-density subsection in `BCH_DOMINANCE_ANALYSIS.md`.
Seven exhaustive checks validate that transfer, two validate the new
integer subdivisions, and joint-gradient/direct-witness checks pass.
Three selected individual types from weak boxes have direct-refined
positive margins, demonstrating search slack at those types. This is
not a complete t128 result. `verify_bch_full_reference_v3.py` is prepared
for the new transfer but has not run. After the dense job, compute the
t128/s20 Q5..256 interval with `close_bch_sparse_tail_v2.py`, then replay
the full reference. Preserve the existing source-bound checkpoints.

Both BCH blocks now have complete audited references at K=2^20, t64/s20:
BCH-64 has 10.305616126 full margin bits (0.081331731 bits lost to higher
occupations); BCH-128 has 32.769596142 (0.000001995 bits lost).
`verify_bch_full_reference_v2.py` replays every selected Q2..4 composition,
all higher intervals, complete disjoint type coverage, and selected
90-digit dense witnesses. These are binary64 diagnostics.

The exact integer kernel scan has 46 positive first-moment lower bounds,
all replayed at 90 digits. For both blocks at K=2^20, it obstructs t64/s7..8,
t128/s8..16, and t256/s9..20. The updated 130-row engineering report has
two full Q1-dominant references, 46 obstructions, and 82 sparse-only rows.
Read `landscape_db/BCH_DOMINANCE_ANALYSIS.md` and run
`landscape_db/report_bch_evidence_v2.py` for the joined plot and evidence.
Next test t128/s17..20, the t64 knees, and selected K extremes. The goal
remains active; no full-margin surface over the unresolved range is claimed.

Latest engineering study: `landscape_db/CONSTITUENT_ENGINEERING_SURFACES.md`
extends the uniform-refresh Q1 analysis to RM and random constituents.
It adds K/s interaction slices, family comparisons, and a long-region
model with the actual cancellation probability 1/(2^s-1). Exact RM stops
at length 512; random references stop at 512 for this focused study.
Random ensemble averages and conditional spectrum caps are distinct.
The next question is whether selected higher occupations follow these
Q1 trends. BCH-256 work remains with its owning thread.
The completed study contains 1,292 family rows over 387 geometries. All
85 tests pass, including the new kernel and model checks. A 90-digit replay
checks 1,155 coefficients. At K=2^20, the long-region model predicts the
44 selected family/state cases within 0.105 bits without fitted parameters.
This is measured agreement, not an outward error guarantee.

The following update records the preceding BCH-only stage.

Current research update: `landscape_db/BCH_GROWTH_ANALYSIS.md` supersedes
the earlier four-size BCH fit and implementation-first priorities below.
A new four-state Q1 transfer removes repeated density loss after exact
uniform refresh. Its 212-tuple study shows a nearly one-bit cost per
doubling of k and explains the high-state BCH-64/128 plateaus with a
first-activation model. All 77 tests pass; an independent 90-digit replay
checks 236 coefficients at four reference points. These are Q1 diagnostics,
not new full-distance certificates. Next, tighten selected higher
occupations before using the size trend for full-code extrapolation.
The BCH workstream reports its completed BCH-256 outward closure; references
below to that closure being open are historical. Its bounded spectrum
remains outside our primary exact-spectrum study.

GitHub now carries source, documentation, and catalog configuration only.
All generated result paths below refer to local ignored outputs. The
existing data in this worktree is preserved. Do not push the earlier
experiment history or create a data release. Check new Git objects with
`landscape_db/check_source_only_git.py` before pushing.

Latest grid result: every one of the 3,108 native tuples has full evaluated
occupation coverage, comprising Q1..64 and a typed Q65..L cover, with
additional composition bounds. Most dense bounds remain loose. Read
`landscape_db/GRID_FINDINGS.md` for the final parameter comparisons and
`landscape_db/complete_landscape_audit.json` for the verified snapshot.
`landscape_db/NEXT_GRID_REFINEMENTS.md` records the next implementation
work. Exact-region composition boxes give stronger selected RM full
diagnostics; three random [512,256] full diagnostics are near the shared
60-bit setup budget. No new outward certificate is claimed.

The completed finite-grid scope is in `landscape_db/COMPLETE_GRID_PLAN.md`:
four exact BCH spectra through length 128, four exact RM spectra through
length 512, and random references through length 1024. BCH-256 and
partial-spectrum RM(5,11) are excluded. The primary exact-family Q1
extrapolations, holdout errors, and t/s cost comparisons remain distinct
from full-distance claims. Schema version 4 explicitly records the setup
event for conditional random-spectrum bounds and authenticates event
containment when combining receipts. Numerical producers are resumable
and must run sequentially.

Historical parameter-study update: read
`landscape_db/SMALL_STATE_AND_Q2_STUDY.md` for the smaller-state extension
and exact-spectrum Q2 comparison. The database has 1,012 preferred Q1
observations and 36 current Q2 screens; higher occupations for these maps
remain open. The original tranche is described in
`landscape_db/ACTIVATION_PARAMETER_STUDY.md`, which records the initial
715-observation activation-aware Q1 grid. It prioritizes exact-spectrum BCH/RM inputs and
keeps BCH-256 as a secondary bounded-spectrum anchor. The database now
exposes historical activation-review status and excludes the old RM receipt
from its current-certificate view. The separate BCH task has since closed
the t64_s20 full target below 2^-50; the partial-coverage status below is
historical. Its full receipt remains in the owning ba80 worktree and is not
yet imported here.

This document is the current entry point for the finite-parameter lane. It
supersedes `NEXT_CHAT_HANDOFF.md` for the grid search, parameter model, and
BCH--RM2Sub work. The older handoff remains useful as a history of rejected
constructions.

## Historical workstream context

The practical target is a rate-one-half Structured SPIN code at message length

\[
  k=2^{20},\qquad N=2^{21},
\]

with relative distance at least 10 percent and setup-failure probability below
\(2^{-40}\). A distance near 10.9 or 11 percent is a bonus. The current search
does not discard results below 40 bits; signed margins are used to understand
scaling.

The user requires one fixed outer constituent, reused in every outer row.
Independent resampling of the constituent in each row is not an acceptable
replacement. Setup may still sample routing permutations and RM2Sub
multipliers.

Three pieces are now in place.

1. A local SQLite database indexes 414 finite RM2Sub results with explicit
   evidence labels and source hashes.
2. A narrow diagnostic model extrapolates the BCH occupation-one curves and
   constructs admissible RM2Sub capacity staircases.
3. A separate active task is attempting a complete BCH-derived
   \([256,128,\ge38]\) plus RM2Sub certificate. The original
   \((t,s)=(128,15)\) route has a proved first-moment obstruction. The active
   proof-first route now uses \((t,s)=(64,20)\).

The \((64,20)\) route has an outward Q1 margin of approximately 50.439 bits.
It also has outward coverage of every occupation in

\[
  \{1\}\cup[512,1024]\cup[2048,4096].
\]

The union of these covered classes is below \(2^{-50}\). The full theorem is
not closed. The missing occupations are

\[
  [2,511],\qquad[1025,2047],\qquad[4097,8192].
\]

The immediate task is to close these gaps without changing the construction.
Only after a full positive proof should state-size and performance optimization
resume.

## Repository and ownership

The canonical GitHub checkout for this lane is:

```text
C:\Users\peter\.codex\worktrees\3061\permute_conv_github_publish
```

It tracks `main` at `https://github.com/ladnir/permute_conv.git`. The relevant
commits are:

```text
5e41724 Add finite parameter extrapolation to SPIN paper
36d270d Add RM2Sub finite-result landscape database
f5fc176 Initial curated permute_conv snapshot
```

At this handoff, the canonical checkout is clean and commit `5e41724` is on
GitHub.

The active BCH--RM2Sub proof task is separate:

```text
task:     bch spectrum
thread:   01a06a73-305c-7571-9b2b-a6a75c5639ab
worktree: C:\Users\peter\.codex\worktrees\ba80\permute_conv
branch:   codex/bch-m22-proof-notes
folder:   workstreams/bch_rm2sub_bridge
```

Inspect that worktree read-only while its task is active. Do not edit its
files, restart its sessions, or run a competing large certificate job. Its
generated files are not yet committed to the canonical repository.

The older checkout

```text
C:\Users\peter\.codex\worktrees\3061\permute_conv
```

contains historical local changes. Do not use it as the clean publication
checkout.

## Required startup reading

A new task should read these files before changing the proof or paper:

1. `collaboration/README.md`;
2. `collaboration/INTEGRATION_CONTRACT.md`;
3. `collaboration/TASK_ASSIGNMENTS.md`;
4. `SPIN_NAMING.md`;
5. `RIFFLE_NEXT_WORK_ROADMAP.md`;
6. `workstreams/finite_asymptotic_theory/BRIEF.md`;
7. this document;
8. `workstreams/finite_asymptotic_theory/small_k_replay/TEN_PERCENT_PARAMETER_STUDY.md`;
9. `workstreams/finite_asymptotic_theory/small_k_replay/RM2SUB_PERSISTENCE_UNCONFOUNDED.md`;
10. `workstreams/finite_asymptotic_theory/landscape_db/README.md`;
11. `paper/scaling_complexity.tex`, especially the finite-parameter subsection;
12. the active task's `workstreams/bch_rm2sub_bridge/REVIEW_HANDOFF.md`;
13. the active task's `workstreams/bch_rm2sub_bridge/LARGER_STATE_PROOF_FIRST.md`;
14. the active task's `workstreams/bch_rm2sub_bridge/ZERO_STATE_FIRST_MOMENT_OBSTRUCTION.md`.

Use the Controlled Writing for Cryptography skill for theorem statements,
proofs, or paper revisions. Keep proved claims separate from diagnostic and
estimated claims.

## Design decisions that must remain fixed

### One constituent is reused

Fix one binary linear constituent

\[
  C:\mathbb F_2^{K_B}\longrightarrow\mathbb F_2^B.
\]

The outer encoder divides the message into \(L=k/K_B\) rows and applies this
same map to every row. It does not sample independent row codes. Correlation
between equal or related row messages is therefore part of the proof problem.

### The main distance is 10 percent

The baseline search uses relative distance 0.10. The default failure margin is
40 bits. A result below 40 bits remains useful during the landscape scan. A
10.9- or 11-percent claim should not force an otherwise expensive design.

### RM2Sub is the inner under study

The current database and parameter model contain RM2Sub inners only. A row
labeled `random` uses a random outer-spectrum reference with RM2Sub. It does
not use RandomStepConv.

RandomStepConv, Toeplitz maps, and other random convolutions remain conceptual
comparators. Their transfer bounds do not automatically apply to RM2Sub.

### Performance is a later selection constraint

The desired implementation should remain near or below 11 ms on the Peach
Ryzen 9 7950X reference machine. Do not benchmark a proof candidate until its
proof status makes the result useful. Never run two benchmarks or large
certificate jobs simultaneously.

### The outer at length 256 is BCH-derived

The power-of-two candidate is a deterministic
\([256,128,\ge38]\) subcode of an extended BCH code. It is often called
`BCH [256,128]` in conversation. It is not a pure dimension-128 BCH
constituent. The distinction matters when citing spectrum information.

## Construction and probability space

The active finite construction fixes the BCH-derived code

\[
  C:\mathbb F_2^{128}\longrightarrow\mathbb F_2^{256}
\]

from the BCH-spectrum task. The proof uses its authenticated algebraic and
spectrum constraints; it does not assume a complete weight enumerator.

The message has \(L=8192\) rows. The outer length and final length are

\[
  k=128L=2^{20},\qquad N=256L=2^{21}.
\]

The setup performs the following operations.

1. Encode every row with the same fixed map \(C\).
2. Sample an independent coordinate permutation for each encoded row.
3. Transpose the resulting \(8192\)-by-\(256\) bit array into 256 regions.
4. Sample an independent permutation of the \(8192\) positions in each region.
5. Serialize the regions and apply the RM2Sub recurrence.

For an epoch length \(t\) and state size \(s\), fix maps

\[
  A:\mathbb F_2^s\to\mathbb F_2^t,
  \qquad
  B:\mathbb F_2^t\to\mathbb F_2^s,
  \qquad BA=0.
\]

The current maps use \(B=A^{\mathsf T}\). Identify the state space with a
fixed basis of \(\mathbb F_{2^s}\). For successive input epochs \(X_i\), set

\[
  q_0:=0,\qquad
  Y_i:=X_i+Aq_i,\qquad
  q_{i+1}:=\alpha_iq_i+BX_i.
\]

Each \(\alpha_i\) is sampled independently from
\(\mathbb F_{2^s}^{\times}\). State continues across region boundaries.
There is no final flush. The setup samples all routing permutations and
multipliers once, and every message uses that same setup.

The active cutoff is

\[
  H=209716,
\]

which is one larger than \(\lfloor0.1N\rfloor\). For occupation \(Q\), let
\(Z_Q\) count nonzero messages with exactly \(Q\) nonzero outer rows whose
output weight is at most \(H\). The positive first-moment route seeks outward
bounds \(U_Q\) such that

\[
  \mathbb E_{\mathrm{setup}}[Z_Q]\le U_Q,
  \qquad
  \sum_{Q=1}^{8192}U_Q<2^{-40}.
\]

This inequality gives the desired failure bound by Markov's inequality. It
must cover every occupation with no gaps.

## Evidence vocabulary

Use these terms consistently.

- **Certified:** An authenticated outward bound with a replay and complete
  coverage for the claim it states.
- **Diagnostic:** A reproducible numerical calculation that uses binary64,
  an incomplete witness grid, or incomplete occupation coverage.
- **Estimated:** A calculation using an estimated outer spectrum or another
  explicitly stated empirical input.
- **Reference:** An ensemble calculation used only for comparison.

A certified occupation interval is not a full distance certificate. A full
certificate requires all occupations, all setup terms, and outward aggregation.

## The finite landscape database

The database is under:

```text
workstreams/finite_asymptotic_theory/landscape_db/
```

Important files are:

```text
schema.sql
catalog.json
build_landscape_db.py
query_landscape.py
spin_landscape.sqlite3
landscape_export.csv
test_landscape_db.py
```

The current snapshot contains:

- 414 observations;
- 17 authenticated source receipts;
- 16 outer models;
- 23 RM2Sub configurations;
- 126 BCH diagnostic rows;
- 126 Reed--Muller diagnostic rows;
- 155 random-outer reference rows; and
- 7 rows from the prior outward RM(4,9) calculation.

The deterministic database SHA-256 is:

```text
2386af5e595de6e9faf1ad839020657ca6442c3843fcdaa8a1e9c0dad1a08128
```

Rebuild and query it from the database directory:

```text
python build_landscape_db.py
python query_landscape.py summary
python query_landscape.py curves --family bch
python query_landscape.py curves --family rm
python query_landscape.py curves --family random
python -m unittest -v test_landscape_db.py test_extrapolate_parameters.py
```

The tests currently pass. The source receipt remains authoritative. Database
rows do not upgrade diagnostics or estimates into theorems.

## Grid search completed so far

### Outer families

The scan compares rate-one-half BCH-derived, Reed--Muller, and random-outer
constituents. One fixed structured constituent is repeated in all outer rows.
The random-outer rows use the analytic expected spectrum only as a reference.

The available exact or authenticated structured spectra are:

- extended BCH \([8,4,4]\);
- extended BCH \([32,16,8]\);
- shortened XBCH \([64,32,12]\);
- extended BCH \([128,64,22]\);
- RM\((1,3)=[8,4,4]\);
- RM\((2,5)=[32,16,8]\);
- RM\((3,7)=[128,64,16]\); and
- RM\((4,9)=[512,256,32]\).

The random reference ladder includes block lengths through 1024. Its purpose
is to show the scale available from a random-like outer, not to define the
desired repeated deterministic constituent.

### Persistence-matched family scan

Write

\[
  e:=\log_2 k,
  \qquad
  p:=s+\log_2t.
\]

The neutral scan fixes \(t=64\) and uses

\[
  s(e)=\max(7,e-4).
\]

Thus \(p=e+2\) for \(e\ge11\). This schedule matches the reset exponent of
the earlier RandomStepConv calibration. It is an empirical coordinate, not a
theorem that margin depends only on \(p\).

At \(k=2^{20}\), the exact-spectrum BCH occupation-one margins are:

| Constituent | Block size | Minimum distance | Q1 margin |
|---|---:|---:|---:|
| extended BCH \([32,16,8]\) | 32 | 8 | 2.846 bits |
| shortened XBCH \([64,32,12]\) | 64 | 12 | 10.089 bits |
| extended BCH \([128,64,22]\) | 128 | 22 | 32.338 bits |

The RM\((4,9)\) matched-persistence curve peaks near 41.462 bits at
\(k=2^{16}\) and falls to 38.941 bits at \(k=2^{20}\). The random-outer
\([256,128]\) reference has 59.817 bits at \(k=2^{20}\). These are Q1
diagnostics or references, not full certificates.

### Fixed-state family scan

A second scan fixes \((t,s)=(64,16)\). This separates state growth from
message-length growth. For extended BCH \([128,64,22]\), the Q1 margin falls
from 38.716 bits at \(k=2^{12}\) to 32.338 bits at \(k=2^{20}\). The same
configuration appears in the matched-persistence scan only at \(k=2^{20}\).

This experiment shows why fixed state is not a neutral scaling schedule. It
also shows that a 40-bit line is a selection target, not a useful filter for
the exploratory database.

### Epoch-length and state calibration

At \(k=2^{16}\) and fixed persistence \(p=20\), the corrected RM\((4,9)\)
screens are:

| RM2Sub | Epochs per region | Q1 margin | Q2 screen |
|---|---:|---:|---:|
| \(t=64,s=14\) | 4 | 42.583 | 87.424 |
| \(t=128,s=13\) | 2 | 42.147 | 85.947 |
| \(t=256,s=12\) | 1 | 41.406 | 83.721 |

The shorter epoch is strongest at every tested occupation in this
equal-persistence comparison. Longer epochs may still reduce implementation
cost because they perform fewer state updates.

An earlier Q2 run falsely reported a large negative margin for
\((t,s)=(128,13)\). It split the tilt witnesses across runs and compared the
aggregate results. The correct proof takes the pointwise minimum over the
union of witnesses for each outer-weight pair before summation. The corrected
Q2 margin is 85.947 bits.

### Fixed-RM occupation ladder

The RM\((4,9)\), \((t,s)=(64,14)\) ladder was extended through all 256
occupations. Binary64 diagnostics gave 42.577982 bits for the complete union,
with Q1 as the bottleneck. A later outward bundle was imported into the
database as certified.

The BCH--RM2Sub integration subsequently found that the old two-state transfer
invariant does not cover the first live state after zero-state activation.
The imported RM\((4,9)\) receipt has not yet been reconciled with that finding.
Treat its `certified` database label as **under re-audit for publication**.
A separate moment-domination argument may preserve its numerical bound, but no
such repair has been recorded in the canonical lane.

This issue does not invalidate the RM2Sub construction. It identifies a gap
in one proof invariant.

## Diagnostic parameter model

The reproducible model is:

```text
workstreams/finite_asymptotic_theory/landscape_db/extrapolate_parameters.py
workstreams/finite_asymptotic_theory/landscape_db/parameter_extrapolation.json
workstreams/finite_asymptotic_theory/landscape_db/test_extrapolate_parameters.py
```

It uses only 15 occupation-one binary64 rows. The rows have exact BCH-derived
spectra, \(B\in\{32,64,128\}\), \(16\le e\le20\), \(t=64\), and
\(p=e+2\). The fitted model is

\[
  \widehat\Lambda_1
  =-14.754+2.142d(C)-0.471(e-20).
\]

The in-sample root-mean-square residual is 0.343 bits. The maximum residual is
0.862 bits. Leave-one-message-exponent-out validation has 0.417-bit RMS error
and 1.072-bit maximum error. Leave-one-block-size-out validation has 1.137-bit
RMS error and 1.941-bit maximum error.

The observed distances of the BCH-derived rate-one-half ladder at block sizes
32, 64, 128, and 256 are 8, 12, 22, and at least 38. The size heuristic

\[
  d_{\mathrm{BCH\text{-}derived}}(B)
  \approx1.188\frac{B}{\log_2B}
\]

fits those four values. Combining the two empirical formulas projects a
40-bit Q1 threshold near \(B=157\) at \(k=2^{20}\). The next power-of-two
candidate is therefore \(B=256\).

This model does not use a BCH-derived \([256,128]\) spectrum or RM2Sub transfer
receipt. It extrapolates Q1 only. The active bridge's 49--50-bit Q1 results
support the choice of block size qualitatively, but they use different state
sizes and certified spectrum constraints. They are not validation points for
the fitted equation.

The model's projections through \(k=2^{40}\) keep the continuous Q1 threshold
below 192 bits. This long extrapolation is useful only for experiment ordering.
It must not be presented as evidence that one fixed \(B=256\) constituent
works asymptotically.

## RM2Sub capacity staircases

The selected \(A\) map is an \(s\)-dimensional subcode of RM\((2,m)\) for
\(t=2^m\). Therefore

\[
  s\le1+m+\binom m2.
\]

The model chooses the smallest tested epoch exponent \(m\ge6\) that satisfies
this capacity condition.

The aggressive diagnostic track uses \(p=e+2\):

| \(e=\log_2k\) | 16 | 20 | 24 | 28 | 32 | 36 | 40 |
|---|---:|---:|---:|---:|---:|---:|---:|
| \(t\) | 64 | 64 | 64 | 128 | 128 | 256 | 256 |
| \(s\) | 12 | 16 | 20 | 23 | 27 | 30 | 34 |

The prior RM\((4,9)\) outward point motivated a conservative track with
\(p=e+4\):

| \(e=\log_2k\) | 16 | 20 | 24 | 28 | 32 | 36 | 40 |
|---|---:|---:|---:|---:|---:|---:|---:|
| \(t\) | 64 | 64 | 64 | 128 | 128 | 256 | 256 |
| \(s\) | 14 | 18 | 22 | 25 | 29 | 32 | 36 |

The active proof-first choice \((t,s)=(64,20)\) at \(e=20\) has offset six.
It intentionally lies above both tuning tracks. The user authorized this
larger state to obtain a positive proof before optimization.

Capacity alone does not certify a new map. Every staircase point needs a
fixed map, an exact \(A\)-spectrum, an exact \(\ker B\)-spectrum, algebraic
checks, and a complete transfer calculation.

## Estimating a missing BCH spectrum

Exact spectra are preferable, but an estimate can guide the experimental
frontier. For a fixed code \(C\subseteq\mathbb F_2^B\), sample a uniform
\(r\)-subset \(S\subseteq[B]\). Define the shortened code

\[
  C_S:=\{c\in C:\operatorname{supp}(c)\subseteq S\}.
\]

Let \(X_w(S)\) count its weight-\(w\) words. For \(w\le r\),

\[
  \mathbb E_S[X_w(S)]
  =A_w(C)\frac{\binom rw}{\binom Bw}.
\]

Thus

\[
  \widehat A_w
  :=\overline X_w\frac{\binom Bw}{\binom rw}
\]

is unbiased over independent uniform choices of \(S\). This is only an
estimator. Without a concentration theorem, its output must remain labeled
`estimated`.

Before using this estimator for the 256-bit constituent:

1. validate it against the exact extended BCH \([128,64,22]\) spectrum;
2. hold out samples and compare the resulting SPIN functional, not only each
   coefficient;
3. repeat several shortening sizes and random seeds;
4. record the generator hash and every seed; and
5. keep certified shell caps separate from statistical point estimates.

The active BCH task has developed deterministic shell constraints. Prefer
those constraints whenever they suffice. Statistical spectrum estimation is
an experiment-planning fallback, not a replacement for a proved cap.

## BCH-derived \([256,128]\) plus RM2Sub history

### First bridge: \((t,s)=(128,15)\)

The first activation-aware bridge corrected the zero-to-live state transition.
Its Q1 margins at \(k=2^{20}\) were:

| RM2Sub | Q1 margin |
|---|---:|
| \(t=64,s=16\) | 49.8589 bits |
| \(t=128,s=15\) | 49.4813 bits |
| \(t=256,s=14\) | 48.8336 bits |

The \((128,15)\) route then obtained outward bounds for Q1 through Q1655 and
several dense partial classes. Those partial results remain valid for that
specific map.

A later lower bound proved a genuine obstruction to the unconditional
first-moment route. At Q2620, the actual expected number of bad messages is
larger than \(2^{19600}\). The argument uses many low-weight BCH words whose
permuted epoch inputs all fall in \(\ker B\). The RM2Sub state then remains
zero and the inner emits its input unchanged.

This lower bound proves that no valid unconditional first-moment upper bound
can close the full \((128,15)\) construction. It does **not** prove that setup
failure is likely. A large expectation can be caused by rare setups with many
bad messages. Turning the obstruction into a failure-probability statement
would require concentration, such as a useful second-moment bound.

The second-moment investigation produced exact dependency identities,
Johnson-scheme constraints, and many verified pair-type bounds. It did not
produce a full second moment. The user then authorized a larger state and
prioritized a positive proof.

### Active bridge: \((t,s)=(64,20)\)

The active task snapshots two larger-state candidates:

| Configuration | Minimum \(A\) weight | Minimum kernel weight |
|---|---:|---:|
| \(t=128,s=19\) | 48 | 6 |
| \(t=64,s=20\) | 16 | 8 |

The map audit checks full rank, \(B=A^{\mathsf T}\), \(BA=0\), distinct
nonzero columns of \(B\), algebraic degree at most two, the complete
\(A\)-spectrum, and the complete kernel spectrum.

The current proof-first candidate is \((64,20)\). Its certified status is:

| Occupations | Status |
|---|---|
| Q1 | outward and independently replayed; about 50.439 bits |
| Q512--Q1024 | outward range receipts and 512-bit replay |
| Q2048--Q4096 | outward range receipts and 512-bit replay |
| covered union | exact sum below \(2^{-50}\) |
| Q2--Q511 | open |
| Q1025--Q2047 | open |
| Q4097--Q8192 | open |

The exact all-one band uses Bernoulli parameter one and counting cost one.
This fixes an earlier surrogate loss for that individual band. At Q8192, all
13 pure-band screens pass, but their adaptive combination still fails. The
remaining loss may require a split by the number of all-one rows or another
fixed-composition refinement.

The current coverage ledger is:

```text
workstreams/bch_rm2sub_bridge/generated/larger_coverage_dense_first.json
```

Its `full_target_proved` field is `false`. Do not infer full coverage from a
successful partial sum.

## Activation-state issue

When the state is zero and the input is a unit vector \(e_J\), the update is

\[
  q'=B(e_J).
\]

The fresh multiplier does not randomize this first nonzero state. The
activation law is supported on only \(t\) states, not on all
\(2^s-1\) nonzero states. A proof cannot immediately replace this law by a
near-uniform live-state distribution.

The active BCH bridge uses an explicit activation-aware state class. Historical
two-state screens without `_activation_` in their filenames are not current
certificates. The discovery also creates a review obligation for the imported
RM\((4,9)\) certificate and any other proof that used the same invariant.

Do not infer that every old numerical value is false. A separate moment
domination could validate the same value. The invariant itself is not a valid
justification for activation.

## Paper integration

Commit `5e41724` added `Finite parameter extrapolation` to
`paper/scaling_complexity.tex`. The section defines the first-moment margin,
the probability space, evidence classes, the BCH diagnostic fit, the RM2Sub
staircases, and the shortening estimator.

The paper builds successfully as a 35-page PDF. Its remaining LaTeX warnings
are pre-existing missing bibliography entries. Repository hygiene and all
eight landscape/model tests pass.

Before publication, revise the sentence that treats the RM\((4,9)\) result as
a complete outward anchor unless the activation-aware audit repairs it. The
paper already labels the BCH-256 projection as diagnostic rather than
certified.

## What to do next

### 1. Finish the active positive proof

Do not restart from the \((128,15)\) first-moment route. Continue with
\((64,20)\) and fill the three missing ranges:

\[
  [2,511],\quad[1025,2047],\quad[4097,8192].
\]

Use occupation-specific witnesses when one shared witness fails. Retain every
failed row as a screen, but include only passing rows in the coverage ledger.
For the dense tail, split the exact all-one band or fixed compositions if the
adaptive maximum remains too loose.

Every new range must have:

1. an outward producer;
2. an independent higher-precision replay;
3. hashes of the producer, inputs, and receipt;
4. an explicit interval of occupations;
5. a passing-row ledger; and
6. exact final aggregation with no coverage gaps.

The full positive target is

\[
  \sum_{Q=1}^{8192}U_Q<2^{-40}.
\]

Do not optimize \(s\), \(t\), or XOR count before this target either closes or
meets a proved obstruction.

### 2. Import the result into the landscape database

After the active task reaches a stable checkpoint, add a new catalog entry.
Use `fixed_certified_constraints` for the BCH-derived outer. Record each
occupation interval separately and add a `full_distance` row only after the
complete union passes.

Do not copy multi-megabyte discovery grids into GitHub. Import compact receipts,
source hashes, verification scripts, and the final coverage ledger. The
repository policy excludes bulk experiments and files larger than 5 MiB.

### 3. Re-audit the RM\((4,9)\) anchor

Replay its Q1 and sparse transfer with the activation-aware classes. Determine
whether a moment-domination lemma preserves the old bound. Until then, expose
the review status in the database and paper rather than silently relying on
the old `certified` label.

### 4. Replace the BCH projection with measured data

Once the BCH-derived \([256,128]\) RM2Sub receipt is imported, rerun
`extrapolate_parameters.py`. Add the measured point without fitting to a
looser full-range cap as if it were an exact spectrum. Keep Q1 scaling and
full-certificate scaling separate.

The first tuning comparison at \(k=2^{20}\) should include:

- \((t,s)=(64,16)\), the offset-two performance candidate;
- \((64,18)\), the offset-four candidate;
- \((64,20)\), the proof-first candidate; and
- \((128,19)\), the larger-epoch safety candidate.

Compare complete proof margin, XOR count, state-update count, and measured
runtime. Benchmark only one implementation at a time.

### 5. Extend the size landscape

After the 256-bit point is stable, evaluate several message exponents under
one fixed scheduling rule. Useful values are

\[
  e\in\{16,18,20,22,24\}.
\]

For each point, retain Q1, Q2, a sparse intermediate occupation, a dense
occupation, and the complete union when available. This prevents a Q1-only
model from hiding a changing middle-occupation bottleneck.

### 6. Complete application obligations

A finite distance proof still does not prove the production application. The
remaining application obligations are:

- an arbitrary-length wrapper with its rate and distance loss;
- a seeded-randomness argument or an explicit fresh-randomness requirement;
- equivalence between the analyzed and implemented maps;
- exact XOR counts for the chosen \(A\), \(B\), and field multiplication;
- ordinary and transposed encoder benchmarks; and
- a final manifest binding code, parameters, receipts, and source hashes.

## Common mistakes to avoid

1. Do not sample a new outer constituent for each row.
2. Do not replace powers of a realized spectrum by powers of its expectation.
3. Do not call a random-outer reference a certificate for a fixed code.
4. Do not transfer a RandomStepConv or Toeplitz bound to RM2Sub without proof.
5. Do not use the old near-uniform live-state invariant at activation.
6. Do not infer high failure probability from a large first moment.
7. Do not interpret a vacuous upper bound as a distance counterexample.
8. Do not infer interval coverage from a receipt's filename.
9. Do not merge binary64 diagnostics into an outward certificate.
10. Do not call the \([256,128,\ge38]\) subcode a pure BCH code in formal text.
11. Do not run two benchmarks or large certificate jobs simultaneously.
12. Do not upload large discovery grids or solver work directories to GitHub.

## Best current assessment

The block-size model selected \(B=256\) before the activation-aware bridge was
available. The bridge's 49--50-bit Q1 margins support that choice. The open
question is no longer whether the one-row BCH spectrum is strong enough. The
open question is whether the RM2Sub transfer can suppress every middle and
dense occupation under one reused constituent.

The \((128,15)\) experiment answered that question negatively for the
unconditional first-moment method at those parameters. It did not refute the
code. The larger-state \((64,20)\) route has no analogous proved obstruction
and already covers 2,563 of 8,192 occupations with a union below \(2^{-50}\).
It is therefore the best current route to a first complete BCH-based finite
certificate.

If the remaining dense tail closes, the next phase is optimization downward
from \(s=20\). If it does not close, determine whether the failure is an actual
first-moment obstruction or only an adaptive-envelope loss before changing the
construction.
