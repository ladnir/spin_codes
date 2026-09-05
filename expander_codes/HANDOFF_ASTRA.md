# Astra handoff: expander-code enumerator project

Status date: 2026-09-05

This document is the working context for an independent review and completion
of the expander-code paper and its companion libOTe implementation.  Treat the
paper proofs and numerical certificates as the primary deliverable.  The
implementation is related, but one implementation mode is deliberately
heuristic and is not covered by the paper's distance proof.

## 1. Objective and project boundary

The project gives new minimum-distance proofs for Expand--Accumulate (EA) and
Expand--Convolute (EC) codes.  The published EA and EC analyses replace the
recursive output process by comparatively loose concentration or Markov-chain
bounds.  This paper retains the exact input--output weight enumerator of the
recursive map.  The result is a sharper analysis of the original ensembles,
new regular-expander ensembles, and finite certificates with parameters near
the Gilbert--Varshamov (GV) benchmark.

The paper has two application-facing cases:

1. Binary codes model the bit choices in Silent OT.
2. Finite-field codes model field-valued correlations in Silent VOLE, with an
   emphasis on fields of order near `2^128`.

The paper proves statements about the minimum distance of a sampled generator.
It does not prove a decoder, a deterministic construction, or a complete
Silent OT or Silent VOLE security theorem.  Proposition 9.1 shows how the
code-sampling failure probability enters a later protocol theorem as an
additive term.

Attack work on undocumented aggressive parameter choices was intentionally
dropped.  The present paper focuses on specified ensembles and provable
parameters.  Do not revive the old implementation attack as a paper claim
unless the project scope changes explicitly.

## 2. Authoritative locations and Git state

### Paper repository

- Working tree: `C:\Users\peter\.codex\worktrees\30ff\permute_conv`
- Project root: `C:\Users\peter\.codex\worktrees\30ff\permute_conv\expander_codes`
- Active local branch: `codex/expander-codes-paper`
- Complete-package commit on that branch: `c7698d6`
- Paper source: `expander_codes/paper/main.tex`
- Built PDF: `expander_codes/output/pdf/exact-enumerator-bounds-for-expander-codes.pdf`
- Submission checklist: `expander_codes/SUBMISSION.md`
- CWC audit: `expander_codes/paper/CWC_REVIEW.md`
- Certificate manifest: `expander_codes/certificate_manifest.json`

The Overleaf Git remote is:

```text
https://git.overleaf.com/69e3be6a47477e151481835e
```

Overleaf `master` contains the complete package at commit `3576f47`.  The
separate synchronization worktree is:

```text
C:\Users\peter\.codex\worktrees\30ff\permute_conv_overleaf_sync
```

It uses branch `codex/overleaf-expander-sync`.  Use that worktree to update
Overleaf because the active research branch and the Overleaf history have
different ancestry.  Fetch before pushing and never force-push.

The current paper working tree was clean when this handoff was started.

### Preserved source material

The old enumerator-paper sources were not destroyed.  They are parked under:

```text
expander_codes/notes/previous_draft/
```

The original enumerator is in
`expander_codes/notes/previous_draft/Expander.tex`.  Additional research notes
are in `expander_codes/notes/`.  The top-level legacy `enumerator_paper/`
directory was restored and is not the active manuscript.  Do not edit the
unified paper by moving the old draft back over it.

Related but non-authoritative repositories include:

- `C:\Users\peter\repo\gen-BAA`, for the generalized Block--Accumulate and
  large-field proof ideas.
- The SPIN work in this repository, for random convolutions and regular-noise
  implementation ideas.

The Toeplitz-prefix tangent in
`workstreams/finite_asymptotic_theory/TOEPLITZ_PREFIX_ARGUMENT.md` was reviewed
and deliberately ignored for this paper.  It does not change the current
proofs.

## 3. Terminology and constructions

The generator has the form `G = B T`.  The sparse matrix `B` is the expander
stage.  The invertible rate-one map `T` is either an accumulator or a recursive
convolution.  The paper uses "expander" as the established construction name;
it does not claim a separate graph-expansion theorem.

### Binary EA

EA applies a binary sparse matrix and then the accumulator

```text
y[t] = u[t] + y[t-1].
```

The accumulator has one memory word and a two-state weight process.

### Wrapped binary EC

The main binary EC construction uses memory `m`, an all-zero initial state,
and recurrence

```text
y[t] = u[t] + y[t-m] + sum_{j=1}^{m-1} b[t,j] y[t-j].
```

The coefficients `b[t,j]` are sampled binary coefficients.  "Wrapped" means
that the oldest feedback coefficient, on `y[t-m]`, is fixed to one.  Time does
not wrap around and out-of-range terms are omitted.  This wording caused
confusion before and is now explicit in both the paper and implementation.

The nonwrapping alternative samples the oldest coefficient as well.  Appendix
B gives its exact transfer matrix and first moment.  The finite certificates
in the paper target the wrapped construction; the paper makes no finite
nonwrapping claim.

### Regular expander terminology

"Left regular" means every message coordinate has exactly one edge in each
output region.  "Two-sided regular" additionally fixes every output
coordinate's degree.  At rate one half, `d_L = 2 d_R`.

For binary two-sided regular codes, `d_R` must be odd.  If `d_R` is even, the
all-one message maps to zero before the recursive stage.  This rules out the
apparently natural degree `8/4` point.  The useful rate-one-half left degrees
are therefore congruent to two modulo four.

Two-sided regularity appears to help several schemes, but the paper does not
claim a universal result.  In the present analysis it improves binary EC
substantially.  It does not remove binary EA's weight-two bottleneck.

### Finite-field EC

Let `p` denote a prime-power field order.  The notation `F_p` therefore
includes extension fields such as `F_(2^128)`.

The proved finite-field ensemble differs materially from a binary code viewed
over a larger field:

- Every expander edge has an independently sampled nonzero field label.
- The convolution uses coefficients sampled from the full field.
- Projective counting identifies nonzero scalar multiples of a message.
- Constraint traces count only the independent equations that remain after
  expander collisions and convolution-state dependencies.

Scalar extension preserves the Hamming distance of a binary generator, but
it acts as independent binary maps on the field components.  It does not give
the component mixing required by field-valued VOLE.

## 4. Paper structure and proof map

The manuscript title is *Exact Enumerator Bounds for Expand--Accumulate and
Expand--Convolute Codes*.  The current PDF has 45 letter-sized pages.

- Section 1 motivates the exact-enumerator method, states the contributions,
  separates binary Silent OT from field-valued Silent VOLE, and reviews
  related work.
- Section 2 defines distance, the GV benchmark, the code architecture, all
  expander distributions, and the regular-ensemble sampling cost.
- Section 3 defines expected input--output weight enumerators, exchangeable
  serial composition, transfer matrices, positive markers, and the generic
  certified support-partition lemma.
- Section 4 gives the exact and asymptotic binary Bernoulli EA analysis, then
  finite same-ensemble certificates.
- Section 5 gives the exact wrapped-EC transfer matrix, the endpoint analysis,
  the asymptotic distance theorem, and finite same-ensemble comparisons.
- Section 6 introduces left-regular and two-sided regular incidence.  It gives
  the regional shell law, point-mass bounds, the binary parity obstruction,
  central-block bounds, the rate-one-half frontier, and scalar extension.
- Section 7 gives the labeled finite-field construction, projective union,
  constraint traces, two-sided regular traces, finite-field EC, and the
  singleton refinement.
- Section 8 states the verifier model, proof-to-verifier correspondence, and
  all headline finite certificates.
- Section 9 interprets convolution, regularity, field size, the implementation
  boundary, application composition, and remaining open questions.
- Appendix A gives fixed-column and fixed-row binary ensembles.  The fixed-row
  finite table is a theorem, not a diagnostic.
- Appendix B gives nonwrapping EC without a finite claim.
- Appendix C closes the wrapped-EC endpoint perturbation argument.
- Appendix D gives exact artifact commands and trust boundaries.

The highest-value proof labels for an independent audit are:

- `lem:certified-support-partition`
- `lem:ec-endpoints` and `lem:ec-two-class-perturbation`
- `thm:binary-biregular-first-moment`
- `lem:degree-five-central-block`
- `lem:degree-three-central-block`
- `lem:odd-degree-central-block`
- `thm:regular-binary-certificates`
- `thm:p-ea-trace` and `thm:p-ea-biregular-trace`
- `thm:p-ec-trace`
- `lem:singleton-free-regions`
- `prop:certificate-soundness`
- `prop:application-composition`

## 5. Main certified results

Every probability below is over the sampling of the generator.  It is not an
adversarial search probability.  The code seed is sampled rather than chosen
by an adversary.  The user judged roughly 40 bits sufficient for one sampled
code and roughly 80 bits sufficient when many codes are sampled.  The paper
reports the actual certified exponents and does not impose either convention
as a theorem assumption.

### Binary rate-one-half frontier

All EC rows reach approximately `99.993%` of the binary rate-one-half GV
distance.

| recursive map | expander | `d_L/d_R` | memory | failure bits |
|---|---|---:|---:|---:|
| accumulator | left regular | `64/-` | 1 | 20.0864 |
| wrapped EC | two-sided regular | `18/9` | 4 | 21.6746 |
| wrapped EC | two-sided regular | `14/7` | 6 | 26.0875 |
| wrapped EC | two-sided regular | `10/5` | 15 | 32.4990 |
| wrapped EC | two-sided regular | `6/3` | 79 | 21.4222 |

The `14/7`, memory-6 point minimizes the paper's simple proxy `d_L + m`.
The `10/5`, memory-15 point is the main binary implementation profile.  The
`6/3`, memory-79 point minimizes degree but has a costly recursive stage and a
slow certificate.

The paper also sharpens the original Bernoulli ensembles without changing
their construction.  At rate `1/5` and relative distance `0.05`, the binary EA
density threshold improves from `C > 4.938` to `C > 1.773`.  For wrapped EC
with memory 21, it improves from `C > 10.4344` to `C > 1.111623`.  The dense
regime reaches the corresponding random-linear-code/GV exponent.

### Fields near 128 bits

At rate one half over `F_(2^127-1)`, the paper gives this EC frontier at the
floored field-relative GV cutoff:

| `d_L/d_R` | memory | failure bits |
|---:|---:|---:|
| `28/14` | 3 | 27.6486 |
| `26/13` | 4 | 71.4724 |
| `24/12` | 6 | 39.5003 |
| `22/11` | 12 | 27.2951 |

The `22/11`, memory-12 result reaches the full floored `p`-ary GV cutoff at
length 2,097,150.

For field-valued Silent VOLE over `F_(2^128)`, the paper certifies `98.415%`
of the rate-one-half GV distance:

| `d_L/d_R` | memory | failure bits |
|---:|---:|---:|
| `24/12` | 5 | 45.031739 |
| `26/13` | 4 | 135.201876 |

The 135-bit margin is much larger than the intended sampling margin.  We kept
the profile because reducing the requested margin did not produce a meaningful
implementation speedup.  Do not describe 135 bits as a protocol security
requirement.

With degree `24/12`, memory 5, and a field-relative cutoff, the certified
failure exponents are 45.84, 45.03, 44.24, and 38.78 for field orders
`2^127`, `2^128`, `2^129`, and `2^136`.  The decline is real for this moving
target.  It comes from the increasing GV-relative distance cutoff, not from
weaker field mixing at a fixed output-weight cutoff.  Other rows through
`2^144` are explicitly marked floating-point diagnostics.

## 6. Enumerator and union-bound ideas

The common first-moment argument counts nonzero messages that encode below a
cutoff.  Markov's inequality is used only after computing or bounding this
expected count.  The improvement comes from retaining the exact recursive
state rather than applying a loose concentration bound to it.

The finite-field union over message values needs extra structure.  The proof
uses the same broad idea as recent large-field RAA and Block--Accumulate work:

1. Fix the support geometry.
2. Identify scalar-multiple messages projectively.
3. Count the independent cancellation equations.
4. Group these constraints before applying the final union bound.

For EA, expander collisions determine the constraint rank.  For EC, the
convolution state also determines which equations are fresh.  The trace
transfer records both output weight and fresh constraints.

The finite verifiers partition every support `1..k` into exact, marker, and
dense regimes.  A valid certificate must cover every support without gaps or
overlap.  It is not enough to sample representative weights.

## 7. Certificate package and current validation

The machine-readable map is `expander_codes/certificate_manifest.json`.  It
contains 12 claim groups and pins 18 frozen artifacts.  It also records the
tested source-tree digest, arithmetic precision, commands, and expected
failure exponents.

The tested environment is:

```text
Python 3.14.4
numpy 2.4.4
python-flint 0.9.0
scipy 1.18.0
```

From `expander_codes/`, run:

```powershell
python -m pip install -r requirements.txt
python -m unittest discover -s scripts -p "test_*.py"
python scripts/verify_manifest.py --strict-versions --list
python scripts/verify_manifest.py --run binary-biregular-10-5
python scripts/verify_manifest.py --run headline
python scripts/verify_manifest.py --run all
```

The manifest driver executes certificate commands serially.  Preserve that
property.  Never run two benchmarks at the same time, and do not run two
large certificate jobs concurrently.

The final package audit established:

- 205 unit tests passed.
- Strict dependency, source-hash, and artifact-hash checks passed.
- The complete `10/5`, memory-15 certificate replayed successfully at
  32.4990025 bits.
- Earlier individual audits replayed the large-field EA result at 29.9715
  bits and both `F_(2^128)` EC results at 45.0317 and 135.2019 bits.
- The LaTeX build completed with no warnings, undefined references, overfull
  boxes, or underfull boxes.
- All 45 PDF pages were rendered and visually inspected.
- The final PDF SHA-256 digest was
  `8056F7C1221B2711572562F842C24DDD312BEAC74155CE05DE49184CCF44491F`.

The degree-`6/3`, memory-79 verifier is unusually slow.  Its structural check
passed during development, and its frozen certificate supplies the displayed
result.  A late full replay was stopped after 16 of 17 long stages rather than
allowed to monopolize the machine.  Astra should make a fresh, uninterrupted,
serial replay of this certificate a release-blocking check.  Do not report the
complete artifact set as independently reproduced until that run finishes.

The specialized Python verifiers remain part of the trusted base.  Only the
two binary two-sided formats have independent standard-library structural
checkers.  The interval monotonicity and convexity arguments are mathematical
claims; the numerical verifier evaluates their outward-rounded endpoints but
does not prove symbolic shape facts.  A second implementation would reduce
this trust but is not currently required by the manuscript.

## 8. Highest-priority review risks

The manuscript is complete, but it has not received an independent expert
proof review.  Review these items before cosmetic work:

1. **Wrapped-EC endpoint perturbation.**  Lemma C.1 was closed late.  It uses a
   two-class Schur complement at a semisimple double eigenvalue and claims a
   uniform quadratic remainder and Perron-eigenvector ratio.  Verify the
   compressed generator, branch selection, uniformity in the marker parameter,
   and passage to the logarithm.
2. **Binary central support.**  Check the conditional point-mass arguments for
   right degrees 3, 5, and general odd degree.  Verify the parity conditioning,
   variance lower bounds, endpoint monotonicity, and complementary-support
   map used by the interval certificates.
3. **Finite-field trace accounting.**  Check that projective grouping removes
   exactly one scalar degree of freedom, that occupied coordinates and fresh
   equations are counted consistently, and that the singleton refinement is
   valid in both prime and characteristic-two fields.
4. **Claim-to-verifier correspondence.**  For each row of Table 6, trace every
   numerical branch to its named equation or lemma.  Confirm that every table
   and theorem is present in the manifest.
5. **Slow degree-6/3 certificate.**  Complete the full replay described above.
6. **Literature and attribution.**  Check the related-work section against the
   final versions of the 2026 RAA and Block--Accumulate papers.  The Khabbazian
   ePrint title was updated to *Sharp Minimum-Distance Lower Tails for RAA
   Codes* during the final audit.

Do not weaken a formal claim silently if a proof issue appears.  Record the
counterexample or missing step, identify the dependent theorem and certificate,
and then decide whether to repair the proof or narrow the claim.

## 9. Companion libOTe implementation

The implementation lives in a separate repository and worktree:

```text
C:\Users\peter\.codex\worktrees\expander-codes\libOTe
```

- Branch: `codex/regular-ec-gf128`
- Commit: `c71e9d1`
- GitHub PR: <https://github.com/osu-crypto/libOTe/pull/187>
- PR state at handoff: open, mergeable, and green on Ubuntu, macOS, Windows,
  and the Curve25519 sodium-fallback job.

The PR adds three distinct encoders:

1. `RegularEcCode<10,15>` is the binary rate-one-half profile for Silent OT.
2. `RegularEcFieldCode` is a materialized reference implementation of the
   independently labeled finite-field paper ensemble.
3. `RegularEcStreamingFieldCode<26,4,G>` is the compact field-valued encoder
   integrated into Silent VOLE.

The streaming field encoder is heuristic.  It regenerates a compact seeded
schedule on demand, uses structured permutations, and uses unit or sign labels
instead of independent nonzero labels.  In characteristic two, the sign labels
collapse to one.  Its convolution still generates full-field coefficients and
therefore mixes field components, but the paper's finite-field distance proof
does not cover its expander distribution.

The PR intentionally exposes only two user-facing choices:

- one paper-aligned materialized reference;
- one practical streaming heuristic.

Avoid adding a collection of small ablation variants to the public API.  The
ablation scripts and notes remain evidence, not products.

Silent OT accepts only the binary `10/5`, memory-15 profile.  Silent VOLE
accepts only the field-valued `26/13`, memory-4 profile and rejects the binary
profile.  The implementation uses `CoeffCtx` generically and includes
GF(`2^128`) and Goldilocks tests.  Parameter adapters pad requested dimensions
to a compatible rate-one-half instance.

The PR is proof-motivated, not a performance improvement.  At `2^20`
GF(`2^128`) correlations, the field streaming profile measured 0.792 seconds
median, or 60.6 MiB/s.  The existing `BlkAcc3x32` default measured 0.210
seconds, or 228.5 MiB/s.  The new profile is therefore opt-in and does not
replace the default.

The implementation includes algebraic equivalence tests, configuration and
rank checks, Silent OT integration tests, GF(`2^128`) and Goldilocks Silent
VOLE tests under regular and stationary noise, and end-to-end correlation
checks after each timed VOLE benchmark.  The frontend follows libOTe's existing
single-main dispatch; no standalone test-library `main` was added.

Relevant implementation files are:

```text
libOTe/Tools/ExConvCode/RegularEcCode.h
libOTe/Tools/ExConvCode/RegularEcFieldCode.h
libOTe/Tools/ExConvCode/RegularEcStreamingFieldCode.h
libOTe/Tools/ExConvCode/RegularEcFastSchedule.h
libOTe/TwoChooseOne/ConfigureCode.h
frontend/RegularEcBench.cpp
libOTe_Tests/RegularEcCode_Tests.cpp
```

The evidence behind the heuristic is recorded in:

```text
expander_codes/notes/STREAMING_DISTANCE_ABLATION.md
expander_codes/notes/STREAMING_FIELD_HEURISTIC_AUDIT.md
expander_codes/notes/STREAMING_LABEL_ABLATION.md
expander_codes/scripts/streaming_ec_*.*
```

## 10. Decisions that should remain stable

- Keep the original permute-convolute paper restored.  The unified manuscript
  belongs under `expander_codes/`.
- Use "left regular" when only the message-side degree is fixed.  Use
  "two-sided regular" when both degrees are fixed.
- Keep "wrapped EC" for the fixed-oldest-tap binary construction and explain
  the term at first use.
- Separate binary Silent OT from field-valued Silent VOLE.  Do not present a
  binary matrix over `F_(2^128)` as field mixing.
- In the proved finite-field ensemble, retain independent nonzero expander
  labels and full-field convolution coefficients.
- Keep the streaming field implementation labeled heuristic until a proof
  covers its structured schedule and labels.
- Treat failure bits as code-sampling error.  Do not imply that an adversary
  searches over seeds unless a protocol explicitly gives that capability.
- Keep the implementation opt-in.  Current performance does not justify
  replacing libOTe's default code.
- Preserve serial benchmark and certificate execution.

## 11. Recommended completion sequence

1. Read `paper/main.tex`, `paper/CWC_REVIEW.md`, and this handoff before making
   changes.
2. Perform an independent proof audit in the priority order from Section 8.
3. Run the full degree-6/3 verifier without interruption.  Then run
   `verify_manifest.py --run all` serially if resources permit.
4. Compare the final verifier outputs with every numerical table and the
   introduction summary table.
5. Review the related-work citations and theorem wording.  Keep theorem
   probability spaces explicit under the CWC rules.
6. Review libOTe PR 187 independently.  Preserve the proof/heuristic boundary
   in both comments and user-facing configuration.
7. Choose the target venue.  Apply its class file, page limit, anonymity rule,
   author affiliations, contact metadata, conflict form, and artifact form.
8. Rebuild the final PDF, render every page, and repeat the manifest and unit
   checks after any substantive edit.
9. Update the Overleaf synchronization worktree, fetch first, and push without
   rewriting remote history.
10. Tag or archive the exact paper and artifact commit used for submission.

## 12. Definition of done

The project is ready for submission when all of the following hold:

- An independent reviewer accepts the high-risk proof steps or records and
  resolves every defect.
- The slow degree-6/3 certificate and the remaining manifest entries replay
  successfully in serial execution.
- Every finite number in the paper matches a frozen, hashed certificate.
- The final PDF builds without warnings and passes a complete visual review.
- The venue metadata and formatting are complete.
- The paper and libOTe documentation make no unsupported claim that the
  streaming field heuristic samples the proved ensemble.
- The submitted source, artifact archive, and recorded commit identifiers are
  mutually consistent.

At handoff, the venue-independent manuscript and artifact package are complete.
The remaining work is independent proof review, one expensive full certificate
replay, independent PR review, and venue-specific submission preparation.
