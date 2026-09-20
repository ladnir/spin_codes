# Paper architecture merge summary

## Decision

Restart the paper around the intended SPIN families. Do not incrementally edit
the current construction story into the new one.

The new paper has one encoder interface,

\[
E_N=I_N\circ\Pi_N\circ O_N,
\]

and three principal instantiations:

| Family | Outer | Interleaver | Inner | Scope |
| --- | --- | --- | --- | --- |
| Accumulator SPIN | independent random blocks | uniform | accumulator | asymptotic |
| Random SPIN | independent random blocks | uniform | random recursive | asymptotic |
| Structured SPIN | one-sampled known-spectrum constituent | structured | RM2Sub-S19 | scalable rate-`1/2`, distance-`0.11`, linear-time theorem; finite specializations separate |

The outer, interleaver, and inner option space remains useful for experiments.
The paper should not attempt to name or analyze every combination.

## Central proof interface

Keep the existing first-moment idea, but present it after the encoder and bad
event are concrete. For interface class `tau`, the exact finite identity is

\[
\mathbb E[Z_d]
=\sum_\tau A_\tau^{\rm out}Q_\tau(d).
\]

For a uniform interleaver, `tau=w` gives the familiar weight-indexed sum

\[
\mathbb E[Z_d]
=\sum_w A_w^{\rm out}p_w^{\rm in}(d).
\]

This is enough common framework. Each family should introduce its own
envelopes and parameter schedule locally.

## What can be written now

The following sections do not require new research results:

1. the SPIN interface and the three-family overview;
2. the exact finite first-moment theorem;
3. the uniform-interleaver specialization;
4. the accumulator definition and exact enumerator;
5. the mathematical definition of the candidate random recursive inner;
6. the Structured SPIN family interface with an explicit parameter selector;
7. the parameterized BA-3 constituent interface and the currently certified
   one-sampled Golay--BA-3/RM2Sub-S19 instantiation;
8. the frozen `N=2^21` reference-member definition as an appendix specialization;
9. the exact ParityFanout weight-transition law;
10. the finite-certificate schema and evidence-status rules;
11. the frozen implementation and performance description.

These parts are now drafted. The scalable certificate bundle and the tightened
`39/4` schedule are integrated. The focused proof appendix now exposes the
analytic reductions, transfer matrices, occupation split, and accepted margins;
the exhaustive rational boxes and machine receipts remain in the artifact
bundle. Remaining work concerns finite-specialization proof, implementation
cleanup, and an optional formal checker.

## What remains open

The Accumulator SPIN and Random SPIN asymptotic theorems are complete for the
adopted independent-injection outer. Their concrete rate-half points are
interval certified.

### Structured SPIN

The main scalable result is now integrated. At native lengths
`L_m=128m`, `N_m=L_m b_m`, with `b_m` the least multiple of `24` satisfying
`b_m>=(39/4) log_2 N_m`, the one-sampled Golay--BA-3/RM2Sub-S19 ensemble has rate
`1/2`, relative distance `0.11` with high probability, and `O(N_m)` ordinary
and transposed encoding work. Zero extension of the largest preceding native
member gives every sufficiently large requested output length with
`o(1)` loss.

The full certificate bundle is imported under
`workstreams/paper_architecture/certificates/single_sampled_ba_rm2sub/` and
its 31 manifest entries pass hash verification. Tightening the block constant
beyond `39/4` would not change the construction or paper spine.
The frozen ParityFanout member remains a separate finite specialization whose
distance ledger still needs a certified actual spectrum and outward arithmetic.

## ParityFanout classification

ParityFanout belongs to the known-spectrum outer branch. It applies a tractable
randomized transformation to the base constituent. Uniform disjoint-set
sampling gives an exact conditional transition from input weight to output
weight.

The current expected transformed spectrum is conditional on a modeled base
spectrum. It is not the exact enumerator of every fixed sampled fanout.
Accordingly:

- the transformation law can be presented as mathematics;
- the current transformed-spectrum receipt can be presented as diagnostic;
- the final finite-distance theorem waits for a certified base spectrum and
  outward arithmetic.

## Old manuscript material

Retain the following as proof sources:

- the bad-codeword first moment and uniform-slice lemma;
- the accumulator enumerator;
- the local random-recursion lemmas;
- range-splitting and exponent techniques from the dense integration;
- finite-ledger and outward-rounding techniques from earlier certificates.

Do not carry forward these old claims as SPIN theorems:

- the sliding-dense outer integrations;
- the `0.109` dense+dense result;
- the RM/EBCH full-split construction;
- BCH projection heuristics.

They concern different component combinations.

## Main-text and appendix split

Main text:

- motivation and construction progression;
- encoder definitions;
- exact first-moment theorem;
- one principal transfer statement per family;
- main distance statements;
- concise finite-certificate partition;
- asymptotic and complexity conclusions;
- benchmark summary.

Appendix and artifacts:

- detailed combinatorial estimates;
- parameter-range casework;
- complete numerical ledgers;
- outward-rounding implementation;
- manifests and hashes;
- implementation-correspondence details;
- historical comparisons.

## Recommended first drafting wave

Create a new manuscript skeleton with section stubs matching `OUTLINE.md`.
Then fill only the interface, first-moment theorem, accumulator enumerator,
ParityFanout transition, and frozen reference-member definition. Mark every
missing integration theorem with a short statement of its required inputs.

This produces a readable partial paper without pretending that the open outer,
routing, or arithmetic work is complete.

## Manuscript integration wave

The user and integration owner expanded this task to create the new manuscript
after the architecture restart.

### Archive

The exact legacy root-level draft was moved to
old/paper_draft_pre_spin_2026-08-31/. The directory contains the 21 TeX files
and the `llncs.cls` dependency present on Overleaf immediately before the new
draft was pushed. ARCHIVE_MANIFEST.md records their byte lengths and SHA-256
hashes after the move. No construction directory, frozen source, receipt,
subdirectory paper, or planning document was moved.

### New entry point

The new manuscript entry point is paper/main.tex. The draft currently
contains:

- a family-level introduction;
- the outer--interleaver--inner SPIN interface;
- an exact class-indexed first-moment theorem and proof;
- the uniform-interleaver specialization;
- the independent uniform-injection outer definition and generating function;
- explicit credit to Block-Accumulate codes for the truncated
  block--interleaver--accumulator architecture;
- explicit related-work credit to Chosen-Block BAA for reusing one
  constituent code across all outer positions;
- credit to Silver for the earlier banded recursive encoder and to
  Expand-Convolute for the closer iid time-varying random-convolution
  architecture;
- the BA local random-block enumerator and independent-block global enumerator;
- the accumulator definition, exact enumerator, sharp contraction, and a
  distance theorem with a certified rate-half point
  `B>=3 log_2 N`, `delta=0.0037`;
- the all-random-tap Random SPIN recursion, exact transfer matrix, asymptotic
  exponent and sparse-range lemmas, and a distance theorem
  with certified rate-half point `B>=17 log_2 N`,
  `m=ceil((51/50)log_2 N)`, `delta=0.11002`;
- the Structured SPIN design narrative, complete modular encoder definition,
  constituent enumerator interfaces, exact routing factorization, structured
  inner recurrence and bijectivity lemma;
- a construction-at-a-glance handoff that gives the full outer--route--inner
  dataflow, parameter roles, and setup law before the component definitions;
- a scalable BA-3 interface parameterized by a base constituent and base
  length, followed by the currently certified Golay--BA-3/RM2Sub-S19 schedule,
  distance theorem, requested-length wrapper, and linear-time analysis;
- a finite-length specialization boundary whose base constituent, parameters,
  certificate, and benchmark receipt remain explicitly open;
- finite-certificate semantics;
- scaling and complexity sections;
- an implementation and performance section;
- explicit red TODOs for unresolved claims.

### Build status

TeX Live 2026 built the manuscript successfully with:

    pdflatex -interaction=nonstopmode -halt-on-error main.tex

BibTeX plus repeated LaTeX passes produced a 34-page PDF with resolved citations,
resolved cross-references, no overfull or underfull boxes, and no document or
package warnings on the final pass. A visual review of the rendered Structured
SPIN section, scalable certificate appendix, and adjacent transitions found no
clipping, overlap, broken tables, or unreadable displays. The packaged output
is `output/pdf/spin_codes_draft.pdf` (SHA-256
`42867ca73f43dbf9b79abbab6abb2bbff6037327f43da2c3cca255724a992fa1`).

### Claims and evidence

The draft proves the general first-moment theorem, its uniform-slice
specialization, and two asymptotic distance theorems for independently sampled
uniform injections. It proves that every realized Structured SPIN encoder is
injective once its local maps are injective. It also states the certified
one-sampled Golay--BA-3/RM2Sub-S19 theorem at rate `1/2`, relative distance
`0.11`, and linear ordinary and transposed work. The block constant `39/4` is
identified as the current certified value. The finite-length section now uses
a parameterized BA-3 outer and makes no distance or performance claim before
the base constituent, finite certificate, and benchmark receipt are fixed.

### Merge risks and dependencies

- The two completed theorems condition each independently sampled block map
  to be injective. The BA enumerator separately covers independent
  unconditioned matrices; the distinction remains explicit.
- The Chosen-Block BAA citation is currently an unpublished 2026 manuscript
  entry. Replace it with its ePrint identifier when that identifier is assigned.
- The proved Random SPIN recurrence samples every time-varying tap vector
  independently. The archived fixed-tap candidate is not a theorem source;
  its audit found a reversed monotonicity step and an invalid independence
  decomposition.
- The scalable certificate bundle is now paper-owned and hash verified. The
  manuscript gives a focused nine-page proof appendix: exact outer spectrum,
  one-sample selection, route domination, three-state and fixed-occupation
  transfers, growing-sparse closure, and theorem completion. Exhaustive box
  covers and outward endpoints remain in the bundle.
- The appendix deliberately separates analytic lemmas from finite rational
  witnesses. Formalization is outside the paper workstream, and the present
  draft makes no formal-proof claim.
- The Structured SPIN exposition pass now defines the selector length,
  probability space, bad-message counts, transfer boundary vectors, spectral
  radius, Collatz witnesses, continuum output fugacity, and weighted norm
  before use. It also gives the $Q=1,2$ certificate its own proposition and
  closes the theorem with a probability union bound.
- The certified block constant is currently `39/4`. A later certificate may
  tighten the schedule equation and endpoints without changing the architecture.
- The finite member needs a selected base constituent and parameters, an
  outward-rounded first-moment certificate, and a manifest-bound benchmark
  receipt before it receives finite-distance or performance claims.
- The implementation section keeps the BA-3/RM2Sub architecture fixed but
  leaves the base constituent open. It still needs the finalized interface,
  correctness receipt, operation ledger, and build and test commands.
- The Silver attribution is intentionally qualified: its banded triangular
  solve is an earlier structured recursive encoder of convolutional form,
  while Expand-Convolute is the closer prior iid random-convolution ensemble.

### Next smallest useful task

Select the finite BA-3 base constituent, then add its parameters and certificate
once the finite proof workstream closes. The performance section should receive
the finalized receipt in the same pass.

### Framework and narrative hardening update

The class-indexed framework now defines the exact conditional failure
probability of a fixed outer word and accepts a certified upper envelope for
each class. Exact equality is recovered when the class determines that
probability. Structured SPIN uses the ordered row-weight profile, and route
domination supplies the required envelope. The framework also states why
conditioning on an outer-only spectrum event preserves independence.

Every sufficiently large requested output length now uses zero extension from
the largest preceding native member. This wrapper preserves injectivity,
absolute minimum distance, and linear work without requiring a nonconstructive
puncturing subspace. The
abstract and introduction now lead with the central rate-one-half,
distance-`0.11`, linear ordinary/transposed result and explain the data-movement
goal. Obsolete reader-facing draft notes and the legacy BCH/ParityFanout finite
narrative were removed; two red TODOs remain because the finite certificate and
implementation receipt are genuinely open.

The related-work subsection now follows the BA paper's three-part taxonomy:
serially concatenated codes, codes for correlation generation, and linear-time
proving systems. It separates lineage from claimed contributions, credits BA
and Chosen-Block BAA for the repeated-block outer construction, and distinguishes
Silver's structured recursive encoder from the iid random-convolution lineage.
The bibliography grew from four entries to eighteen, with the BA citation
updated to its CRYPTO 2026 publication data.

The abstract, introduction opening, and contribution list now lead with the
central theorem rather than the generic concatenation framework. Their main
tension is now intrinsic to the SPIN problem: combine strong binary distance,
linear work in both evaluation directions, and predictable data movement in one
structured family. Silver is not part of the abstract, opening problem
statement, or contribution list; its conjectural distance evidence remains a
precise comparison in the Random SPIN lineage and related work. The four-item
contribution list separates the main code, its modular architecture, the
class-indexed proof framework, and the two explanatory families. It does not
make a finite-performance claim before the benchmark receipt closes.

The opening now gives a concrete cryptographic reason for all three cost
criteria. It identifies the large public-matrix products used by pseudorandom
correlation generators, the row encodings used by code-based proof systems,
and the adjoint identity that makes transposed evaluation operationally
relevant. It also states why long-range memory traffic matters at these vector
lengths. The paragraph is deliberately short so that the Structured SPIN
theorem remains the opening's focal point.

The redundant `Results and scope` subsection was removed. Its headline claims
already appear in the abstract and contribution list, while its schedule
constants and finite-length qualifications remain in the construction,
scaling, and certificate sections where they can be checked against the
corresponding theorems. The introduction now moves directly from the common
first-moment calculation to the organization paragraph.

The preliminaries now follow the useful ordering of the BA paper while keeping
only SPIN-relevant material. They define index, permutation, sampling, and
binary-entropy notation; distinguish a linear encoder from its image code;
define rate, minimum distance, relative distance, and the zero-code bookkeeping
convention; and state that a rate-and-distance claim for a `K`-to-`N` encoder
includes injectivity. They also define the transposed encoder, its adjoint
identity, the bit-operation cost model, sampled-family probability spaces, and
the meaning of `o(1)` and high probability. The former local definition of
binary entropy in Random SPIN was removed.

The Structured SPIN formal exposition now begins with one direct member and a
complete dataflow diagram. It explains the roles of the outer spectrum, the
two-stage route, and the recursive inner before declaring the fixed and sampled
construction data. Family selectors and requested-length scaling no longer
interrupt the transition from the design story to the component definitions.

The scalable outer is now parameterized by an even base length `a`, a rate-half
base map `H_a`, and its certified spectrum. The two BA permutations and
accumulators tile that base map to the logarithmic outer-block length. The
existing theorem remains explicitly attached to the current `a=24` extended
Golay certificate; another base constituent requires a replacement outer
certificate but does not change the route or RM2Sub inner.
