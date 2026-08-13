# Second independent session: close the `g=4` grouped-permutation proof

You are joining an active proof/implementation project for an independent
second opinion.  Please work from the repository, inspect the actual artifacts
and code, challenge our current interpretation, and write your conclusions to

```text
CLAUDE_G4_CLOSURE_SESSION_REPORT.md
```

This is a research audit, not a request to rubber-stamp the current plan.  Label
each important conclusion as **PROVED**, **DIAGNOSTIC**, **PLAUSIBLE**, or
**SPECULATIVE**, and identify precisely what artifact, computation, or lemma
would promote it to the next category.

## Objective

We want an end-to-end first-moment distance certificate for the frozen Riffle
construction with packet-permutation group size `g=4`, target relative distance
`0.09`, and total failure probability at most `2^-40`.  This is the methodology
ladder toward the deployment target `g=8`.

The immediate question is narrower:

> What is the most principled and economical way to eliminate the present
> `g=4` bottleneck, and does the evidence still support the belief that the
> end-to-end `g=4` certificate will close?

We especially want you to distinguish:

1. a genuinely insufficient outer/inner inequality at a concrete profile;
2. optimization failure within an adequate witness family;
3. BSP/cell geometry and lattice-accounting slack; and
4. a soundness defect in the proof architecture.

## What changed since your previous report

First read `CLAUDE_SECOND_OPINION_REPORT.md`, especially Sections 6--10.  Its
main recommendation was an exact per-source-cell dominance BSP, avoiding global
retriangulation.  We implemented that recommendation.  It was highly effective:

- the top-128 BSP frontier now closes diagnostically;
- the top-512 frontier was driven down until every processed cell was below the
  then-unprocessed frontier;
- the source-cell breadth was expanded to the top 1024 cells;
- every split and leaf is represented exactly, but the current runs are still
  discovery/diagnostic runs rather than the final outward theorem artifact.

The current top-1024 replay has:

```text
largest processed source-cell contribution    about +17712.39 bits
largest unprocessed source-cell ceiling        about +12078.98 bits
current leading cell                           s1fr032620
```

The leading cell has one rounded residual integer profile:

```text
[444306, 21517, 37620, 17842, 3003]
```

At that profile, the sequence of dedicated attacks was approximately:

```text
ordinary coordinate tuning miss                  7207.07 bits
cross-profile Powell miss                          746.08 bits
two local continuations miss                       422.45 bits
same frozen inner + exact-graph outer miss          300.25 bits
inner reoptimization under that outer miss          346.62 bits
```

The last regression matters: blindly increasing Powell depth is no longer our
preferred move.  The remaining miss may be an outer-family issue, an
outer/inner tradeoff that our staged optimization does not expose, or a cell
mixture/refinement issue rather than a true point obstruction.

A previous leading cell, `s1fr033687`, had two rounded residual profiles.  Both
eventually closed pointwise, by about `2046.82` and `38.89` bits respectively,
and its cell contribution fell from roughly `+18584.99` to `+13144.95` bits.
This is evidence that small pointwise margins can still be operationally useful,
but it is not evidence that the current profile must close in the same family.

## Current proof architecture

The domain is split into exact positive-support strata.  Each stratum has an
exact integer hull and a mechanically verified rational triangulation.  The
global profile union uses exact-support disjointness and cell-local lattice-count
upper bounds.  Within a source cell, an exact rational BSP may subdivide the
cell without requiring conformity with neighboring source cells.  A leaf is
certified by one fixed witness or a fixed convex mixture, outward-evaluated at
all exact leaf vertices.  Convexity then covers the whole leaf.

All witness bounds share the multinomial normalization term `-log2 Q(a)`; their
remaining profile dependence is affine.  This makes pairwise dominance
boundaries affine and supports exact BSP discovery.  Mixture selection is a
small minimax LP; column generation prices the complete finite eligible atlas,
so a small candidate cap is not being mistaken for a proof.

The graph/puncture outer upgrade conditions on the exact 24-dimensional graph
subcode spectrum and replaces the old adversarial 128-hole tax.  A previous
hard profile was closed outward by this branch with only about `0.0596` bits of
local margin.  The current stubborn profile benefits by about `122.2` bits from
the exact-graph upgrade but remains open.

The diagnostic tuner is accelerated by native C++ kernels for the 65-state
transport recurrence and point-cap construction.  Python remains the reference
and the outward verifier is independent.  Runtime is not the immediate blocker.

## Required audit questions

### 1. Re-audit the BSP architecture against its implementation

Inspect the producer and verifier, not just the prose.  Confirm or refute:

- a per-source-cell BSP may be internally nonconforming without invalidating
  the original exact mesh ledger;
- the split-tree representation really covers each parent exactly;
- leaf vertex enumeration is exhaustive over exact rationals;
- fixed witnesses/mixtures are held unchanged over an entire leaf;
- support eligibility and zero fugacities cannot leak across a stratum;
- shared BSP boundaries and source-cell boundaries are only overcounted, never
  omitted;
- the cell-local count and final log-sum-exp cannot accidentally count a
  discovery sample in place of a certified region.

Name any missing verifier invariant.  Do not infer soundness from agreement
between producer and verifier if both share the same flawed helper.

### 2. Perform a forensic decomposition of `s1fr032620`

Use the stored artifacts to separate, in bits, as far as possible:

- exact outer branch value and its individual relaxations;
- graph/puncture correction;
- inner 65-state bound and the point-cap/representative-fugacity relaxation;
- the best frozen combined witness value at the rounded integer profile;
- the best cell/leaf mixture value at every relevant exact vertex;
- lattice-count charge for the source cell or leaves;
- any rounding gap between a failing rational BSP vertex and its rounded
  integer tuning profile.

We need to know whether “300 bits short” is the actual local proof-content gap,
merely the rounded-profile diagnostic, or neither.

### 3. Audit the optimization parameterization

Check whether the present tuning actually performs the right joint
optimization.  In particular:

- Are the outer charge variables, exact-graph branch parameters, inner
  fugacities, pole, and any representative-fugacity choices optimized jointly
  where joint optimization is mathematically legitimate?
- Does the exact-graph upgrade freeze an outer tilt selected for a different
  relaxation, leaving an obvious saddle-point mismatch?
- Are there gauge freedoms, active constraints, bad coordinate transforms, or
  nonsmooth branch selections that make Powell/coordinate continuation report a
  false plateau?
- Can a fixed convex mixture of existing outer/combined witnesses close the
  relevant leaf even if no individual witness closes the rounded point?
- Would a small deterministic global method, interval branch-and-bound, or
  multi-start in a better coordinate system be justified for this one profile?

Please derive first-order/KKT or dual diagnostics where useful.  Prefer a
falsifiable numerical test over generic advice to “try more seeds.”

### 4. Challenge and extend the outer families

Inspect all currently available outer bounds, including total-weight,
linear-Brascamp--Lieb, total-spectrum, and exact-graph-conditioned variants.
Determine whether the current profile suggests a missing interpolation or a
specific lost correlation.  Candidate questions include:

- Can the exact graph conditioning be combined with a shaped/full spectrum
  rather than only the current total-spectrum parameterization?
- Can two valid outer bounds be combined before Cauchy extraction more sharply
  than choosing their pointwise minimum afterward?
- Is there a support/profile-conditioned outer enumerator that is still cheap
  enough to certify and likely to recover at least 300 bits?
- Does the exact BCH spectrum already contain the needed information but the
  current charge family fails to expose it?
- Is the target profile even realizable by the frozen outer code at the claimed
  multiplicity, and is there a rigorous inexpensive exclusion test worth doing?

For every proposal, state the new lemma required and whether it is a legitimate
upper bound rather than a heuristic.

### 5. Decide the next bounded experiment

Recommend one primary experiment and at most two fallbacks.  For each give:

- exact inputs/artifacts;
- code path to change or invoke;
- expected runtime and parallelism;
- what numerical result would count as success;
- what a negative result would teach us;
- estimated likely bit gain; and
- whether the output is diagnostic or can feed the outward certificate.

The primary experiment should be deliberately bounded.  Do not launch a long
Peach run as part of this audit.  Never run two benchmarks simultaneously.

### 6. Give a stop/go assessment for `g=4` and implications for `g=8`

Answer plainly:

- Is there positive evidence that `g=4` is true, beyond the fact that our search
  has not found a counterexample?
- How many qualitatively distinct proof obstacles remain after
  `s1fr032620`?
- At what observation should we stop tuning this inequality family and redesign
  it?
- Which parts of the `g=4` method scale to `g=8`, and which are likely to
  explode even if `g=4` closes?

Give a ranked plan from the current state to a genuine outward end-to-end
certificate, with explicit gates.

## Required reading and artifacts

Read in this order:

1. `PROOF_STATUS.md`, especially the complete `g=4` ladder section.
2. `CLAUDE_SECOND_OPINION_REPORT.md`, especially Sections 6--10.
3. `explorations/g2_triangular_cover_math.md` and
   `explorations/g4_delaunay_cover_math.md` for the convex-mixture and
   cell-ledger lemmas.
4. `scripts/probe_packet_group_g4_cell_bsp.py`
5. `scripts/probe_packet_group_g4_cell_bsp_batch.py`
6. `scripts/certify_packet_group_g4_cell_bsp.py`
7. `scripts/probe_packet_group_joint_inner_opt.py`
8. `scripts/upgrade_packet_group_exact_graph_outer.py`
9. `scripts/probe_packet_group_exact_graph_puncture_outer.py`
10. The native diagnostic kernels only as needed:
    `scripts/packet_group_shared_drive_native.cpp`,
    `scripts/packet_group_native.py`, and
    `scripts/packet_group_drive_stratified.py`.

Inspect at least these artifacts:

```text
out/g4_cell_bsp_top1024_exact_graph_round1.json
out/g4_cell_bsp_top1024_exact_graph_round1_failures.json
out/g4_cell_bsp_s1fr033687_last_continue3_native.json
out/g4_cell_bsp_s1fr032620_profile.json
out/g4_cell_bsp_s1fr032620_tuned_native.json
out/g4_cell_bsp_s1fr032620_cross_powell_native.json
out/g4_cell_bsp_s1fr032620_continue2_native.json
out/g4_cell_bsp_s1fr032620_continue2_exact_graph.json
```

Some later replay artifacts may exist only on Peach or may not yet have been
copied back.  Treat `PROOF_STATUS.md` as the index of those results, but clearly
label any conclusion that could not be checked against the underlying JSON.

## Deliverable structure

Please write `CLAUDE_G4_CLOSURE_SESSION_REPORT.md` with:

1. headline verdicts;
2. soundness findings, with file/line references;
3. a quantitative decomposition of the stubborn profile/cell;
4. optimization audit;
5. ranked proof refinements;
6. the single best bounded experiment;
7. stop/go criteria and a path to the outward `g=4` certificate;
8. implications for `g=8`; and
9. a short list of missing artifacts or facts that prevented a firmer answer.

Keep theorem claims separate from diagnostic evidence.  If you find a fatal
flaw, lead with it.  If not, tell us exactly why the present approach remains
credible and what result would most efficiently increase or decrease that
confidence.
