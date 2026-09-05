# Paper plan: exact enumerators for EA and EC codes

## Working title

**Exact Enumerator Bounds for Expand--Accumulate and Expand--Convolute Codes**

The title names the two code families and the principal proof technique.  The
abstract and contribution list will advertise regular expanders, prime fields,
and certified near-GV parameters.

## Central claim

The published distance analyses for Expand--Accumulate (EA) and
Expand--Convolute (EC) lose substantial parameters when they bound the state
process by concentration.  Exact input--output enumerators remove this loss.
They also identify low-support messages as the remaining obstruction.  Regular
expander ensembles reduce that obstruction, while random edge labels and
projective counting extend the method to prime fields.

The paper must distinguish three kinds of results:

1. **Analysis-only improvements.**  The binary Bernoulli EA and EC ensembles
   remain unchanged; only their proofs improve.
2. **Construction improvements.**  Region-stratified regular expanders replace
   independent Bernoulli incidence.
3. **Field generalizations.**  Prime-field constructions add random nonzero
   edge labels.  Prime-field EC also samples feedback coefficients from the
   full field.

## Contributions

### C1. Exact enumerator framework

Give an exact composition rule for a sparse random expander followed by a
recursive rate-one linear map.  Instantiate the rule with:

- the two-state accumulator transfer matrix;
- the reduced `(m+1)`-state convolution transfer matrix;
- Bernoulli, fixed-row, region-stratified left-regular, and two-sided regular
  expander ensembles.

The paper still uses a first-moment bound for the minimum-distance event.  It
replaces concentration bounds for the internal state process, not Markov's
inequality in the final first-moment step.

### C2. Sharper proofs for the original binary EA and EC ensembles

For binary EA, prove that the dense-message condition equals the binary GV
condition.  Isolate the additional expander-density condition imposed by
sparse messages.  Compare both conditions directly with the original EA
theorem and give finite certificates for the same ensemble.

For binary wrapped EC, derive the exact reduced-state enumerator and prove the
corresponding asymptotic distance theorem.  Compare against the published EC
conditions while holding the construction and boundary condition fixed.  Give
finite same-ensemble comparisons.

The introduction should use two concrete comparisons:

- EA at relative distance `0.05`: the density threshold decreases from
  `C > 4.938` to `C > 1.773`.
- wrapped EC at relative distance `0.05` and memory `21`: the density threshold
  decreases from `C > 10.4344` to `C > 1.111623`.

### C3. Regular expander ensembles

Define the region-stratified left-regular ensemble and derive its exact shell
law.  Explain why regional independence permits composition with the recursive
map without adding an independent interleaver.

At rate one half, report:

- binary EA: exact left degree `64`, versus Bernoulli expected degree `75`, at
  essentially the binary GV distance and a failure target below `2^-20`;
- binary wrapped EC: exact left degree `28`, versus Bernoulli expected degree
  `42.2127`, with memory `9` and a certified distance within `14.276` output
  symbols of the binary GV shell.

Present regularity as a change to the ensemble, not as a tighter proof of the
Bernoulli ensemble.  Compare equal rates, lengths, distance cutoffs, and total
edge counts whenever possible.

### C4. Prime-field extension

Attach independent labels from `F_p^*` to expander edges.  This makes the law
of a fixed nonzero message depend on its support rather than its values.  Count
projective message lines instead of all nonzero scalar multiples.

For prime-field EA:

- derive the exact two-state zero/nonzero accumulator kernel;
- give Bernoulli and regular first-moment bounds;
- use support grouping and constraint traces to control the large projective
  union;
- state the asymptotic fixed-field theorem separately from large-field finite
  certificates.

For prime-field EC:

- sample each feedback vector from the full field `F_p^m`;
- explain why restricting every tap to `F_p^*` destroys the reduced-state law;
- combine the convolution constraint trace with two-sided regular incidence;
- separate structural low-equation traces from field-suppressed traces;
- use the singleton-free-region lemma for the smallest supports.

### C5. Reproducible finite certificates and the degree--memory frontier

Give a short verifier model in the paper.  Floating-point optimization chooses
candidate markers, but the trusted verifier treats those markers as fixed
decimals and recomputes every bound with outward-rounded Arb arithmetic.

The flagship theorem should state:

- field: `p = 2^127 - 1`;
- rate: `1/2`;
- length: `n = 2,097,150`;
- left and right degrees: `22` and `11`;
- convolution memory: `12`;
- cutoff: the exact floored p-ary GV cutoff;
- code-sampling failure probability: below `2^-27.2951`.

Also report the certified frontier `(m,d_L) = (3,28), (4,26), (6,24),
(12,22)`.  Call `-log_2 U` the **failure exponent**, not security bits.

## Proposed manuscript structure

### 1. Introduction

Start with the tension: EA and EC are sparse, fast linear-code ensembles, but
their published distance guarantees are much weaker than their observed or
intended parameters.

Then give one concrete EA comparison and one concrete EC comparison.  State
the enumerator idea in one paragraph.  List Contributions C1--C5 and give a
single headline-results table.  End with scope and organization.

The introduction must distinguish:

- the original binary ensembles;
- the regular variants;
- the labeled prime-field variants;
- asymptotic theorems versus finite machine-checked certificates.

### 2. Code ensembles and distance events

Define only the objects needed throughout the paper:

- Hamming weight, minimum distance, rate, and the p-ary GV distance;
- the sparse expander matrix `B`;
- the accumulator `A` and convolution `C`;
- the generated code `G = BA` or `G = BC`;
- the probability space for each ensemble;
- the bad event `d_min(G) <= L`.

Give a one-page construction table.  Each row should specify the field,
expander distribution, edge-label distribution, recursive map, and random
choices.  Defer detailed enumerators until they are motivated.

### 3. Enumerator method

Introduce the input--output weight enumerator through the minimum-distance
question.  Then prove:

1. the first-moment reduction from bad codewords to low-weight enumerator
   mass;
2. exchangeable serial composition;
3. transfer-matrix coefficient extraction;
4. Chernoff evaluation with a fixed positive marker;
5. safe aggregation over contiguous message-weight bands.

This section should contain the common proof skeleton used by every later
result.  It should not yet specialize to all expander ensembles.

### 4. Binary Expand--Accumulate

Define the Bernoulli expander activation probability and the exact two-state
accumulator kernel.  Derive the exact first moment before giving any spectral
bound.

Then present, in order:

1. the spectral tail bound;
2. the sparse and dense endpoint exponents;
3. the asymptotic linear-distance theorem;
4. the comparison with the original EA theorem;
5. finite same-ensemble certificates.

End with the conclusion that dense messages attain the GV exponent and sparse
messages determine the required expander density.

### 5. Binary Expand--Convolute

Use wrapped EC as the principal construction because it gives the strongest
and cleanest statement.  Define its boundary condition explicitly.

Present:

1. the exact `(m+1)`-state transfer matrix;
2. the exact EC first moment;
3. sparse and dense asymptotics;
4. the asymptotic distance theorem;
5. a same-ensemble comparison with the published EC analysis;
6. finite rate-one-half and near-GV certificates.

The nonwrapping transfer matrix belongs in an appendix unless a main theorem
requires it.  Do not include the legacy path-counting derivations in the main
paper.

### 6. Regular expander variants

Motivate regularity with the obstruction identified in Sections 4 and 5:
small message supports can place too few useful edges or create too many
cancellations.

Define the region-stratified left-regular ensemble, prove its shell law, and
derive its point-mass and Poisson-comparison bounds.  Apply the same machinery
first to EA and then to wrapped EC.  Finish with a table comparing Bernoulli
and regular profiles at fixed rate and distance.

Do not claim that two-sided regularity improves all schemes.  Record that
extension as future work until it has a proof and an equal-cost comparison.

### 7. Prime-field expander enumerators

Begin with the problem created by field values: support no longer determines
whether an unlabeled coordinate sum vanishes.  Introduce random nonzero edge
labels as the remedy.

Derive the labeled Bernoulli and regular shell laws.  Introduce projective
message lines only when the minimum-distance union requires them.  Then give
support grouping, projective-surplus cancellation, and the shared
constraint-trace method.

This section supplies field machinery but does not yet introduce random
convolution.

### 8. Prime-field EA and EC

First instantiate Section 7 with the prime-field accumulator.  State the
fixed-prime asymptotic theorem and the regular finite result.

Then define the prime-field random convolution.  Explain the full-field tap
distribution before stating its reduced-state kernel.  Derive the two-sided
regular EC trace bound, the singleton-free-region lemma, and the complete
finite certificate.

End with the degree--memory frontier and an interpretation: convolution memory
trades against expander degree in the sparse-support regime, with diminishing
returns beyond the middle of the observed frontier.

### 9. Certified evaluation

Describe the certificate format, verifier inputs, interval coverage checks,
outward rounding, and trusted base.  Give reproducibility commands and point to
the machine-readable artifacts.

Keep optimization heuristics separate from soundness.  A verifier theorem
should state that any accepted certificate implies its displayed probability
bound.

### 10. Discussion and open questions

State the main limitations once:

- the results concern minimum distance and do not supply a decoding theorem;
- most guarantees hold with high probability over a sampled code ensemble;
- the large-field flagship is a finite theorem for one specified prime and
  length, not an asymptotic theorem with a growing field;
- the current comparison does not prove that two-sided regularity helps every
  EA or EC ensemble;
- encoder-cost comparisons must include regional sampling and data layout, not
  only edge count and convolution memory.

Discuss two-sided regularity, smaller fields, lower degree, and deterministic
constructions as future work.

## Appendices

### Appendix A. Additional expander ensembles

Place fixed-column and exact-row formulas here unless a headline theorem uses
them.  Include the positive recurrences needed by verifiers.

### Appendix B. Nonwrapping convolution

Give the exact nonwrapping transfer matrix and any comparison needed for the
historical EC construction.

### Appendix C. Deferred analytic proofs

Move long endpoint calculations, uniformity estimates, convex domination, and
technical coefficient bounds out of the main narrative.

### Appendix D. Certificate specifications

List every certified parameter set, support partition, precision, certificate
file, and verifier command.  Avoid embedding raw optimizer transcripts.

### Appendix E. Artifact manifest

Map each theorem and table row to its script, certificate, result JSON, and
expected digest or summary value.

The aggressive-construction attack notes and legacy combinatorial derivations
remain research notes.  They are not part of the planned appendices.

## Headline-results table for the introduction

The introduction should contain one compact table with columns:

| ensemble | field | rate | expander | recursive map | distance target | failure bound | comparison |
|---|---:|---:|---|---|---|---:|---|
| original EA | 2 | 1/5 | Bernoulli | accumulator | 0.05 | `< 2^-20` | same ensemble, tighter proof |
| original wrapped EC | 2 | 1/5 | Bernoulli | memory 5 | 0.05 | `< 2^-32.08` | about 10 extra failure bits |
| regular EA | 2 | 1/2 | left degree 64 | accumulator | 99.988% of GV | `< 2^-20.08` | 14.67% fewer edges |
| regular wrapped EC | 2 | 1/2 | left degree 28 | memory 9 | 99.994% of GV | `< 2^-22.22` | 33.67% fewer edges |
| regular EC | `2^127-1` | 1/2 | degrees 22/11 | memory 12 | floored p-GV | `< 2^-27.295` | flagship field result |

Verify every rounded value against its certificate before copying this table
into the manuscript.

## Notation policy

The new manuscript must not inherit the old notation conflicts.

- `p`: prime field order only.
- `rho`: Bernoulli edge density.
- `R = k/n`: code rate.
- `d_L, d_R`: left and right expander degrees.
- `m`: convolution memory.
- `a_r`: activation probability for a message support of size `r`.
- `delta`: target relative distance.
- `L = floor(delta n)`: integer distance cutoff.
- `U`: certified upper bound on code-sampling failure.
- `b = -log_2 U`: failure exponent.

Use **left-regular expander** when only left degrees are fixed and
**two-sided regular expander** when both degrees are fixed.  Use **regular
expander** only when the distinction is immaterial or already fixed by the
section.

## Drafting order

Write the paper in this order:

1. construction table and notation ledger;
2. theorem and lemma dependency graph;
3. Sections 2 and 3;
4. binary EA and EC sections;
5. regular section;
6. prime-field sections;
7. certificate section and artifact manifest;
8. introduction and abstract;
9. discussion and appendices.

The next concrete task is to build the construction table and theorem
dependency graph.  Those two artifacts will determine which old lemmas are
reused, merged, renamed, or discarded before new TeX prose begins.
