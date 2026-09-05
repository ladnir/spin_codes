# Controlled-writing review

This checklist records the section-by-section review against Controlled
Writing for Cryptography. It is an editorial artifact, not part of the paper.

## Resolved audit findings

- CWC-03 and CWC-09: quantified the distance target in the wrapped-EC
  asymptotic theorem before it appears in the hypothesis.
- CWC-05: preserved the paper-wide meanings of $p$ (field order) and $R$
  (code rate). Local support fractions now use $\theta$, and central support
  endpoints use $r_{\mathrm{mid}}$.
- CWC-08 and CWC-09: made the expander and convolution sampling experiments
  explicit in the finite wrapped-EC, regular-expander, and nonwrapping
  statements.
- CWC-03 and CWC-06: exposed the parameter of the degree-three conditional
  random variables as $Z_s^{(3)}(\theta)$ and separated its parity index from
  block endpoints.
- CWC-05: normalized the field-size theorem to the global field-order symbol
  and the established notation $\delta_{p,\mathrm{GV}}(R)$.
- CWC-05: separated the failure exponent $b$ from the field-order bit length
  $s$ in the adjacent certificate and discussion sections.
- CWC-15: the field-size discussion distinguishes interval-certified rows
  from frozen-marker diagnostics and identifies the source of the observed
  decline as a moving target rather than weaker mixing at fixed distance.

## Abstract and introduction

- Leads with the construction, the proof gap, and the concrete result.
- Separates the binary case from genuine finite-field mixing near 128-bit
  field orders before presenting either parameter family.
- Defines EA, EC, PCG, the enumerator idea, relative distance, and the role of
  the GV benchmark before relying on them.
- Separates changes to the analysis from changes to the construction.
- States the probability observer and excludes decoding and protocol-security
  claims.
- Compares the proof method with Block--Accumulate and recent large-field RAA
  analyses. The comparison defines RAA before use, identifies the shared
  projective accounting principle, and states why EA and EC need different
  structural traces.

## Code ensembles and distance events

- Defines messages, generators, support, Hamming weight, minimum distance,
  relative distance, rate, Hamming balls, and the GV benchmark.
- Gives the accumulator and convolution as explicit interfaces before using
  their properties.
- States all randomness and independence assumptions locally.
- Defines left and right vertices and both meanings of regularity; it also
  warns that "expander" does not assert a graph-expansion theorem.

## Enumerator method

- Defines the enumerator as an expectation and identifies the probability
  space in the first-moment step.
- States the symmetry required for composition. The finite-field hypothesis is
  monomial invariance; the binary specialization needs only permutations.
- Introduces state, input marker, and output marker before coefficient
  extraction. The positive-marker relaxation follows the exact identity.
- Requires an explicit partition of all nonzero support sizes.  The certified
  support-partition lemma now states the shared endpoint calculation and the
  verifier's exact responsibilities.

## Binary EA

- Keeps the construction fixed and says so before presenting the proof.
- Defines the activation probability and every state of the two-state matrix.
- Places the exact first moment before spectral and asymptotic relaxations.
- Separates sparse, intermediate, and dense supports and distinguishes
  asymptotic comparisons from finite certificates.

## Binary EC

- Defines the wrapping convention and the reduced state before the matrix.
- Explains every exceptional transition, including the fixed oldest tap and
  the all-zero state.
- States which prior comparison changes the analyzed ensemble and which holds
  the ensemble fixed.
- Defers only the lengthy endpoint calculation to a named appendix.

## Regular expander variants

- Presents regular incidence as a construction change, not a proof trick.
- Defines a unit vector and the regional shell recurrence before its use.
- Explains Poisson conditioning and point-mass terminology in place.
- Retains recursive state across region boundaries.
- Defines the binary two-sided probability space before its exact transfer.
- Presents the certified degree-$10/5$ and degree-$6/3$ statements as a
  degree--memory tradeoff, together with the intermediate $14/7$ and $18/9$
  profiles.  It explains why degree $8/4$ is structurally invalid.
- Defines scalar extension and limits its use. The proof identifies the binary
  basis components and thereby explains why scalar extension does not provide
  the cross-component mixing required by field-valued Silent VOLE.
- Introduces the common odd-degree conditional enumerator before using its
  real-rooted factorization in the central-block lemma.
- Places the right-degree parity proposition before all usable binary
  parameters.
- Grounds the local-limit variables before the point-mass claim and states
  the monotonic endpoint facts used by the interval verifier.
- Introduces the degree-three conditional Bernoulli variables before the new
  central-block lemma.  Its proof identifies the conditioning event, the
  observer, the variance comparison, and the complementary-support mapping.

## Finite fields and field-size scaling

- Motivates random edge labels by the exact symmetry they provide.
- Defines projective message lines before counting them.
- Separates the theorem for each fixed field from the finite
  cryptographic-size certificates. The argument is stated for prime-power
  field orders, including characteristic two.
- Defines occupied coordinates and fresh-equation markers before the
  constraint matrices. It explains why full-field feedback sampling is needed.
- States which factors depend on the field order and distinguishes the
  fixed-field asymptotic theorem from finite results near 128 bits.

## Certified finite evaluation

- Defines the failure exponent and the outward-rounded ball representation.
- States the trust boundary: the optimizer and cache are untrusted, while the
  verifier checks domain, coverage, arithmetic, and the final endpoint.
- Identifies every parameter in the headline theorems and states which
  support ranges each verifier covers.
- Gives both $\mathbb{F}_{2^{128}}$ results their own labeled finite-field
  certificates instead of deriving an application claim by scalar extension.
- The binary and $\mathbb{F}_{2^{128}}$ theorem parameters match their
  checked-in certificate files and complete support partitions.
- Separates exact, outer, central, and full-support contributions for the
  binary two-sided regular EC result.
- Grounds the frozen marker file and the independent structural checker before
  claiming that verification does not run the optimizer.
- Gives a branch-to-equation correspondence table.  Every numerical verifier
  validates its schema, marker domains, parameter identities, and support
  partition.  The independent standard-library checker is claimed only for
  the two binary two-sided formats that implement it.

## Discussion and appendices

- Discussion separates distance, decoding, deterministic construction, and
  application security.
- The implementation paragraph assigns the binary profile to Silent OT and
  the field-valued profile to Silent VOLE. It distinguishes the independently
  labeled proof ensemble from the structured streaming heuristic.
- The field-size discussion defines fixed-cutoff and field-relative-cutoff
  comparisons before presenting them.  It separates interval certificates
  from the remaining frozen-marker diagnostics.
- The additional-ensemble appendix defines fixed-column, fixed-row,
  Krawtchouk, conditioning, and dense Fourier terms before conclusions.
- The nonwrapping appendix distinguishes the sampled-oldest-tap construction
  from wrapping and gives its own exact matrix.
- The analytic appendix preserves the order state classes, effective
  generator, optimization, uniform prefactor, then endpoint conclusion.
- The artifact appendix maps every finite result to the checked-in manifest.
  The manifest pins certificate and source hashes and supplies one serial
  runner.

## Final scope checks

- The wrapped-EC appendix proves a uniform quadratic perturbation remainder
  through a two-class Schur-complement argument.
- The construction section states sampling and storage costs.  The discussion
  gives the exact additive rule for composing code-sampling failure with a
  protocol theorem, without claiming a protocol-specific reduction.
- The local field-size profile is interval-certified at orders $2^{127}$,
  $2^{128}$, $2^{129}$, and $2^{136}$.  Other displayed orders are labeled as
  diagnostics in both the paper and manifest.
