# Focused 9% proof attempt

## Objective and probability space

Fix message length (n=2^{20}), output length (N=2^{21}), and threshold
(d=\lfloor0.09N\rfloor=188743).  For each nonzero message, probability is
over setup: the outer-coordinate permutations, routing randomness, and
FieldCheckpoint multipliers.

The calculation uses the existing density envelope for the modeled
([256,128,38])-shaped spectrum.  The envelope replaces every regular active
outer word by 256 independent fair candidate bits at a factor smaller than
two.  It excludes the all-one outer word.

## Exact one-block result

For one active outer block, every transposed region contains one fair
candidate.  MultiBlockStripe restricts that candidate to one of (q) lane
slices.  The lane advances cyclically with the outer coordinate.

The exact two-state transfer gives the same 9% contribution for every tested
(q\in\{1,2,4,8,16,32\}).  The resulting margin is 55.950662585 bits.  Thus
the (q=32) balance does not weaken the layer that dominates the original
low-occupation sum.

## Proof-safe relaxation for reused shifts

For an occupation (a), let (L_k(z)) denote the exact transfer through one
lane slice containing (k) fair candidates.  The script forms an entrywise
upper bound over every composition

\[
 k_0+\cdots+k_{q-1}=a.
\]

It permits a new worst composition in every transposed region.  This
relaxation dominates the reused-shift construction even when region
dependencies are adversarial.

For (q=32), the relaxation retains 55.483806625 bits of aggregate margin
through (a=16).  It first fails the 40-bit target at (a=47).  At (a=64),
the relaxed contribution is (2^{2339.04}), so the relaxation cannot close
the proof.

The failure has a concrete shape.  The relaxation can place nearly all
candidates in one 256-bit epoch in every region without paying for repeated
lane coincidences.  The reused-shift construction can correlate those
coincidences across coordinates, so this loss cannot be repaired by a local
constant-factor estimate.

## Fresh-shift result

MultiBlockStripeFresh samples the group shift independently for every outer
coordinate.  Suppose first that the (a) active blocks belong to distinct
groups.  Their lane choices are then independent in every region.  The exact
averaged region transfer is

\[
 R_a(z)=\frac{a!}{q^a}[u^a]
 \left(\sum_{k=0}^{a}\frac{L_k(z)}{k!}u^k\right)^q.
\]

The calculation raises (R_a(z)) to the 256th power and optimizes the
Chernoff tilt.  For (q=32), every occupation (1\le a\le64) closes at 9%.
The aggregate margin is 55.950662585 bits.  The one-block layer remains
dominant.

This result is exact for active blocks in distinct groups.  It is not a proof
for arbitrary active-block sets.  Blocks in one group share a shift and occupy
distinct lanes.

## Exact shared-group law

Fix active-lane patterns (A_1,\ldots,A_m\subseteq\mathbb Z_q) in (m)
active groups.  In one transposed region, the construction samples independent
shifts (\delta_1,\ldots,\delta_m\in\mathbb Z_q).  Define

\[
 n_r(\boldsymbol\delta)
 :=\sum_{u=1}^{m}\mathbf1[r\in A_u+\delta_u].
\]

If (L_k(z)) is the exact transfer for a lane containing (k) fair
candidates, then the exact region transfer is

\[
 R_{A_1,\ldots,A_m}(z)
 :=q^{-m}\sum_{\boldsymbol\delta\in\mathbb Z_q^m}
   \prod_{r=0}^{q-1}L_{n_r(\boldsymbol\delta)}(z).
\]

Fresh shifts make the 256 region transfers independent.  The complete inner
moment is therefore

\[
 e_0^T R_{A_1,\ldots,A_m}(z)^{256}\mathbf1.
\]

The shared-group analysis script implements this formula directly.  Its cost
is exponential in the number (m) of active groups, so the implementation is
an exact small-profile oracle rather than a complete enumerator.

## Refutation of iid-lane domination

The earlier calculation proposed replacing the exact shared-group transfer by
the iid-lane transfer.  That pointwise domination is false.

Consider one active group with (A_1=\{0,1\}).  At the rational Chernoff
variable (z=63/64), exact rational arithmetic gives

\[
 \left(R_{\{0,1\}}(z)-R^{\mathrm{iid}}_2(z)\right)_{0,1}
 =0.00085191996941532178406595353447313066185435567721623\ldots>0.
\]

Thus (R^{\mathrm{iid}}_2(z)) does not dominate every shared-group matrix
entrywise.  The verifier derives both matrices from exact fractions.  The
positive difference has a 49,409-bit numerator and a 49,419-bit denominator.

The floating-point Chernoff optimization measures the consequence at the 9%
threshold.  The adjacent pair is 0.054577936 bits worse than the iid-lane
model.  After charging only its actual cyclic orbit and group placement, it
still retains 138.014709 bits of modeled margin.  This profile therefore
refutes the compression lemma but does not refute 9% distance.

## Status

The focused attempt supports (q=32) and the balanced 256-by-256 permutation
sizes only as a performance candidate.  Reusing one shift per group leaves a
proof obstruction at moderate occupation.  Fresh per-coordinate shifts close
the distinct-group model, but the exact adjacent-pair calculation refutes the
proposed reduction of shared groups to that model.

The result is not a complete 9% certificate.  Any continuation must retain
group-pattern types or prove a different compression inequality; entrywise
iid domination is unavailable.  A complete certificate would also inherit
the existing obligations for the all-one outer word, the modeled outer
spectrum, the higher occupations, and outward-rounded arithmetic.

The receipts are `receipts/focused_delta09.json`,
`receipts/q32_low64_delta09.json`, and
`receipts/q32_independent_lane_low64_delta09.json`.  The shared-group receipts
are `receipts/shared_group_small_profiles_delta09.json` and
`receipts/shared_group_adjacent_pair_exact_counterexample.json`.
