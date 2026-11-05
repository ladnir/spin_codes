# Component counts and a three-support constraint

2026-10-30. The construction and full45%/40-bit goal are unchanged.
The [partial-state audit](QUARTET_TRANSITION_AUDIT.md) identified harmless
overcounting. This continuation tests a more selective geometric witness.
All statements concern the actual regional route and quartet partitions.
They do not certify whole-code distance.

## What the sampled optimal cuts contain

The existing min-cut diagnostic uses seed20261027 and ONE fixed inactive-row
set of size2100. These are exact fixed-geometry facts, not ensemble tails.
Write F_lambda(S)=lambda*A-C-(37/50)D for its support score.

At lambda=0.25 the optimum is3988.44. Its inclusion-minimal maximizer has
30324 positive interaction components. All but one have at most64 state
vertices. Only827 components have C>0; their total score is447.30,
A=13120,C=1238,D=2155. The largest component has1772 vertices, spans10
updates, and has score149.08. Thus most of the score comes from runs that
cover no inactive inputs, but the positive-coverage components still matter.

At lambda=0.26 the optimum is7832.88. The all-selected support already has
score7659.52. The inclusion-minimal complementary hole set has16266
vertices, uncovers4589 inactive inputs, and has252 births. It forms one
component spanning40 updates, with gain173.36 over the all-selected support.
This differs from the earlier inclusion-minimal SELECTED support because
the optimum can have ties. Both cuts have independently verified objective.

The two slopes are proof parameters, not changes to the encoder. The data
suggest sparse selected components below the inactive density2100/8192,
and sparse holes above it. One sample does not establish that behavior
uniformly over row sets or setups.

## A quota-aware certificate for zero-coverage runs

Fix a row set T of size2100 and arbitrary quartet partitions. A state
coordinate (t,j) is free if its own input is active and, unless t is final,
every input of its quartet at t+1 is active. Exactly the free coordinates
can be selected when C=0. Optimization then separates into512 vertical
sequences of maximal free runs.

Let lambda be0.25 or0.26, tau=0.74, and alpha=3lambda-tau. A nonterminal
free run of length l contributes max(0,lambda*l-tau). If W3,W4 count its
length-three and length-four windows, this contribution equals

    alpha*W3+(lambda-alpha)*W4.

A terminal run contributes lambda*l instead; the difference is at most tau.
Retain only windows starting in updates1 through4092. The omitted windows
and terminal corrections cost at most4lambda*512 in total. Hence

    max_{C(S)=0} F_lambda(S)
      <=4lambda*512 + sum_{t<=4092,j} Z_(t,j),
    Z_(t,j)=alpha*I3+(lambda-alpha)*I4.

The event Il requires4l+1 distinct input positions to be active: one own
input followed by l quartets. Also I4 implies I3. Each input position is
read by at most17 of the Z variables, for every fixed set of partitions.

Put a=1-2100/8192. Under independent active indicators of probability a,
the read17 Holder bound gives

    E[b^(sum Z)] <= [1-a^13+(a^13-a^17)b^(17alpha)
                       +a^17 b^(17lambda)]^(2095104/17).       (1)

For completeness, the read-k inequality follows from relative entropy.
For product measure P and an arbitrary Q, the chain rule and convexity give
sum_i D(Q_{E_i}||P_{E_i})<=k D(Q||P) whenever each coordinate belongs to
at most k sets E_i. Apply the variational formula for log E_P exp(sum f_i),
then bound each marginal term by log E_P exp(k f_i)/k. This proves (1)
without assuming that the Z variables are independent.

The same upper bound applies to the ACTUAL regional quotas. Expanding
b^Z gives a polynomial in the active indicators with nonnegative
coefficients, since0<=alpha<=lambda and b>=1. For any specified d positions
in one region, the probability that they are all active is

    binom(8192-d,2100)/binom(8192,2100) <= a^d.

Products of these inequalities apply across independent regions. Thus the
actual expectation is at most the independent-bit expectation, term by
term. There is no conditioning penalty and no independence claim within
a region. Finally union-bound over binom(8192,2100) row sets.

At lambda=13/50, B=25604, and log2(b)=5/8, directed arithmetic gives

    log2 Pr[exists T,S: |T|=2100, C_T(S)=0, F_lambda(S)>B]
      <= -141594448094614557/93458488360960 < -1515.

The helper uses rational upper roots and directed logs. At lambda=1/4,
B=18171 and log2(b)=9/16, the corresponding bound is positive1804.64,
so that calculation does not close.

This new geometric certificate does NOT add distance coverage by itself:
the existing field/activity argument already excludes C=0 in the hard
row class. Its role is to validate a selective counting method without
the enormous harmless-history sum. No new allowance is charged to the
full distance ledger until it is used in a complete decomposition.

## Complementary holes and the missing constraint

Let H_t be the complement of S_t and set H_0 to all512 coordinates, since
S_0 is empty. Let A_H count holes and B_H count hole births after update1.
Let U_H count uncovered inactive inputs. For t>1 an inactive input (t,j)
is uncovered precisely when j belongs to H_t and its entire preceding
quartet belongs to H_(t-1). At t=1 only membership in H_1 is required.
Then, with256 regions,

    F_lambda(S)=lambda*N-256r+J(H),
    J(H)=U_H-lambda*A_H-tau*B_H.                           (2)

This identity keeps both boundaries correct. The hole score is
supermodular. Its interaction components connect consecutive holes and
the five holes of every completed inactive-input condition. The score is
additive over these components. An inclusion-minimal maximizer has only
strictly positive components. Collections of components must still be
bounded jointly; their probabilities are not asserted independent.

There is a stronger local consequence of maximality. Suppose j is a new
hole at update t>1: j belongs to H_t but not H_(t-1). Its own input cannot
be uncovered, because the preceding quartet contains j, which was not a
hole. If this hole does not help uncover any input at update t+1, deleting
it loses no U_H. Deletion saves lambda if the hole persists at t+1, or
lambda+tau if it does not. Either change strictly improves J.

Therefore EVERY maximizing support satisfies the following gate:

    A new hole at (t,j) must lie in an all-hole quartet at update t+1,
    and that quartet must have an inactive input whose outgoing
    coordinate is also a hole.

There are no new terminal holes, since no following input could justify
them. This gate holds for every maximizer, without a tie-breaking
assumption. It can therefore be combined with the earlier filters for
the inclusion-minimal selected maximizer. It must not be imposed on
arbitrary codeword histories: it is a witness reduction for the geometric
maximum in the sufficient distance argument.

## An exact marked local kernel

The gate couples THREE consecutive selected supports. A sufficient state
is (u,b), where u=|S_t| and b=|S_(t-1) outside S_t|. The second set records
new holes and is disjoint from S_t. Initial state is(0,0); a maximizing
witness has final b=0, with its final u unrestricted.

Within one quartet of update t+1, record

    (a,e,v,d,l,c):
    incoming selected positions, marked new holes, outgoing selected
    positions, outgoing deaths, inactive labels, covered inactive labels.

For a fixed canonical incoming mask of size a, enumerate outgoing and
inactive masks. Its multiplicity is binom(4,a)binom(4-a,e). If e>0,
require a=0 and at least one inactive position outside the outgoing mask.
Outgoing deaths are determined by the two selected masks. The helper
`marked_local()` generates the resulting integer polynomial with265 terms.
An independent enumeration of20736 colored-mask cases verifies every term.

Raise this polynomial to128 quartet factors. Extract incoming degrees u,b
and outgoing degrees v,d, and divide by

    binom(512,u)binom(512-u,b).

This is the fixed-incoming-pair normalization, not a fresh-message law.
The inactive-label marker must remain through16 updates before extracting
the regional quota and dividing by binom(8192,r). Attach the score weights
and the earlier canonical constraints inside the same polynomial.
There are131841 possible (u,b) states. No full transfer has been run.

A simpler exact gate audit illustrates its potential strength. Ignoring
the inactive-child requirement, the probability that every marked new
hole lies in an all-hole quartet is

    [x^u z^b] ((1+x)^4+(1+z)^4-1)^128
      / [binom(512,u)binom(512-u,b)].

At u=256 and b=43,44,64 its directed log2 upper bounds are approximately
-98.996,-100.914,-137.263. At(u,b)=(128,32),(384,32) they are-38.119 and
-110.841. Exact integer coefficients verify these values. These factors
are NOT independent of the old kernel and must not be multiplied into it.
The marked polynomial is the proposed way to retain the dependence.

## Checks and next action

The new scalar helper checks8192 exact free-run identities, exact quota
monomial inequalities,45 exhaustive small colored partitions, and20736
marked local cases. The prior audit additionally checks65536 complementary
score identities. An exact16-state toy DP checks16 optimal triples,
including two nonvacuous internal births and the terminal condition.
The full-size slope0.25 replay verifies30631 birth-bearing quartet gates.
Both full-size cuts retain flow feasibility and independent objective checks.

All numerical work was serial and remote. C++ builds were temporary and
removed; no raw data was retained. Stored graph and audit arrays are bounded
below512MiB, with a768MiB process limit. Scalar tests have a30s watchdog.
The previous frozen proof helper remains unchanged.

Next test the marked kernel on the known failing partial-state family,
retaining the actual quota and inherited births jointly. First establish
whether it removes that family; do not start a large131841-state transfer
or multiply marginal gate probabilities. Full r=2100 closure and the
remaining row counts are still open.
