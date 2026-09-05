# Goal 04: boundary-prefix certificate

**Status:** Complete. The exact certificate closes 1,677 first-node strata and
moves the cutoff to node 20,976. See
`proof/GOAL_04_BOUNDARY_PREFIX_PROOF.md` and
`receipts/goal04_boundary_prefix_audit.json`.

## Objective

Resolve the first three nonsuffix first-node strata for authenticated
support-33 words. These strata have first occupied node

\[
r\in\{22{,}650,22{,}651,22{,}652\}.
\]

The goal first tests a deterministic strengthening of the autonomous-prefix
certificate. If that strengthening fails, the goal must expose the exact
conditional PacketMul law and produce a finite probability obligation.

## Probability space

Fix one authenticated support-33 outer word \(x\). Setup samples independent
multipliers

\[
A_i\gets\mathbb F_{16}^{\times}
\]

for its active packet coordinates. Setup independently samples the unchanged
uniform packet permutation \(\Pi\). Define

\[
y:=\Pi(D_Ax).
\]

Let \(L(y)\) denote the first-lap terminal state. Conditioning on
\(L(y)\ne0\) may correlate the active multiplier values. No proof step may
reuse their unconditional independence after this conditioning.

## Deterministic route

For zero-input nodes in the second lap, consecutive outputs satisfy

\[
o_{j+1}=U(o_j),
\qquad
U(o):=\operatorname{Acc}(P(o)).
\]

The first route is to certify a four-node lower bound strong enough to replace
the three-node bound near the suffix cutoff. A four-node lower bound of 34
would imply

\[
\operatorname{wt}(o_0,\ldots,o_{22{,}647})
\ge 5{,}662\cdot34
=192{,}508
>188{,}766.
\]

This would close every nonzero-terminal stratum with \(r\ge22{,}648\),
including the three target strata, without conditioning on PacketMul values.

## Work packages

1. Construct an exact four-node autonomous-window search.
2. Enumerate every nonzero output of weight at most eight and all four window
   alignments that contain it.
3. Either prove that every nonzero four-node window has weight at least 36 or
   record an exact counterexample.
4. Reconstruct the search independently with a distinct evaluation path.
5. Insert the resulting cutoff into the Goal 03 event ledger.
6. If the deterministic route fails, state the exact conditional probability
   reduction for the three boundary strata.

## Completion criteria

The goal is complete when one of the following outcomes is independently
verified.

- A deterministic certificate closes all three target strata.
- An authenticated counterexample refutes deterministic closure.
- A local counterexample refutes the proposed four-node lemma, and the report
  identifies the smaller remaining tail event.

The result must distinguish a local autonomous counterexample from an
authenticated construction counterexample.

## Excluded work

This goal does not close the entire support-33 complement. It does not prove
the full construction. It does not run a performance benchmark.

## Outcome

Primary and independent enumerations prove that every nonzero autonomous
four-node window has weight at least 36. The target three strata are closed.
The first open support-33 row now has nonzero terminal state and first occupied
node at most 20,975.
