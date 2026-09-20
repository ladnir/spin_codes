# Riffle PacketMul-WrapMul-2Lap g=4: proof plan

## Current question

The parent candidate leaves the nonzero first-lap terminal state correlated
with the driven second lap. WrapMul should remove that dependence at the state
interface. The first task is to state this benefit in the correct setup
probability space and expose the remaining finite list bound.

## Goal 01: complete

1. Freeze the terminal-state field and setup chronology.
2. Prove conditional uniformity for every fixed nonzero terminal state.
3. Audit transfer of the terminal-zero and autonomous-prefix results.
4. Convert the remaining support-33 probability budget into an integer cap on
   bad wrapped states.
5. Define the first list-size lemma without assuming independence across outer
   words.

The exact uniform cap is 316,505 bad states per fixed driven word when the 26
support-33 words consume the entire remaining ledger budget.

## Goal 02: active large checkpoint

Goal 02 targets the complete support-33 nonzero-terminal row. It first extends
the boundary argument to every first occupied node \(r\ge7060\) through a
tail-code list bound. It then treats \(r\le7059\) through the low-weight
spectrum of the full retained-state response code.

The proof and refutation tracks run together. The proof track targets
\(A_{\le377{,}530}\le316{,}605\). The refutation track searches for a reachable
affine list above the active ledger cap and quantifies its setup probability.

## Proof direction after Goal 01

For a fixed driven word \(y\), count the nonzero states \(z\) for which

\[
\operatorname{wt}(F(F(y))\mathbin\oplus J(z))\le188{,}765.
\]

A sufficiently small uniform list closes the nonzero-terminal support-33 row.
Placement strata and autonomous-prefix bounds can reduce this full-output
list problem to smaller tail lists.

## Refutation direction

Search for driven offsets whose bad-state list is too large. A single offset
with a large list does not by itself refute the setup claim; the offset must
occur with enough probability under the inherited packet experiment.

## Scope

The full construction remains open. Results from the parent transfer only
after the Goal 01 audit records the required invariance.
