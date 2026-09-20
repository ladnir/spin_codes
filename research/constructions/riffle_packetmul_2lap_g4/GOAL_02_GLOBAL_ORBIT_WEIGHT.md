# Goal 02: prove or refute the global orbit-weight lemma

**Status:** `COMPLETE_PROVED` (2026-08-21)

The proof is recorded in `proof/GOAL_02_GLOBAL_ORBIT_WEIGHT_PROOF.md`.
Primary and independent exhaustive certificates prove that every consecutive
four-state window has nibble weight at least 24. Since four divides 32,772,
this proves the target global orbit bound and discharges the Goal 01
mixed-character lemma.

## Objective

Prove or refute the remaining mixed-character lemma from Goal 01.

Let \(T\) be the 64-bit zero-input state map from the frozen construction.
Let

\[
S:=T^{\mathsf T},
\qquad N:=32{,}772.
\]

Write \(\operatorname{wt}_4(x)\) for the number of nonzero four-bit nibbles
in \(x\in\mathbb F_2^{64}\). The target statement is

\[
\sum_{t=1}^{N}\operatorname{wt}_4(S^t\chi)\ge6N
\qquad(\chi\ne0).
\tag{OW}
\]

Goal 01 shows that (OW) implies \(q_\chi\le5/8\). This cap closes the
authenticated support-33 terminal-zero row.

## Proof track

1. Reconstruct the transpose recurrence and prove the exact identity between
   zero coefficient symbols and zero nibbles in the transpose orbit.
2. Compute the exact two-state symbol branch number

   \[
   \min_{x\ne0}\bigl(\operatorname{wt}_4(x)+
   \operatorname{wt}_4(Sx)\bigr).
   \]

   If the value is at least 12, pair the \(N\) orbit states and prove (OW).
3. If two-state pairing fails, test increasing block lengths that divide
   \(N\). For a block length \(w\), the sufficient local statement is

   \[
   \min_{x\ne0}\sum_{j=0}^{w-1}\operatorname{wt}_4(S^jx)\ge6w.
   \]

4. If every practical local statement fails, retain boundary state and derive
   an amortized transition inequality. The inequality must telescope over the
   complete \(N\)-state interval.
5. Use the irreducible-component decomposition only when it reduces the
   certificate. Do not enumerate component-support cases without an explicit
   complexity bound.

Every claimed minimum or transition bound requires an independently
checkable certificate. A numerical mixed-integer optimizer can find a target,
but its status output alone is not a proof.

## Refutation track

1. Minimize the exact full-orbit weight from pure-component characters, every
   XOR combination of pure minimizers, and known mixed witnesses.
2. Run reproducible coordinate, packet-coordinate, and component-coordinate
   descents while retaining a nonzero high-component projection.
3. Replay every best witness independently.
4. If a character violates (OW), record its full component support, exact
   orbit weight, zero-symbol count, and coefficient histogram.

The refutation track continues while the proof track develops. Failure to
find a counterexample is diagnostic evidence only.

## Completion criteria

This goal ends with one of the following outcomes.

1. A rigorous proof of (OW), with primary and independent certificate checks.
2. An exact nonzero character that violates (OW), with independent replay.
3. A rigorously demonstrated obstruction to the attempted certificate route,
   together with one smaller finite lemma whose proof implies (OW). The new
   lemma must reduce either the state dimension, the block length, or the
   number of cases by a quantified amount.

## Excluded work

This goal does not prove the full construction. It does not cover other outer
supports or placement strata. It does not benchmark PacketMul or optimize the
encoder.
