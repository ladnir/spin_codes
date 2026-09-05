# Riffle ParallelAcc g=4: diagnostic construction

## Purpose

Riffle ParallelAcc g=4 tests whether the proof can separate a recursive constituent
from the deterministic structure of its input. It is a diagnostic construction,
not the active deployment candidate.

## Frozen outer and permutation

The construction retains the outer stage and packet law of Riffle DP g=4.
The outer stage produces (524{,}352) packets in \(\mathbb F_2^4\). Setup
samples one uniform bijection of these packet positions. No packet multiplier,
second permutation, or fresh per-codeword randomness is added.

## Accumulator constituent

Let \(x_1,\ldots,x_n\in\mathbb F_2^4\), where \(n=524{,}352\), denote the
permuted packets. The accumulator starts in state \(s_0=0\) and computes

\[
s_t:=s_{t-1}+x_t,
\qquad
y_t:=s_t
\]

for \(1\le t\le n\). Addition is bitwise XOR. The binary output is the
concatenation of \(y_1,\ldots,y_n\).

This definition gives four lane-parallel binary accumulators. It does not
require a termination constraint. The final state is

\[
s_n=\sum_{t=1}^n x_t.
\]

## Construction delta

Relative to Riffle DP g=4, this diagnostic replaces the complete deterministic
inner map by the four-bit accumulator above. The outer binary matrix and the
uniform packet permutation remain unchanged. The replacement therefore defines
a distinct binary matrix ensemble.

## Cost scope

The forward binary map needs one four-bit XOR per packet. In the transposed PCG
evaluation, the four lanes require about one field-element XOR per binary
coordinate. This diagnostic is proof-facing and has no performance claim.

## Claim scope

No distance claim is attached to Riffle ParallelAcc g=4. Goal 01 studies its exact
input-output spectrum on reduced instances and identifies the deterministic
riffle lemma that a full proof would require.
