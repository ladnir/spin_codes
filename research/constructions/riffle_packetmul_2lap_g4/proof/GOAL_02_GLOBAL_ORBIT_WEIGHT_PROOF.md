# Goal 02: proof of the global orbit-weight lemma

## Result

The remaining mixed-character lemma from Goal 01 holds.

Let \(T\) be the zero-input state map of Riffle PacketMul-2Lap g=4. Define

\[
S:=T^{\mathsf T},
\qquad N:=32{,}772.
\]

For every nonzero \(\chi\in\mathbb F_2^{64}\),

\[
\sum_{t=1}^{N}\operatorname{wt}_4(S^t\chi)
\ge 6N
=196{,}632.
\tag{1}
\]

Here \(\operatorname{wt}_4\) counts nonzero four-bit nibbles. Consequently,
the zero-symbol fraction from Goal 01 satisfies

\[
q_\chi\le\frac58.
\]

This result proves the mixed-character cap needed for the authenticated
support-33 terminal-zero row. It does not prove the complete construction.

## Orbit interpretation of coefficient symbols

Let \(J_s(y)\in\mathbb F_2^{64}\) place the field element
\(y\in\mathbb F_{16}\) in packet slot \(s\). The terminal contribution of
cell \((s,t)\) is

\[
z_y(s,t)=T^tJ_s(y).
\]

Fix a terminal character \(\chi\). Transposition gives

\[
\langle\chi,z_y(s,t)\rangle
=\langle S^t\chi,J_s(y)\rangle.
\tag{2}
\]

The right side is the binary pairing between \(y\) and nibble \(s\) of
\(S^t\chi\). The trace pairing on \(\mathbb F_{16}\) is nondegenerate.
Therefore,

\[
a_\chi(s,t)=0
\quad\Longleftrightarrow\quad
\text{nibble }s\text{ of }S^t\chi\text{ is zero}.
\tag{3}
\]

There are \(M=16N=524{,}352\) packet cells. Equations (2) and (3) give

\[
q_\chi
=1-\frac1M\sum_{t=1}^{N}\operatorname{wt}_4(S^t\chi).
\tag{4}
\]

Thus (1) is equivalent to \(q_\chi\le5/8\).

## Exact four-state lemma

The proof uses one finite local statement.

> **Four-state lemma.** For every nonzero \(x\in\mathbb F_2^{64}\),
> \[
> \sum_{j=0}^{3}\operatorname{wt}_4(S^jx)\ge24.
> \tag{5}
> \]

The bound is exact. The state

\[
x=\mathtt{0x96000624}
\]

produces four nibble weights

\[
(5,8,4,7),
\]

whose sum is 24.

## Exhaustive certificate for the four-state lemma

Suppose that a four-state window has total weight at most 23. At least one
state in that window has nibble weight at most five. Re-anchor the window at
such a state. Its position in the window is one of four possible phases.

The primary certificate enumerates every nonzero 64-bit state of nibble
weight at most five. The exact number of states is

\[
\sum_{r=1}^{5}\binom{16}{r}15^r
=3{,}411{,}004{,}740.
\]

For each state, the certificate evaluates all four anchor phases. It finds no
window of weight below 24. It also finds the weight-24 witness above.

The independent certificate uses a different cover. It enumerates each set of
five nibbles and every nonzero 20-bit word on that set. It therefore checks

\[
\binom{16}{5}(2^{20}-1)
=4{,}580{,}175{,}600
\]

words. States supported on fewer than five nibbles occur more than once. Every
state required by the low-anchor argument occurs at least once.

The implementations also reconstruct the recurrence differently. The primary
implementation computes the inverse of the BCH generator polynomial. The
independent implementation obtains the systematic BCH half by binary Gaussian
elimination. It describes the accumulator on each basis bit directly. Both
implementations find the same minimum window and the same witness.

These two exhaustive checks prove (5).

## Global orbit proof

Fix a nonzero character \(\chi\). The map \(S\) is invertible, so every state
\(S^t\chi\) is nonzero. Also,

\[
N=32{,}772=4\cdot8{,}193.
\]

Partition the orbit interval \(t=1,\ldots,N\) into 8,193 consecutive blocks
of four states. Apply (5) to the first state of each block. Summing the block
bounds gives

\[
\sum_{t=1}^{N}\operatorname{wt}_4(S^t\chi)
\ge 8{,}193\cdot24
=196{,}632
=6N.
\]

Equation (4) now gives

\[
q_\chi
\le1-\frac{196{,}632}{524{,}352}
=\frac58.
\]

This proves the global orbit-weight lemma.

## Why shorter local bounds fail

The analogous two-state and three-state bounds are false. Exact recurrence
replay gives a two-state window with weights

\[
(3,8),
\]

whose total is 11 rather than 12. A separate low-anchor search gives a
three-state window with weights

\[
(4,8,5),
\]

whose total is 17 rather than 18. The next states in these trajectories are
dense. The four-state certificate captures that amortization.

## Terminal-zero consequence

Goal 01 proves

\[
\beta_\chi=\frac{16q_\chi-1}{15}.
\]

The bound \(q_\chi\le5/8\) implies

\[
-\frac1{15}\le\beta_\chi\le\frac35,
\qquad
|\beta_\chi|\le\frac35.
\]

Substitution into the exact cross-value Parseval and distinct-cell ledger from
Goal 01 bounds the aggregate terminal-zero probability of the 26 authenticated
support-33 outer words by a value in

\[
[2^{-40.980718082062},2^{-40.980718082061}].
\]

The authenticated support-33 terminal-zero row is therefore closed below
\(2^{-40}\), with more than 0.980718082061 bits of margin.

## Scope

The four-state lemma, global orbit bound, zero-symbol cap, and support-33
terminal-zero consequence are exact. Other outer supports, placement strata,
and output-weight events remain outside this checkpoint. No performance claim
is made.

## Receipts and sources

- Primary certificate: `../receipts/goal02_window4_branch_primary.json`
- Independent certificate: `../receipts/goal02_window4_branch_independent.json`
- Goal 02 audit: `../receipts/goal02_audit.json`
- Primary source: `../../../scripts/certify_riffle_packetmul_2lap_g4_triple_branch.cpp`
- Independent source: `../../../scripts/verify_riffle_packetmul_2lap_g4_window4_branch.cpp`

