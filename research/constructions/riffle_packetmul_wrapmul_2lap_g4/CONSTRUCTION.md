# Riffle PacketMul-WrapMul-2Lap g=4: construction

## Candidate boundary

Riffle PacketMul-WrapMul-2Lap g=4 modifies
Riffle PacketMul-2Lap g=4 at the wrap between laps. Setup samples one
additional nonzero 64-bit field scalar. The encoder multiplies the first-lap
terminal state by this scalar before the retained-state injection.

The outer stage, packet multipliers, packet permutation, convolution map, and
lap count remain unchanged.

## Terminal-state field

Write a terminal state as

\[
s=\sum_{i=0}^{63}s_i\beta^i,
\qquad s_i\in\mathbb F_2.
\]

The bit at position \(i\) is the coefficient of \(\beta^i\). Arithmetic uses

\[
\mathbb F_{2^{64}}
:=\mathbb F_2[X]/(X^{64}+X^4+X^3+X+1).
\]

The implicit-top reduction constant is 0x1b. Goal 01 audits irreducibility of
the modulus.

## Setup

The inherited setup samples the packet permutation \(\Pi\) and the packet
multiplier schedule \(A\). Setup also samples

\[
b\gets\mathbb F_{2^{64}}^\times
\]

independently of \((\Pi,A)\). Rejection of the all-zero 64-bit string gives the
uniform law on \(\mathbb F_{2^{64}}^\times\).

The realized scalar \(b\) is public encoder data. Setup samples it once and
reuses it for every codeword. Encoding does not sample fresh randomness
between laps.

## Encoding map

For an outer packet vector \(x\), define

\[
y:=\Pi(D_Ax).
\]

Let \(F\) be the unchanged one-lap output map. Let \(L(y)\) be its terminal
state, and let \(J(s)\) be the output caused by entering the second lap in
state \(s\) with zero drive. Define field multiplication by

\[
M_b(s):=bs.
\]

The new retained-state two-lap map is

\[
G_b(y):=F(F(y))\mathbin\oplus J(M_b(L(y))).
\]

The encoder computes

\[
\operatorname{Enc}_{\Pi,A,b}(x):=G_b(\Pi(D_Ax)).
\]

For every fixed \((\Pi,A,b)\), this map is binary linear. The map \(M_b\) is
invertible and preserves whether the terminal state is zero.

## Randomness scope

The same scalar \(b\) acts on every outer word. Wrapped states from different
outer words are therefore correlated. The distance proof may use a union
bound across outer words, which does not require independence between them.

The multiplier \(b\) is not a secret and provides no pseudorandomness claim.
Its role is only to randomize the nonzero terminal-state orbit in the setup
experiment.

## Scope

This specification makes no distance or performance claim. In particular, it
does not claim that a uniform wrapped state makes the second-lap output
uniform. Goal 01 reduces that output question to a finite bad-state list.
