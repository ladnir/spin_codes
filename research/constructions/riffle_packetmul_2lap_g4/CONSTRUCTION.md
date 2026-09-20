# Riffle PacketMul-2Lap g=4: construction specification

## Candidate boundary

Riffle PacketMul-2Lap g=4 modifies Riffle DP-2Lap g=4 at one interface. It
randomizes every four-bit input packet before the existing packet permutation
and two-lap convolution. All other mechanisms remain unchanged until a later
version explicitly changes them.

## Packet field

Interpret a packet

\[
(x_0,x_1,x_2,x_3)\in\mathbb F_2^4
\]

as the field element

\[
x_0+x_1\alpha+x_2\alpha^2+x_3\alpha^3\in\mathbb F_{16},
\]

where

\[
\alpha^4=\alpha+1.
\]

Thus the polynomial basis uses the irreducible polynomial
\(X^4+X+1\), encoded as `0x13`.

## Setup randomness

Let \(I\) be the set of input packet coordinates. Setup samples

\[
A_i\gets\mathbb F_{16}^{\times}
\qquad(i\in I)
\]

independently. Setup uses rejection sampling to obtain the uniform law on the
15 nonzero field elements. The multiplier schedule is independent of the
packet permutation and every other randomized schedule.

A realized multiplier schedule is public encoder data. It is not secret key
material. The same schedule is used for both laps because multiplication
occurs once, before the first lap.

## Encoding map

Let \(x=(x_i)_{i\in I}\) be the packet vector produced by the unchanged outer
stage. Define the packet-local map

\[
(D_Ax)_i:=A_i x_i.
\]

Let \(\Pi\) be the unchanged packet permutation from Riffle DP-2Lap g=4. Let

\[
G(y):=F(F(y))\mathbin\oplus J(L(y))
\]

be the unchanged retained-state two-lap map. Riffle PacketMul-2Lap g=4
computes

\[
\operatorname{Enc}_{\Pi,A}(x):=G(\Pi(D_Ax)).
\]

The encoder samples no fresh multiplier, packet permutation, or state
permutation between laps.

For every fixed \((\Pi,A)\), the map
\(\operatorname{Enc}_{\Pi,A}\) is binary linear. Each packet multiplier is
invertible and preserves whether a packet is zero. Therefore, \(D_A\)
preserves packet support but not binary Hamming weight inside a nonzero
packet.

## Distributional effect

Fix \(x_i\ne0\). Since \(A_i\) is uniform in
\(\mathbb F_{16}^{\times}\), the product \(A_ix_i\) is uniform over the
15 nonzero packet values. Products from distinct packet coordinates are
independent before conditioning on any other event.

This distributional statement is the intended proof benefit. One shared
multiplier for several packets is not part of this candidate.

## Scope

This specification freezes the candidate identity. It does not assert a
distance bound, a terminal-zero bound, or a performance result. The outer
code, packet permutation law, convolution matrix, retained-state rule, and
lap count require separate inherited-mechanism audits.
