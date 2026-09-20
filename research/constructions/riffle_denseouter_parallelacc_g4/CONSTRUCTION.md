# Riffle DenseOuter-ParallelAcc g=4: construction

## Parameters

Fix a message length \(k\), outer memory \(M\), and packet width \(g=4\).
Set \(x_t=0\) for \(t\le 0\) and \(t>k\).

## Terminated dense outer

For each \(t\in[k+M]\), sample an independent vector

\[
a_t\gets\mathbb F_2^M,
\qquad
g_t:=(1,a_t).
\]

Define

\[
p_t:=\left\langle
g_t,(x_t,x_{t-1},\ldots,x_{t-M})
\right\rangle.
\]

The outer output contains the \(k\) systematic bits and all \(k+M\) parity
bits:

\[
C_{\mathrm{out}}(x):=(x_1,\ldots,x_k,p_1,\ldots,p_{k+M}).
\]

Its unpadded length is \(L_0=2k+M\). Append fewer than four zero coordinates
to obtain length \(L=4N\), and divide the result into \(N\) consecutive
four-bit packets.

Emitting the final \(M\) parity bits is essential. It ensures that every
connected support cluster pays a complete tail of \(M\) active parity
positions. A long zero gap can split one message span into several clusters.

## Packet permutation and inner map

At code construction, sample one uniform permutation \(\Pi\gets S_N\).
Apply \(\Pi\) to the packet positions. No within-packet permutation or packet
multiplier is used.

For permuted packets \(u_1,\ldots,u_N\in\mathbb F_2^4\), set

\[
s_0:=0,
\qquad
s_t:=s_{t-1}+u_t,
\qquad
y_t:=s_t
\]

for \(t\in[N]\). Addition is bitwise XOR. The binary codeword is the
concatenation of \(y_1,\ldots,y_N\).

The outer taps and \(\Pi\) are setup randomness. Encoding a message after
setup is deterministic.

## Candidate boundary

Changing the termination parities, packet width, packet-permutation law, or
accumulator recurrence defines another candidate.
