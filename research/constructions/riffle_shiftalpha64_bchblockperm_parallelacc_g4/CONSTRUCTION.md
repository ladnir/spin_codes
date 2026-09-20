# Riffle ShiftAlpha64-BCHBlockPerm-ParallelAcc g=4

## Outer code

Let

\[
\mathbb F:=\mathbb F_2[X]/(X^{64}+X^4+X^3+X+1),
\]

and let \(\gamma\) denote the residue class of \(X\). For data symbols
\(m_0,\ldots,m_{16383}\in\mathbb F\), define

\[
p_0:=\sum_{i=0}^{16383}m_i,
\qquad
p_1:=\sum_{i=0}^{16383}\gamma^{64+i}m_i.
\]

The coefficients are distinct and nonzero. The systematic field code
therefore has block distance three.

The encoder computes

\[
q:=\sum_{i=0}^{16383}\gamma^i m_i
\]

with the existing Horner recurrence. It then computes

\[
p_1:=\gamma^{64}q.
\]

The final multiplication uses 64 applications of the sparse shift-and-reduce
map. It occurs once per complete outer word, not once per data symbol.

## BCH blocks and permutations

Encode every data symbol and both parity symbols with the extended binary BCH
code \([128,64,22]\). Independently permute all 128 bits inside every BCH
block. Divide each permuted BCH block into 32 consecutive four-bit packets.

Apply one uniform global permutation to all 524,352 packets. Setup samples
all permutations once. Encoding is linear and deterministic after setup.

## Parallel accumulator

Let \(u_1,\ldots,u_N\in\mathbb F_2^4\) be the globally permuted packets. Set

\[
s_0:=0,
\qquad
s_t:=s_{t-1}+u_t,
\qquad
y_t:=s_t.
\]

Addition is bitwise XOR. The binary output is the concatenation of the
\(y_t\).

## Candidate boundary

This candidate differs from **Riffle BCHBlockPerm-ParallelAcc g=4** only in
the exponent offset of the double-parity coefficients. Any other offset is a
different schedule and must be identified explicitly.
