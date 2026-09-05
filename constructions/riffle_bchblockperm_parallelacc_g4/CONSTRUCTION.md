# Riffle BCHBlockPerm-ParallelAcc g=4

## Parameters

The message contains 16,384 symbols in \(\mathbb F_{2^{64}}\). The inherited
double-parity outer appends two field symbols. It therefore produces 16,386
symbols.

Let

\[
B:\mathbb F_2^{64}\to\mathbb F_2^{128}
\]

be the frozen extended primitive narrow-sense BCH \([128,64,22]\) encoder.
Applying \(B\) to every field symbol gives binary length

\[
L=16{,}386\cdot128=2{,}097{,}408.
\]

## Permutations and packetization

For each BCH block index \(i\in[16{,}386]\), setup samples an independent
permutation

\[
\sigma_i\gets S_{128}.
\]

The encoder applies \(\sigma_i\) to all 128 coordinates of BCH block \(i\).
It then divides that block into 32 consecutive four-bit packets. Packets do
not cross BCH block boundaries.

After packetization, setup samples an independent global permutation

\[
\Pi\gets S_{524{,}352}.
\]

The encoder applies \(\Pi\) to the complete packet sequence. The permutations
are fixed after setup. Encoding is linear and deterministic after setup.

This candidate does not reuse one block permutation across different BCH
blocks. Permutation reuse defines a different ensemble.

## Parallel accumulator

Let \(u_1,\ldots,u_N\in\mathbb F_2^4\) be the globally permuted packets, where
\(N=524{,}352\). Define

\[
s_0:=0,
\qquad
s_t:=s_{t-1}+u_t,
\qquad
y_t:=s_t.
\]

Addition is bitwise XOR. The binary output is the concatenation of
\(y_1,\ldots,y_N\).

## Candidate boundary

Changing the BCH code, block-permutation law, packet width, global
packet-permutation law, accumulator recurrence, or termination rule defines a
different candidate.

