# Construction

Let \(g=4\), \(\sigma=18\), and \(N=524352\). The message contains 4096
groups of 256 bits. Fix binary bases for \(\mathbb F_{2^{256}}\) and
\(\mathbb F_{2^{22}}\).

## Data constituents

For each data-group position \(i\), setup samples

\[
(a_i,b_i)\gets
\mathbb F_{2^{256}}^2\setminus\{(0,0)\}.
\]

For a group input \(m_i\in\mathbb F_{2^{256}}\), define

\[
E_i(m_i):=(a_i m_i,b_i m_i)\in\mathbb F_2^{512}.
\]

For every fixed nonzero \(m_i\), the value \(E_i(m_i)\) is uniform over
\(\mathbb F_2^{512}\setminus\{0\}\). The map \(E_i\) is injective for every
allowed setup value. Different data groups use independent pairs.

The field outer retains two fixed MDS parity equations over
\(\mathbb F_{2^{64}}\). Each parity symbol is encoded by the extended BCH
\([128,64,22]\) code. Setup samples one independent coordinate permutation
for each parity codeword. The outer encoder forms four-bit packets and applies
one uniform global permutation to all \(N\) packets.

## Inner convolution

Write each input-state pair as

\[
v_j:=(x_j,s_j)\in\mathbb F_2^{22},
\qquad
s_j\in\mathbb F_2^{18}.
\]

For each packet position \(j\), setup samples

\[
c_j\gets\mathbb F_{2^{22}}.
\]

Initialize \(s_1:=0^{18}\). Compute

\[
(y_j,s_{j+1}):=c_jv_j,
\]

where the fixed binary basis splits the product into four output bits and
18 next-state bits. The encoder discards \(s_{N+1}\).

For every fixed nonzero \(v_j\), multiplication by a uniform \(c_j\) produces
a uniform element of \(\mathbb F_2^{22}\). Therefore this scalar family gives
the same live-step transition law as a uniform random \(22\times22\) binary
matrix. Independence across positions preserves the existing transfer proof.

## Randomness scope

All randomness belongs to setup. Encoding is deterministic after setup. The
candidate does not use a computational pseudorandom generator.
