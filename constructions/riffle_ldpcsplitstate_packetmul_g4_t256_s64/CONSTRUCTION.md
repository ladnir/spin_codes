# Construction

`Riffle LDPCSplitState PacketMul g=4 t=256 s=64` modifies
`Riffle LDPCSplitState g=4 t=256 s=64` by adding one packet-local linear map.
The outer code, fixed four-block packet groups, packet permutations, and
LDPCSplitState inner are unchanged.

For every transposed region `c` and logical four-bit packet `j`, setup samples

\[
\beta_{c,j}\gets\operatorname{GF}(16)^*.
\]

Before applying the packet-position permutation, the encoder replaces packet
`X_{c,j}` by

\[
\widetilde X_{c,j}:=\beta_{c,j}X_{c,j}.
\]

Setup fixes all multipliers.  The resulting encoder remains binary linear.
The multiplier schedule is public and makes no secrecy or pseudorandomness
claim.

There are `256*2048=524288` multipliers.  A four-bit representation per
multiplier occupies 256 KiB before schedule packing.  The intended evaluator
fuses multiplication into packet routing; no implementation is selected yet.

## Packet law

Fix a packet group containing `r` Bernoulli(`p`) candidate bits.  The
unmultiplied packet is zero with probability `(1-p)^r`.  Conditioned on being
nonzero, multiplication by an independent uniform element of `GF(16)^*`
makes its value uniform over the 15 nonzero nibbles.  Thus its binary-weight
generating function is

\[
G_r(u)=(1-p)^r+\left(1-(1-p)^r\right)
\frac{(1+u)^4-1}{15}.
\]

This law is averaged over setup randomness in the distance argument.  A final
fixed-code claim must use the resulting ensemble bound or exhibit a fixed
multiplier schedule that satisfies the required certificate.
