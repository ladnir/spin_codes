# Proof reduction

Fix a shared-group profile `n=(n1,n2,n3,n4)` and the Bernoulli envelope
parameter `p`.  Define

\[
q_r:=1-(1-p)^r.
\]

A width-`r` packet survives as a nonzero packet with probability `q_r`.
Every surviving packet is uniform over `GF(16)^*`, independently of its
original width.  Therefore

\[
K_n:=\sum_{r=1}^4\operatorname{Binomial}(n_r,q_r)
\]

is the number of nonzero packets in one region, with independent summands.

Let `U_k(z)` be the exact two-state region matrix conditioned on placing `k`
independent uniform nonzero nibbles in a uniform `k`-subset of the 2048 packet
positions.  The complete profile matrix factors as

\[
R_n(p,z)=\sum_{k=0}^{|n|}\Pr[K_n=k]U_k(z).
\]

The factorization is exact.  It replaces the four-variable matrix coefficient
calculation by:

1. one table `U_0,...,U_m`, where `m` is the largest packet-group count; and
2. one scalar Poisson-binomial distribution for each profile.

Direct comparison with the original typed calculation agrees to about
`1.6e-15` relative error on the tested profiles.

## Current result

At occupation 128, `p=1/2`, and log-surprisal
`-3.546242271905631`, all 16,335 profiles have positive 9% margin.  The
all-quad profile `(0,0,0,32)` is worst with 7,105.481 bits.
Summing the bounds for all profiles changes this margin by less than float64
resolution because the second-worst profile is 94.6 bits safer.

PacketMul does not improve that endpoint.  A width-four fair input is already
uniform over all 16 nibbles, so multiplication preserves its law.  The
all-quad termination inflation remains 1,007.389 bits.  The corresponding
tilted expected termination count remains about 123.234 by the same argument.

PacketMul also does not produce entrywise domination by the fully split
profile.  Its benefit is the survivor-count reduction, not a new matrix
order.

The next proof target is a scalable bound on

\[
\sum_n N(n)\,\mathbf e_Z^T R_n(p,z)^{256}\mathbf 1,
\]

where `N(n)` is the exact active-block multiplicity.  The survivor
factorization should permit a one-dimensional generating-function or saddle
bound without enumerating every profile at larger occupations.

## Active-packet collapse audit

The inner transfer can indeed be expressed entirely through the surviving
packet count `K`.  The fixed four-block profile remains correlated across the
256 regions, so it cannot be averaged independently inside each matrix
factor.

`ACTIVE_PACKET_COLLAPSE.md` tests a profile-free exponential envelope.  The
envelope reduces the complete outer sum to univariate occupation
coefficients, but it is numerically too weak: the four-transition version has
`-12034.886` bits of margin at occupation 128, compared with `+7105.481` bits
for the exact profile sum.  The loss is concentrated in a path that remains
in the zero state for 255 regions and activates in the last region.

The next reduction should retain an exact or piecewise survivor-count law
during zero-state runs.  After the first activation, the coarser count
envelope may be sufficient.
