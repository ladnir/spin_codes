# Deterministic gamma feedback with state dimension 31

The fixed-gamma argument also permits state dimension 31. The earlier public
whitelist stopped at 30; even state dimension is not a proof requirement.
This note records the extension as a fallback construction, not a distance
certificate. It uses the joint return argument in
[the gamma-feedback note](HD_FIELD128_GAMMA_RETURN.md).

## Maps and recurrence

Let E=GF(256), let F=GF(2^128), and fix the existing embedding E into F.
Use the existing byte labels for elements of E. Let A be the 64-by-31
Lagrange evaluation matrix from nodes 128..158 to nodes 0..63, and set C=A^T.
Thus each column of A evaluates a polynomial of degree at most 30.

Let `H31` be the 31-by-31 Lagrange evaluation matrix from nodes 0..30 to
nodes 31..61. The two node sets are distinct and disjoint. Lagrange's formula
writes `H31` as a Cauchy matrix multiplied by invertible diagonal matrices.
Every square minor is therefore nonsingular. This argument does not require
either 31-node set to be an additive subspace.

Fix gamma in F of degree 16 over E. In the existing sparse tower

    F = E[Y]/(Y^16 + Y^5 + Y^3 + 0x32),

gamma=Y is such an element; `0x32` is the existing AES-field byte constant.
This choice uses the previously checked irreducibility of that polynomial.
It adds no random setup choice. Multiplication by gamma can use the existing
isolated multiply-by-Y primitive, but there is no state 31 encoder here.

Keep the existing full-field Toeplitz outer ensemble and paired route.
For each input packet X_t in F^64 and state Q_t in F^31, define

    Y_t = X_t + A Q_t,
    Q_(t+1) = D2_t H31 D1_t (Q_t + gamma C X_t).

The diagonals `D1_t` and `D2_t` have independent entries uniform in F minus zero,
independently across packets and from the outer setup. The initial state is
zero. State persists across regions, and the final state is unrestricted.

## Deterministic principal nonsingularity

The 64 evaluation nodes form the existing additive binary subspace. Its power
sums vanish through degree 62. Products of columns of A have degree at most 60,
so A^T A=0. Each Lagrange row sums to 1; characteristic 2 therefore gives
every diagonal entry of A A^T equal to 1.

Fix any coordinate set D of size c. Put

    M_D = A_D A_D^T,       p_D(t) = det(I_c + t M_D).

The symmetric matrix M_D has unit diagonal. In the determinant expansion,
cycles of length at least 3 cancel with their reversals. The remaining terms
are fixed points and disjoint transpositions. Consequently

    p_D(t) = R_D(t)^2                 if c is even,
    p_D(t) = (1+t) R_D(t)^2           if c is odd,

with R_D(0)=1. Sylvester's identity also gives

    p_D(t) = det(I_31 + t A_D^T A_D),

so deg p_D<=31. In the even case its degree is even and at most 30; in the
odd case removing 1+t leaves degree at most 30. Both cases give deg R_D<=15.
The degree 16 minimal polynomial of gamma cannot divide R_D, and gamma!=1.
Thus p_D(gamma)!=0 for every D. This is deterministic, with no bad-setup event.

## Local caps and limitations

The MDS rank and joint return proofs do not require even state dimension.
For occupancy j, c active zeros, and e inactive zeros, put u=min(31,j-c).
The returning nonzero-syndrome probability is bounded by

    (q^min(u,31-e)-1) q^(j-c-u) / q^j,       e<31,

and is zero for e>=31, where q=|F|. The local transfer matrix uses
d=q^31-1 and refresh density B=d/(q-1)^31. Its entries retain the form
K, B(T-K), J_gamma/d, and B V. The packet, route, and outer dimensions
remain unchanged.

State 31 changes the inner construction. The emission code is RS[64,31],
with minimum distance 34 instead of 35 for state 30. The refresh uses 62
general-field diagonal products per packet instead of 60. The larger return
denominator can help, but the emission tradeoff prevents a monotonicity claim.

No state 30 numerical ceiling transfers to this construction. State 31 requires
fresh local coefficients, support sums, exact witness metadata, and directed
coverage. This note supplies no encoder implementation, timing measurement,
or directed numerical certificate for state 31.
