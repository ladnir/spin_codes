# A systematic RS inner: certified 40% distance at about 6.3 ms

This construction replaces the layered LRS inner maps with a single RS
emission map. At K=2^20 and rate 1/2, its complete certificate gives at
least 40% GF65536-symbol distance with more than 50.06 bits of setup-failure
margin. The native forward encoder takes about 6.3 ms for 128-bit records,
each carrying eight parallel GF65536 symbols under the same setup.
Setup is precomputed and excluded from that encoding time.

## Construction and setup

Let F=GF65536, q=65536, R=8192, K=128R, and N=256R. Retain the
[certified outer and route](HD_SUMRANK_FEEDBACK_CERTIFICATE.md): independent
RS[256,128,129] rows over the embedded AES-GF256, independent uniform
nonzero F scalars on all outer output coordinates, and an independent
uniform permutation of R positions within each of 256 coordinate regions.
There is one regional routing stage. Each region supplies 128 packets of
64 symbols. All setup choices are fixed and reused for every message.

Let n_j be the normalized additive-polynomial basis on the ordered binary
AES-GF256 basis. In particular deg(n_j)=j. Define the 64-by-16 matrix A
by A[i,j]=n_j(i), for i=0,...,63 and j=0,...,15, and put C=A^T. These
are fixed F-linear maps. Every restriction of A to h rows has rank
min(16,h), because this is evaluation of all polynomials of degree below
16 at distinct points. Also A^T A=0: the additive 64-point subspace has
a subspace polynomial with constant nonzero derivative, so the sum of a
polynomial of degree below 63 over that subspace is zero. Each product
n_j n_k has degree at most 30.

For each packet independently sample M=I+u v^T, with u uniform nonzero
in F^16 and v uniform conditioned on v^T u != 1. Independently sample
T=D2 H D1, where both diagonals have independent uniform nonzero F
entries and H is the fixed sixteen-point Cauchy interpolation map from
nodes 0,...,15 to nodes 16,...,31. These coins are independent of the
outer and route. Encode each packet by

    Y = X + A Q,
    Q_next = M Q + T C X,       Q_initial = 0.

The final state is unrestricted. Given the previous state and Y, recover
X=Y+A Q, then compute the next state. Thus the inner is invertible on F^N
and the full code has dimension K. The failure event concerns the
existence of any nonzero message after the entire setup has been revealed.

## Exact packet moments

Fix a potential input set U of size j, and an output-zero set Z. Write
c=|U intersect Z| and d=|Z minus U|. The identity feedthrough fixes c input
variables. The remaining quiet constraints have rank min(16,j-c), while
the remaining constraints on the incoming state have rank min(16,d).
These ranks follow directly from the MDS coordinate restrictions of A.
No sum-rank theorem or tile-occupancy profile is required.

Let X be uniform on F^U, including artificial zero coordinates, and let
0<z<1. Write K_j for E[1[CX=0] z^wt(X)], T_j for E[z^wt(X)], and V_j
for E[z^wt(X+A Q)] with Q independent uniform nonzero in F^16. Expanding
the Hamming mark over positive output-zero constraints gives

    K_j = sum_(c=0)^j binom(j,c) z^(j-c)(1-z)^c
                         q^(-c-min(16,j-c)),
    T_j = (z+(1-z)/q)^j,
    V_j = T_j sum_(d=0)^min(15,64-j) binom(64-j,d)
                z^(64-j-d)(1-z)^d (q^(16-d)-1)/(q^16-1).

For nonzero Q, A Q has weight at least 49. Changing j input coordinates
leaves at least max(0,49-j) nonzero outputs. The fixed-state quiet moment
is therefore bounded by the minimum of K_j,
q^(-min(16,j)) z^max(0,49-j), and (q^16-1)V_j. These exact rational
moments enter the same positive two-component feedback-refresh operator
as in the preceding proof. There are now only 65 packet classes, with
multiplicities binom(64,j).

Pair consecutive packets before applying the fractional moment bound.
The existing regional coefficient extraction enforces exactly r potential
positions in every region. The outer counting-measure comparison and
projective message count are unchanged. The inverse-output dense tail
still applies because the inner is invertible.

First analyze an ideal independent transitive T. The
[diagonal-Cauchy density comparison](HD_DIAGONAL_CAUCHY_REFRESH.md) then
costs 32768 log2((q^16-1)/(q-1)^16) < 11.541649 bits for the implemented
distribution. This is an explicit loss, not an assertion that T is uniform.

## Complete certificate and validation

The [pinned replay](theory/hd_systematic_mds_witnesses.py) contains 78
exact rational witnesses, plus the dense tail. Its 256-bit and independent
384-bit runs both cover every r=1,...,8192 and give

    log2(failure) < -50.06028067205724.

This bounds the probability that any nonzero message has output weight
at most 838860. Outside that event the distance is at least 838861/N,
which exceeds 40%. The shared setup suffices for eight SIMD field lanes;
it does not require eight independent certificates.

The [exact-map tests](theory/hd_systematic_mds_tests.py) check
self-orthogonality, concrete global mask ranks, every output-zero mask
of a smaller construction, and agreement between rational and floating
moments. The theorem uses the universal MDS rank argument above, not
sampled ranks. The native F,A,C maps match the independent Python digest
11358755321456897600. Native basis, random-state, routed-offset, and
full aliasing tests agree with an independent dense scalar oracle.

## Implementation checkpoint

State stays in polynomial-coefficient coordinates. Emission evaluates four
sixteen-point FFTs, and feedback uses their transposes. This removes the
Moore maps and tile-weight transforms of the layered construction. The
outer retains fused scalar multiplication and sequential output planes;
the inner retains next-region L2 prefetch and scalar routed record loads.

On the established Peach CPU2/4.5-GHz harness, the first 101-call median
is 6.188 ms and a separate 301-call median is 6.234 ms (p10/p90
6.209/6.272 ms). The checksum is 976174626021906176. Separate component
controls measure 2.411 ms outer, 3.651 ms routed inner, and 3.128 ms
pre-routed inner; these isolated measurements are not additive phase
measurements of one call. Reproduce with the isolated target
`spin_char2_hd_sumrank_forward --bench 301 full 20 mds`.

This is a research checkpoint, not a production-library change. The
older certified LRS version and all its proof inputs remain preserved.
The fixed-GFNI-affine implementation has identical output and the same
certificate. Two 301-call medians are 6.261 and 6.297 ms around a
6.343-ms ordinary-multiplication control. An earlier comparison was
6.323 ms versus 6.597/6.474-ms controls. The advantage is small and there
is run-to-run variation; use about 6.3 ms as the current checkpoint, not
an improvement claim relative to the earlier single 6.188-ms sample.
Reproduce this variant with `--bench 301 full 20 mds-affine`.

Next, reduce the fixed arithmetic and refresh costs against this anchor.

## Packet-local sorting

The proof permits a deterministic order inside each sampled packet. First
sample each region's uniform permutation as before. Then sort the 64 source
indices within each packet, without moving any index between packets.
This defines a different routing distribution and generally changes the
encoded word; it is not a claim of identical output.

Fix a set of r potential outer rows. Sorting preserves every packet's
potential count j. Its joint count distribution is therefore unchanged:
within a region it has mass product_t binom(64,j_t)/binom(8192,r), subject
to sum_t j_t=r. Conditional input positions inside a packet need not be
uniform after sorting. The rank argument above, however, gives the same
K_j,T_j,V_j and quiet cap for every set of j positions. The proof uses
only these counts and bounds, so every pinned witness and the final
50.06-bit margin transfers without loss. All other setup choices remain
independent and unchanged. Sorting across packet boundaries is not
covered by this argument.
