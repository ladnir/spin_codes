# Exact regional coefficients for the randomized-outer quiet bound

Exact regional coefficients extend the certified quiet classes for the same
randomized-outer construction. The new bounds cover all nonzero quiet words
with 954 through 962 active rows, and quiet words of weight at most 838860
with 5865 through 8192 active rows. Their respective union bounds are below
2^-226.2046 and 2^-355.6458.

Together with the existing sparse result, quiet occupancies 1 through 962
and 5865 through 8192 are covered. Quiet occupancies 963 through 5864 remain
open. This note does not change the construction, cover the remaining
nonquiet words, or establish a whole-code distance certificate.

The numerical comparison separates two losses. Optimizing the coefficient
marker improves the earlier coarse marker grid. Extracting the regional
coefficient then gains another 1610 to 1980 bits at the seven initial
occupancies. The interior witnesses still fail by substantially more than
this coefficient correction.

## The probability bound and its coefficient

Use the setup and symbol-weight convention in
[the randomized-outer proof](HD_RANDOMIZED_OUTER_PROOF.md).
There are R=8192 outer message rows, 256 independently permuted coordinate
regions, and P=32768 chronological packets. A region contains 128 packets,
each of width 64. The scalar field has q=65536 elements. Put D=838860 and

    beta = q^256/(q-1)^128.

The probability is over the sampled outer scalars and regional routes.
The bad event quantifies over every nonzero message in the stated class.
Mixers have no effect along a quiet trajectory: the incoming state is zero
before every emitted packet. The final outgoing state is unrestricted.

Fix an active-row set I of size r. Under the artificial uniform measure
used in the outer comparison, its routed potential input positions form
a uniform r-subset of each region, independently between regions. A
potential position carries a uniform field element, which may be zero.

For a packet subset S, define

    K_S(z) = E_X[1[CX=0] z^wt(FX)],
    T_S(z) = E_X[z^wt(FX)],

where X is uniform on the coordinates in S and zero elsewhere. Let
K_j(z) and T_j(z) average these quantities over uniform j-subsets S.
Write Khat_j(z) for the upper bound on K_j from the
[joint quiet/output calculation](HD_RANDOMIZED_OUTER_JOINT_QUIET.md).
That calculation evaluates T_j exactly. Here, 0<z<=1 and 0<alpha<=1.

For a fixed route, first average the amplitudes and apply the clipped
first-moment bound for I. The inequality min(1,x)<=x^alpha then leaves
products of K_S(z)^alpha, with T_S(z)^alpha in the final packet.
Conditional on all packet occupancy counts, the within-packet subsets are
independent uniform subsets. Concavity gives

    E_S[K_S(z)^alpha] <= K_j(z)^alpha <= Khat_j(z)^alpha,
    E_S[T_S(z)^alpha] <= T_j(z)^alpha.

Thus the following polynomials retain the exact regional occupancy law:

    U(x) = sum_(j=0)^64 binom(64,j) Khat_j(z)^alpha x^j,
    V(x) = sum_(j=0)^64 binom(64,j) T_j(z)^alpha x^j.

Define

    A_r = [x^r] U(x)^128 / binom(R,r),
    H_r = [x^r] U(x)^127 V(x) / binom(R,r).

The normalized coefficient A_r bounds one ordinary region. H_r bounds
the final region, including its free last packet. Independence between
regional permutations and the union over active-row sets give

    B_r = binom(R,r) [beta^r z^(-D)/(q-1)]^alpha A_r^255 H_r.

This is an upper bound for the entire r-row quiet class. The active-row
binomial remains outside alpha; the projective witness division remains
inside alpha. The word "exact" refers to extracting coefficients of these
upper polynomials. It does not assert that Khat_j equals K_j.

## Stable extraction and the precise marker loss

Fix a rational scaling marker 0<p<1. Set

    u_j = binom(64,j) p^j(1-p)^(64-j) Khat_j(z)^alpha,
    v_j = binom(64,j) p^j(1-p)^(64-j) T_j(z)^alpha,
    S_K = sum_j u_j,             S_T = sum_j v_j,
    c_r = binom(R,r) p^r(1-p)^(R-r).

Let pi_K be the coefficient of x^r in (sum_j u_j x^j/S_K)^128.
Let pi_H be the same coefficient with 127 ordinary factors and one
normalized factor sum_j v_j x^j/S_T. Then

    A_r = pi_K S_K^128/c_r,
    H_r = pi_H S_K^127 S_T/c_r.

These identities hold for every marker. The marker rescales the polynomial;
it does not change the exact coefficients A_r and H_r. The helper forms
the two normalized positive polynomials and extracts their coefficients
with 192-bit Arb polynomial arithmetic. It rejects a computation if the
coefficient enclosure does not establish positivity.

The earlier marker bound drops pi_K and pi_H. Consequently its exact
excess over log2 B_r is

    -255 log2(pi_K) - log2(pi_H).

This identity isolates coefficient extraction from any improvement to the
local moments. The marker is selected near a saddle point so that the
required coefficient is numerically well scaled. The stationary equation
includes the final free factor. A floating bisection proposes the marker;
the resulting rational marker is replayed with Arb before acceptance.

## A bounded comparison in the former middle range

The following table compares the optimized marker with exact extraction
at the same z and alpha. All bound values were replayed at 192-bit Arb
precision; the displayed decimals are rounded for comparison.

| Active rows r | z | alpha | Optimized marker log2 bound | Exact coefficient log2 bound | Coefficient gain |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 954 | 1 | 3/32 | -120.7308 | -1788.7050 | 1667.9742 |
| 1024 | 1 | 1/32 | 7007.1852 | 5397.0467 | 1610.1385 |
| 1536 | 1/4 | 1/32 | 35168.8773 | 33495.7205 | 1673.1568 |
| 2048 | 1/4 | 1/32 | 54231.0437 | 52528.6599 | 1702.3838 |
| 4096 | 1/4096 | 1/32 | 76784.4999 | 75014.7913 | 1769.7086 |
| 5120 | 1/4096 | 1/32 | 43828.8814 | 42088.1566 | 1740.7247 |
| 5939 | 1/16384 | 3/8 | -42798.9686 | -44778.9596 | 1979.9910 |

Witness selection used z=2^-h for h in {0,1,2,4,6,8,10,12,14,16},
plus 2/[3(q-1)] and 133/[128(q-1)]. It used
alpha in {1,2,3,4,6,8,12,16,24,32}/32. Every marker was optimized by
the stationary equation, rather than selected from the earlier coarse grid.

This finite comparison does not optimize over every z and alpha. In
particular, alpha approaching zero gives the trivial active-row subset
count. Positive entries in the table establish failure of those witnesses,
not failure of the construction or of every possible witness.

## The sparse extension from one coefficient

At z=1, the local quiet moment is exactly q^-min(j,4), and T_j(1)=1.
Fix alpha=3/32. Couple a uniform r-subset of a region with a uniform
(r+1)-subset by adding a uniformly selected missing position. Every packet
count either stays fixed or increases. Each local factor
q^(-alpha min(j,4)) is nonincreasing, while the final free factor is one.

Applying this coupling independently in every region proves that

    G_r := A_r^255 H_r

is nonincreasing in r. For 954<=r<=962, both the binomial factor and
beta^(alpha r) increase. Hence

    B_r <= B_954 * binom(R,962)/binom(R,954) * beta^(8 alpha).

Only B_954 requires a polynomial coefficient calculation. Directed
evaluation gives

    log2 B_954 <= -1788.7050398873592452304481136452968752561766,
    log2 max_(954<=r<=962) B_r <= -229.3745503597183845401043120308651536413558,
    log2 sum_(r=954)^962 B_r <= -226.2046253582760721771968341429695206238361.

Because z=1, this extension covers all nonzero quiet words, regardless of
weight. The proof covers every integer in the interval by monotonicity;
it does not interpolate between successful sampled occupancies.

## The dense extension from four coefficients and two endpoints

Use z=1/16384 and alpha=3/8. Exact regional coefficients give

| Active rows r | Directed log2 upper bound, rounded upward |
| ---: | ---: |
| 5865 | -355.6458487 |
| 5866 | -962.2176384 |
| 5867 | -1568.6193425 |
| 5868 | -2174.8509053 |

For the remaining interval, fix the rational marker

    p = 12860472643/12983074412.

The marker bound is affine in r minus 255 log2 binom(R,r). Binomial
log-concavity makes this bound convex. Two directed endpoints therefore
cover every occupancy from 5869 through 8192:

    log2 B_marker,5869 <= -674.4750457990390728358568773584721715291965,
    log2 B_marker,8192 <= -413276.8331093622093392767812068640423988153309.

After unioning this interval, its failure probability is below
2^-663.2926. Adding the four exact-coefficient endpoints proves

    Pr[exists a quiet word with 5865<=r<=8192 and weight<=838860]
        < 2^-355.6458.

At r=5864, the same z and alpha give a positive exact-coefficient log
bound of about 251.096. This is only the boundary of the stated witness.

## What stronger residual ranks can recover

There is a useful limitation at each fixed z and alpha. In the positive
output-zero expansion, let rho be the existing residual-rank lower bound,
and let d be the actual combined residual feedback rank. Then

    rho <= d <= rho+2.

For rho=0, every residual image is zero. For rho=1, exactly one residual
tile image is nonzero and it is a line, so d=1. For rho>=2, the ambient
feedback space has dimension four, giving d-rho<=2.

Let J_j be the uncapped existing joint-moment bound. The same positive
expansion with actual ranks is the exact K_j. Term by term,

    K_j >= J_j/q^2 >= Khat_j/q^2.

The second inequality remains valid after taking the independent syndrome
cap Khat_j=min(J_j,q^-min(j,4)). Thus replacing every residual rank by
its exact value, while retaining the same averaged local moment and
fractional-power composition, can improve log2 B_r by at most

    2 alpha (P-1) log2 q.

At alpha=1/32 this maximum is 32767 bits. The four interior witnesses
below remain above one even after this largest possible improvement:

| Active rows r | Lower bound on log2 of the best rank-only envelope at this witness |
| ---: | ---: |
| 1536 | 728.7204 |
| 2048 | 19761.6598 |
| 4096 | 42247.7912 |
| 5120 | 9321.1566 |

At z=1 the syndrome-rank moment is already exact, so residual-rank
refinement cannot improve the positive r=1024 witness at all.

This limitation applies to the specified z and alpha and the same
packet-averaged composition. It does not rule out different parameters,
a different order of route averaging and fractional powers, additional
outer-row information, or another proof method.

## Reproduction and next step

The new helper is
[theory/hd_randomized_outer_quiet_exact.py](theory/hd_randomized_outer_quiet_exact.py).
Its eight tests include exhaustive fixed-subset coefficients, rational
fractional powers, marker invariance, a closed-form saddle point, exact
two-packet counts, all residual-rank profiles, and the two whole-interval
receipts. From the theory directory:

    python -B -m unittest hd_randomized_outer_quiet_exact_tests -v
    python -B hd_randomized_outer_quiet_exact.py
    python -B hd_randomized_outer_quiet_exact.py --extensions-only

The last command replays only the certified extensions. No native code,
benchmark, remote execution, or existing file was changed for this work.

The next useful step is to preserve information discarded by the constant
outer density comparison, especially the active outer row's minimum
support, or retain route-dependent joint moments before the fractional
power. The coefficient and fixed-witness rank limits above make those
changes more promising than a wider repetition of the same marker search.
