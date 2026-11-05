# Actual-state profiles with the same exact occupancy moment

2026-10-31. The previously verified run/occupancy moment also counts
actual state-support histories. Their field weights have a stronger
death factor than the availability-run envelope. This note derives
several exact branches using that same moment, proves separate
actual-state census exclusions, and records two exact point certificates.
The full 45%-distance, 40-bit goal remains open. The construction,
accepted setup distribution, route and free terminal are unchanged.

The moment and routing certificate are in
[the run/occupancy checkpoint](QUARTET_RUN_OCCUPANCY_CERTIFICATE.md).
The starting field count is the projective joint-rank argument in
[the rank reset](QUARTET_GLOBAL_RANK_RESET.md), with the per-support
weighted count stated in
[the grouped-footprint reduction](QUARTET_GROUPED_FOOTPRINT.md).
Here q=2^32, N=2^21, K=2^20, W=943718, width n=512 and T=4096.
There are 8192 outer rows and 256 independently permuted regions.
The accepted-setup density correction remains gamma<=1024/989<2.

## Define actual-state statistics separately

Let S_t be the support of the actual state Q_t, with S_0 empty. For the
update-t quartet partition, Gamma_t(U) is the union of quartets meeting
U. Define

    V_t = Gamma_t(S_(t-1)) union S_t,
    B_t = S_t \ Gamma_t(S_(t-1)),
    Ls = sum_t |S_t|,       factual = sum_t |B_t|,
    F = sum_t |V_t|,        G = sum_t |Gamma_t(S_(t-1))|
                              = F-factual.                 (1)

These are actual-state statistics. In particular, Ls is not the
availability count L, and factual is not the availability forced-birth
count. A numerical tuple with the same entries in the two conventions
does not define the same event.

Write D for the number of coordinate deaths, a=|S_T| for the terminal
support size, b for the total number of coordinate births, and i for the
number of coordinate survivals. Telescoping from S_0 empty gives

    b = D+a,       i = Ls-b,       factual <= b.             (2)

Only births outside a whole old-hit quartet contribute to factual.
Other coordinate births still contribute to b and to the field weight.
There is no death after T and no expansion of S_T into a hypothetical
next update.

Let H be the outer rows with at least 129 covered positions, h=|H|,
and let R be the retained outer parity rank. The footprint definitions
give the same row identities as before:

    C_H=F-R+128h,       C_out=R-128h,
    129h<=C_H<=256h,    0<=C_out<=128(8192-h).               (3)

RS shortening forces every input outside H to be zero. At a position
of B_t the incoming quartet is zero and the new state is nonzero, so
the input there is nonzero. Hence B lies in H.

For an output of weight at most W, every coordinate birth and death
forces a nonzero output. These events are disjoint, so

    2D+a = b+D <= W.

Also S_(t-1)\S_t is contained in Gamma_t(S_(t-1))\S_t.
Consequently every relevant profile satisfies

    d_min := max(0,factual-512) <= D
      <= d_max := min(F-Ls,floor(W/2)).                    (4)

Further necessary conditions used by the helper are

    1<=factual<=Ls<=F<=N,       factual<=C_H,
    G in 4Z,   0<=G<=N-512,   G/4<=Ls<=G+512,
    2 factual-512 <= W.                                    (5)

These follow from actual supports. Availability-only conditions are
not imported. For example, in two updates the support history
S_1=[n], S_2=empty has F=2n, Ls=n and factual=n. It violates the
availability inequality 2F-L-f<=N-n when L and f are replaced by
Ls and factual. This is a counterexample within support-history
enumeration, not a claim that the toy history is a codeword.

## The actual-state field weight

Choose 0<x<=1. For a fixed actual support history, the weighted projective field count
is bounded by

    (gamma/q) q^(-R) x^(-W) (qx)^b x^D (1+qx)^i.           (6)

It already sums the compatible output supports and amplitudes.
Put lambda_x=1+qx and t_x=qx^2/lambda_x. Substituting (2) gives

    (qx)^b x^D lambda_x^i
      = lambda_x^Ls t_x^D (qx/lambda_x)^a.

The last factor is at most one. Thus, conditional on route and
partitions, the probability associated with this fixed S is at most

    min(1, B_x t_x^D),
    B_x=(gamma/q) q^(-R) x^(-W) lambda_x^Ls.                (7)

One must not replace b by factual in (6). The latter statistic controls
the routing restriction, while b determines the birth/death algebra.

## Reuse the exact moment and routing bound

Fix the entire candidate S history before averaging over partitions.
The positive 513-state matrix from the run/occupancy checkpoint counts
arbitrary binary masks, so the same identity holds with S in place of A:

    M(s,d0,z) = sum_S E_partitions[s^Ls d0^D z^G].           (8)

The frozen parameters are

    s=75/64,   z^4=1/2,
    d0=(q+65536)/(256(q+257)) >= 2^-8.

The existing integer witness proves M<=2^22323, including the free
terminal. The actual-state helper recomputes and checks that witness
before its exact profile evaluations. Reusing a binary-mask moment does
not reuse an availability profile constraint.

For fixed S and partitions, V and B are independent of routing. The
same forced-birth slot count therefore gives

    Q = binom(F-factual,C_H-factual)
        binom(N-F,256h-C_H) / binom(8192,h)^255.             (9)

The routing probability is at most min(1,Q), after the single global H
union. It is not raised to a fractional field exponent. The expected
history count in (8) is not capped at one.

## Six branches share the moment

For theta in [0,1], min(1,u)<=u^theta. Apply this inequality to the
conditional field bound (7), before summing histories. Each branch below
bounds that conditional probability by 2^B d0^D. Let [d_-,d_+] be an
integer interval containing all D in the event being counted; initially
use [d_min,d_max] from (4).

| Branch | Residual exponent B |
|---|---|
| x=2^-8, theta=1 | -30+8W-32R+24Ls |
| x=2^-16, theta=1/2 | 8+8W-16R+8Ls |
| x=2^-16, theta=1 | 16+16W-32R+16Ls-8d_- |
| x=2^-32, theta=8/33 | ceil(8[-31+32(W-R)+Ls]/33) |
| x=2^-32, theta=1 | -31+32(W-R)+Ls-25d_- |
| Geometry only | 8d_+ |

The constants follow from finite rational inequalities. At x=2^-8,
lambda_x=2^24(1+2^-24) and

    (1+2^-24)^Ls <= 1/(1-N2^-24)=8/7<2,
    t_x < 2^-8 <= d0.

At x=2^-16, lambda_x=2^16(1+2^-16) and t_x=1/65537<2^-16.
For m=65536, the binomial series gives
(1+1/m)^m<=sum_j 1/j!. Its terms through j=4 sum to 65/24;
the tail starting at j=5 is at most 1/100, because its successive
ratios are at most 1/6. Therefore

    (1+2^-16)^N < (1631/600)^32 < (11/4)^32 < 2^47,

where the last comparison is the integer inequality 11^32<2^111.
Thus t_x^(1/2)<=d0 and t_x/d0<=2^-8.

At x=2^-32, lambda_x=2 and t_x=2^-33 exactly. Hence
t_x^(8/33)=2^-8<=d0 and t_x/d0<=2^-25. Finally,
d0^(-D)<=2^(8d_+) proves the geometry-only branch. No floating
logarithm establishes any of these constants.

Let b_Q be the exact integer ceiling of log2 Q and let Mbits be the
freshly verified moment exponent. Since 75^1024>2^6378,
the common integer profile bound is

    E = Mbits-floor(117Ls/512)+(F-factual)/4
          +min(0,b_Q)+min_(six branches) B.                 (10)

The quarter is an integer for actual profiles. The minimum is taken
inside each fixed profile; it incurs no union over branches. The moment
has already summed D and every terminal state.

There are at most 8192*K*N^3=2^96 nonzero profiles
(h,R,Ls,factual,F). Thus the whole class satisfying E<=-172 has total
failure below 2^-76. With the physical interval (4), this is a
whole-profile class theorem. With an interval refined by the exclusions
below, it is a theorem for the profile remainder outside those classes.
Their failure allowances must then be added separately.

## Prove actual-state census exclusions independently

For j>=0, a one-coordinate mask history has transfer

    H_j = [[1,1],[2^(-2j),1]].

The positive vector v=(1,2^-j)^tr satisfies
H_j v=(1+2^-j)v and 1<=2^j v entrywise. Over 512 coordinates,
including every free endpoint,

    sum_S 2^(-2jD) <= 2^(512j)(1+2^-j)^N <= 2^c_j.         (11)

Use c_0=N, c_3=361984, and otherwise
c_j=512j+ceil(3N/2^(j+1)). The general bound follows from ln2>2/3,
which one obtains by integrating 1/t>=4/(1+t)^2 on [1,2], with strict
inequality for t>1. For j=3,
9^64<2^203 gives c_3=512*3+11N/64. This census counts actual S
masks directly; availability recovery is not part of its proof.

At x=1, (7) gives a per-S bound strictly below 2^-31 q^(Ls-R), since

    gamma(1+1/q)^N <= (1024/989)/(1-N/q) < 2.

Combining it with (11), and separately using the x=2^-16 bound with
t_x<2^-16, proves the following two whole-history classes, EACH with
failure below 2^-76:

    32(R-Ls)-2jD >= c_j+45,                               (12)
    32R-16Ls+(16-2j)D >=16W+c_j+92.                        (13)

There is no 96-bit profile union in these two proofs: (11) already sums
all histories. The helper uses j=0,...,16, so the 34 cut classes together
cost at most 34*2^-76<2^-70. A separate direct 2^N mask union proves
R-Ls>=65537 below 2^-63. Its boundary slice is not implied by (12) at
j=0: the latter starts at integer gap 65538. Keep the additional
2^-63 allowance explicit.

For the remainder, intersect (4) with the integer complements of
(12)--(13). For j>0, the rank complement gives

    D >= ceil((32(R-Ls)-c_j-44)/(2j)).

For the closure complement, put B0=32R-16Ls and k=16-2j. Then

    B0+kD <=16W+c_j+91.                                  (14)

Positive k gives an upper endpoint; negative k gives a lower endpoint.
At j=8, k=0, so reject the remainder precisely when
B0>16W+c_8+91. The code handles this case explicitly. If the intersection
is empty, every D lies in an independently bounded prior class; this
does not mean the original geometric profile is empty. Otherwise use
the refined lower endpoint in the two plain closure-gain branches and
the refined upper endpoint in the geometry branch of (10).

## Checked points, controls and remaining scope

Before adding the census cuts, the exact helper recomputed all 513
moment rows and obtained Mbits=22323. Its whole-profile results were:

| (h,R,Ls,factual,F) | Exact exponent | Selected branch |
|---|---:|---|
| (2048,900000,900000,20000,1150000) | -34041 | x=2^-8 |
| (7168,917504,1310720,33280,1720320) | -134165 | x=2^-32, theta=8/33 |

These are actual-state classes. They neither assert that these profiles
are realizable nor certify availability profiles with equal numeric
labels. Exact checks passed 2307089 physical support histories, 1917
rational field-weight identities and 513 signed ceiling comparisons.
Independent mathematical and source audits agree with the underlying
branches. The final root-run checks also passed 60667 signed census
complement comparisons, including the zero coefficient at j=8.

The earlier coarse scan, before the actual-state census cuts and the
optional optimized geometry control, evaluated 448811 sampled profiles:
319895 had floating exponents at or below -172 and 128916 remained open.
Those are historical diagnostic counts, not a coverage certificate or a
description of the final remainder. A finite grid cannot prove coverage.

The refined scout, with the independently proved actual census cuts and
the optional optimized geometry control, sampled the same 448811 physical
profiles. It assigned 290674 to the actual rank-gap class and 40475 to
the actual census classes. Of the 117662 remaining sampled profiles,
58884 had floating bounds at or below -172 and 58778 remained open.
These counts describe only the sampled grid; they are not percentages
of proved coverage. The largest reported remaining value was about
+977144 at (h,R,Ls,factual,F)=(4096,786432,786432,4,1048576).

The unresolved diagnostics now emphasize small factual and R near Ls.
Two bounded searches that also mark forced births remained positive:
about +910296 at (2048,819200,819200,4,950272), with x=2^-16 and
theta=1/2, and +2166328 in the reported h=4096 test with x=2^-16 and
theta=1. These are failures of the tested upper bounds, not codewords
or counterexamples to the construction.

The implementation is
[quartet_actual_state_profiles.py](quartet_actual_state_profiles.py).
Its exact evaluator and floating scout use the same refined actual
death interval by default and distinguish physical emptiness from
prior-class coverage. The --whole-profile control omits those cuts and
reproduces the scope of the two recorded certificates. --profile takes
the five actual-state coordinates in the order shown above. --certify
freshly verifies the moment before evaluating exact routing binomials.

The optional floating control from
[quartet_occupancy_geometry.py](quartet_occupancy_geometry.py) optimizes
the closure-free count

    (1+s)^512 [1+z^4((1+s)^4-1)]^(128(T-1)).

This identity holds for arbitrary S masks for every fixed partition
sequence. The scout may take the smaller of this optimized count times
the routing bound and its fixed-moment diagnostic. The exact certificate
does not consume floating results. The geometry helper's exact tests
passed 1653 physical coefficient comparisons, 35 rational moments and
53 endpoint identities; 89 separate floating optimizer checks passed.

The scout also retains the previously proved exclusions depending only
on F or R. It imports no availability rank-gap or census label. All
class allowances remain separate when assembling a complete proof.
The next step is to sharpen the global bound on the remaining small-birth
profiles with R near Ls, retaining the route quotas and free endpoint,
then certify a covering family with a controlled total failure budget.
The full goal is still open; these point certificates and finite scans
establish no obstruction.
