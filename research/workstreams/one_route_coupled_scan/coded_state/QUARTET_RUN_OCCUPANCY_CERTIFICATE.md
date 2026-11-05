# An exact run, quartet-occupancy and routing certificate

2026-10-31. This checkpoint proves a new whole-history class bound for
the fixed quartet construction. It combines the availability-run field
bound, the concentration of available positions inside random quartets,
and the regional routing quota. All closure counts are summed in one
moment. The selected point and its half-width-512 profile box satisfy the
new class criterion by exact integer verification. The full 45%-distance,
40-bit certificate remains open.

The subsequent [actual-state profile checkpoint](QUARTET_ACTUAL_STATE_PROFILES.md)
reuses this binary-mask moment with a different field weight and distinct
statistics Ls and factual. Its actual-state class labels must not be
identified with availability labels merely because their values coincide.

The construction, accepted setup distribution, outer code, route and free
terminal are unchanged. Local seven-type refinements remain parked. The
probability space and projective field count are those of the
[global rank reset](QUARTET_GLOBAL_RANK_RESET.md) and
[grouped-footprint reduction](QUARTET_GROUPED_FOOTPRINT.md). The recovery
identity comes from the [availability census](QUARTET_AVAILABILITY_CENSUS.md),
and the field estimate below uses the
[all-run closure lemma](QUARTET_ALL_RUN_CLOSURES.md).

## Objects and probability space

Work over J=GF(2^32), with q=2^32, K=2^20, N=2^21 and W=943718.
The event to be bounded is the existence of a nonzero J-codeword with
output weight at most W. There are 8192 RS[256,128] outer rows and 256
regions. Every region independently permutes one coordinate from each
outer row into its 8192 physical slots. Each of the T=4096 updates uses
512 state coordinates and an independent uniform partition into quartets.
Local field maps are sampled once, subject to the existing acceptance
rule; the global density correction is gamma<=1024/989<2. All this setup
is fixed for every message.

For an actual state-support history S, let V_t be its covered footprint
at update t. Write Gamma_t(U) for the union of update-t quartets meeting
U. Define A_t to be the positions in V_t whose next quartet is wholly
covered by V_(t+1), with A_0 empty and A_T=V_T. The canonical identities are

    V_t = Gamma_t(A_(t-1)) union A_t,
    B_t = A_t \ Gamma_t(A_(t-1)).                       (1)

Here B_t is the canonical forced-birth set. Introduce the statistics

    F = sum_t |V_t|,            L = sum_t |A_t|,
    f = sum_t |B_t|,            D = sum_t |A_(t-1) \ A_t|,
    G = sum_t |Gamma_t(A_(t-1))| = F-f.                (2)

Thus D counts closed availability runs. There is no closure after T and
no additional quartet expansion of A_T. In particular, the terminal
state remains free. G is divisible by four.

Let H be the outer rows containing at least 129 covered positions and
h=|H|. Write C_H and C_out for covered positions inside and outside H.
The retained outer parity rank is

    R = 128h+C_out,       C_H = F-R+128h.              (3)

RS shortening forces every input outside H to be zero. In particular,
every forced birth lies in H. A nonzero profile has h,R,L,f,F all positive.

## Correct the three earlier empty boxes

Each eligible row contains between 129 and 256 covered positions.
Consequently

    129h <= C_H <= 256h,       h <= F-R <= 128h.       (4)

Three earlier diagnostic boxes had h=2048, half-width 1024 in every
listed coordinate, and centers

| Center (R,L,f,F) |
|---|
| (900000,1000000,40000,1200000) |
| (900000,1100000,20000,1200000) |
| (900000,900000,20000,1200000) |

All three are empty: throughout each box,

    F-R >= 1198976-901024 = 297952 > 262144 = 128h.

Their earlier floating bounds therefore provide no evidence of a
remaining low-birth obstruction. The corresponding correction is also
recorded in the [old box checkpoint](QUARTET_ALL_RUN_CLOSURES.md).
Replacing F=1200000 by F=1150000 passes this aggregate screen. Passing
necessary aggregate inequalities does not establish that a profile is
realizable. The theorem below bounds every realizable history in the
replacement profile and box, without asserting their existence.

## Sum runs and quartet concentration together

Fix a complete candidate availability history A before averaging over
partitions. Formula (1) reconstructs its candidate footprint. Apply the
canonical field estimate only when this reconstruction is canonical;
then drop the canonical indicator in a nonnegative sum. This order
avoids conditioning a fresh partition on an availability mask that was
itself recovered from that partition.

For a fixed old mask of size u, the probability that a uniform quartet
partition has exactly k hit quartets is

    P_(u,k) = binom(128,k)
      [X^u] ((1+X)^4-1)^k / binom(512,u).             (5)

One may obtain (5) by uniformly assigning the u selected positions to
512 slots in a fixed quartet partition. Choose the k hit quartets and
choose a nonempty subset in each. For positive markers s,t,z, put

    M_u(z) = sum_k P_(u,k) z^(4k),
    T_(u,v)(s,t) = s^v [y^v](t+y)^u(1+y)^(512-u),
    K(s,t,z) = diag(M_u(z)) T(s,t).                    (6)

The coefficient in T sums all new masks of size v. Each old position
that is absent from the new mask receives t; every new position receives
s. Partition averaging depends only on u. The exact unrestricted moment is

    M(s,t,z) := e_0 K(s,t,z)^4096 1
      = sum_(A_1,...,A_T) E_partitions[s^L t^D z^G].  (7)

The final vector is all ones. The sum includes arbitrary availability
masks, so it bounds the canonical subset while retaining both run
closures and quartet concentration. It is an expected count, not a
probability, and must not be capped at one. No union over partitions is
introduced.

## Keep the eligible-row quota and forced births

Fix A and the partitions, hence V and B, and fix an eligible row set H
of size h. Within each region the H-labelled slots form a uniform
h-subset of 8192 slots. The 256 choices are independent, giving exactly
binom(8192,h)^256 equally likely slot sets.

Put S_H=256h and C=C_H=F-R+128h. A favorable slot set contains every
forced birth, contains C covered positions, and contains S_H-C uncovered
positions. Dropping its individual regional quotas leaves at most

    binom(F-f,C-f) binom(N-F,S_H-C)

choices. Union over the single global H, of which there are
binom(8192,h) choices, gives the routing bound

    Q(h,R,f,F) =
      binom(F-f,C-f) binom(N-F,S_H-C)
      / binom(8192,h)^255.                            (8)

This counts all eligible row identities once for the complete code.
Imposing the actual row-degree classification could only reduce the
count. Thus the routing probability is at most min(1,Q). This cap applies
after the H union; it does not cap the expected number of availability
histories in (7). Formula (8) holds uniformly for each reconstructed
canonical footprint with the stated profile, so it can be combined with
the partition expectation in (7).

## Match the field closure weight exactly

For any 0<x<=1, define c=(q+1)x, lambda=1+c and

    d(x) = max(min(1,lambda/(q+1)), (1+qx^2)/lambda).

The all-run lemma bounds the weighted state sum by
c^f lambda^(L-f) d(x)^D. Therefore, conditional on route and partitions,
the probability of a bad output with a fixed canonical footprint is at
most

    (gamma/q) q^(-R) x^(-W)
      c^f lambda^(L-f) d(x)^D.                       (9)

All possible state and output amplitudes in that footprint have already
been summed. There is no additional state-support union. Closed runs
receive their factor once; no separate forced-closure bonus is added.

Freeze the parameters

    x = 1/256,             s = 75/64,       z^4 = 1/2,
    lambda = (q+257)/256,
    t = d(x) = (q+65536)/(256(q+257)).                (10)

For a fixed profile (h,R,L,f,F), equations (7)--(10) give the bound

    (gamma/q) q^(-R) x^(-W) c^f lambda^(L-f)
      s^(-L) z^(-(F-f)) min(1,Q) M(s,t,z).            (11)

The factor t^D in the moment is exactly the field closure factor. Thus
(11) sums every D, including every terminal endpoint, without a D grid
or a D union allowance.

## Exact moment certificate

At (10), every entry of K is rational: z appears only as z^(4k)=2^(-k).
The verifier proposes a positive integer vector v and checks every one
of the 513 inequalities

    K v <= 32v.                                     (12)

The checked free-endpoint comparison is

    v_0 / min_u v_u <= 2^1843.

Positivity now gives

    M(s,t,z) <= 32^4096 v_0/min_u v_u <= 2^22323.     (13)

The proposal uses floating power iteration, but the proposal itself
certifies nothing. Every row of (12), the endpoint ratio, and the final
exponent are verified using integers. The certificate is recomputed by
the profile helper; the value 22323 is not an unchecked input constant.

For completeness, the exact matrix action avoids constructing a dense
rational matrix. Write t=a/b and s=c_0/d_0. For a positive integer vector v,
first compute

    N_j = sum_l binom(512-j,l) c_0^l d_0^(512-j-l) v_(j+l).

Then (Tv)_u is the integer numerator

    sum_j binom(u,j) a^(u-j) [c_0(b-a)]^j N_j

divided by d_0^512 b^u. This follows by expanding
(t+y)^u=[t(1+y)+(1-t)y]^u. Every summand is nonnegative.
The occupancy factor M_u has denominator binom(512,u)2^128 and integer
numerator

    sum_k binom(128,k) 2^(128-k)
      [X^u] ((1+X)^4-1)^k.

Cross multiplication therefore verifies (12) with exact integers.

## An integer criterion for an entire class

The remaining constants also admit integer checks. Since L<=N and c<=lambda,

    c^f lambda^(L-f) <= lambda^L
      = 2^(24L)(1+257/q)^L < 2^(24L+1).

Indeed, the binomial series gives
(1+257/q)^N <= 1/(1-257N/q)=2048/1791<2.
Also 75^1024>2^6378 proves s^1024>2^234, hence

    s^(-L) <= 2^(-floor(117L/512)).

Define b_Q to be the smallest integer with Q<=2^b_Q. It is computed
from exact integer binomials by bit lengths and one shifted comparison,
without logarithms. Combining gamma<2, x^(-W)=2^(8W), (11) and (13)
gives a strict per-profile bound 2^E, where

    E(h,R,L,f,F) = 22323-30+8W-32R
      +24L-floor(117L/512)+(F-f)/4+min(0,b_Q).         (14)

There are at most

    8192 * K * N^3 = 2^(13+20+3*21) = 2^96

nonzero integer profiles (h,R,L,f,F). Consequently the entire class

    E(h,R,L,f,F) <= -172                             (15)

has total failure probability strictly below 2^-76. This union covers
all five profile coordinates, including h. D has already been summed in
(7). No additional D, R-range, state-support or eligible-set union is
needed. Other independently budgeted class theorems keep their own
failure allowances when combined with this one.

This is a J-field class theorem. Its profile label need not be preserved
by extension-field basis projection. A proof covering every J-word would
transfer the distance guarantee in the established way; (15) alone does
not cover every remaining profile.

## A uniform box certificate

A finite box can be covered by a uniform per-profile upper bound, then
included in the same 2^96 profile union. This is different from bounding
one routing event that unions a range of R: such an event requires
summing its possible C values. Here each profile has its own bound (14),
and the global profile union already performs that summation.

Fix h, write a_0=128h and S_H=256h, and let each box coordinate range over
a closed integer interval. First set f=f_lo: Q is nonincreasing in f,
and z^f is also nonincreasing. The dropped factor (c/lambda)^f would
only strengthen this reduction. Set L=L_hi, since
24L-floor(117L/512) is increasing.

For the remaining (R,F) rectangle, verify binomial support throughout:

    C_min = F_lo-R_hi+a_0 >= f_lo,
    C_max = F_hi-R_lo+a_0 <= S_H,
    R_lo >= a_0,       R_hi-a_0 <= N-S_H,
    0 <= F_lo <= F_hi <= N.                          (16)

With C=F-R+a_0 and f=f_lo, adjacent raw routing ratios are

    Q(F,R+1,f)/Q(F,R,f)
      = (C-f)(N-F-S_H+C) / ((F-C+1)(S_H-C+1)),        (17)

    Q(F+1,R,f)/Q(F,R,f)
      = (F+1-f)(S_H-C) / ((C+1-f)(N-F)).              (18)

The maximum of (17) over the rectangle is attained at (F_hi,R_lo).
Checking that ratio <=q proves that q^(-R)Q decreases with R. The maximum
of (18) is attained at (F_lo,R_hi). If its numerator and denominator are
n_F and d_F, the exact check 2 n_F^4<=d_F^4 proves that the ratio is at
most z; hence z^(-F)Q decreases with F. A varying F interval has F_lo<N,
so this denominator is positive. For a fixed F no F-increment check is
needed. The reduction in f must precede these ratio checks.

These inequalities place the upper bound at
(R_lo,L_hi,f_lo,F_lo). The routing cap needs care: reduction in R remains
valid through a capped interval, but reduction in F through a plateau
Q>=1 need not. If the corner Q<=1, decreasing F reaches the largest raw
Q after the R reduction, so the cap causes no problem. Otherwise the
helper retains the uncapped corner bound as a safe fallback. For a
noncanonical corner it rounds (F-f)/4 upward only after applying the
continuous monotonicity bound. Canonical profiles themselves have
F-f divisible by four.

## Verified results and scope

The exact checks and full certificate were run serially on Peach.
The recomputed moment has 513 verified rows, multiplier 32,
terminal allowance 1843 bits and total moment allowance 22323 bits.
The resulting profile certificates are:

| Profile or box, h=2048 | Exact integer upper exponent |
|---|---:|
| (R,L,f,F)=(900000,900000,20000,1150000) | -34041 |
| Every profile within 512 of those four coordinates | -5285 |

For the point, C_H=512144, S_H-C_H=12144 and b_Q=-482914.
For the box, the checked upper corner is

    (R,L,f,F)=(899488,900512,19488,1149488).

The full support rectangle and both ratio tests pass, and the corner's
raw routing factor is at most one. Both results are far below the -172
threshold in (15). They include every D and all free terminal states.
They prove that any realizable history in these ranges belongs to the
new certified class; they do not establish realizability.

The implementation and verification sources are:

- [quartet_run_occupancy.py](quartet_run_occupancy.py): the exact moment
  definition, positive numerical assembly, and physical small-instance checks.
- [quartet_run_certificate.py](quartet_run_certificate.py): integer
  occupancy coefficients, the two exact positive transforms, all-row
  verification and the free-terminal comparison.
- [quartet_profile_certificate.py](quartet_profile_certificate.py): exact
  constants, forced-birth routing, profile and box exponents, and toy
  checks of the ratio and corner inequalities.
- [quartet_profile_entropy.py](quartet_profile_entropy.py): earlier
  floating diagnostics and independent finite routing/census checks;
  its floating results are not inputs to (14).

The guarded runner injects the occupancy and run-certificate definitions
before invoking the profile helper with --checks or --certify. The
default --certify run regenerates the proposed vector and verifies it;
it does not read a saved witness or an earlier numerical result.
Small-instance checks compare physical masks and partitions, exact
matrix actions, free endpoints, quota counts with forced births,
integer exponent rounding, binomial recurrences and box corners.
Independent mathematical and source audits agree with the certificate.
The numerical jobs used one thread, a 175-second guard, a 2-GiB address
space limit and a 20-GiB free-space check. They wrote no raw data files.

The next step is to combine (15) with the existing class exclusions and
screen the remaining (h,R,L,f,F) ranges before evaluating new moments.
A complete proof must cover that entire remainder with a controlled
total failure budget. This checkpoint supplies a new exact class theorem
and resolves the replacement test box; it does not close the full
45%/40-bit goal or identify a genuine obstruction.
