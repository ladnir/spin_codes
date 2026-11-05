# An exact geometry-only profile class

2026-10-31. A direct coefficient bound now certifies another global
profile class without an eigenvector or field-weight moment. The selected
profile has exact integer exponent -72111 and satisfies the class
criterion. The construction is unchanged, and the full 45%-distance,
40-bit certificate remains open.

This bound applies separately to actual state-support histories and to
canonical availability histories. Their statistics must retain their
respective meanings. It uses only conditions common to both descriptions;
no availability-specific separation or closure condition is imposed.

## Count the geometry, including the terminal mask

Let n=512, T=4096, N=nT, and J=n(T-1)/4=524160. Fix any sequence of
quartet partitions. For masks S_0 empty, S_1,...,S_T, define

    L = sum_i |S_i|,
    G = sum_i |Gamma_i(S_(i-1))|.

Each nonterminal mask S_i is counted through the next partition. A
quartet contributes 1 when its mask is empty and
z^4((1+s)^4-1) otherwise. The terminal mask contributes (1+s)^n, with
no further expansion. Thus the exact counting polynomial, for every fixed
partition sequence, is

    (1+s)^n [1+z^4((1+s)^4-1)]^J.                         (1)

For G=4k, the number of histories with statistics (L,G) equals

    C(J,k) [s^L] ((1+s)^4-1)^k (1+s)^n.                (2)

Its support is 0<=k<=J and k<=L<=G+n. For any positive integers a,b,
evaluation at s=a/b bounds (2) by

    C(J,k) ((a+b)^4-b^4)^k (a+b)^n
    ---------------------------------.                (3)
              a^L b^(G+n-L)

Equation (3) uses integers only. It sums all coordinate death counts and
retains the free terminal. Since (1) is independent of the fixed
partitions, it also counts histories averaged over their random choice.

## Combine with routing on the same profile

For a profile (h,R,L,f,F), use G=F-f and

    C_H = F-R+128h.

The [existing forced-birth routing bound](QUARTET_RUN_OCCUPANCY_CERTIFICATE.md)
is the following ratio, capped at one:

    Q = C(F-f,C_H-f) C(N-F,256h-C_H) / C(8192,h)^255.    (4)

Here f counts actual forced births for actual-state profiles and
canonical forced births for availability profiles. The history count
relaxes the specified birth count and any other consistency conditions;
(4) still uses the specified f for the same profile. The conditional
field-event probability is bounded by one, so no field-density correction
or recomputed eigenvector witness enters this class criterion.

Let g(L,G;a,b) be the smallest integer upper bound on the base-two
logarithm of (3), and let r be the corresponding integer bound on (4).
The implementation chooses

    E(h,R,L,f,F) = min_(a,b) g(L,F-f;a,b) + min(0,r),    (5)

over the fixed markers s in {1/2,1,3/2,2,4,8}. Each marker bounds the same
profile, so this minimum costs no family union. Exact powers, binomial
coefficients and integer comparisons establish every ceiling in (5).

There are strictly fewer than

    8192 * 2^20 * (2^21)^3 = 2^96

nonzero physical profiles: h<=8192, R<=2^20 and 1<=f<=L<=F<=2^21.
The restriction f<=L already makes the count strict. Consequently all
profiles satisfying

    E(h,R,L,f,F) <= -172                               (6)

form a class with total failure probability below 2^-76. Actual-state and
availability applications give separate valid class bounds. Combining
them with other exclusions must retain the project's probability
accounting and the correct profile definitions.

## Exact checks and the selected profile

[quartet_geometry_profile_certificate.py](quartet_geometry_profile_certificate.py)
reuses the generic integer routing and binomial helpers from
quartet_profile_certificate.py. It needs no run-moment witness. All work
was run serially under the existing 175-second, 2-GiB and
20-GiB-free-space guards.

The small checks passed 87 coefficient comparisons against independent
labeled-mask censuses, 522 rational coefficient bounds, 522 minimal
integer logarithm ceilings, two free-terminal checks, and six rational
marker-rescaling checks. An independent read-only audit passed the
coefficient identity, integer ceilings, routing cap, marker minimum,
common support screens, global union count, and the uniform consequence
below.

For

    (h,R,L,f,F) = (3584,524288,526950,4,954368),
    G = 954364,    k = 238591,

the exact certificate chose s=1 and returned

    geometry count exponent <=  1453770,
    routing exponent        <= -1525881,
    combined exponent       <=   -72111.

The run took 1.73 seconds. This profile satisfies (6). Passing the
necessary support screens does not assert that the profile is realizable;
the result is a valid bound in either case.

The same certificate covers all admissible L and f at the fixed triple

    (h,R,F) = (3584,524288,954368),
    f >= 4,   f congruent to F modulo 4.                 (7)

To see this, retain the winning marker s=1. Its geometry bound is
C(J,k)15^k2^512, independent of L. Increasing f by four decreases k by one
and multiplies that bound by

    k / (15(J-k+1)).                                    (8)

At f=4, k=238591 and J=524160, so k<=15(J-k+1). This inequality persists
as f increases. At fixed h,R,F, C_H is also fixed. Along the same step,
the uncapped routing ratio is multiplied by

    product_(i=0)^3 (C_H-f-i)/(F-f-i) <= 1.             (9)

All factors in (9) are defined for a feasible next profile; capping the
routing ratio at one preserves this monotonicity. Hence both integer
ceilings are nonincreasing, and every admissible profile in (7) has
combined exponent at most -72111, for every admissible L. Conditions
include f<=L<=F and (F-f)/4<=L<=F-f+512, as well as the common routing
support conditions. Because F is divisible by four, every nonzero
admissible f at this fixed triple is already covered by (7).

This uniform family is a subset of the same global profile class (6).
The existing 96-bit union already counts every L and f; no second union
penalty is added for the consequence (7).

Next apply this exact geometry class alongside the existing field-weight
classes when covering the common remaining profile ranges. The selected
point does not establish coverage of that complement or complete the
full-code goal.
