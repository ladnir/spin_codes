# Exact death slices of the hard actual-state profile

2026-10-31. Two tails of the fixed hard profile admit rigorous bounds:
all D<=32768 have probability below 2^-14587, and all D>=196608 have
probability below 2^-26713. The low slice can be extended to D<=34569
with exponent -179. These bounds already sum their death intervals.
They leave the central interval open and do not certify the whole
profile or the full 45%/40-bit goal. The construction is unchanged.

The profile is

    (h,R,Ls,f,F) = (4096,786432,786432,4,1048576).

All statistics are actual-state statistics as defined in
[the actual-state note](QUARTET_ACTUAL_STATE_PROFILES.md). Write L=Ls,
N=2^21, q=2^32, W=943718, and

    G=F-f=1048572,       C_H=786432,       C_out=262144.

The physical death interval is 0<=D<=F-L=262144. The terminal state
is free throughout. An independent audit passed the derivation below.

## The existing census complements remove nothing here

The implemented rank exclusions have left side
32(R-L)-2jD=-2jD<=0, below their positive threshold c_j+45.
For the closure exclusions, maximize

    32R-16L+(16-2j)D-16W

over the physical interval. For j=0,1,2,3 its maxima and thresholds are

| j | Maximum left side | c_j+92 |
|---|---:|---:|
| 0 | 1677728 | 2097244 |
| 1 | 1153440 | 1573468 |
| 2 | 629152 | 787548 |
| 3 | 104864 | 362076 |

For 4<=j<8 the maximum left side is already negative. For j>=8 it is
maximized at D=0 and is again negative. Thus all 34 implemented census
exclusions fail throughout this profile. The separate rank-gap exclusion
also fails because R-L=0. The coverage complement F<=1736703 is automatic.
These complements remain useful globally, but inserting their multipliers
cannot sharpen this exact profile by excluding any of its death values.

## A uniform exact routing bound

The forced-birth routing ratio, including the eligible-set union, is

    Q = C(F-f,C_H-f) C(N-F,256h-C_H) / C(8192,h)^255.

Here N-F=F, C_H=3F/4, and 256h-C_H=F/4. Dropping the requirement that
the four forced births lie in H only increases the numerator, so

    Q <= C(F,F/4)^2 / C(8192,4096)^255.                  (1)

The binomial theorem at 1/3 gives
C(F,F/4)<=4^F/3^(3F/4). Since 3^12>2^19 and F is divisible by eight,
the numerator in (1) is below 2^(13F/8)=2^1703936.
The central binomial coefficient exceeds 2^8192/8193>2^8178, so its
255th power exceeds 2^2085390. Consequently

    Q < 2^-381454.                                      (2)

This bound is uniform in D and in the compatible partition geometry.
It is conservative; no floating binomial estimate is used.

## Low deaths: reuse the verified moment

At x=1 the fixed-support rank bound is

    P_S < 2^-31 q^(L-R) = 2^-31.                       (3)

For example, its correction is bounded using
gamma/(1-N/q)<=(1024/989)(2048/2047)<2.

Invoke the existing [integer-verified run/occupancy theorem](QUARTET_RUN_OCCUPANCY_CERTIFICATE.md):

    M(s,t,z) <= 2^22323,
    s=75/64,
    t=(q+65536)/(256(q+257)) >= 2^-8,
    z^4=1/2.

It counts arbitrary state masks, including the free terminal, so it
applies to actual states. It sums s^L t^D z^G over all histories and
partitions. Since 75^1024>2^6378 and L/1024=768,

    s^L > 2^179712,       z^-G=2^262143.

For any D<=d, t^D>=2^-8d. The expected number of histories in this
entire low-death slice is therefore at most

    2^(22323-179712+262143+8d) = 2^(104754+8d).

The specified f is relaxed only in this history count; (2) still bounds
routing on the same fixed profile. Combining (2)-(3) gives

    Pr[bad output in this profile with D<=d]
      < 2^(8d-276731).                                 (4)

At d=32768 the exponent is -14587. At d=34569 it is -179, still below
the usual -172 profile-class threshold; d=34570 would give -171 under
these particular conservative constants. No union over the individual
death counts is needed.

## High deaths: exact output coefficients and a geometry count

The [optional-output coefficient bound](QUARTET_MDS_REFINEMENT_LIMIT.md)
gives, uniformly over every terminal state,

    P_S <= min(1,2^-31 q^(L-R-Delta) C(L-D,Delta)),
    Delta=max(0,L+D-W).                                 (5)

The geometric-series correction gamma/(1-N/q)<2 is included in (5).
For this profile D_0=W-L=157286. If D>=196608 then Delta>=39322 and
n=L-D<=20Delta. The binomial theorem and 17^10<2*16^10 imply

    C(n,Delta) <= (17/16)^n 16^Delta < 2^(6Delta).

Thus the field probability throughout the high-death slice is strictly
below

    2^(-31-26Delta) <= 2^-1022403.                       (6)

To count the histories, the exact geometry polynomial with no death
marker is

    (1+s)^512 [1+z^4((1+s)^4-1)]^J,
    J=(N-512)/4=524160.

For G=4k with k=262143, its fixed-(L,G) coefficient is bounded at s=3 by

    C(J,k) 255^k 4^512 / 3^L.

Use C(J,k)<=2^J, 255<2^8, and 3^L>2^(19L/12). This gives the uniform
history count

    count(L,G) < 2^(J+8k+1024-19L/12) = 2^1377144.     (7)

The count in (7) includes all death counts and all f; it therefore also
bounds the specified high-death slice. Multiplying (2), (6) and (7),

    Pr[bad output in this profile with D>=196608]
      < 2^(1377144-381454-1022403) = 2^-26713.            (8)

Again there is no death-count union: one whole-slice count is multiplied
by a uniform field bound.

## Scope, verification and the remaining interval

Equations (4) and (8) concern this fixed profile. They are not numerical
coverage claims for neighboring profiles. A global class can use the
usual 96-bit union over (h,R,L,f,F) only after proving the corresponding
criterion for every included profile. For bookkeeping, even the larger
low slice (-179) and the high slice (-26713), summed with a factor-two
allowance and then charged 96 bits, give an exponent below -82. This
does not assert that the same numeric inequalities hold elsewhere.

The small self-contained helper
[quartet_hard_death_slices.py](quartet_hard_death_slices.py) checks every
finite integer/rational comparison above and all 34 census maxima.
Its default --checks mode invokes the previously proved moment theorem;
it does not claim to reverify the 513-row witness. The optional --certify
mode, with the existing run/occupancy and certificate definitions
injected, proposes a vector and freshly verifies every row with integers.
It returns a witness hash without writing the witness to disk.

The supporting source files at this checkpoint have SHA256 values:

| Source | SHA256 |
|---|---|
| quartet_run_occupancy.py | 3e066747f42a56fa25bda16b42c8ff53677419d3ecaf550a9985b5d46d6c285e |
| quartet_run_certificate.py | ae08f37541905bc68eaf63e824f116ed7c546f686a3c9cb747044c0a706c5dff |

These identify the supporting source, not a saved numerical witness.
The prior theorem verified all 513 rows with multiplier 32 and included
a terminal allowance of 1843 bits, giving 4096*5+1843=22323.

The root agent subsequently ran --certify serially. All 34 census
comparisons and the fresh 513-row integer certificate passed in
3.29 seconds, reproducing the 22323-bit bound and every slice exponent
above. The fresh witness SHA256 was

    8b406fcb326fe9faf0d57d50e8e541940cd2ebdcc4671c967e78c5b9525c24e1

The hash encodes the comma-separated positive decimal witness entries,
with no trailing comma. The helper returned only compact metadata and
the hash; it wrote no witness file.

The remaining interval after the larger low slice is

    34570 <= D <= 196607.

For D<=157286 inside this interval, the output cutoff does not bind:
every output compatible with S has weight at most L+D<=W. Sharpening
only the output-weight envelope cannot resolve that part of the support
count. The next step must retain the joint (F,L,D,f) count and routing
constraints, or prove additional field-rank information. No codeword
witness or failure of the construction is implied by the open interval.
