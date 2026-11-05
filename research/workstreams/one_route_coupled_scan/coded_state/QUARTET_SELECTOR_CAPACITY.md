# Global rank from the capacity of the available quartets

2026-10-31. The hard profile that defeated the previous state-history
moments is excluded by a simpler global routing event. The argument
bounds how many eligible outer positions ANY collection of live
quartets can contain. It does not enumerate state histories. At the
profile (h,R,Ls,f,F)=(4096,786432,786432,4,1048576), the resulting
failure bound is strictly below 2^-4094, for every death count and
every free terminal state. The full 45%/40-bit certificate remains open.

The encoder, once-sampled setup, outer code and routing are unchanged.
Use the actual-state definitions from
[the global rank proof](QUARTET_GLOBAL_RANK_RESET.md) and
[the actual-state profile note](QUARTET_ACTUAL_STATE_PROFILES.md).
There are M=8192 outer rows, 256 independently permuted coordinate
regions, N=2^21 positions, and J=N/4=524288 quartets. These are the
quartets in actual updates, with no added terminal update.

## Why this avoids the previous overcount

Fix a candidate eligible-row set H of size h. Condition on all quartet
partitions before exposing the route. Every region contains a uniform
h-subset of H-labelled positions. Let j_Q be the number of these labels
in quartet Q. For r in {0,1,2,3}, define the excess capacity

    Z_(H,r) = sum_Q max(j_Q-r,0).

For EVERY collection A of k quartets, simultaneously,

    sum_(Q in A) j_Q <= r k + Z_(H,r).                 (1)

This is just j<=r+(j-r)_+ summed over selected quartets, followed by
including the nonnegative excess from the unselected quartets. The
collection may depend on the entire sampled setup or on a message.
No union over collections is required.

For a genuine state history, its live quartets number k=G/4, where
G=F-f. Forced births outside these quartets lie in H, because outer
rows outside H shorten to zero. Write C_out=R-128h. The number of
H-labelled positions in live quartets is exactly G-C_out. Thus (1)
implies the deterministic condition

    z_r := (1-r/4)G-C_out <= Z_(H,r).                  (2)

The same inequality would follow if some declared births were outside
H, by upper-bounding their eligible contribution by f. Equation (2)
depends on G,R,h, but not on a particular state history, its death count,
or its field amplitudes. This is the reason the large history-counting
factor disappears.

For example, at h=4096, half the slots carry H labels. Selecting half
the quartets preferentially selects those with many H labels. But this
selection has a finite capacity: it cannot arbitrarily place three
quarters of all covered slots inside H. The hard profile requires an
unusually large total excess Z_(H,2).

## Use the exact regional quotas without a conditioning penalty

For a single quartet in a region, j_Q has distribution

    Pr[j_Q=j] = C(h,j) C(M-h,4-j) / C(M,4).             (3)

Although the quartet counts in that region are dependent, products of
nonnegative increasing functions of disjoint quartet counts satisfy

    E product_Q g_Q(j_Q) <= product_Q E g_Q(j_Q).       (4)

Here is a direct proof in the probability space we use. Condition on
the count j in the first quartet. The remaining slots contain a uniform
(h-j)-subset. Generate these subsets as nested prefixes of one uniform
ordering. The conditional expectation of the product of the remaining
increasing factors decreases with j, whereas the first factor increases.
Opposite-monotone functions of one random variable have nonpositive
covariance; apply this to separate the first factor. The remaining
expectation is unconditional under the original uniform h-subset of M
slots. In the next separation, previously unused quartet slots stay in
that sampling population. Repeating proves (4) with the SAME marginal
distribution (3) in every factor.

Choose beta>1 and let

    B_(h,r)(beta) = sum_(j=0)^4
        C(h,j) C(M-h,4-j) beta^max(j-r,0) / C(M,4).

The function beta^max(j-r,0) is nonnegative and increasing. Applying
(4) within each region and independence between regions gives

    E beta^Z_(H,r) <= B_(h,r)(beta)^J.

Union over the C(M,h) candidate sets H only. For every integer z,

    Pr[exists H, |H|=h, Z_(H,r)>=z]
      <= C(M,h) B_(h,r)(beta)^J beta^(-z).              (5)

There is NO regional quota divisor, state-history union, death-count
union, or profile union in (5). The actual quotas have already been
used in (3)-(4). The inequality holds conditional on arbitrary quartet
partitions and hence after averaging them. Setup acceptance concerns
the independent field maps and does not change this geometry law.

## An exact certificate for the whole hard profile

Set h=4096, r=2 and beta=3/2. Its hypergeometric moment is smaller than
77/64. For an explicit check, write M=2m and D=(2m-1)(2m-3). Then

    B_(m,2)(3/2) = 1+(m-2)(13m-15)/(16D),
    77/64-B_(m,2)(3/2) = (60m-81)/(64D) > 0.           (6)

The hard profile requires

    z_2 = (1048576-4)/2-(786432-128*4096)
        = 262142 = J/2-2.

Using C(8192,4096)<=2^8192, equation (5) is strictly below

    2^8192 (3/2)^2 (5929/6144)^(J/2).

The integer inequality 8*5929^64<6144^64 implies that the last power
is below 2^-12288, since J/2=262144 is a multiple of 64. Also (3/2)^2<4.
The full bound is consequently strictly below

    2^(8192+2-12288) = 2^-4094.                         (7)

One elementary derivation of that integer inequality is
30*5929<29*6144 and

    sum_(i=0)^4 C(64,i)/29^i
      = 5807265/707281 > 8.

Thus (30/29)^64>8 and 5929/6144<29/30. All constants can be checked
with small integer or rational arithmetic; no floating eigenvalue
or sampled spectrum is needed.

More generally, the same event excludes EVERY h=4096 profile satisfying
G/2-C_out>=262142, simultaneously, with the same failure bound. The
formerly unresolved middle death range is included. This supersedes
the need to close that profile using the
[death-slice bounds](QUARTET_HARD_PROFILE_DEATH_SLICES.md) or the failed
exact-birth and regional-bin moments. Those results remain preserved.

## Turn this into a rank bound for all remaining profiles

For each 1<=h<M and r in {1,2,3}, choose an integer threshold z_(h,r)
and a beta>1 that certify the right side of (5) at most 2^-92. There
are fewer than 2^15 such pairs, so all their bad geometry events have
total probability below 2^-77. Choosing the best beta deterministically
for one threshold does not introduce an additional union.

On the complement of these events, every genuine profile satisfies

    R >= 128h + max_(r=1,2,3)
                     ((1-r/4)G-z_(h,r)+1).             (8)

This is a simultaneous piecewise-linear lower bound on the retained
outer rank. It can be intersected with the previously proved rank and
field-count classes. There is no union over Ls,D,f,F or individual
messages at this step. Thresholds above J(4-r) give empty bad events;
h=M is deterministic and can be handled directly.

The exact complement identity

    Z_(H,r) = 256h-rJ+Z_(H^c,4-r)

can also reflect thresholds between h and M-h. It follows pointwise
from (j-r)_+ = j-r+(r-j)_+ and the fixed total 256h labels.

## A sharper bound retaining the quota in the moment

Negative association gives a short proof, but retaining the complete
regional quota can improve the numerical rank planes. Define

    Q_r(x,beta)=sum_(j=0)^4 C(4,j) x^j beta^max(j-r,0).

For one region the exact moment is

    E beta^Z = [x^h] Q_r(x,beta)^2048 / C(M,h).

For every positive rational x, positivity of the coefficients bounds
this coefficient by Q_r(x,beta)^2048/x^h. Independence between regions
and the single union over H therefore give

    Pr[exists H, |H|=h, Z_(H,r)>=z]
      <= Q_r(x,beta)^J /
         (x^(256h) C(M,h)^255 beta^z).                  (9)

This is still an all-history setup event. The minimum of the bounds
(5) and (9) is valid for the same threshold, without another union.
The implementation chooses rational x by bisection of x Q'/Q=4h/M;
every candidate x gives a valid bound even if it is not optimal.

The checker uses exact rational arithmetic and directed integer
intervals for logarithms. An exact binomial coefficient, with a
directed logarithm of its leading 64 bits, handles C(M,h). It never
forms a power of a polynomial value to exponent J. The beta family is
9/8, 5/4, 11/8, 3/2, 2, 4, 8.

## Verification and coverage limits

Two independent mathematical and source audits passed (1)-(9),
including the negative-association induction, quota coefficient,
all support/message quantifiers, and the free endpoint. A serial,
file-free Peach replay passed 379 directed logarithm checks, 551
entropy checks, 756 exact negative-association toy comparisons, 441
exact quota coefficient comparisons, and 12 threshold rows. The
standalone hard-profile comparison verifies the strict exponent -4094.

A separate all-h replay constructed all 24573 threshold rows and
verified their integer tail inequalities. There are 19397 rows using
the quota bound, 5161 using negative association, and 15 deterministic
impossibility rows. The ordered (h,r,threshold,kind) table has SHA-256
`81aa4241dde1a1751093d2001ef76a6764d7aedd0029a518655abfa0d53081fb`.
The full collection of rank planes has total exceptional probability
below 2^-77. The replay took 15.6 seconds and wrote no remote files.

The existing bounded profile scout was rerun with these rank planes.
Of 117662 profiles surviving its earlier exact class screens, 82351
are now excluded by selector capacity. This removes 34422 of the
58778 formerly FLOAT-open points; 24356 remain FLOAT-open. The grid
contains 448811 profile samples, plus 21 separate footprint screens.
These counts describe a diagnostic grid, not coverage of its complement.

The remaining worst sampled tuple is
(h,R,Ls,f,F)=(3584,819200,845414,4,1048576). Scalar capacity alone
does not exclude it. The next question is whether H-rich selected
quartets can also support the required state overlap across consecutive
independent partitions. The full 45%/40-bit certificate remains open.
