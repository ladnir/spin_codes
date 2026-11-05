# Why the filtered geometry count cannot close

2026-10-30. The construction remains fixed. This note answers the coefficient
question left by [the marked-region test](QUARTET_MARKED_TRANSITION.md).
Even exact coefficients cannot close that filtered first-moment bound.
This is an obstruction to one proof calculation, not to the code's distance.

Fix one inactive-row set T of size2100. Each of256 regions contains exactly
2100 inactive labels, uniformly placed among8192 positions. The quartet
partitions are independent at each update. For a selected state history S,
let A count selected coordinates, C count covered inactive inputs, and D
count deaths. The geometric score is

    F(S)=A/4-C-(37/50)D.

The sufficient geometric cutoff is18171. We construct a subfamily of
histories exceeding this cutoff that pass every local filter currently
implemented in the marked kernel. Their expected number exceeds2^10134.
The expectation is over the actual route and quartet partitions. Histories
are not asserted to be realized by codewords or to maximize the score.

## A directed regional subcount

Use the fully filtered153-term polynomial from the marked-region test.
A local term records incoming selected count a, inherited hole births e,
outgoing selected count v, deaths d, inactive labels l, and coverage c.
Fix the following numbers of quartet types (a,e) at one update:

| (a,e) | 131-label update | 132-label update |
|---|---:|---:|
| (0,0) | 10 | 10 |
| (0,1) | 21 | 19 |
| (0,2) | 10 | 11 |
| (0,3) | 1 | 1 |
| (1,0) | 2 | 3 |
| (2,0) | 21 | 21 |
| (3,0) | 40 | 41 |
| (4,0) | 23 | 22 |

Both columns contain128 quartets,256 incoming selected positions, and44
inherited births. Sum ALL outgoing assignments allowed by these types,
subject to v=256,d=44,c=27 and the stated label count. This is less
restrictive than fixing all six local statistics quartet by quartet.

The incoming-type weight is

    128! / product_(a,e) n_(a,e)!
      * product_(a,e) [binom(4,a)binom(4-a,e)]^n_(a,e)
      / [binom(512,256)binom(256,44)].

For each incoming type, divide the local polynomial by its incoming
binomial multiplicity before summing outgoing assignments. Use the outgoing
birth count g=v-a+d. Live quartets have e=0 and l=c; empty quartets have
d=c=0. Thus the coefficient splits into a three-dimensional live calculation
and a two-dimensional empty calculation, joined by total outgoing births44.

`quartet_grouped_type.py` computes directed lower and upper bounds using
positive integer fixed-point recurrences. Rational tilts change only the
normalization, not the counted polynomial. The frozen replay does not use
a floating-point proposal. Its per-step log2 intervals are

    131 labels: [465072964179459/1099511627776,
                 58134155925175/137438953472],
    132 labels: [463950926216745/1099511627776,
                 463951211949341/1099511627776].

Take12 specified updates of the first type and4 of the second. Divide by
binom(8192,2100) ONCE for the regional label quota. The region has

    A=4096, C=432, D=704, F=1776/25=71.04.

The log2 expected continuation count lies in

    [47394462823793/1099511627776,
     23699502208735/549755813888]
      = approximately [43.1050128316,43.1091433871].

This continuation count is identical for every incoming pair consisting of
256 selected positions and44 disjoint inherited births. Coordinate symmetry
therefore permits concatenation without assuming independence between
successive supports. Region coins are independent; the quota is retained
within each region.

## Complete the initial and terminal boundaries

Use one prefix region,254 middle regions of the preceding form, and one
suffix region. The state starts empty and has no imposed terminal reset.

In the prefix, keep the state empty for14 updates. Choose256 selected
coordinates at update15, then replace44 of them at update16. Require all
2100 inactive labels to occur in the first14 updates. This passes all local
filters and has A=512,C=0,D=44. Its expected count is

    P=binom(512,256)binom(256,44)^2
        * binom(7168,2100)/binom(8192,2100).

In the first suffix update, segregate the256 selected coordinates and their
256-coordinate complement into homogeneous quartets. Require exactly the
complementary positions to be inactive at this update, and keep the support
unchanged. Inherited births then vanish. At the second suffix update,
require each quartet to have two incoming selected positions, and select
all512 outgoing positions. Keep the state full thereafter. Place the
remaining1844 inactive labels anywhere in those last15 updates.

The conditional probability of this suffix is

    Q=binom(128,64)*6^128*binom(7680,1844)
        / [binom(512,256)^2*binom(8192,2100)].

The suffix has A=7936,C=1844,D=0. A full outgoing quartet passes the
canonical deletion filter for every inactive mask when its incoming weight
is two or four. Segregated empty quartets also satisfy the inherited-birth
budget. Hence the complete family passes the same local filters at both
boundaries, including the no-terminal-birth condition.

The complete histories have

    A=1048832, C=111572, D=178860,
    F=91398/5=18279.6 >18171.

Let Z count these histories. Directed arithmetic gives

    log2 E[Z] >=11143504282125779/1099511627776
              >10134.95.

The counted geometries also satisfy the full-width burst bound for this T.
Middle updates contain131 or132 inactive labels, so a full-state update
contributes at most128-l<=0 before death costs. In either boundary region,
the sum of positive contributions is at most

    16*128-4*128-52=1484,

using l<=512 and sum l=2100. Thus every full-width burst has score at most
2968<3702. This assertion is for the counted T; it does not establish
burst-goodness simultaneously for every other row set.

## What this rules out

For every b>=1, the filtered witness sum used by the geometric Markov bound
satisfies

    E sum_(locally filtered S) b^(F(S)-18171)
      >=2^10134.95 * b^108.6.

Improving coefficient extraction cannot make this expression small.
The local filters admit too many witnesses. This does not lower-bound the
probability that any witness exists: many witnesses can occur on the same
rare geometry. It also does not rule out counting canonical maxima,
imposing stronger global conditions, or retaining field constraints.

In particular, these histories have C+D=290432. The existing field argument
excludes low-output codewords with C+D>255463 in this row class. Consequently
the subcount does NOT demonstrate bad code distance. Conversely, the field
cutoff cannot simply be imposed on a geometric maximizing witness: replacing
a codeword support by an unconstrained maximizer need not preserve C+D.

## Reproducibility and next test

The helper checks150 exact coefficient cells,401 bridge-mask cases, and35
small partitions. Coefficients use24-bit probability bounds and38-bit state
bounds; assertions keep each product below64-bit overflow. Final convolutions
use Python integers. Frozen replay takes under one second on Peach, with a
90-second watchdog, under64MiB numeric storage, and no files written.
`quartet_marked_type.py` preserves the bounded proposal method separately.

The companion [whole-column bound](QUARTET_COLUMN_RUN_BOUND.md) controls
zero-coverage components without summing arbitrary supports. The next
geometric task is its positive-coverage remainder. Another valid direction
is to keep the field cutoff through a penalized objective, rather than
silently restricting the current canonical count. No new whole-code
failure allowance is charged by this obstruction result.
