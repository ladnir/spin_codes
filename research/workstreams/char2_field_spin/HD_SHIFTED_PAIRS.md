# Alternating adjacent-pair routes

This is a proof/implementation candidate, not a certified replacement.
It retains the certified grouped Toeplitz outer and h4 inner. Its purpose
is to break permanent source-packet equivalence classes without giving up
two-record loads. Setup and allocation remain outside timed encoding.

## Sampled construction

Index the 8192 outer groups cyclically. Each of the 256 regions contains
one field symbol from every group. There are two perfect matchings:

    pairing 0: (0,1), (2,3), ..., (8190,8191)
    pairing 1: (1,2), (3,4), ..., (8191,0).

In each region, sample a uniform permutation of its 4096 pairs. Preserve
the two lanes within each pair. Consecutive pairs fill the 64-record inner
steps with the same striped coordinate assignment as the unshifted pair
experiment. The fixed h4 maps, independent inner refreshes, and outer
setup are unchanged.

Two setup schedules are implemented separately. Mode 1 independently
samples either matching with probability one half in each region. Mode 2
uses pairing `region % 2`, deterministically. Mode 2 still samples a fresh
uniform pair permutation in every region. It avoids an unnecessary random
choice and forces both matchings to occur equally often.

For contiguous loads, the implementation copies each region's first
record into one padding slot after its last record. A pair ending at the
cyclic boundary can then use the same unaligned 256-bit load as every
other pair. All other padding stays unused. This adds 256 record copies
per full encoding, not a new pass over the data. Odd pairs can cross cache
lines, so their cost must still be measured.

## What changes in the support count

Fix an active outer-group set U. Let r be its cardinality and let a and b
be its numbers of fully occupied pairs under matchings 0 and 1. The two
matching histograms are respectively

    (singleton pairs, full pairs) = (r-2a, a), (r-2b, b).

For a nonempty proper U, its number c of cyclic runs satisfies

    c = r-a-b.

Thus a support cannot remain almost entirely in full pairs under both
matchings unless it has few runs. Those supports have much lower
multiplicity. The number of indexed cyclic binary strings of length L,
weight r, and c runs is exactly

    L/c * binom(r-1,c-1) * binom(L-r-1,c-1).

Strings are not identified under rotation. A valid upper bound for the
number with a specified triple (r,a,b) is the minimum of this run count
and either matching's support count

    binom(L/2,a) * binom(L/2-a,r-2a) * 2^(r-2a).

There is a stronger exact formula for the original obstruction, where
all active pairs are full in one matching. If r=2a, contract those pairs
to a cycle of length L/2. The exact count for b=a-c is

    (L/2)/c * binom(a-1,c-1) * binom(L/2-a-1,c-1).

For example, at r=1024 and a=512, the old full-pair class has about
2220.71 bits of support multiplicity. Subclasses with c=128 and c=256
have about 1200.82 and 1831.59 bits, respectively. These savings are in
the support union, not in the encoder's claimed security margin.

The exact joint count also has a small transfer representation. For
bits u,v in {0,1}, put T_z(u,v)=x^v z^(uv). Then the coefficient of
x^r y^a z^b in

    trace((T_y T_z)^(L/2))

counts precisely the indexed supports with that triple. The implementation
uses an equivalent integer dynamic program for small validation instances.
Three tests compare it with exhaustive strings, marginal binomial counts,
run counts, pairing symmetry, and the proposed support-count upper bounds.

## Remaining proof obligation

First dominate the expected outer counting measure for a fixed U by
beta^r times the product-uniform measure on its 256r potential coordinates,
where beta=256^256/255^128. This is the established Toeplitz outer bound.
The actual encoded coordinates across regions are not independent. The
dominating comparison coordinates are independent, as are the independently
sampled route permutations and inner refreshes. Apply regional composition
under this comparison measure, retaining its beta^r density factor in the
whole-message bound.

For this fixed U, let R0(U) and R1(U) be sound regional positive operators
for its two matching histograms. In mode 1 the regional average is
(R0+R1)/2, composed across all 256 independent regions. In mode 2 it is
(R0 R1)^128. The latter product must retain its matrix order; the two
operators need not commute.

For h4, a histogram alone does not specify the orientations of singleton
pairs among local emission blocks. A proof must retain those orientations
or dominate all feasible assignments. It must also retain the dependence
between both histograms of the same U. Counting supports by the triple
does not by itself solve the inner moment bound.

A complete result additionally needs sound regional
operators, all feasible support classes, an aggregate union bound, and
directed arithmetic. The combinatorial observations above do not establish
35% distance or any whole-code margin for this route.
