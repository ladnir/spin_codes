# Closing the aligned six-defect frontier branch

For the unchanged three-shear model, the complete no-full-reset branch
h=6,a=1,c=7 now has a restricted 46% distance bound below 2^-78. This includes
all total live lengths T, early episode lengths, noncrossing episodes,
sharing patterns, and active-row counts. It is not the whole six-defect
certificate or a whole-code certificate.

The earlier joint-region verifier left exactly 63 values T=1878,...,1940
in this branch. Its worst was T=1938, with upper exponent +6491. The new
count closes all 63, and the unchanged verifier covers the remaining T.

## Geometry forced by this branch

Here h counts defective incoming batches, a indicates Q_1 is active, and
c counts nonzero incoming-state intervals. The c-a=6 later activations
already consume all six defective batches. There are no internal births
or one-sided pairs. The six closed early episodes are clean.

In the residual length range, some state is full by the retained proper-
state bound. Since the branch forbids a later reset from a full state,
the final live interval extends through Q_(L-1). Every earlier interval
is proper and has length l_i in {1,...,8}. The first activation is input
batch zero, so its l_1+1 input batches lie in one region.

Let E=sum_i l_i and d=max_i l_i. Then d+5<=E<=6d. Let p be the number of
crossing early intervals. Each other early interval lies within one region,
and 0<=p<=5 because the first is aligned. Before merging shared regions,
the early intervals contribute m=6+p fragments. The terminal interval
touches r=ceil((T-E+1)/16) regions. It is safe to give these terminal
regions full occupancy in the counting bound.

There are R_v=m+r virtual regions. Let g count duplicate regional copies
removed when fragments from different intervals share an actual region.
The actual number of touched regions is R=R_v-g. RS minimum weight forces
R>=129. Thus k_v=R_v-128 and 0<=g<=k_v-1. On the residual T range, k_v<=4.
All other count combinations are either excluded or conservatively included.

## One-region activation factor

For a clean episode contained in one region, fix its nonzero row count n.
If its initial support has size w, it consists of w/2 complete matching
pairs. All later inputs use the same coordinates. Consequently its joint
route/reset contribution is at most

    sum_(even w>=2^l)
      binom(b/2,w/2) binom(lw,n-w) / binom(M,n) * (l/q)^w.

Use b=512, M=8192, q=2^32. Invalid terms are zero. This applies both to
the aligned first episode and to a noncrossing episode elsewhere; translation
within a region changes no subset count. The matching phase can change,
but every phase partitions the same b coordinates into b/2 pairs.

For crossing intervals use the previous two-region factor. The verifier
takes a coefficientwise maximum over lengths 1,...,d separately for each
one-region or two-region factor. This relaxes the exact sum E in the local
product. It retains the actual E in the route cap and terminal-region count.
The relaxation is an upper bound: every actual factor is dominated, and
the actual maximum length still forces s>=2^d active rows.

## Actual route cap with arbitrary sharing

The sum of early fragment fractions is (E+6)/16. A row eligible for any
nonzero compatible RS word must occupy at least 129-r early fragments
in the virtual comparison. With m fragments and mean fraction
mu=(E+6)/(16m), a positive factorial-moment upper bound is

    beta = min(1, min_j [binom(m,j)/binom(129-r,j)] * mu^j),

for 1<=j<=129-r. If m<129-r, eligibility is impossible. The elementary-
symmetric bound justifies replacing the actual fractions by their mean.
Splitting merged regional fractions only increases these symmetric
polynomials, after applying the same-set moment bound in each actual region.
It does not assert independent row routing into different fragments.

For s specified rows, the joint bound is beta^s. Choose the row cap C so
binom(M,C+1)*beta^(C+1)<=2^-282, using directed rational logarithms. There
are at most (4097)^14<2^182 endpoint tuples for seven intervals. The route
exception for all relevant geometries is therefore at most 2^-100.

A real region may contain more than two fragments when noncrossing episodes
are present. For r_0 fragments, their independent-denominator comparison
is bounded by r_0^s*(M/(M-s+1))^s. There are at most (s+1)^(r_0-1) input
allocations to sum. Since r_0<=2^(r_0-1), multiplying over all merged
regions is bounded by the existing correction

    [(s+1)*2^s*(M/(M-s+1))^s]^g.

This explicitly pays for dependence and allocations even when several
episodes occupy one region. Treating those regional routes as independent
without the correction would be invalid.

## RS zeros and the joint count

There are six episode factors: p crossing and 6-p single-region. Their
order is irrelevant to the upper envelopes, not to the actual code. The
global endpoint tuple still fixes the actual order and matching phases.
The dynamic program retains total virtual missing coordinates. For s
active RS rows and Z actual zeros among R touched regions, its budget is
Z+gs. The message count is at most binom(Rs,Z)*q^((R-128)s-Z).

Sum over Z=0,...,(R-129)s and s=2^d,...,C. Include the shared-region
correction, row choices binom(M,s), and one global density factor two.
Each endpoint tuple determines T,E,d,p,g. A maximum over those count
classes followed by the endpoint union covers all of them; adding a
separate union over every count class would count the same partition twice.

## Integration with the retained branch verifier

At T=1938, the screen covers 1048 count classes, of which 656 are excluded
by the eligible-row cap. The new good-route exponent is -127. The complete
integration checks all 63 previously residual T values, covering 23763
count classes. Their worst good-route exponent remains -127, first at
T=1926. Adding the common new route exception gives a bound at most2^-99
for this entire previously residual range.

For every other T in the same no-full-reset branch, the unchanged verifier's
good-route exponent is at most -3806. Its retained route exception is
at most2^-80. Combining the four contributions conservatively gives

    2^-3806 + 2^-80 + 2^-127 + 2^-100 <= 2^-78.

The retained verifier uses the 46% output-weight threshold. The new
activation/route exclusion itself needs no output-weight restriction.
This proves the stated complete branch bound at46%, without altering the
previous five-defect certificate or the certified uniform-feedback control.

## Checks, reproducibility, and remaining work

`diagnose_pair_aligned_frontier.py --test` checks 66 exact single-region
subset counts, 4746 multifragment denominator ratios, and 105 allocation
counts. `--legacy-branch` recomputes the old residual list rather than
assuming a saved list. `--integrate` then checks every residual value and
combines both proof contributions. The integrated run took about4.7s on
Peach under256MiB and a45s process timeout. Local tests began with more
than22GiB free. All jobs were serial and generated no data files or caches.

Other h=6 branches, including a=0,c=6 and branches with an internal
defective-batch budget, remain open. The full-reset branch is also not
covered by this integration. The next direct target is a=0,c=6: five
closed early episodes, with no forced aligned first interval. Reuse the
single/crossing factors and first recompute that branch's residual T list.
After the clean activation-maximal branches, use the all-charge defect
envelope in branches with internal defects. Whole-code and extension-field
closure remain unproved; no structural reason for redesign has emerged.
