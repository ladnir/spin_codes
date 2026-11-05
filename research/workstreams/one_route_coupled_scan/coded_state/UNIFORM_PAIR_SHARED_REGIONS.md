# Shared regions and the actual RS zero count

The fixed three-shear model now excludes the following shared-region
family except with probability at most 2^-99. There are five clean proper
closed episodes, each of length one through four and crossing one regional
boundary, followed by a live tail. Exactly one region is shared between
intervals, and the total number of touched regions is 131. The tail's
internal activity is unrestricted. This is not whole-code distance closure.

The preceding proof multiplied factors from disjoint regions. That product
does not apply directly here: two intervals draw their positions from one
regional permutation. The new argument retains that dependence and the
number of actual zero coordinates in the RS words.

## Splitting a region only in the counting argument

Fix a nonzero message and all interval endpoints. In a shared region,
let n be the number of nonzero outer-row coordinates, with n<=s for s
active rows. The regional permutation maps these coordinates to a uniform
n-subset of M=8192 positions. If n_1 and n_2 land in its two interval
fragments, then n_1+n_2=n.

For a fixed allocation and fixed activation supports, the two fragments
occupy disjoint positions. Their favorable subset counts multiply, but
the denominator is binom(M,n), not binom(M,n_1)binom(M,n_2). Comparing
with two virtual independent regions introduces the ratio

    binom(M,n_1)binom(M,n_2)/binom(M,n)
      <= 2^s * (M/(M-s+1))^s.

To see this, extract binom(n,n_1) and write the remaining ratio with falling
factorials. Bound its numerator by M^n and its denominator by (M-n+1)^n.
Then use n<=s. This upper bound is uniform over the activation supports
and allocations. It also permits relaxing a tail fragment to all positions.

The allocation is random, not fixed by the message. Sum over at most s+1
allocations for each shared region. For g shared regions the correction is

    [(s+1) * 2^s * (M/(M-s+1))^s]^g.

There are at most two intervals in each region in this crossing-episode
family: an interval must cross its left or right boundary, and disjoint
intervals cannot both cross the same boundary. The terminal tail does not
invalidate that count. The result here uses g=1; no larger-g coverage is
claimed from the general formula alone.

## Virtual zeros are not message degrees

Write R=128+k for the actual touched-region count. Splitting g regions
creates R+g virtual regions but does not change the outer code. Its message
count must still use k, not k+g.

More importantly, let Z be the total number of actual zero coordinates
among the R positions of the s nonzero RS rows. Then 0<=Z<=(k-1)s.
For each exact zero pattern, RS interpolation bounds the number of messages
by q^(ks-Z). Summing over patterns gives the upper bound

    binom(Rs,Z) * q^(ks-Z).

Patterns assigning too many zeros to one row have no nonzero witness;
including their number only increases this upper bound on the valid patterns.
Here q=2^32, and the forced zeros outside the R regions are already included
in the shortened dimension k.

Splitting one region preserves its total nonzero count and adds exactly
s virtual zeros. Therefore the early-fragment deficit obeys

    sum_i (2s-n_(u_i)-n_(v_i)) <= Z+gs.

The earlier relaxed bound used the maximum (k+g-1)s while charging every
message q^(ks). It allowed a message to avoid early input constraints
without paying for the RS zeros used to do so. Retaining Z repairs this loss.

For each s and Z, use the activation-route dynamic program with deficit
budget Z+gs. Its episode factors retain complete matched activation pairs,
containment of all later inputs in the initial support, and reset weights.
Multiply by the sharing correction, binom(M,s), and the message count
above. Sum over Z, then over s, and include the one global factor two for
accepted-map density. No independence of the actual shared routes is used.

## Eligible-row cap and endpoint union

The source-count cutoff still follows from an event uniform over all messages.
For a real shared region, its allowed-position fraction is the sum of its
fragment fractions, or a smaller value if a fragment was relaxed. Replacing
one weight a+b by two weights a,b increases every elementary symmetric
polynomial. This comparison is made after the actual regional factorial
moment bound, where each region uses its single allowed-position union.
It does not assume negative dependence for different destination sets.

Thus the prior eligibility majorant with virtual dimension k+g is an
upper bound for actual eligibility. Outside an event of probability at
most 2^-100 after the endpoint union, it gives a cap C on active rows.
The first moment still chooses rows from M; the cap only restricts the
range of s. It does not supply a second independent route discount.

At most 2^156 endpoint choices cover the six input intervals. Their
endpoints determine the location of the shared region, interval ordering,
splits, and lengths. No further union over those choices is needed after
taking a uniform maximum. The verifier computes 56 length multisets for
the stated length range one through four.

## Results and checks

Without the actual-Z count, the length-four profile gave a positive upper
exponent of 694 after endpoints. With Z retained, all 56 length multisets
close. The worst becomes five length-one episodes, with s=2 and Z=0:

- Good-route contribution: at most 2^-106.
- Route exception: at most 2^-100.
- Combined restricted-family bound: at most 2^-99.

All source counts allowed by the eligibility cap are included. No output
weight restriction is used. The older disjoint-region certificates remain
unchanged. This result does not cover all lengths, multiple shared regions,
noncrossing episodes, internal defects, or the whole code.

The helper `diagnose_pair_shared_regions.py` tests 3860 finite-population
ratios, 225 exact allocation/image counts, and 995 RS-zero/splitting identities.
The RS test uses the exact nonzero GF4 RS[4,2] spectrum across up to six
active rows. The complete corrected 56-profile run took about 7.3 seconds
on Peach. A local all-length sweep and a remote all-length sweep of the
earlier relaxed bound both hit their 40-second cap; neither establishes
full coverage. The source now also has a per-profile time check.

Runs were serial. Remote jobs used a 256 MiB address-space limit and
45-second process timeout. Only source and compact notes were retained.
Local space initially exceeded 23 GiB; a later check returned about 0.1 GiB.
That final check and tiny test had been placed in one shell call, so the
test ran before the low-space result was inspected. No further local
numerical runs are allowed until the safety floor is revalidated.

## Next discriminating step

Continue proof work, not redesign. The new gain comes from a constraint
already present in the outer: zeros consume message degrees. Extend the
joint Z/activation count to longer clean episodes in bounded profile groups.
Then incorporate internal births and one-sided pairs, where input support
can leave the initial activation support. That is a genuinely different
interface, and must not be covered by the clean-episode formula unchanged.
