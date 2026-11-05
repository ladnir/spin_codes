# Reset bounds at widths four and eight

This note extends the [width-two obstruction](WIDTH_TWO_OBSTRUCTION.md)
without treating occupied width-w batches as a uniform subset. At larger
widths, several occupied pair slots may share a batch. The argument below
works directly with the uniformly placed pair slots and rounds their empty
gaps inward to complete batches.

The main finite conclusion is negative: **width four cannot give 45%
distance with setup-failure probability at most 2^-40 under the initial
outer and paired route.** Its failure probability at 45% exceeds 0.018.
This does not rule out the 35% target at width four. The present argument
does not rule out either 35% or a near-50% target at width eight.

No simulations, numerical searches, or benchmarks were used. The formulas
and conservative rational inequalities below are explicit checking targets.

## Geometry in pair slots, before batching

Write M=4096 for the number of pair slots in each region. The region has
2M records. For an even width w, put h=w/2, so a complete input batch
contains h consecutive pair slots. Here h divides M.

Choose 2p consecutive outer groups, beginning with group zero. Their message
space has dimension D=256p over F. They occupy exactly p source pairs in
each even region and p+1 source pairs in each odd region. This remains true
even if some actual output coordinates of the selected groups vanish.

For a region with b occupied source pairs, let their destination slots be

    0 <= s_1 < ... < s_b < M.

They form a uniform b-subset of the M pair slots. The empty pair-slot gaps,
including the initial gap, form a uniform weak composition of M-b into
b+1 parts. Regions are independent. Swaps inside the pairs do not change
this geometry.

The occupied *batch* set is not assumed uniform, and its size is not fixed.
The proof never replaces the pair-slot distribution by such a distribution.

## Rounding an empty pair-slot gap to complete batches

Consider the empty gap after an occupied slot s and before the next occupied
slot t, or before t=M for the final gap. Its length is z=t-s-1 pair slots.
Keep only complete width-w batches contained in this gap, and let ell be
their number. The batch containing s is not kept.

The exact count is

    ell = max(0, floor(t/h) - floor(s/h) - 1).

For a final gap, t=M and this formula still holds. In all cases,

    h ell >= z - 2(h-1).                                    (1)

To see this, put r=s mod h. Before truncating at zero, the count is
floor((r+z+1)/h)-1. Multiplying by h gives at least
r+z+1-(h-1)-h >= z-2(h-1). Truncation at zero only improves the bound.
The final gap has a smaller rounding loss, but that improvement is unused.

Choose l largest post-occupied-slot gaps in a region, excluding its initial
gap. Impose state zero immediately after the batch containing the occupied
slot that precedes each chosen gap. There are at most wl linear equations.
Every kept complete batch then emits zero. If several constraints refer to
the same state, their rank can only decrease.

The kept intervals are disjoint. If G is the total length of the selected
pair-slot gaps, and Z is the fraction of the region's records in the kept
complete batches, (1) gives

    Z >= G/M - 2(h-1)l/M,
    0 <= Z <= 1.                                            (2)

This is a deterministic statement for each routed setup. It accounts for
collisions of occupied pair slots inside the same width-w batch.

## A finite expectation and failure bound

For l<=b, define

    mu(b,l) = ((M-b)/M) * (l/(b+1)) * (1+H_b-H_l),
    nu_w(b,l) = mu(b,l) - 2(h-1)l/M,                         (3)

where H_j=sum_(i=1)^j 1/i. The Dirichlet-multinomial and conditional Jensen
argument in the width-two note gives E[G/M]>=mu(b,l). Consequently

    E[Z] >= nu_w(b,l).                                      (4)

The value nu_w is a lower bound on the mean fraction we can erase, not an
upper bound on the amount that some better reset strategy might erase.

Choose an integer l_r for every region r, with l_r no larger than its
occupied-pair count b_r. Suppose

    w sum_r l_r < 256p.                                    (5)

Then a nonzero message survives all state-reset constraints. With Z_r as
above, its output weight is at most

    N (1 - average_r Z_r).                                 (6)

Let nu_bar be the average of the 256 rational bounds nu_w(b_r,l_r).
If nu_bar>1-delta, regional independence and Hoeffding imply

    Pr_setup[d_min <= floor(delta N)]
      >= 1 - exp(-2 * 256 * (nu_bar-(1-delta))^2).            (7)

More precisely, Hoeffding bounds the event average Z_r<=1-delta. Outside
that event, (6) is strictly below delta N, which suffices for the displayed
integer bound. The randomness in (7) is only the route. All other setup
choices may be fixed arbitrarily first.

The selected message may depend on every setup choice. No rank-independence
assumption is made about the reset equations. Initial regional gaps are
unused, and no final-state constraint is added.

## Finite K20 rejection of 45% at width four

Set w=4, h=2 and p=64. Thus we select 128 consecutive groups, with message
dimension D=16384. There are 64 occupied pairs in every even region and 65
in every odd region.

Select l_r=16 gaps in 255 regions and l_0=15 gaps in region zero. The rank
budget is

    4 * (255*16+15) = 16380 < 16384.                         (8)

For the ordinary regions, the following rational lower bounds suffice:

    H_64-H_16 > 34/25,
    H_65-H_16 > 137/100,

    nu_4(64,16)
      > (4032/4096)*(16/65)*(59/25) - 32/4096 > 557/1000,

    nu_4(65,16)
      > (4031/4096)*(16/66)*(237/100) - 32/4096 > 557/1000.

For the exceptional even region, H_64-H_15>71/50 gives

    nu_4(64,15)
      > (4032/4096)*(15/65)*(121/50) - 30/4096 > 53/100.

Therefore

    nu_bar > (255*(557/1000) + 53/100)/256 > 139/250 = 0.556.

At delta=0.45, the needed deleted fraction is 0.55. Formula (7) yields

    Pr_setup[d_min <= floor(0.45 N)]
      > 1 - exp(-2*256*(0.556-0.55)^2)
       = 1 - exp(-1152/62500)
       > 0.018.                                            (9)

This is already incompatible with failure probability at most 2^-40.
It also excludes every larger proposed distance at that failure target.
The obstruction is unchanged by pair-antipodal state wiring, since it
requires only the state dimension and linearity.

At delta=0.35, these bounds do not erase the required fraction 0.65. Thus
(9) does not reject width four for the initial 35% target. It also does not
prove that width four passes every reset obstruction or has any positive
distance guarantee.

## Width eight: a valid bound, but no rejection at the intended distances

The same choice p=64 can use l_r=8 in 255 regions and l_0=7 in one region.
Its constraint count is

    8*(255*8+7) = 16376 < 16384.

Its exact lower bound on mean deleted fraction is

    [127 nu_8(64,8) + 128 nu_8(65,8) + nu_8(64,7)] / 256.    (10)

The rounding loss in nu_8 is 6l/M, because h=4. Equation (10) remains a
sound finite expression for the actual route, with all batch collisions
accounted for. No failure bound at 35%, 45%, or near 50% is claimed from
this expression.

Passing this particular group-cluster strategy is not a certificate.
Other group subsets, jointly selected cross-region gaps, partial-state
constraints, or properties of the fixed outer may give stronger attacks.

## Asymptotic interpretation of this particular attack

There is also a clean asymptotic statement. Keep the 256-region structure
and fixed even width w. Let M grow, and choose p with p tending to infinity
and p/M tending to zero. Select 2p adjacent groups and use

    l = floor((p-1)/w)

gaps in every region. Then (5) holds. The occupied-pair counts are p and p+1,
and the batch-rounding loss in (3) tends to zero.

For uniform simplex spacings, the sum of the largest fraction 1/w converges
in probability to

    (1+log w)/w,

where log is natural. This follows by representing the spacings as normalized
independent exponentials. The upper 1/w fraction starts at log w, and its
mean mass is the integral of x exp(-x) from log w to infinity.

The discrete route has the same limit. Conditional on its simplex vector,
the gaps are multinomial with total M-b. The expected total variation
error between normalized counts and that vector is at most the square root
of (b+1)/(M-b), which tends to zero. Selecting the largest l coordinates is
Lipschitz in their l1 distance. The unused initial gap also has vanishing
mass. These observations justify the limit without assuming independent
exponential gaps for the finite route.

For every fixed epsilon>0, the reset construction therefore gives, with
probability tending to one,

    d_min/N <= 1 - (1+log w)/w + epsilon.                    (11)

The limits in (11) are about 15.34%, 40.34%, and 61.51% at widths two, four,
and eight. These are asymptotic *upper bounds* from one attack, not finite
K20 certificates or predictions of typical distance. The width-eight value
is weaker than the rate-half Singleton bound and supplies no useful new
ceiling near 50%.

## Next discriminating question

For width four, the finite result rejects 45% but leaves 35% open. The next
proof question is whether a continuation bound can control low-weight words
below the clustered-reset ceiling, including all group unions.

Width eight avoids the demonstrated sub-50% ceiling, but costs and partial
state patterns remain unresolved. Compare its rank charges and fixed state
wiring before selecting an implementation. More coefficient randomness
cannot repair the state-dimension obstructions proved here.
