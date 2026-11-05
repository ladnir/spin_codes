# Two-update counting must retain the regional counter

Checkpoint, 2026-10-30. The construction and45%/40-bit goal are unchanged.
The full certificate remains open. This continues the
[one-step envelope audit](QUARTET_CLOSED_GLOBAL_ENVELOPE.md), using the
same seven pending types and closed paid sets.

There are two separate dependencies to retain. Ending a selected run
creates a blocked coordinate that constrains the next partition. Also,
each region has exactly2100 original-inactive inputs. A better temporal
calculation cannot compensate for replacing those regional quotas by one
global count.

## Exact two-update rows

Let a=b^lambda, d=b^(-tau), T=1+xy, and let h weight the destination types.
Use the local destination sum K_m and non-full polynomial N(h) from the
preceding checkpoint. Suppose all512 incoming coordinates have type L2.
The first update gives each coordinate either type B with weight dT or
type L2 with weight aT. These weighted choices factor across coordinates.
In the second partition, a quartet containing both types is forbidden.
Therefore the exact two-update product-weighted row is

    { T^4 * [d^4*N(h) + a^4*T^4*(d*h_B+a*h_L2)^4] }^128
        / h_L2^512.

This includes every intermediate death count and every final population.
It does not choose the most favorable intermediate count separately.

Other homogeneous sources admit a similar reduction. For source Z, the
first-step weighted coordinate choices are

    Z: x,   V0: 1,   L0: a*T.

They factor independently. Substituting these choices into the second
quartet kernel and taking its128th power gives the exact two-update row.
Sources L0 and L1 have deterministic intermediate types. Source L2 has the
two choices above.

For source B, each first-step quartet must contain at least one Z after
the update. For sources V0 and V1, a full first-step quartet instead
produces four V1 or four B coordinates, respectively. These correlations
must not be replaced by independent coordinate choices.

Let w be the width, g=w/4, t the intermediate count of this full-quartet
target, and z the intermediate Z count. For V0 or V1, t must equal4f for
an integer f. Conditional on the intermediate coordinate labels, the
exact probability that the first partition is compatible is

    binom(g,f)/binom(w,t)
      * [v^z]((1+v)^4-1)^(g-f) / binom(w-t,z).

It is zero for incompatible counts. For source B, use t=f=0. The first
factor requires all target coordinates to occupy complete quartets. The
second requires every remaining quartet to contain a Z.

Construct the second-quartet polynomial using the first-step coordinate
weights, with markers for t and z. Raise it to power g, then multiply
each coefficient by the compatibility probability above. This gives an
exact two-update row for every homogeneous source type. The computation
uses positive coefficients; there is no alternating inclusion-exclusion.

For an incoming mixture of B and L2, the first partition is possible only
when the B count is a multiple of4. Its probability is
binom(128,B/4)/binom(512,B). A separate positive upper bound keeps this
factor and the second partition's cohort count, but relaxes the requirement
that each first-step B quartet produce a Z. It is not the exact mixed row.

## Checks and numerical result

`quartet_closed_pair_population.py` compares the formulas with direct
composition of the original kernel at widths4 and8. All22 exact row
comparisons and mixed-row upper comparisons pass. Another182 rational
evaluations check the original-inactive generating polynomials. Each
comparison uses degree+1 distinct points, so it checks every coefficient
for the stated small profile, not just a single marker value.

A coefficient-adjoint calculation proposes product weights without
rerunning the long positive recurrence for each derivative. Its49 small
finite-difference checks pass. Proposal derivatives are FLOAT; no certificate
depends on optimizer convergence.

At the preceding weights, the exact L2 pair row has log2 value1080.96335,
versus1289.76130 from twice its one-step row. The other homogeneous rows
show smaller improvements, or require their age weights to be rebalanced.
A bounded homogeneous-source optimization takes71 seconds and leaves a
positive diagnostic of about638052 in the full-stream expression. This
does not cover all source populations and is not a margin.

The mixed B/L2 upper bounds include a negative repeated-row diagnostic
at B=128 and at B=256. These are not whole-code certificates: a valid
trajectory need not remain in those source populations, and the marker
calculation below still loses the regional quotas.

## A closed suffix defeats the global-count relaxation

The older [regional-quota audit](QUARTET_REGIONAL_QUOTAS.md) established a
related obstruction for all paid sets. Closing the paid sets changes their
multiplicity, so that result cannot simply be reused. The following count
verifies the obstruction directly for closed sets.

Let N=2097152, R=537600, and choose a suffix of L positions, aligned to
complete512-coordinate updates. Choose c paid inactive positions inside
the suffix and R-c unpaid inactive positions outside it. The suffix is
fully allowed. Its entire state support is feasible, with no terminal
death charge, so the canonical score is at least lambda*L.

Every chosen payment is touched by this support. Conversely an unpaid
position cannot be touched by a feasible support. Thus every resulting
paid set is closed. The exact number of pairs(I,P) is

    binom(L,c) * binom(N-L,R-c).

This count uses one prescribed suffix and does not count overlapping
suffix events separately. It holds for every realization of the quartet
partitions. All3388 tested small suffix witnesses are closed and attain
the stated score lower bound.

If one extracts only the global original-inactive count R, the normalized
r=2100 expression, including the row-set union, is therefore at least

    binom(L,c)*binom(N-L,R-c) / binom(8192,2100)^255
      * b^(lambda*L-(6/5)*c-31669.6).

This is a lower bound on the relaxed moment expression, not on code
failure. It also lower-bounds any positive common-marker evaluation of
that same global coefficient.

Two fixed witnesses have directed log2 lower endpoints:

| L | c | range of log2(b) | trusted lower endpoint |
| ---: | ---: | --- | --- |
| 910336 | 121999 | at least3/2 | 48118430667851/2147483648 >22406 |
| 1077760 | 94501 | at least5/2 | 73290348952119/536870912 >136513 |

Both counts respect the surviving payment cap255463. The exponent of b
is positive in each row, which proves the stated ranges from the endpoint
calculations. The current log2(b)=5/2 is included.

Both events are impossible under the actual regional distribution. Every
complete region in a fully allowed suffix must have all2100 inactive
positions paid. The two suffixes contain at least111 and131 complete
regions, requiring at least233100 and275100 payments. These exceed their
specified c values. Thus exact temporal counting alone cannot close the
current global-quota relaxation. Smaller moment parameters are not ruled
out by these two witnesses, but do not repair the missing regional law.

## Carry the pair operator through real region boundaries

Write T_i for one update with exactly i original-inactive inputs, retaining
both endpoint populations. The two-update coefficient is

    P_j = sum_(i=0)^j T_i T_(j-i).

The intermediate pending types are contracted, not reset or terminally
tested. A real region consists of eight such pairs. Its raw operator is
the coefficient of x^2100 in (sum_j P_j*x^j)^8; divide by
binom(8192,2100), then carry its outgoing pending types into the next region.
There is one terminal test after the256th region. Scalar markers used to
bound a coefficient must not silently replace this counter by a single
whole-stream count.

`quartet_closed_pair_regions.py` checks this composition on the exhaustive
width-four, two-region instance. It reproduces all11 paid-count/score
coefficients of the direct closed-set enumeration, with peak69 states.
Two negative controls fail as intended: replacing the regional quotas by
one total, and adding a terminal test at the internal region boundary.

The next test is a quota-conditioned two-update bound for mixed source
populations and their exits. Keep the intermediate blocked/selected
interaction and the original-inactive counter together. Do not spend a
longer-run computation merely sharpening the known-failing common-marker
expression. Full-size coverage, other row counts, directed upper bounds,
and the final failure-budget audit are still required.

All experiments were serial, remote, source-only, and file-free, with
explicit180-second deadlines and small arrays. No construction, library,
paper, or frozen-helper changes; no commits, native benchmarks, or dumps.
