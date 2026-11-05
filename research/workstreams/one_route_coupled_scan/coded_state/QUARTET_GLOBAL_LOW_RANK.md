# Global rank: grouped outputs and the low-rank remainder

2026-10-31. The construction is unchanged. The active goal is the full
K=2^20, rate-1/2, 45%-distance, 40-bit certificate. Local seven-type
refinements are parked. This note builds on
[the inverse-rank reset](QUARTET_GLOBAL_RANK_RESET.md).

The new result bounds all bad outputs whose covered footprint has size
at least 53N/64. Their total probability is strictly below 2^-6359.
This is an all-history class theorem, not the full minimum-distance claim.
The remaining histories have smaller footprints and fail the earlier
profile-sensitive rank threshold. They still need a complete bound.

## Sum optional outputs before taking the moment

Retain q=2^32, N=2^21, K=2^20, W=943718, M=8192, and T=4096 updates.
Setup is sampled once. Its global acceptance-density factor is gamma<2.

Fix a nonzero state-support history S and the quartet partitions, before
exposing the route. Let b count state births, d state deaths, and i
surviving nonzero states. These are totals over all coordinates and
updates. A birth or death forces a nonzero output; a surviving state
permits either output choice. The terminal state is unrestricted.

As in the reset note, L contains every coordinate of a live incoming
quartet. Let B0 contain new state coordinates outside L, and put V=L union B0.
Unlike b, B0 counts only births in incoming-zero quartets. V depends on S
and the partitions, but not on the optional output choices at survivors.
Write R(V)=sum_r min(128,|V intersect outer row r|).

For each optional-output support, normalize its first nonzero output
coordinate to 1. Scaling preserves all support conditions and the
setup event. The joint death/RS rank argument then saves one field
dimension. Conditional on the route and partitions,

    Pr_field[there is an output with state support S and weight <= W]
      <= min(1, (gamma/q) q^(b-R(V))
             sum_{e=0}^{min(i,W-b-d)} binom(i,e) q^e).          (1)

The sum is zero when b+d>W. This counts all messages in this state class.
It is not a fixed-message estimate.

Choose 0<theta<=1 and 0<z<=1. Define

    a=q^theta,  m=a z,  c=(1+m^(1/theta))^theta.

Apply a cumulative-weight marker h=z^(1/theta)<=1 to (1), then take
the theta power. The resulting per-state-history bound is

    (gamma^theta/a) z^(-W)
      E_geometry[a^(-R(V)) m^b z^d c^i].                     (2)

This order matters. Optional output supports are summed within one S
before taking the fractional moment. Summing them afterwards would
replace c by the larger 1+m. We never take a fractional power after
summing all state histories.

## A rational all-history coverage certificate

We prove the stated class theorem with exact rational inequalities.
Fix a history and partitions with |V|>=F0, where F0=53N/64. Select any
deterministic F0-subset U of V before exposing the route. Then R(V)>=R(U).
The selected positions in region j map to a uniform subset of outer rows.

Condition independent row/region Bernoulli incidences on their exact
regional counts. Each conditioning event is a binomial mode, so their
joint probability is at least (M+1)^(-256). For t>=1 and rho>0, the
row function t^((128-k)_+) rho^k is convex in k: it is the maximum of
two exponentials. Equalizing regional probabilities therefore increases
the row expectation, while preserving total coverage F0. With p=F0/N,

    E_route[a^(-R(V))]
      <= a^(-K) (M+1)^256 rho^(-F0)
         [E_{Z~Bin(256,p)} a^((128-Z)_+) rho^Z]^M.             (3)

This bound is uniform over all such histories and partitions. It retains
the actual regional quotas through conditioning; it does not replace
the route by independently sampled row labels.

Set theta=1/4, a=256, z=1/256, p=53/64, and rho=176/53. The row moment
in (3) is bounded using the two positive exponential terms:

    E a^((128-Z)_+) rho^Z
      <= a^128(1-p+p rho/a)^256 + (1-p+p rho)^256
       = 2 (187/64)^256.                                    (4)

After (3), discard the coverage restriction when summing state histories.
No partition variable remains in their weights. Each of the 512 state
coordinates has the two-state transition matrix

    H = [[1, 1], [1/256, 2^(1/4)]].

Rows denote the previous state, columns the next state. Start at zero
and accept either terminal state. For the positive vector v=(1,107/512),

    H v <= (619/512) v.

Indeed, 2^(1/4)<609/512, and direct rational arithmetic verifies the
second row. Since 1/(107/512)<8, the complete history sum is at most
2^1536 (619/512)^N. Including the zero history only enlarges this bound.

Also gamma^(1/4)/256<2^-7 and (M+1)^256<2^3584. Combining (2)--(4)
gives a failure bound strictly below

    2^[-7+8(W-K)+M+1536+3584]
      [(115753/32768)^64 (53/176)^53]^(N/64).

The exact integer inequality

    115753^64 * 53^53 < 2^985 * 176^53

bounds the bracket by 2^25. The final exponent is

    -7+8(W-K)+8192+1536+3584+25N/64 = -6359.

This proves the all-history coverage claim over J=GF(2^32), with the
actual accepted setup distribution. It includes every output amplitude
and the free terminal state. As with the earlier class theorem, an
extension-field output need not retain the same coverage class in each
J-basis component; do not transfer this restricted class by itself.

Independent reviews checked both specializations, including the first
F0=7N/8 bound below 2^-56023. The stronger threshold above uses the same
argument and sharper rational constants. Both integer verifications are
retained in `quartet_global_coverage_probe.py`.

## What remains and why total coverage alone loses information

The remaining event has

    |V| < 1736704,
    R(V)+d < j+g(j,d),   1<=j<=943718,

with g from the reset note. Allocating the two proved class bounds uses
less than 2^-62 of the final budget. No claim is made for the remainder.

A total-coverage-only scalar envelope is still loose. In floating
diagnostics, grouping the output choices reduces its log2 upper bound at
|V|=3N/4 from about 1.53 million to 326700. This remains positive.
Those evaluations are diagnostic, not directed certificates or optimized
lower limits on the method.

Temporal concentration explains one important loss. A full-state burst
whose footprint occupies at most 128 regions cannot contain a nonzero
RS outer word. Equalizing its regional coverage creates artificial rows
with more than 128 covered coordinates and discards that exact exclusion.
A dense burst with a sparse fringe needs a bound retaining that fringe.

There is a second necessary control. Fix 228 state coordinates and 3648
outer rows. Suppose every quartet partition keeps those coordinates a
union of quartets, and every region routes those rows into their slots.
The resulting nonzero subcode is supported on 4096*228=933888 outputs,
below W, for every accepted choice of the local field maps. This is a
real, rare bad-setup family. For those fixed sets, its probability is

    [binom(128,57)/binom(512,228)]^4095
      * binom(8192,3648)^(-256).

A raw first moment rewards this event by q^(theta*466944) and can become
enormous. The conditional probability must remain capped at one, or use
valid fractional moments separately for different state/coverage classes.
This example does not contradict a 40-bit ensemble certificate.

## Two tools for the lower-coverage remainder

### Exact quartet transfer

The grouped bound admits an exact kernel with only the incoming state
weight and footprint counters. Let y mark outgoing weight, r incoming
weight, xL live-quartet coverage, and xB births in empty quartets. Put

    F(r,y)=(1+m y xB)^4
       + xL^4 sum_{u=1}^4 binom(4,u) r^u
                          (z+c y)^u (1+m y)^(4-u).

For incoming weight u and outgoing weight v, the transition coefficient is

    [r^u y^v] F(r,y)^128 / binom(512,u).

This sums every outgoing support and averages the actual uniform quartet
partition. It does not truncate destinations or introduce run ages.
Setting xL=xB retains total V; keeping both permits forced-birth tests.
For the earlier ungrouped moment, use c=1+m instead. Exact enumeration
checks 212 all-destination coefficients at widths four and eight: 106
for the ungrouped kernel and 106 for a grouped rational specialization
q=16, theta=1/2, z=3/16, m=3/4, c=5/4.

### Covered regions and the necessary fringe

Fix S and partitions. Select c0<=128 regions C containing a region r0
with b0>0 births from B0. Let Rcore=sum_{r in C}v_r, where v_r is regional
coverage, and k=129-c0. Deterministically R(V)>=Rcore: these columns
contain at most 128 coordinates per outer row.

The b0 births in r0 map to distinct rows. Each must have at least k
covered coordinates outside C. Write p_r=v_r/M. For each birth row,
choose k fringe regions witnessing these hits. If a region receives d_r
distinct required rows, its exact hit probability is
(v_r)_(d_r)/(M)_(d_r) <= p_r^d_r. Summing the witnesses proves

    E_route[a^(-R(V)) 1_compatible]
      <= a^(-Rcore) min(1, e_k((p_r)_{r outside C})^b0),       (5)

where e_k is the elementary symmetric polynomial. This argument uses
sampling without replacement and no row-independence assumption.
In particular (5) vanishes for footprints in at most 128 regions.

Accepted superregular maps give additional exact filters. Put
Hrows={r:|V intersect row r|>=129}. Inputs outside Hrows are zero by RS
shortening. In a live quartet with incoming support size u, next support
Snext, and output support J0, compatibility requires

    u+|Snext union Hslots|>=5,
    u+|J0 union Hslots|>=5,

where Hslots contains that quartet's positions routed from Hrows. These
follow from superregularity of A and A+delta I. For a singleton incoming
state, every omitted next-state coordinate and every zero-output
coordinate has nonzero input. Together with B0 these are fixed required
hits in Hrows. They prevent the raw field calculation from charging
witnesses that exist only under rejected maps.

## Verification and next action

The selected bounds, 212 exact transfer comparisons, and 1536 exact
core/fringe checks passed in
bounded, file-free Peach runs. The route diagnostics and scalar searches
use floating arithmetic and remain explicitly non-certifying. Runs were
sequential, single-threaded, and below seven seconds of numerical work
each. No benchmarks, construction changes, paper edits, or production
changes were made.

Next retain regional concentration in the grouped global sum: use the
quartet transfer for sustained sparse states, and (5) for a dense core
with a sparse fringe. A candidate compression keeps total coverage and
the coverage of the 128 fullest regions. The latter is a deterministic
rank lower bound; equalization within each half can retain it. This is
a proposed next bound, not a proved coverage result. Test the
invariant-subcode family and the repeated
full-state burst controls before attempting a complete numerical bound.
Do not resume the parked seven-type refinements or treat the positive
scalar envelope as an obstruction to the fixed construction.
