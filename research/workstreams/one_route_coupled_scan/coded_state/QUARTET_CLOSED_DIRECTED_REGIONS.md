# Directed regional bounds with moving populations

Checkpoint, 2026-10-30. The construction and the full45%/40-bit objective
are unchanged. This extends the
[seven-type interface](QUARTET_CLOSED_POPULATIONS.md) from conditional
profiles to restricted regional transfers. It does not cover all populations
or provide the final minimum-distance certificate.

## Extract more than one count

The earlier conditional bounds evaluated most population variables at
positive markers. This loses the requirement that several counts attain
their targets simultaneously. The exact seven-type quartet polynomial
allows additional positive coefficient extraction.

For fixed incoming and outgoing populations n,n', let i count original
inactive inputs. The paid count is i-n'_Z. Thus either i or the paid count
can be extracted, along with selected population coordinates. Any remaining
coefficients are bounded by positive evaluation. This does not assume
independence between the retained counts.

`quartet_closed_extract.py` implements the positive recurrence. Thirty-nine
integer multivariate tests and61 full quartet extractions pass. For the
two previously positive conditional profiles, extracting paid count,
incoming L0,L1, and outgoing L0 gives FLOAT regional values163.27 and177.23,
instead of331.49 and322.29. These particular comparisons still fix twelve
inactive counts131 and four counts132; they are not full regional coverage.

## A directed implementation

At b=2^(5/2), each exact quartet atom has the form

    integer * 2^((707*activity-1905*deaths-3072*paid)/1024).

`quartet_closed_directed.py` uses floating optimization only to choose
markers. It rounds each logarithmic marker to an integer multiple of1/1024.
Any such finite choice is valid; convergence is not assumed.

An integer1024th-root calculation gives an upper dyadic endpoint for
2^(1/1024). Upward fixed-point powering bounds every needed fractional
power. The local upper polynomial therefore has rational coefficients.
Normalize that polynomial exactly, round its coefficients upward on a
2^-31 grid, and propagate positive integer arrays. Each stage rounds upward;
an explicit maximum proves that uint64 multiplication and accumulation
cannot overflow. Truncation above a requested nonnegative degree is exact.
Final logarithms use the existing rational interval helper, with the lower
endpoint selected when subtracting a logarithm.

Eight power checks and54 rational coefficient comparisons pass. The cheaper
three-count extraction gives directed conditional regional bounds

    61842669800051 / 274877906944 = approximately224.982322,
    63452025936835 / 274877906944 = approximately230.837126.

These are upper endpoints; the displayed decimals are not the trusted values.
The directed run takes under one second after loading the sources.

## Sum every allocation of the regional quota

Retain the entire paid-count polynomial for one step, while fixing the
incoming and outgoing populations. The original quota2100 then requires
total paid count2100-16*n_Z when n_Z stays constant. A second positive
coefficient extraction sums every allocation across the16 steps. No
restriction to labels131/132 remains.

For the two fixed test populations, the directed regional upper endpoints
are approximately268.207363 and274.473345. Both values exclude the global
row union, bad-event threshold, initialization, and terminal acceptance.

## Allow the population to move

Let the center be (Z,V0,V1,B,L0,L1,L2). For integers a,c in[-radius,radius],
include the population

    (Z,V0-a-c,V1,B,L0+a,L1+c,L2).

The transfer retains every edge consistent with the exact age identities.
In particular, outgoing L1 equals incoming L0. A polynomial matrix sums
all paths that remain in this set during the region. Its coefficients retain
the original-inactive quota exactly. The regional result is an operator
bound with a positive population potential; that potential must be carried
across boundaries, not discarded or reset.

`quartet_closed_neighbors.py` chooses the potential numerically, converts
it to positive rationals, and verifies the resulting row bound exactly.
Twenty-eight rational polynomial-matrix tests pass. There is no dependence
on convergence of the potential proposal.

For sharper bounds, the one-step coefficient recurrence also retains the
incoming L1 count modulo4 or8. This is a necessary condition for the exact
count, so it gives an upper bound. The recurrence uses positive cyclic
integer convolution, not a Fourier estimate. Fifty-three small modular
coefficient tests pass. The no-overflow bound still applies: the wrapped
slices for one term are disjoint.

Selected directed regional operator bounds are:

| Center | Radius | States / edges | L1 modulus | log2 upper |
| --- | ---: | ---: | ---: | ---: |
| (104,156,60,64,32,32,64) | 1 | 9 / 27 | 4 | 261.880048408 |
| (104,176,72,32,16,16,96) | 1 | 9 / 27 | 4 | 278.324274397 |
| (104,176,72,32,16,16,96) | 2 | 25 / 125 | 8 | 291.230704331 |
| (104,120,64,32,16,16,160) | 2 | 25 / 125 | none | -28.661259882 |

Trusted rational endpoints, in the same order, are

    287940158306143 / 1099511627776,
    153010387995851 / 549755813888,
    320211545776577 / 1099511627776,
    -984793390859 / 34359738368.

Every row sums all internal paths and all allocations of quota2100 within
its stated population set. Transitions leaving that set are not included.
The last row is a regional contraction for the weighted closed-set kernel;
it is not a full-code failure probability. No terminal test occurs inside
these regional operators.

The maximum ratio between the largest and smallest potential, in log2,
is approximately10.25,8.92,17.94, and14.15 respectively. These boundary
costs must remain explicit when connecting different covered sets.

## The global score telescopes

The exact age-flow identities give an additional way to identify the
trajectories that matter. Let S_t be the selected count and L0_t the new
selected count. Then

    D_t = S_(t-1)+L0_t-S_t.

Also paid_t = original_inactive_t-Z_t. For r=2100, define the integer
state score

    e(n) = 707*S -1905*L0 +3072*Z -403200.

Using the exact quota in all256 regions, the current closed-set bad-event
condition is equivalent to

    sum_(t=1)^4096 e(n_t) +1905*S_4096 > 81074176.

The initial selected count is zero. The terminal correction is retained;
its maximum is975360. This identity holds for the chosen field multipliers
and threshold, without any stationarity assumption.

`quartet_closed_energy.py` checks2245 exact local instances. In the
25-state neighborhood centered at(104,176,72,32,16,16,96), the largest
state score is-19886. A trajectory confined to that neighborhood cannot
satisfy the global bad-event condition, even with the maximum terminal
credit. Such regions can still occur inside a bad trajectory, compensated
by higher-score regions elsewhere. They must not be deleted from a
global transfer without accounting for that compensation.

This explains why a large unconditional moment bound for a low-score
neighborhood need not identify a likely distance failure. It also changes
the next priority: cover transitions between different score ranges,
rather than continuing to tune isolated low-score profiles.

Six selected stationary profiles with score above the global threshold
per update have FLOAT conditional regional bounds between-0.35 and-2402.59.
These are location diagnostics only. The25-state directed contraction in
the table is one subsequent rigorous regional check in that part of the
state space; it does not establish a uniform bound over the entire range.

## Remaining work and reproducibility

Next retain the telescoping score while composing covered population
sets. Bound transitions leaving these sets, then cover the other row
counts and combine with the existing failure allocation. Initialization
and the one true terminal acceptance remain part of that audit.

The source helpers are the four files named above. The exact seven-type
kernel is unchanged. All runs were serial on Peach, with180-second time
bounds and at most256MiB array bounds. The longest selected run took63
seconds; the wider25-state modular run took28 seconds. They wrote no files.
The first in-memory launch exceeded Windows command length before starting
the remote process; trimming unused module bodies resolved it. No
benchmark, construction, production, or paper changes were made.

From this directory on a host meeting the resource guard, reproduce the
regional cases serially with `python3 -B quartet_closed_neighbors.py
--case first`, replacing `first` by `second`, `wide`, or `score` for the
other table rows. `python3 -B quartet_closed_directed.py --regions` runs
the fixed-population all-quota checks. `quartet_closed_energy.py` checks
the exact score identities and the selected FLOAT location probes.
