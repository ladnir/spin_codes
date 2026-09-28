# Uniform coordinate shuffles with GFNI: better counting, about 6.49 ms

Follow-up: VECTOR_ROUTE.md measures a code-equivalent implementation at
about 5.40 ms, or 1.76--1.78x against matched production controls. The timings
below document the preceding implementation and its proof/implementation tradeoff.

2026-09-24. This iteration recovers a stronger permutation distribution at
about 9% extra encoding time relative to canonical blocks. It also proves
a sharper one-group delayed-activation bound. A subsequent output-weight
calculation bounds all rank-one messages in one group with 45.24 bits of
margin at 10% distance. The 2x objective and the full
10%-distance certificate remain open. Production code and the paper are
unchanged.

## Distribution and layout

Use fixed groups of four outer rows. Independently for each group, sample
one uniform permutation of all 256 coordinates and share it among those
rows. Divide the permuted columns into 64 macroregions of four columns.
Each macroregion independently permutes the 2048 row groups. Within each
column, independently permute the four row labels. This is the rectangular
g=4,c=4 distribution from SECOND_ITERATION.md, not the canonical-block
distribution. The mathematical bounds below concern these ideal uniform
draws. Timing seeds select concrete instances of the implementation.

The encoder writes each 16-element bundle contiguously in shuffled column
order. A four-row group occupies 1024 elements plus one cache line of
padding. A table maps each canonical BCH column to its shuffled offset.
All routing and permutation tables are prepared before timing.

GFNI preparation can read each dense BCH input once through that table.
However, the first implementation is slow. Direct scattered reads from
the large intermediate buffer give about 8.5 ms for the full encoder,
versus about 6.0 ms for canonical blocks. A phase diagnostic attributes
roughly 3.3 ms to inner/routing in both cases. The corresponding BCH phase
is about 5.2 ms with direct shuffled reads, versus 2.6--2.7 ms for canonical
blocks. These instrumented means are diagnostics, not replacement timings.

The winning implementation first copies one complete 16 KiB group
sequentially into an aligned local buffer. GFNI then reads shuffled
coordinates from that buffer. This adds a copy but removes much of the
scattered-read penalty. It keeps the same code and permutation; no
probabilistic assumption changes. The buffer is stack-owned in a
non-coroutine helper and reused across all groups. There are no new
encode-time heap allocations or per-row runtime function-pointer calls.

## Performance evidence

The host and protocol match GFNI_BCH.md: Peach Ryzen 7950X, GCC 15.2,
CPU 15, Release/znver4, K=2^20, 128-bit elements, half-rate BCH[256,128],
IMT(128,19), mask seed 2. All three shared benchmark locks are held.
Each process runs three warmups and 101 in-place encodes, without resetting
input. Setup, allocations, initialization, validation, and checksums are
excluded. No benchmarks run concurrently.

A seed-17 screen motivated the sequential copy:

| Full encoder, uniform shared coordinate shuffle | Median (ms) |
|---|---:|
| Direct shuffled GFNI loads | 8.497175 |
| Additionally prepare systematic inputs | 8.130090 |
| Prefetch next group to L2 | 7.973657 |
| Prepare systematic inputs and prefetch | 7.598878 |
| Sequential local copy, then prepare systematic inputs | 6.768126 |
| Same copy with next-group prefetch | 6.984912 |

An earlier separate screen of full canonical-order repacking gave
7.938532 ms. Prefetching alone is insufficient. The copy result supports
a cache-locality explanation, but no hardware-counter analysis isolates
the individual causes.

The final comparison adds a version that reads systematic inputs directly
from the local copy. Three processes per mode and seed, with mode order
reversed in the middle repetition, give these medians of process medians:

| Full encoder | Route seed 1 (ms) | Route seed 17 (ms) |
|---|---:|---:|
| Production control | 9.621754 | 9.625170 |
| Canonical blocks with GFNI | 5.944528 | 5.924991 |
| Uniform shuffle, local copy plus systematic preparation | 6.666707 | 6.649965 |
| Uniform shuffle, local copy and direct systematic reads | **6.484326** | **6.494094** |

The last row is 1.4838x and 1.4821x faster than the production controls.
Its process medians range from 6.480108 to 6.496579 ms for seed 1, and from
6.492791 to 6.539810 ms for seed 17. It costs 9.08% and 9.61% more time
than canonical blocks in this comparison. Canonical performance remains
about 1.6x; small changes from the preceding run are not claimed as a new
optimization gain.

## Exact delayed-activation bound

The benefit of the uniform coordinate shuffle is a simpler support law.
Fix a nonzero tuple of four BCH words. Let h be its binary rank and u the
size of its union of bit supports. The shuffled support is a uniform
u-subset of 256 positions, regardless of the original coordinate pattern.
Let T_h(u) count rank-h tuples with union support at most u. We can use the
existing bit-support CDF caps from bch_joint_support.py directly; no fixed
canonical-block enumerator is needed.

Consider one active four-row group, with all other messages zero. Let J
be the number of nonempty macroregions in a prefix of ell macroregions.
The local feedback census gives q=1/9216 as an upper bound on the
zero-feedback probability for every nonempty column-weight shape. For
rank one use q=0: all nonzero columns have the same weight, and neither
zero-feedback shape has that property.

Condition on the unordered sets of original coordinates assigned to each
macroregion. Within each macroregion, their order remains uniform. Row
permutations and group positions are independent across macroregions.
The group position gives a uniform choice among the eight aligned
16-bit windows of an IMT epoch. Thus the probability that the state
remains zero throughout the prefix is at most q^J. Empty epochs preserve
zero under any transvection, so this statement does not require averaging
over the transvections.

For a fixed union weight u, averaging over the uniform support gives

\[
f_{\ell,q}(u)=\frac{[z^u]\bigl(1+q((1+z)^4-1)\bigr)^\ell
                    (1+z)^{256-4\ell}}{{256\choose u}}.
\]

Here [z^u] extracts the coefficient of z^u. An empty prefix block contributes
one; a nonempty block contributes q times its support-counting polynomial.
We use 0^0=1 when q=0. The function decreases with u: adding a support
position cannot decrease the number of occupied prefix blocks.

Summation by parts combines this decreasing function with the upper CDF
caps using only nonnegative coefficients. Union-bound over all four ranks
and the 2048 possible active groups. No independence between messages or
between these group events is needed. Exact arithmetic gives:

| Prefix length | log2 one-group zero-prefix union bound |
|---:|---:|
| 48/64 | -13.037 |
| 49/64 | -23.631 |
| 50/64 | -34.332 |
| 51/64 | **-45.114** |
| 52/64 | -55.884 |
| 53/64 | -66.709 |

At 51 macroregions the code verifies the strict inequality against 2^-45
using exact rational arithmetic. The rank-one contribution has log2 bound
-51.167; rank four dominates the total. The corresponding canonical-block
analysis did not cross 40 bits until 57 macroregions.

This is not a distance certificate. It controls a completely zero state
prefix for messages contained in one group. Later cancellation, output
weight, and messages spanning multiple groups remain open. The roughly
20% remaining length after 51 macroregions also leaves little room for
a 10%-distance argument without a full output-weight calculation.

`random_block_support.py` checks the occupancy polynomial against exhaustive
subsets at small lengths, including q=0 and q=1. It verifies 67 rational
shortened-code witnesses and authenticates the same 163 BCH dependency
files as the earlier support calculation. It does not replay every
historical BCH proof. Run it with Python and `-B`; it writes no data files.

## Validation and reproduction

The two new GFNI preparation variants match the original BCH routine on
all 131072 input basis vectors, separately for coordinate-permutation
seeds 1 and 17. The six prior GFNI variants are retained in the same
check. Every timed encoder compares its complete in-place buffer against
the materialized-route reference, including the unchanged suffix.

ASan/UBSan checks pass at K=2^14 for all eight shuffled-route implementation
variants, for seeds 1 and 17. The two-column copy variant also passes at
both seeds. The four-column route passes dense-transpose and adjoint checks,
and all four existing library tests pass. The remote log is
`/tmp/spin-locality-JjY7yR/sanitize-mapped-gfni.log`. This does not constitute
a full-size sanitizer campaign.

Use `run_mapped_gfni.sh` for the basis checks and direct-load/repack screen,
`run_warm_gfni.sh` for cache experiments, and `confirm_random_gfni.sh` for
the final comparison. All use the existing C++ experiment and serial locks.
Raw measurements stay under `/tmp/spin-locality-JjY7yR/measurements` on
Peach, outside version control.

The confirmed release executable SHA-256 is
`c8fab538e7c98adf86e24b1c9002ffabe992c5d12cc9b0be623484e4334ed33d`.
The unchanged production library SHA-256 is
`ecbebbfa16d15b4b5460c8811053720d1dd534f2aa81ea5b82757bebd153f0d7`.

## Rank-one output-weight bound

`random_rank_one.py` now bounds the probability that any rank-one message
confined to one group has output weight at most 209715. Here K=2^20 and
N=2^21. The bound includes all 2048 possible groups, every nonempty subset
of their four rows, and all nonzero BCH words. It averages over the ideal
uniform route above and the independent IMT transvections. It does not
cover higher-rank messages or messages spanning multiple groups.

A rank-one message repeats one BCH word in a nonempty subset of rows.
Let a be the number of selected rows. If a macroregion contains b occupied
columns, each occupied column has a ones. Independent lane permutations
and a uniform aligned window give

    n_ab = 8 binomial(4,b) binomial(4,a)^b

equally likely input patterns. The script enumerates all sixteen pairs
(a,b), for 1 <= a,b <= 4. None of these patterns has zero feedback Bx.
It also reconstructs exact output-weight and cancellation counts using
the current E and B maps. Thus the calculation includes returns to zero
after activation; it does not assume persistent nonzero state.

Use the seven-coordinate weighted-measure envelope from GROUP_DISTANCE.md:
zero mass, arbitrary nonzero mass, and five upper density bounds indexed
by the expansion weights 48,56,64,72,80. The density coordinates do not
assert that a conditioned state is uniform. For each a,b, the exact local
counts replace the earlier sixteen-row orbit counts. The arbitrary-state
output moment is bounded by z^(48-ab), for 0 < z < 1. The lazy-cancellation
bound uses the largest syndrome multiplicity and the least corresponding
output weight. Refresh and empty-epoch transfers are unchanged.

Write W_0(z) for an empty epoch and W_ab(z) for an active epoch. A macroregion
has 256 epochs; the active group's epoch is uniform among them, independently
of its window within that epoch. Its transfer bounds are

    R_0 = W_0^256,
    R_ab = (1/256) sum_{j=0}^{255} W_0^j W_ab W_0^(255-j).

For a fixed BCH word of weight w, condition on the first occupied macroregion.
Let l be the number of following macroregions, and b its occupancy. These
events have probabilities

    binomial(4,b) binomial(4l,w-b) / binomial(256,w).

Conditional on this event, the remaining support is a uniform (w-b)-subset
of 4l positions. Define the matrix polynomial

    R_a(x) = R_0 + sum_{j=1}^4 binomial(4,j) x^j R_aj.

The formal variable x counts occupied coordinates; z weights emitted bits.
Starting at zero state, the conditional output moment is bounded by

    M_albw(z) = e_0 R_ab [x^(w-b)] R_a(x)^l 1 / binomial(4l,w-b).

Coefficient extraction selects supports of the required size. The final
vector of ones sums all envelope coordinates. For each a,l,b,w, minimize
z^(-209715) M_albw(z) over the seven fixed tilts in the script and cap the
result at one. Call this upper bound p_albw.

The first occupied macroregion is shared by all row subsets for this BCH
word. Union-bound over the subsets *conditional on this event*, then cap
that union at one. With A_w the authenticated BCH shell cap, the final bound is

    2048 sum_w A_w sum_{l=0}^{63} sum_{b=1}^4
        [binomial(4,b) binomial(4l,w-b) / binomial(256,w)]
        min(1, sum_{a=1}^4 binomial(4,a) p_albw).

Only terms with 0 <= w-b <= 4l contribute. Neither this union bound nor
the subset union requires independence between the corresponding messages.

Separate outward-rounded Arb runs at 192 and 384 bits give

    failure upper < 2.410310e-14 < 2^-45,
    margin approximately 45.237775 bits.

The weight-38 shell contributes about 2.32987e-14; weight 40 contributes
about 7.48106e-16. These are upper-bound contributions, not measured failure
probabilities. The script checks the final 40-bit inequality directly.

Self-tests compare orbit intersection counts with direct mask enumeration,
region averaging with explicit placements, and the coefficient recurrence
with exhaustive small supports. The first-occupied-macroregion counts sum
exactly to binomial(256,w) for every w. Both full runs authenticate the 163
BCH dependency files and rebuild the local census. As before, authentication
does not replay every historical BCH proof.

Run from the repository root:

    python -B research/workstreams/permutation_locality/random_rank_one.py --test-only
    python -B research/workstreams/permutation_locality/random_rank_one.py --precision 192
    python -B research/workstreams/permutation_locality/random_rank_one.py --precision 384

## Next step

Retain both candidates. The canonical route is faster; the uniform shared
shuffle has a stronger counting bound at a modest additional cost.
The subsequent calculation in RANK_TWO.md covers all rank-two tuples within
one group with 53.38 bits. For this four-row distribution, keeping union
support size exact makes a worst-column-shape envelope sufficient for rank
two. RANK_THREE.md subsequently retains a distinguished rank-two subspace
and verifies 40.643 bits for rank three. Together, ranks one through three
give a 40.584-bit partial union. Rank four and multiple active groups remain
open; these results do not establish the full ensemble's distance.
Further performance work must still remove roughly one quarter of the
uniform-shuffle candidate's runtime to reach the matched 2x target.
