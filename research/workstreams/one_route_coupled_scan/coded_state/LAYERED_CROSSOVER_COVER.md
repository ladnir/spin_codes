# Closing the crossover and extending into suffix 193

Follow-up: [actual matching moments](LAYERED_MATCHING_MOMENT_COVER.md) widen
the suffix-193 regions and close neighboring count regions at suffixes 194
and 195. This note preserves the previous checkpoint and its bounds.

This continues the [count-region certificates](LAYERED_COUNT_REGIONS.md).
The outer, regional permutation distribution, inner, and setup distribution
are unchanged. The output-weight threshold remains W=838860, the exploratory
40% threshold for N=2^21. These results do not certify the whole layered code
or near-half distance. The GF32 sandwich remains the whole-code control.
The history counts here are over GF(2^32). The intended scalar-extension
distance guarantee still requires complete base-field coverage.

Two proof improvements suffice for the next step. Fractional interpolation
closes all previously failing selected neighbors. A sharper evaluation of
the existing matching moment bound then closes the suffix-193 region.

## Coverage obtained

As before, x=(z,r,f,e,B) counts empty pairs, 01 pairs, full pairs, exceptional
full pairs, and internal births. In the table, each of r,f,e,B varies
independently around (24568,589632,368520,49136). Every integer z in the
stated interval is included.

| R | z interval | Radius in each other count | Matching envelope | Boxes | log2 good-route contribution upper |
| ---: | --- | ---: | --- | ---: | ---: |
| 192 | 0..49264 | 128 | previous | 29 | -111 |
| 193 | 0..49264 | 128 | refined | 19 | -261 |
| 192 | 0..50000 | 256 | refined | 6 | -758 |

Each bound includes all active-row counts s=1..8192, all choices of active
rows, and all integer count vectors in the region. The local transfer sums
support locations, exceptional branches, components, and amplitudes. The
box union is included. Binary subdivisions are checked for disjointness and
total volume; a failed leaf prevents a complete-coverage result.

The earlier 243-point screen now passes all 243 points. Under the previous
route envelope, its weakest per-profile upper bound is -564; union over
those points gives -556. This screen is not the justification for the table:
the affine box bounds and convex active-row bounds cover the interiors.

The table bounds expected low-weight nonzero-word counts intersected with
the indicated good-route event. They are not whole-code margins. The
refined route event implies the previous route event, so the previous
certificate remains valid when restricted to the refined event. The route
failure budget is added once, not once per count box or marker choice.

## Fractional interpolation, with exact rational local weights

The previous implementation restricted the interpolation exponent m to
integers. That restriction was computational, not a premise of the proof.
At the worst neighbor and s=8004, optimized integer choices m=24 and m=25
give upper bounds +307 and +327. Neither closes that slice.

Write q=2^32. Choose any rational tau in [1,q] and define

    m = log2(tau),       theta = m/32.

The same interpolation inequality applies because 0<=theta<=1. The exact
local weights can be expressed using tau without irrational arithmetic:

    2^((32-m)*gain-m*zeros) = q^gain / tau^(gain+zeros),
    beta = 1-p+p*tau/q,
    ambiguous-input factor = p+beta/tau.

Terminal component weights are (q/tau)^gain. Every local row entry and row
maximum is therefore rational. No rounded exponent is inserted into these
operators.

The global log bound still contains m. Directed rational-log arithmetic
provides lower and upper scaled bounds m_lo,m_hi. The implementation uses
m_hi in the nonnegative b, P, f, and route-dimension contributions. It uses
m_lo in the subtracted z and (R-128)s contributions. All counts are
nonnegative. Rounding m_hi upward to an even scaled integer also permits
the existing D/2 dimension cap. Each coefficient change is outward.

The retained weight tau=6369051672525773/268435456 is near 2^24.5. At the
same s=8004 slice its verified upper bound is -1257. Together with the
refined m=25 weights and previous bank, it closes the entire active-row
range for all 243 neighbors. The proof weights are not encoder parameters.

The search is bounded and warm-started from retained weights. Numerical
convergence only produces a proposal. Recomputed rational local maxima and
directed global bounds determine acceptance. The historical integer-weight
replay remains available and returns its previous result.

## Why suffix 193 needed a different refinement

At the central counts, suffix 193 did not close through interpolation alone.
A bounded continuous-m search at s=8061 returned m about 24.48975 and a
verified upper bound +5685. This is a failed upper bound, not a low-weight
word or evidence that the construction fails.

The limiting route dimension came from the internal-edge envelope near
u=2411 outer rows. The old envelope permits 92718 internal matching edges
there. Its Hoeffding/Bernstein simplification discards part of the already
proved moment information.

For a fixed u-row subset in one uniform matching, the earlier factorial-
moment argument proves a binomial moment bound. Across R independent
regions, define

    n = R*floor(u/2),
    p = (u-1)/(M-1) if u is even, and u/(M-1) otherwise.

If E is the number of internal edges, then for every rational t>1,

    E[t^E] <= (1+p*(t-1))^n,
    Pr[E>=h] <= (1+p*(t-1))^n / t^h.

Let c_u be the bit length of binom(M,u). A sufficient threshold is

    h >= (n*log2(1+p*(t-1))+c_u+77) / log2(t).

We compute a ceiling using an upper numerator and a positive lower
denominator, both from directed rational logarithms. The retained tilt
set is {5/4,3/2,7/4,2}. For each (R,u), choose the smallest threshold among
these and the previous envelope. This choice depends only on R,u, not on
the sampled route. It needs no union over tilt choices: the selected
threshold itself has the required tail bound. A cap at n is deterministic.

Use complement symmetry above u=M/2. The same union over row subsets,
cardinalities, and suffixes gives matching-event failure at most 2^-56.
The head envelope and its budget are unchanged. The refined matching
event replaces the old matching event; its smaller envelope implies the
old one, so no extra matching-event failure budget is needed.

At R=193,u=2411, the new edge cap is 91116. After maximizing over all u,
the central route-dimension cap drops from 156715 to 155878. The maximum
can move, so the dimension gain is not the entire edge-cap reduction.
The central full-profile log2 upper bound becomes -10674. The box proof
then extends this to the suffix-193 region in the table.

## Current limits and next test

With the refined route envelope and retained weights, the central profiles
at R=194 and R=195 still have upper bounds +9195 and +29245. At R=193,
two bounded region attempts also remain incomplete: increasing the four
radii to 256, or extending z to 49500 while retaining radius 128. Their
failed cells are retained, not silently removed. Neither attempt is a
counterexample.

The next step is to use the refined route envelope throughout the existing
count cover, then target R=194 and the uncovered boundary of R=193. First
test whether weights optimized for the new dimension caps close those
regions. If they stall, retain additional nonzero-input relations jointly
with outer constraints. Do not change the construction merely because an
older marker set fails on a new suffix.

## Reproduction and validation

Run serially from the repository root:

    python -B research/workstreams/one_route_coupled_scan/coded_state/test_layered_fractional_weights.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/test_layered_matching_chernoff.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/cover_layered_location_regions.py --crossover

The fractional tests compare 14 integer-weight realizations exactly and
check 684 rational local rows against an independent floating expression.
Another 160 direct-expression checks exercise both route envelopes.
The matching tests enumerate small perfect matchings, verify 23 factorial
moments and 40 generating-function bounds, and check 240 rational thresholds.
The K20 envelope is checked never to exceed the old envelope and to satisfy
the complement identity. The algebra above, not toy sampling, supplies the
distribution-wide argument.

The original integer-weight and local-history regressions pass. The old
single-profile verifier still returns -5153, and both frozen certificate
helper hashes are unchanged.
Compact rational proposals are in layered_crossover_markers.py. Search
commands print proposals; they write no raw files. No benchmark, library
change, paper update, commit, push, or promotion is part of this work.
