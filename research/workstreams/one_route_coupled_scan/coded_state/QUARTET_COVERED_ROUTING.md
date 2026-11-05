# Covered coordinates and the remaining regional quota loss

2026-10-27. Continue the unchanged quartet construction and the whole-code
45%/40-bit goal. The preceding [full-weight audit](QUARTET_FULL_WEIGHT_FLOW.md)
was progress: it exposed losses hidden by the restricted-band certificates.
This round changes proof bounds and diagnostic helpers, not the encoder.
All large numerical results below are FLOAT screens.

## Integrate genuinely free labels

Fix a regional support history and the quartet partitions. A quartet with
nonzero incoming state covers all four input coordinates. An empty incoming
quartet covers only its outgoing nonzero coordinates. Labels outside this
covered set do not enter its local feasibility or selected-zero polynomial.

Let M=8192 and fix an inactive-row set of size r. Under an auxiliary
Bernoulli-p labeling, write E for the number of uncovered positions in the
region and c for the number of inactive labels on covered positions. The
remaining quota has probability Pr[Bin(E,p)=r-c]. Therefore any nonnegative
integrand F that depends only on covered labels satisfies

    E[F | total inactive labels=r]
      <= kappa_E(p,r) E[F],

where

    kappa_E(p,r)
      = max_j binom(E,j) p^j (1-p)^(E-j) / pi_r(p),
    max(0,r-(M-E)) <= j <= min(E,r),
    pi_r(p) = binom(M,r) p^r (1-p)^(M-r).

The maximum is attained at a binomial mode, clipped to the displayed
interval. In particular kappa_M=1. This argument applies to each fixed
history before summing histories; it does not assume that the eventual
message was chosen before setup.

For a positive free-position marker b, define

    a_b = max_(0<=E<=M) [log2 kappa_E(p,r) - E log2 b].

Then kappa_E<=2^a_b b^E. Insert b^(4-v) into an empty incoming quartet
with outgoing weight v, and insert one into a nonempty incoming quartet.
After 16 updates, multiply by 2^a_b. The full-weight positive operator now
tracks the actual free-position count through its marker, without adding a
third coefficient array. A family of b values gives alternate upper bounds
on the same regional continuation. Their coordinatewise minimum is valid.

This improves the older cap based only on the number of nonquiet updates.
It still maximizes over c instead of enforcing the realized value. That
remaining relaxation is substantial.

## Full-weight numerical diagnostic

`quartet_covered_route.py` implements the cap and its affine envelopes.
At r=2100, theta=7/40, ell=0.00377, and p=112627/167924, the markers
b=0.6 and b=0.9 give a full-weight log2 bound of approximately +622276.17.
The earlier six-iteration screen gave +685155.11. Neither result provides
a nontrivial probability guarantee.

The two affine envelopes touch their exact combinatorial caps at E=4654
and E=3311, respectively. These are envelope parameters, not measured
coverage in actual encoded messages.

To locate the overcount, marked polynomial derivatives retain conditional
means of free positions, deaths, and active covered positions in the local
count. A left/right positive-eigenvector calculation gives the artificial
stationary flow of one selected transfer kernel. For b=0.9:

| Quantity per update | Artificial mean |
| --- | ---: |
| State weight | 488.54909 |
| Deaths | 14.89137 |
| Uncovered positions | 0.00829 |
| Active covered positions | 491.36130 |
| Inactive covered positions, by subtraction | 20.63041 |

The stationarity residual is below 4.5e-15 in FLOAT. The greatest mass is
at weights 487--491. These are diagnostics of an upper-count kernel, not
statistics of the encoder or a counterexample.

A real region must contain exactly 2100 inactive labels, averaging 131.25
per update. Since this artificial flow leaves almost no coordinates free,
its small inactive-covered count cannot persist throughout a real region.
The cap maximizes the remaining-label probability separately from that
count and therefore admits this mismatch. This identifies a concrete loss
in this envelope, not the sole cause of every previous failed bound.

## Keep all regional endpoint weights

`quartet_regional_matrix.py` forms the complete 513-by-513 regional transfer
before taking minima over p. For each p, raise the positive one-update
kernel to its 16th power and divide by pi_r(p). Each matrix bounds the same
fixed-quota continuation, entry by entry. Taking their entrywise minimum
is therefore valid; composing 256 such regions preserves free endpoints.
This removes the previous coarse bins, but still relaxes the quota inside
each region by the full conditioning factor.

Positive matrix products use row/column scaling in small blocks. Entries
with tiny scaled products are recomputed in log space. No numerical zero
from underflow is silently accepted. This implementation is still FLOAT,
not directed arithmetic.

At r=2100 and coverage marker one:

| theta | ell | p choices | log2 upper-bound screen |
| --- | --- | --- | ---: |
| 0.175 | 0.00377 | rho, 112627/167924, 0.94 | +686335.42 |
| 0.175 | 0.00377 | nine odds-spaced choices | +572899.41 |
| 0.5 | 0.00377 | nine odds-spaced choices | +3238777.99 |
| 0.175 | 0.001 | nine odds-spaced choices | +739302.01 |

Here rho=r/M. The odds-spaced grid uses
p_i/(1-p_i)=[rho/(1-rho)]*q^(theta*i/8), for i=0,...,8.
Its endpoints interpolate between no covered-label field penalty and the
fully covered regime. Larger theta and a smaller output marker do not
repair the loss in these bounded comparisons. No parameter optimum is
claimed. Only one inactive-row count was tested.

## Checks and execution

- 275 free-count cap values agree with exact rational enumeration at small
  region sizes. Each tested affine envelope bounds every free count there.
- Empty-quartet label sums and marked activation/reset identities pass.
  Marking incoming weight recovers the exact extracted incoming count.
- 567 matrix-product entries agree with direct log-space sums, including
  110 deliberately exercised tiny-product fallbacks.
- Each full regional matrix agrees with 16 independent vector applications
  to within 3e-11 bits. No tiny-product fallback was needed in those runs.

Commands relative to this directory:

    python -B quartet_covered_route.py --tests-only
    python -B quartet_covered_route.py --markers 0.9 --flow
    python -B quartet_regional_matrix.py --p-levels 9
    python -B quartet_regional_matrix.py --theta 0.5 --p-levels 9
    python -B quartet_regional_matrix.py --ell 0.001 --p-levels 9

All runs are serial with a 60-second watchdog and one BLAS/OMP thread.
The free-coordinate helper has a 64 MiB numeric-array budget with optional
marked diagnostics; the regional-matrix helper has a 32 MiB budget.
Neither creates generated files.

The local disk fell below the 20 GiB guard after the first coverage run.
The next local run stopped before tests or numerics. Remaining tests ran
on the configured Peach host from a compressed in-memory source bundle,
with no uploaded files or caches. The original `choices` function was
loaded by AST from its source, avoiding an unused SciPy import on Peach.
Other helper modules were loaded unchanged. Nine-p matrix screens took
about 31--32 seconds there. These are proof-computation times, not encoder
benchmarks. No files were deleted to bypass the local guard.

## Next discriminating step

Retain the inactive-label total inside the regional state-path calculation,
rather than maximizing its missing quota after the fact. A useful first
test should reject the dense artificial flow above. Exact or directed
coefficient extraction is one route; a conditional bound retaining the
residual quota is another. Keep rowwise RS feasibility as a separate
constraint, rather than interpreting the aggregate coverage floor as its
replacement. Whole-code closure still requires all inactive-row counts,
all histories, and a directed final union.

## Continuation: count intervals, exact profiles, and an activity cutoff

`quartet_quota_bins.py` retains inactive-count intervals during regional
composition. Write K_p(u,v)=sum_t B_t(u,v)p^t(1-p)^(512-t), with B_t>=0.
For a bin I and reference rho=r/8192,

    sum_(t in I) B_t(u,v)rho^t(1-rho)^(512-t)
      <= K_p(u,v) max_(t in I) [(rho/p)^t ((1-rho)/(1-p))^(512-t)].

The maximum is at an endpoint. Take minima over p before composing these
bin-indexed kernels. If each bin has width b, a 16-update history whose bin
indices sum to J can have exactly r labels only if

    ceil((r-16(b-1))/b) <= J <= floor(r/b).

Sum only those J, then divide by pi_r(rho). Truncating the other indices
is safe because they cannot meet the quota. The last bin may be shorter;
using the displayed condition only enlarges the admitted set.

The method excludes the previously identified persistently low-label
paths, but its interval and coefficient bounds introduce substantial slack:

| State-bin width | Count-bin width | Best three-iteration log2 screen |
| ---: | ---: | ---: |
| 8 | 16 | +785175.04 |
| 1 (all 513 weights retained) | 16 | +708280.34 |
| 2 | 8 | +647300.91 |

These use r=2100 and the preceding nine-p grid and markers. None improves
the earlier +572899 full regional screen. Thirty components of exact small
polynomial examples validate the interval upper bound and quota extraction.
The finer runs took about61--63 seconds on Peach, under a 120-second
watchdog and 128 MiB numeric-array budget. They wrote no files.

A cleaner diagnostic is available for nearly full fixed-weight profiles.
`quartet_exact_profile_quota.py` retains incoming weight, outgoing weight,
and inactive count directly. Using complement supports, weight489 requires
only24*24*513 cells. The regional quota is extracted as a coefficient,
without a count interval or a conditioning maximum.

For state weight489 after the first update, arbitrary turnover, r=2100,
and first activation at update1, its FLOAT log2 bound is -377912.25.
The best matched nine-p Cauchy bound for this same profile is -376199.63.
Thus exact quota extraction saves about1712.62bits here; it is not a new
whole-code certificate, nor an enormous improvement over a correctly
retuned fixed-profile bound. Eight small diagonal comparisons and eight
artificial-p invariance checks pass. The full diagnostic takes about15s,
with a 32 MiB numeric-array budget and no generated files.

The more consequential new result is the
[uniform low-activity certificate](QUARTET_LOW_ACTIVITY.md). It bounds
arbitrary low-total-activity histories across every active-row count below
2^-44 using exact integers and rationals. This includes many sparse
histories that the common output marker priced poorly. Next impose its
activity threshold in the remaining distance count, rather than keep
refining interval widths while still counting those histories.
