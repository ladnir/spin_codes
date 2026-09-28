# Average random windows before bounding mature states

The nine-coordinate envelope used a triangle bound for every input window.
That permits every window to overlap an incoming state's support maximally.
The actual window is uniform. Averaging it first gives stronger bounds
without changing the one-column distribution or the encoder.

Fix a nonzero column weight b in {1,2,3,4}. Let X_b be the set of all
128-bit vectors supported in one of the 32 four-bit windows, with weight b.
The input x is uniform in X_b. Let E and B be the authenticated expansion
and feedback maps, and let z=exp(-lambda). Define

    A_b = max_{q!=0} average_{x in X_b} z^wt(Eq+x),
    D_b = max_{y!=0} average_{x in X_b} z^wt(E(y+Bx)+x),
    Z_b = average_{x in X_b} z^wt(EBx+x).

These definitions average over the input before maximizing over the state.
They do not assert that the incoming state is uniform. Conditioning on epoch
assignments leaves the input window independent of the entering state, as
in STATE_MEMORY.md.

For the mature component, M bounds mass and C bounds pointwise density.
The one-input lazy branch has outgoing mass at most A_b M/2. Its outgoing
nonzero density is at most D_b C/2, and its zero mass is at most Z_b C/2.
For the density claim, fix the outgoing state y and substitute q=y+Bx.
Bounding each incoming mass by C leaves exactly the average defining D_b.
Including source-zero terms only enlarges this upper bound.

The refresh branch has zero mass at most A_b M/(2m), where m=2^19-1.
Its class-v density coefficient is at most A_b M N_v/(2m). The existing
class-v lazy-density entry can also be intersected with D_b/(2N_v), since
that source component is bounded pointwise by 1/N_v. Other entries remain
unchanged. The all-one-column factor rho^[b=4] is applied before maximizing
over b.

## Exact finite enumeration

`occupancy_window_average.py` enumerates all 2^19 states and all 480 nonzero
window inputs. It reconstructs E and B from the selected maps. For each
possible output weight w, it upper-rounds exp(-lambda w) to an integer
multiple of 2^-44, using outward arithmetic. All subsequent sums and maxima
use uint64 integers. Each sum has at most 192 terms and cannot overflow.

Python-integer calculations independently replay both maximizing states and
the zero row for each b. Additional outward checks replay selected states.
The maxima are exhaustive, not statistical estimates. Quantization adds
less than 2^-44 to an average, apart from the outward power evaluation;
the verifier uses the actual rounded upper values throughout.

For example, at lambda=.0064 and b=4, the mass, density, and zero averages
are bounded by approximately .7310367720, .6841876347, and .6601308640.
The old common pointwise factor was exp(-.0064*(48-4)), about .7546.

## Results and scope

Combined with fresh-state averages and all-one-column counts, the sixteen-
group class with each support exactly 84 passes 192- and 384-bit outward
replays. Its upper is below 3.618789e-15, or **47.9734146515 bits**.
It includes all ranks and group locations in that class, not all support
vectors. This is separate from complete occupancy coverage.

The corresponding selected-point binary64 scores, after minimizing over
the tested witnesses, are:

| Common support | Log2 contribution upper |
|---|---:|
| 80 | -73.60173 |
| 84 | -47.97341 |
| 96 | +3.35916 |
| 128 | +194.60289 |
| 176 | +39.50305 |

Positive entries remain vacuous. The encoder's distance is not refuted.

A separate ten-coordinate experiment, `occupancy_pair_memory.py`, retains
the nonzero part of two-feedback convolution distributions. Its exact
moment/cancellation component checks pass, but its tested scores are worse
than the simpler fresh-average bound. Its pointwise domination coefficient
loses some of the benefit of an averaged mass. It is not used by the new
certificates and should not be expanded merely because it carries history.

Reproduce the useful refinement:

    python -B research/workstreams/permutation_locality/occupancy_allones.py --fresh --window-average --tilts .0064 .008 .01 --penalties .5 .625 .75 .9 --supports 80 84 96 128 176
    python -B research/workstreams/permutation_locality/occupancy_allones_verify.py --support 84 --fresh --window-average
    python -B research/workstreams/permutation_locality/occupancy_allones_verify.py --support 84 --fresh --window-average --precision 384

The fresh-average checks cover 1,560 exact-map inequalities, with both unit
and one-half all-one penalties. Reference evaluations use 128 more precision
bits than the tested upper coefficients to avoid inconclusive interval
overlap at tight entries. No inequality was weakened to resolve such overlap.

Full-support results combining these operators with improved interval
counting are recorded in [CDF_COVER.md](CDF_COVER.md).

## Several occupied windows

`occupancy_multi_average.py` improves mass and refresh bounds when two or
more windows are active. Fix the entering state, its expansion weight v,
and the column weights a_1,...,a_j. The occupied windows are uniform and
distinct; their masks are independent and uniform at the specified weights.
For z=exp(-lambda), the output-weight moment is at most

    z^v product_i ((1-v/128) z^a_i + (v/128) z^-a_i).

Clipping this expression at one is also valid. To derive it, let b_k count
the expansion bits in window k, for k=1,...,32. For a weight-a mask, let R
be its overlap with those bits. The function phi_a(b)=E[z^(-2R)] is positive
and increasing in b: nested sets of b marked positions give a coupling
under which R increases.

For positive increasing functions f_i and ordered samples X_i without
replacement from a scalar population, E[product_i f_i(X_i)] is at most
product_i E[f_i(X_i)]. Here is an induction proof. Conditional on X_1=x,
let h(x) be the expected product of the remaining factors. Removing a larger
population element leaves a smaller remaining population. Couple the two
remaining samples by matching their common elements and the one replacement,
then sample the same ordered distinct indices. Every coordinate in the
smaller population is no larger, so h decreases. Thus f_1 and h have
nonpositive covariance. Apply induction to the unconditional remaining
sample, which is again an ordered sample without replacement from the
original population. This proves the product bound, including repeated
population values and unequal functions.

Apply that bound to f_i=phi_{a_i}. The convex chord on R in [0,a] gives
z^(-2R) <= 1+(R/a)(z^(-2a)-1). A uniform window and mask have mean overlap
a v/128. Multiplying the resulting factors by z^(v+sum_i a_i) yields the
displayed bound. The code checks 12,460 exact rational inequalities on
small populations, including unequal column weights.

The implementation averages over the fresh-state distributions and applies
the bound separately to each uniform expansion-weight class. It intersects
the resulting entries with the existing envelope. Density entries remain
unchanged: fixing an outgoing state makes the incoming state depend on Bx,
so the same independence argument does not apply.

Used alone, this refinement did not extend the certified occupancy range.
At twelve active groups, a 500-split binary64 full-support screen improved
from log2 upper -17.93 to -22.73. Selected sixteen-group, equal-support
scores improve to -80.55, -55.35, -5.10, +177.19, and +6.75 at supports
80, 84, 96, 128, and 176. Those selected-point screens use the separately
checked dual-distance premise for outer counts; they are not full occupancy
certificates. Combining the multi-window bound with joint witnesses and
parent retention subsequently closed twelve groups at 57.19597 bits in
192- and 384-bit outward runs; see CDF_COVER.md. Larger occupancies remain
open. The selected sixteen-group points do not cover those occupancies.

Reproduce the new inequality checks and the full-support screen:

    python -B research/workstreams/permutation_locality/occupancy_multi_average.py
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 12 --max-splits 500 --fresh --window-average --multi-average --penalties .5 .75 1 --tilts .0032 .004 .005 .0064 .008 .01
