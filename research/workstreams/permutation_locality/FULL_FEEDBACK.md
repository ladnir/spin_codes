# Full feedback distributions

The grouped route places j active four-bit windows into distinct locations
among the 32 windows of an IMT epoch. Earlier bounds used exact zero counts
but coarse maximum probabilities for other feedback values. Full Walsh
inversion now gives every feedback probability, for each input-weight
multiset through eight active windows. This improves the proof bound but
does not close the intermediate-occupancy gap.

## Distribution and exact calculation

Fix a multiset of nonzero window weights, with n_b occurrences of weight b
for b in {1,2,3,4}. Choose distinct windows uniformly, assign those weights
uniformly, and choose a uniform lane mask of the assigned weight in each
window. Let X be their 128-bit union and S=BX its 19-bit feedback.

For a character c in F_2^19, let r_i(c) count the feedback columns in
window i on which c has odd dot product. The signed sum over weight-b
masks in this window is the four-bit Krawtchouk value K_b(r_i(c)).
Consequently the coefficient of product_b z_b^n_b in

    product_{i=1}^{32} (1 + sum_{b=1}^4 K_b(r_i(c)) z_b)

is the signed count of inputs of this shape. Its unsigned denominator is

    D = binomial(32,j) j!/product_b n_b! * product_b binomial(4,b)^n_b.

The character polynomials depend only on the histogram of the 32 values
r_i(c). There are 9218 such histograms among the 2^19 characters. We compute
their coefficients once, expand through the character-to-histogram index,
and apply an integer Walsh transform. Division by 2^19 gives the exact
number of inputs with each value of S.

Every polynomial coefficient is a signed count of at most D assignments.
Before each inverse transform, the implementation computes the exact sum
of absolute character coefficients using Python integers. That sum bounds
every Walsh intermediate. It must fit signed 64-bit arithmetic before the
vectorized transform is permitted. The implementation also checks exact
division, nonnegative recovered counts, total count D, the zero count,
and selected forward character identities.

Every one- and two-window distribution matches an independent direct
enumeration. The three-/four-window zero counts match the independent
rank census on all 320 ordered weight shapes. Expansion-weight class
counts and peaks also match the earlier pair census.

All feedback columns have odd weight. Thus wt(S) mod 2 equals wt(X) mod 2.
Each fixed shape is supported on one parity class, not uniformly on all
2^19 states. The census explicitly checks this constraint. Neither the
construction nor these bounds assume uniform feedback.

## Largest nonzero probabilities

Each maximum below is over every input-weight shape and every nonzero
feedback value at that occupancy. The maximizing shape has all weights 4
in these tests.

| Active windows | Exact maximum |
|---|---:|
| 1 | 1/32 |
| 2 | 1/496 |
| 3 | 1/2480 |
| 4 | 3/35960 |
| 5 | 1/28768 |
| 6 | 13/906192 |
| 7 | 5/560976 |
| 8 | 11/1753050 |

The scripts permit requests through ten windows, but only through eight
has been executed here. Every requested run must pass its integer-width
checks; support in the interface is not evidence of a completed census.

    python -B research/workstreams/permutation_locality/full_feedback_census.py --maximum 8

## Use in the existing state envelope

For a fixed shape, let w=wt(X), p_0=Pr[S=0], and p_* be the largest
nonzero probability. Let p_all=max(p_0,p_*). Let
c_v=Pr[wt(ES)=v], for v in {48,56,64,72,80}; zero is excluded.
The two-update kernel has lazy probability alpha=1/4. The other part
refreshes uniformly over the nonzero states. The tilt is lambda>=0.

Starting from zero, the weighted zero mass, nonzero mass, and nonzero
pointwise density are respectively bounded by

    exp(-lambda*w) p_0,
    exp(-lambda*w) (1-p_0),
    exp(-lambda*w) p_*.

For a lazy return to zero, the entering state must equal S. Then the
output is ES+X, of weight at least |v-w| if wt(ES)=v. A mature component
with pointwise density bounded by C therefore contributes at most

    alpha C sum_v c_v exp(-lambda*|v-w|).

For the uniform expansion-weight-v generator, whose class has size N_v,
the corresponding lazy coefficient is at most

    alpha c_v exp(-lambda*|v-w|)/N_v.

Convolution bounds its outgoing pointwise density by

    alpha exp(-lambda*|v-w|) min(p_all,1/N_v).

Fresh-state bounds use p_* for cancellation and p_all for convolution,
with the pointwise output factor exp(-lambda*max(0,48-w)). The exact
c_v values also bound the low-expansion-weight tails created from zero.
All-one penalties and total-input-weight tilts are applied per shape
before maxima are taken.

Fresh and uniform sources also have a refresh contribution to zero.
The implementation retains an upper for that contribution from the
existing uniform-class columns; it does not subtract a loose lazy bound.
All affected coefficients refer to the same existing coordinates and
components. Taking the minimum with an old coefficient is therefore valid.

The optional --full-feedback flag applies these refinements after the
requested transvection-round transformation. Defaults, the fourteen-group
certificate path, production code, and the paper are unchanged.

## Diagnostic effect and remaining gap

For 64 active groups, equal support 128, tilt .032, and all-one penalty
.75, the cached binary64 screen gives:

| Refine through this local occupancy | Log2 contribution upper |
|---|---:|
| Original envelope | +3809.0603 |
| 2 | +3809.0603 |
| 3 | +3548.8634 |
| 4 | +3418.9178 |
| 5 | +3389.9765 |
| 6 | +3380.7392 |
| 7 | +3378.2765 |
| 8 | +3377.7169 |

An independent reconstruction of the outward local operators through six
windows reproduces the binary64 support-128 score within 2e-9 bits. Its
final selected-point score is still binary64, not an outward certificate.
The new outer containment bounds do not change that point. With those
outer bounds and the same fixed tilt/penalty, supports 72, 80, 96, and 160
give +4.8932, +801.6342, +1847.0632, and +4863.8231. These fixed-witness
values do not replace better values from other tilts. In particular, the
earlier support-72 screen remains better.

    python -B research/workstreams/permutation_locality/full_feedback_refinement.py --maximum 8
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 64 --mixing-rounds 2 --shortening-moments --full-feedback 6 --fresh --window-average --multi-average --prefix-flags --fresh-collision --zero-moment --mature-tail 64 --pair-tail --joint-witness --penalties .75 --tilts .032 --probe-supports 72 80 96 128 160

The total gain at the tested point is about 431 bits. Extending from six
to eight windows adds only about three bits there. This suggests retaining
the refinement without extending the census solely for that point.
No full occupancy certificate or 40-bit code guarantee follows.

The next useful test is to couple total outer weight and all-one-column
counts in the same bound. The current envelope can favor different input
shapes in different coefficients, while the outer count remembers only
limited information about those shapes. A joint weighting would test that
loss directly; it is not yet implemented or known to close the gap.
