# Retaining the mean overlap in a cancellation bound

The feedback census records the expansion weight, but the triangle bound
allows every input bit to overlap a one in that expansion. Exact overlap
means constrain how often this extreme cancellation can occur. The
following local bound uses that constraint without enumerating the full
joint output distribution. It does not change the encoder.

Fix a packet-weight shape h. Let X be uniform over all inputs of that
shape, including distinct-window placement and uniform masks within each
window. Write D for the number of these inputs and W for their common
weight. The fixed linear maps are B: F_2^128 -> F_2^19 and
E: F_2^19 -> F_2^128. The authenticated E is injective.

## Overlap and its character sum

For a target state t, define

    V_t = weight(E(t+BX)),
    H_t = |support(X) intersect support(E(t+BX))|.

Then weight(E(t+BX)+X) = V_t+W-2H_t and 0 <= H_t <= W.
Let a_i be the state character that gives output coordinate i of E.
Define the signed integer

    S_i = sum_(x of shape h) x_i (-1)^(a_i dot Bx).

Because every x has weight W,

    E[H_0] = W/2 - sum_i S_i/(2D),
    E[H_t] = W/2 - sum_i (-1)^(E(t)_i) S_i/(2D)
           <= W/2 + sum_i |S_i|/(2D).

No independence between X and BX is used. These identities concern the
unweighted input distribution; the output tilt is applied only afterward.

Each S_i can be computed by a small multivariate character polynomial.
In a four-bit window w, let r_w count the feedback columns on which a_i
is one. The factor for that window is

    1 + sum_(b=1)^4 K_b(r_w) y_b,

where K_b(r) is the signed sum over all weight-b masks in that window.
Equivalently, K_b(r) is the coefficient of z^b in
(1-z)^r (1+z)^(4-r). For the window containing coordinate i,
replace that factor by

    sum_(b=1)^4 (-1)^d_i K_(b-1)^(3)(r_w-d_i) y_b,

where d_i = a_i dot B(e_i) and the superscript means a three-bit window.
This replacement restricts the chosen mask to x_i=1. Extracting the
coefficient for the shape multiplicities gives S_i. All arithmetic here
is integral; no sampled states or floating polynomial coefficients enter.

## A positive output-moment bound

Fix 0 < z <= 1. Convexity on 0 <= H <= W gives

    z^(v+W-2H)
      <= (1-H/W) z^(v+W) + (H/W) z^(v-W).

Let p_v = Pr[V_t=v] for v>0, and let h_v = E[H_t 1_(V_t=v)].
The source-zero case has V_t=0 and H_t=0, so omitting it loses no
overlap budget. If mu upper-bounds E[H_t], the tilted moment restricted
to nonzero source states is at most

    sum_(v>0) [(p_v-h_v/W) z^(v+W) + (h_v/W) z^(v-W)],

maximized subject to 0 <= h_v <= W p_v and sum_v h_v <= mu.
The coefficient of h_v is positive and decreases with v. Therefore the
maximum is obtained by allocating the budget to increasing v, up to
each capacity W p_v. This allocation is the same for every z in (0,1].
It yields nonnegative coefficients of the powers z^(v-W) and z^(v+W).
Those positive coefficients permit straightforward outward evaluation.

For t=0, the existing exact feedback census supplies every p_v, and the
character calculation supplies the exact mu=E[H_0]. Multiplication by
alpha=1/4 gives the lazy C-to-Z coefficient. Here C bounds the incoming
pointwise density. The argument includes every nonzero source state;
it does not assume that the incoming state is uniform. Apply rho^n4 to
each shape, where n4 counts its weight-four packets, before maximizing
over shapes of the same occupancy. The result may tighten the existing
scalar C-to-Z bound before any coupled mass-column replacement.

For the current cutoff of ten windows, W<=40 and every nonzero expansion
weight is at least 48. In this range the chord upper is no larger than
the corresponding triangle upper, which replaces every H_t by W.
This comparison alone does not establish an improvement over all other
existing local bounds or over a chosen coupled-column mixture.

## Complete local calculation and scope

`candidates/overlap_mean/mean.py` computes both overlap quantities for
all 1000 nonempty shapes through ten windows. Its character polynomials
use unbounded integers. For all-four shapes, selected values are:

| Occupied windows | Input weight W | Exact mean at target zero | Universal translated upper |
|---:|---:|---:|---:|
| 5 | 20 | 503173/50344 | 1009979/100688 |
| 7 | 28 | 5890847/420732 | 23576737/1682928 |
| 10 | 40 | 346845/17342 | 430132033/21504080 |

These means are close to W/2, not W. The translated upper at ten windows
is less than 20.002346. All values average over the stated packet law;
they are not assertions about each input or a sampled encoder's distance.

Exhaustive tests on three small pairs of linear maps check every shape
through three windows and every target state. They check 648 exact
rational moment inequalities; 290 strictly improve their triangle
comparison. The full-map calculation also checks all counting denominators
and coordinatewise signed-count bounds.

```sh
python -B -m unittest discover -s research/workstreams/permutation_locality/independent_rows/candidates/overlap_mean -p 'test_*.py'
python -B research/workstreams/permutation_locality/independent_rows/candidates/overlap_mean/mean.py
python -B research/workstreams/permutation_locality/independent_rows/candidates/overlap_mean/compare.py
```

The C-to-Z adapter is separate from the verifier and the baseline local
family builder. Its selected-point replay and the translated C-to-C
refinement below are now implemented. Complete support covers using
these additional overlap bounds remain to be done.

The full-map local comparison has completed at tilts .068, .072, and .076
with rho=.9. At tilt .072, the ratios of the new shape-maximized C-to-Z
upper to the triangle upper are at most .673311, .658355, .648676,
.642354, .638292, and .635612 for occupancies five through ten.
These are ratios of bounds, not transition probabilities. Other local
bounds and mass-column alternatives may already dominate part of the gain.

The four mean tests pass, including the 648 exhaustive comparisons above,
outward replay against a higher-precision calculation, penalty-before-maximum,
unchanged-entry checks, and rejection of incomplete shape coverage.
Four translated-density tests and four refinement-cache tests also pass.
The isolated selected-point driver applies this refinement after loading
an unchanged exact local memo. It refuses to rebuild a missing memo and
does not store overlap-derived matrices under the baseline cache key.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/overlap_mean/replay.py --epoch-cache tmp/independent-row-local-families --groups 80 --support 200 --tilt .072
```

The completed 192-bit outward replay gives log2 upper +133.142758,
rounded upward, for this selected event. It retunes the coupled
mass-column alternatives so that the old fixed replacement does not
automatically discard the improved C-to-Z entry. The result is still
vacuous as a probability bound and is not full occupancy coverage.

## Extension to the density column

The same argument can use translated class masses. Let n(r) count the
feedback values of this shape. The exact class mass for target t is

    p_v(t) = sum_(s: weight(E(s))=v) n(t+s)/D.

The current density calculation already obtains these translated counts.
Combine them with the universal overlap upper from the character sums,
and perform the same increasing-v allocation separately for each t.
Taking the maximum over nonzero t bounds C-to-C. Unlike the older
triangle estimate, this calculation cannot spend a full overlap W in
every expansion class simultaneously.

Upward-rounded feedback counts also suffice. Suppose n(r)<=u n'(r),
where u is a positive integer and D'=sum_r n'(r). Let G'_v(t) be the
translated class counts for n'. Then p_v(t)<=u G'_v(t)/D. Write mu for
the universal upper on E[H_t] and set J=ceil(mu D/u). At each target,
allocate at most J units across classes in increasing v, with class
capacity W G'_v(t). If j_v(t) denotes this allocation, a valid upper is

    u/(D W) * sum_v [
        (W G'_v(t)-j_v(t)) z^(v+W) + j_v(t) z^(v-W)
    ].

To justify the larger class masses, express each chord contribution as
its positive baseline p_v z^(v+W) plus a nonnegative overlap term.
Increasing the baseline mass and each overlap capacity cannot decrease
the maximum. Rounding the overlap budget upward preserves this order.
The rounded count measure need not be a probability law, and its
normalizing factor u/D must remain in the bound.

For W<=40, all displayed output exponents are nonnegative for the actual
expansion classes. Dyadic upward powers Q_w/2^b therefore satisfy
0<Q_w<=2^b. Every weighted sum is at most W D' 2^b. Choosing b so that
this quantity is below 2^63 permits guarded exact integer evaluation of
every target, followed by the final outward scalar multiplication.
The allocation uses only nonnegative integer arrays and a single integer
budget J; it does not need a floating-point optimization at each target.

`candidates/overlap_mean/density.py` implements this calculation. A bound
for a single uniform expansion class may use the whole overlap budget
in that class; it does not reuse the joint greedy allocation. Each such
bound is divided by the class cardinality. Small-model tests compare all
targets and classes with direct enumeration, including deliberately coarse
upward count rounding.

The full-size census passes all 931 shapes at occupancies five through
ten. Its maximum count inflation is 15335685049/15335681400, below
1.00000024. It improves twelve scalar density entries. The completed
support-200 replay, with both mean refinements and optimized coupled
columns, gives log2 upper +121.548507, rounded upward.

```sh
python -B research/workstreams/permutation_locality/independent_rows/candidates/overlap_mean/replay.py --epoch-cache tmp/independent-row-local-families --groups 80 --support 200 --tilt .072 --density-chord
```

The refined memo has a separate stage key that binds all source files
in the overlap-mean directory and whether translated density is enabled.
It does not overwrite or relabel the baseline family.
