# Couple arbitrary consecutive updates

The full-input pair checkpoint improved the regional bound but did not
close it. This calculation removes the requirement that one input count
equal512. It changes the proof transfer, not the encoder or its setup.

## Sum the intermediate histogram before relaxing it

Use the five ages, potential h, selection weights d, and quartet counts
J_(l,k) defined in [the full-input pair note](QUARTET_FULL_PAIR.md). Fix an
incoming histogram n of total width w. Let l,m be the two allowed-input
counts. Let R,S be positive markers for the final first and second ages.
All probabilities still average the actual once-sampled setup.

Write b3=b^alpha and b4=b^lambda. The following polynomial records the
first selection and the possible second selections from those positions:

    P_n(x,y) = (1+x*(1+d1*S*y))^n0
               (1+x*(1+d2*y))^n1
               (1+x*(1+d3*y))^n2
               (1+x*b3*(1+d4*y))^n3
               (1+x*b4*(1+d4*y))^n4.

The degree in x counts the first selected positions. The degree in y
counts second selected positions. A selected old type2 becomes type3;
its second selection therefore uses d3, not an independently maximized
choice between d3 and d4. The first reward contributes b3 or b4 for the
old types3 and4.

Conditional on first selection size k, there are also w-l inactive and
l-k fresh type1 positions. Define the univariate polynomial

    F_(n,l)(y) = sum_k J_(l,k)/binom(w,k)
                  * (1+d0*R*y)^(w-l) (1+d1*S*y)^(l-k)
                  * [x^k]P_n(x,y)/phi(n).

The exact two-update moment, including final potential divided by initial
potential, is

    sum_j J_(m,j)/binom(w,j) * (h1*R)^(m-j) * [y^j]F_(n,l)(y).

The final first-age count is m-j plus the second selection from the
inactive group. The final second-age count is the second selection from
the type1 group. These are exactly the powers recorded by R and S.
The intermediate potential cancels; no independent histogram is chosen
between the two updates.

## Uniform source-cell bound

For fixed y,R,S, apply the existing coupled coefficient bound in x with
selection weights

    (1+d1*S*y, 1+d2*y, 1+d3*y,
     b3*(1+d4*y), b4*(1+d4*y)).

This bounds each first-selection coefficient uniformly over the incoming
cell. Sum k with its actual quartet multiplicity. Then bound each second
coefficient by y^(-j)F_(n,l)(y), minimizing over positive evaluated y only
after the source maximum. Different j may use different evaluated y.

The second coefficient also admits the earlier Bernoulli point-mass bound.
For each fixed first selection, F is a positive multiple of a product of
w linear factors. At a fixed y, its tilted variance is at least

    (w-l)*v(d0*R*y) + l*min(v(d1*S*y),v(d2*y),v(d3*y),v(d4*y)),
    v(z)=z/(1+z)^2.

This lower bound is uniform over the first selection. Therefore its common
point-mass factor may multiply the positive sum over those selections.
It is not the variance of an arbitrary mixture. Numerical implementation
uses the existing point-mass envelope and remains FLOAT, not directed.

To bound a final bucket, multiply a marked moment by the maximum of
R^(-r)S^(-s) over that bucket. Such bounds may be intersected with existing
edge bounds. A total-row bound must not instead be used to scale every
edge down: those operations have different meanings.

## Checks and initial local results

`quartet_arbitrary_pair.py` checks5979 exact marked moments against direct
two-update enumeration at widths4 and8. It also checks879 FLOAT uniform
moment bounds against exact moments. The tests include zero/full inputs,
mixed ages, both marker directions, and every input-count pair at the toy
widths. They do not establish full-width directed rounding.

The full-width screen holds the ordinary proof potential and width16 age
buckets fixed. Selected total-row results are:

| Source a; n1; n2 | Counts | Separate log2 bound | Coupled moment log2 bound |
| --- | --- | ---: | ---: |
| 312;64..79;16..31 | 312,408 | 844.21635 | 843.00925 |
| 312;64..79;16..31 | 408,312 | 859.38057 | 831.46316 |
| 312;64..79;16..31 | 408,408 | 854.38134 | 852.85302 |
| 408;192..207;96..111 | 408,312 | 862.49746 | 846.86281 |
| 512;0..15;0..15 | 408,312 | 804.28859 | 747.40903 |
| 512;0..15;0..15 | 512,312 | 877.50150 | 767.24703 |

All12 tested rows improve. The effect is small on steady408,408 pairs and
larger on transitions into count312. This suggests that retaining the
intermediate age matters particularly when activity falls. These are local
coefficient reductions, not security margins and not a whole-code result.

The serial one-thread remote run completes in13.01s, under a180s watchdog
and256MiB array limit, without generated files. The previous full-input
pair formula remains sharper on the tested full-input pair and is preserved.

## Regional composition and subsequent obstruction

The marked edge envelope passes6732 exact bucket checks. Evaluating it
over coarser width64 incoming cells permits reuse across the existing
width16 boundary cells; the latter are not merged. On a local test, the
coarse cap helps some falling-activity rows but loses other fine-cell gains.
An underflow assertion caught tiny positive composed terms; the diagnostic
now raises them to a conservative positive floor instead of dropping them.

`quartet_arbitrary_regional.py` applies nine outgoing marker pairs to all
201 coarse source cells, intersecting1,405,231 existing pair edges. The
regional calculation retains2502 boundary states and all paths over the
same eight selected counts. All regional payment allocations at the fixed
r=2100,c=111572 are retained. It completes in729.66s, with a1200s watchdog,
single thread, less than8GiB conservative array storage, and no files.

After three boundary updates the FLOAT upper expression is+255263.80682;
the relaxed balanced lower diagnostic is+212380.35033. This fitted upper
candidate is slightly larger than the previous+255164.00196 candidate.
The new matrix is entrywise no larger, so the old candidate remains a
valid upper bound for it; the fitting procedure is not an exact optimizer.
The explicit counter-floor contribution is over250bits below every output.

More importantly, the subsequent [exact moment obstruction](QUARTET_EXACT_MOMENT_OBSTRUCTION.md)
rules out closing this all-paid-set moment even with exact transfer
calculation and optimized scalar parameters. Preserve these local bounds,
but do not continue their refinement as the route to the full certificate.
The [closed-paid-set interface](QUARTET_CLOSED_PAID_SETS.md) is the next
direction. The encoder, full45%/40-bit goal, and failure allocation remain
unchanged; full coverage and directed verification are still required.
