# From one empty-pair profile to count regions

Follow-up: [crossover refinement](LAYERED_CROSSOVER_COVER.md) closes all six
remaining selected neighbors, certifies the intervening region, and extends
the coverage into suffix 193. This note preserves the earlier bounds.

This continues the [location-aware outer count](LAYERED_OUTER_LOCATIONS.md).
The RS[256,128] outer, regional permutations, and layered pair-shear inner
are unchanged. The field is GF(2^32), with the existing scalar extension to
GF(2^128). We retain forced-input locations, surviving components, and
exceptional updates jointly.

The new result covers four regions of integer history counts at the
exploratory 40% threshold. Each region includes every active-outer-row count,
not just the active-row count used to propose proof weights. Full history
coverage and a near-half distance certificate remain open. The GF32 sandwich
remains the certified whole-code control.

## What is covered

Use N=2^21, M=8192, and output-weight threshold W=838860. For a suffix of R
regions, write x=(z,r,f,e,B). The first four coordinates count empty pairs,
01 pairs, full pairs, and exceptional full pairs. B counts internal births
of active coordinates. These are the counts from the all-pair transfer;
they do not specify locations or amplitudes.

The previous hard profile had R=192 and

    x=(49136,24568,589632,368520,49136).

In the table, z ranges over the entire indicated interval. Each of r,f,e,B
varies independently within the stated radius of (24568,589632,368520,49136).

| R | Empty-pair count z | Radius in each other count | Verified boxes | log2 contribution upper bound |
| ---: | --- | ---: | ---: | ---: |
| 190 | 0..49168 | 128 | 4 | -824 |
| 191 | 0..49168 | 128 | 6 | -3933 |
| 192 | 0..49200 | 64 | 9 | -707 |
| 192 | 0..48000 | 128 | 9 | -417 |

These bounds include the union over integer count vectors in each region,
all active-row sets and sizes, and every subdivision box. For each count
vector, the transfer sums support locations, exceptional branches, component
histories, and amplitudes. Infeasible count vectors contribute no histories.
No coverage claim is inferred from testing corners or neighboring points.

The contribution is the expected number of low-weight nonzero words
intersected with the existing good-route event. Union over the four regions
gives an upper bound of 2^-415; overlap only makes this union bound looser.
The bad-route probability must be added once when assembling a complete
certificate. These exponents are not whole-code security margins.

## Keeping the shared count dependence

A separate maximum for the route dimension loses useful cancellation against
the other count terms. Instead, maximize their sum over each box.

Fix R and a box Q=[lo,hi] of integer counts. Set P=(RM-b)/2 and
D(x)=P-z+r+f-2e. Let H_R(u) and E_R(u) be the existing simultaneous head
and internal-edge envelopes for a u-row subset. Their definitions and route
failure budgets are unchanged. The route dimension bound is

    dbar(x) = min(floor(D(x)/2), max(0,RM-N/2),
                  max(0, max_u [min(D(x), P-z-e+E_R(u),
                                        H_R(u)+r+f-e, Ru)-128u])).

For fixed positive proof weights, the remaining log bound has the form
a.x+c*dbar(x)+terms independent of x, with c>=0. Define

    L_Q(v) = sum_j v_j * (hi_j if v_j>=0, else lo_j).

Thus L_Q(v) is the exact maximum of the affine linear part v.x on Q.
For each u, let F_(u,j)(x) be the four affine expressions inside the
innermost minimum, each minus 128u. Then

    max_(x in Q) [a.x+c*min_j F_(u,j)(x)]
        <= min_j max_(x in Q) [a.x+c*F_(u,j)(x)].

Take the maximum of these bounds over u and the zero branch L_Q(a).
Also take the minimum with the separately valid D/2 and max(0,RM-N/2) caps.
The result upper-bounds a.x+c*dbar(x) throughout Q. This uses
max-min <= min-max, not an interchange of equality. The checker evaluates
every affine maximum using integer coefficients and endpoint choices.

For each marker set, the remaining dependence on active-row count s is
affine(s)+(1-R)*log2 binom(M,s). It is convex, so endpoint upper bounds
cover each complete s interval. The fixed schedule covers 5901..8192;
the earlier two-update count covers 1..5900. A 13-bit allowance safely
unions the s classes, including the already-unioned small-s term.

Each box receives a further ceil(log2(number of integer count vectors))
allowance. Binary subdivision splits only the z interval. Accepted and
failed intervals must together form a gapless partition. The verifier
reports failure if any cell remains uncertified; it never drops such cells.
The region bound then unions the accepted boxes.

## Tight directed logarithms

The old integer-power logarithms were sound, but their rounding loss could
reach several thousand bits after multiplication by history counts. The
new helper bounds log2 of a positive rational using directed integer
fixed-point arithmetic. It does not change a probability model or proof
inequality.

Normalize the rational to 2^e*x with 1<=x<2 and put t=(x-1)/(x+1).
For t<=1/3,

    ln(x) = 2*sum_(j>=0) t^(2j+1)/(2j+1).

The helper evaluates 64 terms at scale 2^96, rounding every multiplication
and division outward. Using t<=1/2, the omitted tail is at most

    4 / (3*(2*64+1)*4^64).

It bounds ln(2) the same way, then divides positive intervals outward.
The returned log2 bounds have scale 2^24; the exponent e is added exactly.
All certificate acceptance uses integer upper bounds. Decimal and floating
evaluations appear only in regression checks or proposals.

## Remaining boundary

The fixed marker schedule passes 237 of 243 selected neighbors whose five
offsets independently belong to {-128,0,128}. This is a diagnostic screen,
separate from the region certificates above. Its worst point is

    (R,z,r,f,e,B)=(192,49264,24440,589504,368648,49008),

with log2 upper bound +2306. Selecting among eight marker sets separately
for every s improves the full-profile bound to +1760, but does not close it.
The remaining positive contributions occur in the s7936/s8016 marker ranges
near their crossover. A new marker closes its assigned 8014..8061 interval;
that local success is not a certificate for the whole profile.

At the original central counts, the fixed-schedule bounds for R=193 and
R=194 are +12704 and +37637. Neither is a low-weight witness. Both are
explicitly retained failures of this bound. We have not covered the full
count domain, all suffix lengths, or a near-half output threshold.

Next sharpen weights near the active-row crossover, then extend boxes into
the six failed neighbors and R=193. If reweighting stalls, retain additional
relations among nonzero fixed inputs in the same outer/component count.
The earlier quotient-rank fixtures demonstrate such unused relations, but
do not yet provide a uniform K20 rank gain.

## Reproduction and checks

Run serially from the repository root, with Python bytecode disabled:

    python -B research/workstreams/one_route_coupled_scan/coded_state/test_rational_log_bounds.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/test_layered_location_regions.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/cover_layered_location_regions.py --tubes
    python -B research/workstreams/one_route_coupled_scan/coded_state/cover_layered_location_regions.py --neighbors
    python -B research/workstreams/one_route_coupled_scan/coded_state/cover_layered_location_regions.py --adaptive --active-step 1

The log tests include 513 exact powers and five rational-series enclosures.
The box tests check 4044 shared-dimension inequalities, including exhaustive
small boxes, and 640 direct evaluations of the original expression. They
also check binomial rounding, replay a gapless region, and require an
uncertified cell to survive an exhausted subdivision budget. These tests
validate implementations; the affine and convex arguments establish
coverage between tested points.

Only source, proof notes, and compact rational weights are retained.
No encoder, setup distribution, production parameter, or paper claim changes.
