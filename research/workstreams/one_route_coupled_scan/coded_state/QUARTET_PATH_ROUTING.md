# Retaining intermediate state weights in the routing bound

2026-10-27. This continuation keeps the field-quartet construction fixed.
The objective remains a whole-code45% distance certificate with40-bit
margin. The operator inequalities below are exact. Their large numerical
evaluations are FLOAT diagnostics, not directed finite certificates.

## The loss being removed

The previous full-weight calculation retains the two endpoints of each
16-update region. It then chooses a Bernoulli proof parameter p for that
entire endpoint class. Different internal state histories can prefer
different p values. Choosing p only after adding them can be substantially
looser than partitioning the histories first.

This is a choice in the proof, not in setup. The actual route still has
exactly r inactive-row positions per region, uniformly placed. The same
sampled code serves every message. No additional permutation or refresh is
introduced.

## The midpoint inequality

Fix an inactive-row count r, the common field moment theta, and every
global marker. Let K_p be the positive513-by-513 one-update kernel from
[the full-weight count](QUARTET_FULL_WEIGHT_FLOW.md). It includes every
incoming and outgoing state weight and every support-death count. Its
local inactive labels are independent Bernoulli-p variables used solely
to evaluate a generating polynomial. Put

    H_p=K_p^8,
    pi_r(p)=binom(8192,r) p^r (1-p)^(8192-r).

Let G_r(u,v) denote the same weighted regional count with exactly r
inactive labels, for incoming weight u and outgoing weight v. Partition
this count by the weight w after eight updates. For each fixed w, positivity
and conditioning the total labels give

    G_r(u,v;w) <= H_p(u,w) H_p(w,v) / pi_r(p).

The proof must use the same p in both factors of this inequality. A finite
set P of proof parameters gives the refined envelope

    G_r(u,v) <= sum_w min_(p in P)
                   H_p(u,w) H_p(w,v) / pi_r(p).                 (1)

The previous envelope places the minimum outside the sum over w. Hence
(1) is no larger, entry by entry, before numerical rounding. All513 values
of w remain. The binomial normalization already present in K_p averages
the uniform quartet partitions; symmetry lets the weight kernels compose.

The global row-set factor, field factor, initial zero state, and free
terminal vector are unchanged. Compose the regional envelope256 times.
No union over a chosen proof parameter is needed: every available p bounds
the same fixed midpoint class.

## Allocate inactive labels before allowing separate half-region bounds

The two halves cannot independently satisfy arbitrary label totals. To
retain this relation, write the half-region polynomial as

    H_p = sum_(s=0)^4096 B_s p^s (1-p)^(4096-s),

where each matrix B_s is entrywise nonnegative. Fix one reference rho in
(0,1). For an integer interval I=[a,b], define

    M_I = sum_(s in I) B_s rho^s (1-rho)^(4096-s).

For any p,

    M_I <= H_p R_I(p),
    R_I(p)=max_(s in I)
        (rho/p)^s ((1-rho)/(1-p))^(4096-s).                     (2)

The logarithm of the ratio is affine in s, so an endpoint attains the
maximum. Take entrywise minima over the proof grid to obtain E_I>=M_I.

Partition the feasible first-half counts into disjoint intervals I. A
first-half count in I requires its second-half count to lie in r-I.
Positivity gives

    G_r <= [sum_I E_I E_(r-I)] / pi_r(rho).                    (3)

Each matrix product in (3) retains its midpoint weight. Enlarging the
second-half count to the whole complementary interval introduces slack
but never omits a count pair summing to r. Only after this allocation may
the two halves use different p values. There is one full-region quota
normalizer, not two independently chosen regional distributions.

One can refine E_I further. Write J_p=K_p^4 and retain the weight x at the
quarter boundary. The valid replacement for (2) is

    E_I(u,v)=sum_x min_(p in P) J_p(u,x) J_p(x,v) R_I(p).       (4)

The implementation takes the minimum of (4) and the previous E_I. Using
(4) in (3) retains all three interior quarter-boundary weights without
storing a four-dimensional state array. It remains an upper envelope,
not an exact extraction of every label coefficient.

A separate diagnostic retains all three quarter-boundary weights for one
fixed endpoint pair (u,v), without count intervals:

    sum_(x,w,y) min_(p in P)
        J_p(u,x) J_p(x,w) J_p(w,y) J_p(y,v) / pi_r(p).          (5)

Each summand is a complete-region bound with one common p. Formula(5)
cannot exceed(1). Computing a few selected entries of(5) does not replace
the entire regional operator or certify its surrounding endpoint classes.

## Diagnostic results

Unless stated otherwise, r=2100, theta=7/40, ell=0.00377, and the grid has
nine p values. There are no weight bins, turnover restrictions, or terminal
restrictions. Entries below are log2 probability upper-bound screens.
Every positive entry is inconclusive for minimum distance.

| Regional calculation | FLOAT log2 screen |
| --- | ---: |
| Endpoints only, matched replay | +572899.41 |
| Retain midpoint before min-p | +404889.64 |
| Midpoint and half-quota intervals of width64 | +404889.64 |
| Midpoint and half-quota intervals of width16 | +404889.64 |
| Midpoint and quarter-refined quota intervals of width256 | +404889.64 |
| Midpoint and global coverage/death marker0.9 | +420216.43 |

The midpoint refinement improves the (272,272) regional entry by879.98bits.
The count-interval alternatives, at widths64 and16, are respectively148.87
and36.43bits worse at that entry; taking their minimum with (1) leaves the
reported global bound unchanged.
Retaining quarter weights with count width256 improves the (272,272)
entry by a further32.00bits, but again leaves the dominant global bound
unchanged. That run took about74s under the stated resource guard.
Formula(5) gives2931.0073bits for the0->489 entry and3217.8255bits for
489->0, improving their midpoint bounds by88.8407 and73.4858bits.
These two entry improvements are much smaller than the roughly1582bits
per region needed to remove the current global gap. This comparison is
diagnostic, not a bound on what a complete finer partition could achieve.
The selected entries were not substituted into the full operator.

An artificial weighted-path diagnostic helps locate the remaining slack.
Under the endpoint-only envelope, the midpoint between region128 and129
has mean weight272.05, with negligible mass below128 or above384. Under
(1), approximately49.61% of its mass lies below128 and50.39% above384;
17.85% is at zero. Its high-weight peak is near489. These distributions
describe the proof envelope, not states observed under the actual encoder.
They identify sparse/dense mixtures as a useful next diagnostic, not a
low-distance codeword or a reason to change the construction.

The regional-edge diagnostic confirms alternation, not just a bimodal
stationary marginal. Its low-to-high and high-to-low masses are0.496064
and0.503936; low-to-low is about6.3e-41 and high-to-high about6.2e-168.
Here low means weight<128 and high means weight>384. These are FLOAT
weighted-envelope statistics at the center of the256-region product.

For the representative edge0->489, the largest midpoint class is348,
carrying38.31% of that edge's bound. For489->0 it is332, carrying35.06%.
Both use proof p=0.70594552. Conditional mean weights in the corresponding
artificial positive-kernel bridges are:

| Update within region | 0->489, midpoint348 | 489->0, midpoint332 |
| ---: | ---: | ---: |
| 0 | 0 | 489 |
| 4 | 0.256 | 486.082 |
| 6 | 29.853 | 466.978 |
| 7 | 151.179 | 426.395 |
| 8 | 348 | 332 |
| 9 | 430.897 | 102.802 |
| 10 | 468.303 | 9.183 |
| 12 | 486.228 | 0.211 |
| 16 | 489 | 0 |

The low-output envelope therefore permits gradual growth and decay across
several updates. It does not merely pay for one abrupt full-state reset.
These conditional means are not a single integral support trajectory.
Any subsequent certificate must bound whole neighborhoods and their
complement, not substitute these displayed mean weights as actual states.

A separate test restores the earlier zero-selection-free interpolation
to the all-weight kernel. With tau=0.1 and the old endpoint envelope, the
screen is +776562.53, worse than tau=0. The free-terminal score potential
telescopes exactly: its total multiplier is
(q*ell)^[theta*tau*(A+D-W)]. This option is retained as a checked bound,
not presented as an improvement.

## Checks and resource limits

`quartet_remaining_activity.py` implements these envelopes. Rational toy
polynomials check27 midpoint components,54 half-quota components, and54
quarter-refined components against exact coefficient extraction. The
selected-edge formula(5) has27 further exact toy comparisons. The
midpoint tests also verify that the refined envelope never exceeds the
old one. Sixteen short paths check the free-terminal score identity.
Earlier exact activity, coverage/death, and routing-quota checks still pass.
Four bridge marginals match exhaustive rational path enumeration.

All numerical jobs run serially, file-free, on Peach. The helper now has
a120s watchdog and128MiB numeric-array budget for retained quarter matrices.
It still requires20GiB free at the execution site and writes no caches or
result dumps. None of these proof-computation timings is an encoder timing.

Whole-code closure still requires the remaining history classes, every
inactive-row count, and directed numerical verification of their final
union. Restricted certificates and the two global exclusions remain intact.

Next analyze the gradually growing/decaying pair of regions with the
remaining coverage/death condition and actual outer-row constraints. The
selected-edge check argues against spending the next round solely on finer
quarter-boundary enumeration. A positive bound on this family is still not
an obstruction in the encoder; a stronger joint count or a genuine rank/
codeword obstruction is needed to distinguish those possibilities.

## Continuation: quarter classes and the new dominant profiles

At theta0.13 and ell0.03, the midpoint envelope's representative alternating
edges are368->12 and12->368. Their largest midpoint classes are111 and120,
using proof parameters p=0.5040675381 and0.4147978567. Marked bridge moments
give the following artificial totals for one complete region:

| Class | A | D | C_T | Inactive labels | Active coverage |
| --- | ---: | ---: | ---: | ---: | ---: |
| 368,111,12 | 650.99 | 565.67 | 70.55 | 3211.78 | 1889.68 |
| 12,120,368 | 3128.34 | 420.18 | 145.18 | 1609.50 | 4516.60 |

Both regions require exactly2100 inactive labels. Their derivatives place
better proof-p values between the two displayed grid points. Rescaling the
two artificial profiles over128 pairs gives A about483755 and A+D about
609943, below the existing activity floor and below W. These are diagnostics
of separate proof envelopes, not a distribution of actual code histories.

The helper `quartet_quarter_classes.py` retains all three quarter-boundary
classes before minimizing over p. Let Q_p=K_p^4. Partition weights into bins
B_i, keeping weight0 separate, and fix any positive potential h. Define

    P_p(i,j)=max_(u in B_i) sum_(v in B_j) Q_p(u,v)h(v)/h(u).

For each fixed sequence of four bin transitions, the product of its P_p
entries bounds the corresponding path sum. Conditioning the full region
therefore permits the bound

    R(i,j)=sum_(a,b,c) min_p
      P_p(i,a)P_p(a,b)P_p(b,c)P_p(c,j)/pi_r(p).

The same p applies to all four factors before the minimum. Sum the bin
paths only afterward. Iterating R for256 regions uses terminal value
max_(u in B_i) 1/[f(u)h(u)] and initial factor h(0). The chosen h comes from
64 positive continuation iterations, but validity does not require an
eigenvector, convergence, or any approximation to actual code states.
Incoming maxima and terminal maxima ensure upper bounds for every weight.

Twelve rational toy cells agree with direct enumeration and dominate exact
quota coefficients; twelve further cells check monotonicity when adding p
choices. All full computations below retain every history, with a free
terminal state, but use FLOAT rather than directed arithmetic.

| Refinement at r=2100, theta0.13, ell0.03 | log2 screen |
| --- | ---: |
| Original nine-p midpoint | +303965.75 |
| Half-quota width256 with internal quarter refinement | +303965.75 |
| Three quarter classes, weight bins16 | +300683.93 |
| Three quarter classes, weight bins8 | +297812.70 |
| Midpoint, eight added p values between0.4148 and0.5041 | +285778.93 |
| Quarter bins16, same refined p set | +289630.74 |
| Midpoint, zero-selection-free interpolation tau0.05, original p grid | +277336.06 |
| Midpoint, tau0.05, eight added p values between0.5041 and0.5931 | +271673.74 |
| Quarter bins16, tau0.05, same refined p set | +275532.67 |

The finer p sets preserve every original choice. Their minima remain
valid and cannot worsen the ungrouped midpoint envelope. Coarsening weights
can offset the benefit of retaining more boundaries, as these results show.
The interpolation uses the previously proved factor
(q*ell)^[theta*tau*(A+D-W)], including its free-terminal correction.

These tests improve the bound but do not approach closure. The dominant
profile after interpolation moves to4->429->4 at region boundaries; fixing
one old profile would miss that change. The more substantive new result is
the [zero-neighbor geometry and activity certificate](QUARTET_SUPPORT_CUT.md#a-proved-geometry-tail-and-stronger-activity-exclusion),
which raises the rigorous activity cutoff at r=2100 to747590. Next extend
that structural exclusion and combine it with the remaining output count,
rather than spend another round only refining quarter bins. All runs were
serial on Peach, with120s/128MiB guards and no retained numeric files.
