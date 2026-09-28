# Sum count bounds inside support intervals

The previous support-box method multiplied the largest count in an interval
by the worst binomial denominator in that interval. Those extrema can occur
at different supports. Their combination caused a large loss even when the
individual support points had ample margin.

This refinement sums within each interval. The code, routing distribution,
K=2^20, IMT(128,19), and bad output threshold 209715 are unchanged.

## A weighted interval bound from upper CDFs

Let a(u) be the unknown number of four-row BCH tuples with support exactly u.
The authenticated increasing function U(u) upper-bounds sum_{v<=u} a(v).
Fix an interval [lo,hi] and positive weights w(u). Set

    W(u) = max_{u<=v<=hi} w(v).

This is a decreasing majorant. Summation by parts gives

    sum_{u=lo}^{hi} a(u) w(u)
      <= U(lo) W(lo) + sum_{u=lo+1}^{hi} (U(u)-U(u-1)) W(u).

The restricted prefix sum on [lo,u] is at most U(u); this is why the first
term uses U(lo), without subtracting U(lo-1). The differences of upper CDFs
are coefficients of this extremal upper calculation, **not claims about
the true shell counts**. All coefficients and terms are nonnegative.

For a support witness p, use w(u)=1/Pr[Bin(256,p)=u]. The resulting folded
count replaces U(hi) divided by the smallest binomial mass in the interval.
The folded count never exceeds that older bound. Its implementation passes
2,430 exact-rational tests with increasing, decreasing, and nonmonotone
weight sequences.

## Combine with the existing moment and cover

Fix q active groups. A support box specifies an interval for each group.
For fixed witness probabilities p_1,...,p_q, the positive region matrix is
the same Poisson-binomial mixture of the occupancy matrices R_0,...,R_q.
The box contribution is bounded by its group-location and label counts,
times the 256-region tilted moment, times the product of the folded counts.
The density coordinate C is excluded from the terminal mass functional.
With `--mature-tail`, its two additional subset bounds are also excluded.

`occupancy_cdf_cover.py` refines the largest estimated box. Its selected boxes
partition the full labeled support domain {38,...,256}^q. The existing
exact split-volume and label-multiplicity checks remain in force. Floating
search proposes rational probabilities; after it passes, every retained
box is replayed outward. The script then checks the total against its
requested failure budget. A search that stops above budget produces no
certificate.

Optional all-one-column penalties use a weighted outer CDF. For each support
prefix, sum the exact-support flag caps by all-one count J, retain the
original total-count cap, and maximize sum_J A(J) rho^-J by integer greedy
allocation. Rounding this rational upper upward gives an increasing integer
CDF. The same folded-count proof then applies to the weighted counts.
Small-code enumeration checks these weighted CDFs at rho=1/2,3/4,1.

## Verified full-support results

Each row includes all ranks, support vectors, and choices of q active groups.
Rows seven through ten use the earlier pair-conditioned nine-coordinate
operators, with no all-one penalty or new window averages.

| Active groups q | Leaves | Margin bits | Outward precision |
|---|---:|---:|---|
| 7 | 85 | 168.89303 | 192 and 384 |
| 8 | 99 | 146.79273 | 192 and 384 |
| 9 | 109 | 110.62991 | 192 |
| 10 | 109 | 82.38367 | 192 and 384 |
| 11 | 309 | 73.76703 | 192 and 384 |
| 12 | 493 | 57.19597 | 192 and 384 |
| 13 | 366 | 75.22360 | 192 and 384 |

The corresponding conservative uppers are 1.439240e-51, 6.471185e-45,
4.978286e-34, 1.585054e-25, 6.221719e-23, 6.057548e-18, and 2.266948e-23. Together with the earlier
one-through-six results, their combined contribution is still below
3.188222e-15, about 48.1562 bits. This is not a full-code certificate:
messages with fourteen or more active groups remain outside this combined
bound.

The eleven-group row adds the all-one-column penalties and fresh/window
averages from WINDOW_AVERAGES.md. Its witnesses use rho in {1/2,3/4,1} and
lambda in {.0032,.004,.005,.0064}. Without those refinements, the tested
eleven-group cover reached only about 21 bits after 1,000 splits. A basic
twelve-group screen remained at log2 +54.50 after 500 splits. Neither basic
run met the requested budget.

The first combined twelve-group attempt remained open. With all-one penalties,
fresh averages, window averages, and six tilts through .01, its 500-split
cover had 1,460 leaves and log2 upper -17.93. It did not enter outward replay.
The largest remaining boxes mix low, intermediate, and high supports, rather
than consisting of a single equal-support point.

Splitting is not necessarily monotone for this relaxation: each child can
reuse a global CDF cap, and the children can select different witnesses.
In that twelve-group run the full cover worsened slightly between 250 and
500 splits. The new optional `--retain-parents` mode retains the better of
a parent bound and the sum of its refined children, while continuing to
explore the children. It selects either the entire parent or all selected
descendants, never both. The final exact-volume check and outward replay
remain unchanged. Synthetic tests check both choices and monotonicity.
The earlier successful replays do not rely on this option.

The optional `--joint-witness` mode optimizes the Bernoulli witness for the
actual mixed-support box. The older proposal optimizes each interval as if
every group used that interval. Joint optimization can instead account for
the other groups' supports. Equal intervals share a probability; the best
existing proposal remains a fallback. Only rational probabilities after
quantization are used in replay. Numerical optimality is not a premise of
the certificate. Combining joint witnesses and parent retention with the
multi-window bound closes twelve groups after 75 splits. Its 192-bit outward
upper is 6.057547359530498638e-18, including every support and rank. The
384-bit rerun gives 6.057547359531440552e-18; both are below the reported
conservative upper. The vectorized floating objective changed the selected
rational witnesses slightly, but each full cover passed its own outward
replay. This run uses only the authenticated BCH spectrum and shortening
caps, not the
optional orthogonal-array count refinement in the selected-point screens.

Reproduction:

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 7
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 7 --precision 384
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 8
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 8 --precision 384
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 9
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 10
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 10 --precision 384
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 11 --max-splits 500 --fresh --window-average --penalties .5 .75 1 --tilts .0032 .004 .005 .0064
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 11 --precision 384 --max-splits 150 --fresh --window-average --penalties .5 .75 1 --tilts .0032 .004 .005 .0064
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 12 --max-splits 500 --fresh --window-average --penalties .5 .75 1 --tilts .0032 .004 .005 .0064 .008 .01
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 12 --max-splits 100 --fresh --window-average --multi-average --retain-parents --joint-witness --penalties .5 .75 1 --tilts .005 .0064 .008
    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 12 --precision 384 --max-splits 100 --fresh --window-average --multi-average --retain-parents --joint-witness --penalties .5 .75 1 --tilts .005 .0064 .008

The initial thirteen-group follow-up used the same refinements, added tilt .01,
and stopped after 150 splits. Its selected 966-leaf cover had log2 upper
+36.92962, so it produced no outward certificate. The largest boxes mostly
have supports between 121 and 175; some extend from 93 to 202. This does not
show a bad codeword. Before growing the split budget, test exact support
points within these boxes to separate interval slack from operator slack.

The subsequent [exact-point investigation](THIRTEEN_GROUPS.md) confirms a
gap even at fixed supports and adds two further refinements. These are
optional and do not change the reproduced certificates above.

[MATURE_TAIL.md](MATURE_TAIL.md) subsequently closes thirteen groups by
retaining two low-expansion-weight subset bounds within mature mass. It
records the new invariant, both outward replays, and reproduction commands.

    python -B research/workstreams/permutation_locality/occupancy_cdf_cover.py --groups 13 --max-splits 150 --fresh --window-average --multi-average --retain-parents --joint-witness --penalties .5 .75 1 --tilts .005 .0064 .008 .01

The count refinement reduces search cost as well as proof slack. Multi-input
window averaging alone reaches log2 upper -22.73 after 500 splits at twelve
groups. The joint-witness cover then closes that occupancy. Full thirteen-
group coverage is now closed by the mature-tail extension. Fourteen groups
remain the next full-coverage target. The measured encoder remains at about
1.9x; no new timing or production change is implied by these proof results.
