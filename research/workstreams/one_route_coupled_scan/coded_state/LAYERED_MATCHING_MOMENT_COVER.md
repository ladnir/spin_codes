# Reweighted suffix regions and actual matching moments

This continues the [crossover refinement](LAYERED_CROSSOVER_COVER.md).
The RS outer, regional permutation distribution, layered inner, and setup
distribution remain fixed. We retain the exploratory threshold W=838860
at N=2^21. The counted histories are over GF(2^32); the intended scalar-
extension distance guarantee still requires complete base-field coverage.

The new certificates widen the suffix-193 region and cover neighboring
regions at suffixes 194 and 195. They are not whole-code certificates or
near-half distance claims. The GF32 sandwich remains the certified control.

## Certified regions

Use x=(z,r,f,e,B) for empty-pair, 01-pair, full-pair, exceptional-full-pair,
and internal-birth counts. In the table, z ranges over the whole interval.
Each other count varies independently around (24568,589632,368520,49136).

| Suffix R | Empty-pair counts z | Radius in each other count | Verified boxes | log2 good-route contribution upper |
| ---: | --- | ---: | ---: | ---: |
| 193 | 0..51000 | 512 | 12 | -398 |
| 194 | 0..50000 | 256 | 5 | -882 |
| 195 | 0..49264 | 128 | 6 | -1084 |

Every box includes all active-row counts s=1..8192 and every active-row set.
The transfer sums support locations, exceptional branches, components, and
amplitudes. The count-vector and box unions are included. Disjointness and
total volume are checked; failed leaves prevent a complete-coverage result.

The bounds control expected low-weight nonzero-word counts intersected with
the good-route event. They are not whole-code failure margins. Union over
these three regions is at most 2^-396 on that event. The route-failure
budget must be added once when assembling a complete certificate.

## What reweighting accomplished

We first kept the previous binomial matching envelope. At the difficult
suffix-193 neighbor

    (z,r,f,e,B)=(49500,24440,589504,368648,49008), s=8060,

a continuous interpolation proposal near m=24.4624 gives verified upper
-3095. The entire region z=0..49500 with radius 128 then closes in nine
boxes, with union upper -368, without sharpening the route bound further.

Suffix 194 remained harder. Reweighting its central profile at s=8115
reduced the slice bound to +1444, at m about 24.4419. The full-profile
bound, after selecting the enlarged marker bank, was +4726. More precise
selection of the previous binomial Chernoff tilt could only modestly lower
the active edge threshold: at u=2415, a denominator-64 tilt screen improved
91789 to 91741. This motivated retaining the actual matching moments.

Both retained proposals are rational proof inputs. They are not construction
parameters. Every acceptance uses recomputed rational local row maxima and
directed exponent bounds, not optimizer convergence.

## The actual uniform-matching moment

Fix an s-vertex subset of an M-vertex set, with M even. Let Z count its
internal edges in a uniform perfect matching. Counting specified collections
of l disjoint internal edges gives the exact binomial moment

    a_l = E[binom(Z,l)]
        = s! / ((s-2l)! * 2^l * l! * product_(j=0)^(l-1)(M-2j-1)).

The earlier proof upper-bounded a_l by a binomial coefficient times p^l.
We now evaluate the positive series before that relaxation. For t>=1,

    g_s(t) = E[t^Z] = sum_(l=0)^floor(s/2) a_l*(t-1)^l.

Writing c_l=a_l*(t-1)^l, its recurrence is

    c_0=1,
    c_(l+1)/c_l = (s-2l)(s-2l-1)(t-1)
                 / (2(l+1)(M-2l-1)).

This retains the dependence among edges in a matching. It uses the same
uniform matching induced by the existing regional route, not a different
permutation family or an independence assumption about matching edges.

Different regions are independent. Thus, for E_R the internal-edge count
over R regions and rational t>1,

    Pr[E_R>=h] <= g_s(t)^R / t^h.

Let c_s be the bit length of binom(M,s). A sufficient threshold is

    h >= (R*log2(g_s(t)) + c_s + 77) / log2(t).

The implementation uses the single fixed tilt t=3/2. It evaluates all terms
of the finite series at scale 2^96, rounding each positive multiplication
and division outward. Induction on l bounds every term and their sum.
The upper sum becomes a rational input to the directed log2 routine.
An upper numerator and positive lower denominator give an integer ceiling.
No floating-point calculation establishes a tail bound.

For each (R,s), take the minimum with the previous envelope. Both thresholds
are deterministic, so the smaller one itself has the required tail bound;
there is no additional union over the two calculations. Use complement
symmetry above s=M/2. The same union over subsets, cardinalities, and suffixes
gives matching-event failure at most 2^-56. The head event is unchanged.

The new good-route event implies the earlier one. It replaces its matching
component with the sharper envelope at the same failure budget. Existing
class certificates remain valid on the stronger event; no separate failure
budget is charged for every box or every proof-weight choice.

## Effect on the dimension cap

The route dimension maximizes over subset sizes, so its improvement is not
the edge-cap improvement at one selected size. The complete maximization
at the central history counts gives:

| R | Previous binomial-moment cap | Actual-moment cap |
| ---: | ---: | ---: |
| 193 | 155878 | 154310 |
| 194 | 159381 | 157740 |
| 195 | 162877 | 161180 |

With the retained marker bank, the central complete-profile upper bounds
become -51731, -35636, and -17034, respectively. The box verifier then
establishes the larger count regions above. No code operation changes.

## Explicit remaining boundary

The attempted suffix-193 region z=0..52000 with radius 512 does not close
under a 256-test subdivision budget. At its difficult corner and s=8024,
continuous reweighting gives verified upper +2698. This is not a low-weight
witness and does not establish optimality of the bound.

At the central counts, retained-weight bounds for R=196,197,198 are
+6294,+32620,+51482. Other count regions, suffix lengths, and higher output
thresholds remain uncovered. None of the region certificates says that
every possible history at its suffix length is covered.

Next use the actual-moment envelope in the broader count-domain cover,
then target the surviving boundary cells rather than repeatedly sampling
the same neighborhood. At the first failing cells, compare dedicated
reweighting with additional relations among nonzero fixed inputs. Preserve
the construction and the certified control.

## Reproduction and checks

Run serially from the repository root:

    python -B research/workstreams/one_route_coupled_scan/coded_state/test_layered_matching_moments.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/test_layered_fractional_weights.py
    python -B research/workstreams/one_route_coupled_scan/coded_state/cover_layered_location_regions.py --suffix-progress

The matching test independently enumerates small perfect matchings and
checks 140 exact generating-function identities, 315 rational tail bounds,
and five 96/128-bit precision comparisons. The K20 envelope never exceeds
the previous envelope and satisfies the complement identity. The fractional
weight tests exercise the original, binomial-refined, and actual-moment
route envelopes against the direct expression.

The final serial replay reproduced all three region bounds and passed
14 integer-weight equivalences, 1368 local fractional-weight checks, and
480 direct-expression majorants. Regression checks also passed: 4044 joint
dimension inequalities, 640 direct location bounds, 128 rank fixtures,
972 colour bounds, 1188 positive local rows, 2207 exhaustive histories,
and 24768 shortening checks. Rational-log and prior binomial-matching
tests passed. The two frozen theory helper hashes remain unchanged.

Bounded searches use refine_layered_location_weights.py with --refined-route
or --matching-moments. The latter selects the actual-moment envelope. Search
outputs are proposals; only compact accepted inputs and proof notes are
retained. No benchmark, raw result file, production edit, paper change,
commit, push, or promotion is part of this work.
