# Extra RS zeros and eligible routed rows

For the unchanged three-shear model, the clean-crossing argument now covers
130 touched regions with a restricted failure bound below 2^-70. This
extends the prior 129-region argument by allowing one extra RS zero per
active row. It does not close the whole code or all six-defect histories.

The family has five clean closed proper episodes, followed by a live tail.
Each closed input interval crosses one regional boundary. Regions touched
by different intervals are disjoint. The tail remains live through the
last internal state; its internal births and one-sided pairs are unrestricted.
As before, no output-weight restriction is needed to exclude this family.

## Why counting only reset cost was insufficient

Let R=128+k be the number of touched regions. An active RS[256,128] row
may now have k-1 zeros among these regions. Charging every region to every
active row is invalid. Conversely, assigning k message degrees to every
row loses useful information: most rows cannot route enough coordinates
into the short intervals to attain that dimension.

Fix the interval endpoints before sampling the route. For outer row i,
let f_i count its coordinate positions inside the allowed input intervals,
including positions where a particular codeword may be zero. Its shortened
RS message space has dimension max(0,f_i-128). Call a row eligible when
f_i>=129. Every nonzero row of a compatible message must be eligible.

Keep a separate upper bound C_j on the number of rows with f_i>=128+j,
for j=1,...,k. These bounds control the message dimension and the available
RS zeros jointly, rather than treating them as unrelated choices.

## A route event uniform over all messages

Write l_i for the length of closed episode i, between one and eight. Its
input interval contains l_i+1 batches. If it splits into u and v batches
across its two private regions, the containment fractions are u/16 and v/16.
Give every tail region full occupancy, which only increases eligibility.

For threshold j, at least r=10-k+j of the ten early-region positions must
lie in their allowed fragments. Let p be their ten containment fractions.
For any d between one and r, the positive factorial-moment majorant is

    beta_j = min(1, min_d e_d(p)/binom(r,d)).

For s specified rows, multiply these positive majorants and expand. Within
each region, uniform sampling without replacement bounds each term by the
corresponding product of containment fractions. Independence between
regions then gives joint probability at most beta_j^s. Thus

    Pr[number of threshold-j rows >= s] <= binom(M,s)*beta_j^s.

This is not an assumption that row routes are independent. Equalizing the
two fractions within one episode increases every elementary symmetric
polynomial: their sum stays fixed and their product increases. Therefore
the verifier may replace both fractions by (l_i+1)/32, simultaneously for
every possible split and every moment order.

There are at most (4097)^12 < 2^156 endpoint choices. The verifier sets
each threshold's failure budget to 2^(-256-ceil(log2 k)). After the
threshold and endpoint unions, the bad-route probability is at most 2^-100.
Endpoints determine the ordered episode lengths. The 792 sorted length
multisets only compute a worst-case bound; they do not require an additional
union factor on top of the endpoint union.

## Charge message dimension and reset constraints together

Fix a route satisfying every threshold cap. Consider s active rows.
The number with shortened dimension at least j is at most min(s,C_j).
The total message dimension is therefore at most sum_j min(s,C_j).

For a clean episode, every input lies in its initial support U_i, and
|U_i|>=2^l_i. If its two regions contain n_u and n_v nonzero inputs, then

    |U_i| >= (n_u+n_v)/(l_i+1).

Put a_i=32-ceil(log2 l_i). The reset bound supplies at least a_i|U_i|
bits of exponent. Assign both regions the integer price

    p_i = floor(a_i/(l_i+1)).

An actual nonzero word of weight 128+j has k-j zero coordinates
among the R touched regions. To lower-bound its reset charge, remove the
k-j largest prices from the list of ten region prices. Let P_j be the
sum of the remaining prices. Separately, backward matching closure gives
the absolute reset charge A=sum_i a_i*2^l_i.

Available shortened dimension does not determine actual word weight. This
distinction was missing in the first version of this note. The following
weighted count repairs that gap; a fresh replay preserves the stated bounds.

A row with f=128+d available positions has at most binom(f,z)*q^(d-z)
nonzero words with exactly z additional zeros there. Each additional zero
saves at most 16 reset bits, but costs one message degree. Its weighted
nonzero enumerator is therefore at most

    q^d * 2^-P_d * sum_z binom(f,z)*2^(-16z)
      <= (256/255) * q^d * 2^-P_d.

The last inequality follows from binom(f,z)<=f^z and f<=256. The exact
unweighted nonzero count is at most q^d. Also (256/255)^128<2. Consequently,
for fixed available dimensions d_1,...,d_s, the message/reset exponent is
bounded by

    32*sum_i d_i - max(A, sum_i P_(d_i)-ceil(s/128)).

Increasing a row dimension by one adds 32 bits of message entropy, but
increases its reset charge by at most 16 bits. The expression is therefore
nondecreasing in every dimension. It is bounded by saturating the nested
caps: set n_j=min(s,C_j), and use n_j-n_(j+1) rows of dimension j.
There are at most binom(C_1,s) choices of active rows on the fixed route.
Also s>=2^(max_i l_i), because an activation batch contains at most one
coordinate per active row.

For each s, a second bound keeps route rarity directly in the first moment:

    binom(M,s) * beta_1^s * q^(ks) * 2^-max(A,P_1*s).

Here the row choices range over M, not C_1; the route cap is not counted
twice. Both expressions bound the same contribution restricted to good
routes, so their minimum is valid for each s. Sum over s and include one
factor two for the accepted-map density relative to the independent uniform
pair-image comparison. The prior reset argument justifies this single
global factor; it is not one factor per episode.

## Exact bounded checks

`diagnose_pair_eligible_rows.py` uses rational route moments and directed
integer logarithmic bounds. It examines all 792 multisets of five episode
lengths. Results include all active-row counts and interval placements:

| Touched regions | Good-route contribution | Route exception | Combined upper |
| --- | ---: | ---: | ---: |
| 129 | 2^-150 | 2^-100 | 2^-99 |
| 130 | 2^-71 | 2^-100 | 2^-70 |
| 131 | bound exceeds one | 2^-100 | unresolved |

The older 129-region proof is stronger and remains unchanged. For 130
regions, 746 of 792 length multisets are excluded because too few eligible
rows exist to activate the required clean support. The remaining multisets
close by the combined dimension/reset and route/reset bounds. The worst is
five length-one episodes. A diagnostic counter named `cells_above_minus80`
reports one cell because its internal threshold is stricter than the final
40-bit target; it does not contradict the 70-bit restricted bound.

In this eligibility-only bound, at 131 regions the worst retained case
has five length-three episodes.
The route caps are (28,18,13), and the per-row reset prices are (56,63,70).
The present bound is positive after the placement union. This is a failure
of the bound, not a witness of a low-weight codeword. The subsequent
[joint activation-route count](UNIFORM_PAIR_ACTIVATION_ROUTE.md) closes
this case and the entire stated 131-region crossing family.

The test mode checks 810 factorial majorants, 512 RS-zero price bounds,
100 fragment-balancing inequalities, and 738 nested-dimension assignments.
The complete three-column diagnostic took about 26 seconds; tests took
under one second. The source enforces a 35-second sweep cap. No numerical
files, Python caches, benchmarks, or concurrent jobs were created.

## Decision and next test

Continue the proof rather than redesign now. The additional RS zero was
handled by a constraint already supplied by the unchanged outer and route.
The next discriminating case is five length-three episodes with 131 touched
regions. Current eligibility counts ignore that each activation support
must be a union of matched pairs, and that later supports must remain
compatible with the subsequent matchings. Test this actual route constraint
before extending the sweep or proposing a new inner operation.

Shared regions and within-episode defects remain outside this result.
Neither this restricted base-field calculation nor its predecessor is a
whole-code certificate over the application's larger alphabet. The frozen
certified control and all previous proof helpers are preserved.
