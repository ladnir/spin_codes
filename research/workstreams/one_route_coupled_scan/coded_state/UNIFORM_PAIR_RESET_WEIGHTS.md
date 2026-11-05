# Joint reset weights and clean crossing episodes

The fixed three-shear construction now has a bound below 2^-142 for the
family specified below. It allows 5--64 clean closed partial episodes,
followed by a live tail whose internal activity is unrestricted. The route
and every active-row count are included. This is not a whole-code result,
nor complete coverage of six-defect histories.

The earlier five-defect 46% certificate is unchanged. This continuation
tests whether the partial-reset argument supplies a reusable interface,
rather than only another isolated count profile. It changes no encoder.

## A weighted count that combines the outer, route, and resets

Use the same q=2^32, M=8192, RS[256,128] outer, regional route, and
512-coordinate inner as in [the prior note](UNIFORM_PAIR_REGIONAL_RANK.md).
Fix the live intervals before exposing the route. A closed interval has
l incoming states, preceded by its activation batch. Its possible input
support therefore occupies l+1 consecutive batches.

For a fixed message and route, its activation state equals the activation
input. If that input has weight w, a reset within l updates costs at most
2*l^w*q^-w. The initial support is known before the future maps; no choice
of those coordinates is added to the bound. Across disjoint episodes use
the product of these weights, but only one total density factor two.

Assign each input position weight zero outside the allowed intervals,
l/q in the activation batch of a closed interval of length l, and one
elsewhere in the allowed input support. In region r let b_r be the mean
of its position weights. For a fixed nonzero RS row, the route-averaged
weight is the product of b_r over its nonzero coordinates.

A nonzero RS word has at most 127 zero coordinates. For an exact zero
set Z of size z, there are at most q^(128-z) such words. Thus the weighted
nonzero enumerator of one row is at most

    F = sum_(z=0)^127 q^(128-z) * e_(256-z)(b_1,...,b_256).

Here e_j denotes the elementary symmetric polynomial: sum the products
over all j-element subsets of the 256 arguments. This bound deliberately
overcounts words, but never uses the zero word as an actual witness.

For s distinct rows, their regional positions are sampled without
replacement. The mean product of s nonnegative position weights is at
most the s-th power of the mean, by the elementary-symmetric inequality.
Expand the nonzero codeword sums with positive coefficients. This gives
a joint bound F^s without assuming independent row routes. Consequently,
the reset contribution for the fixed intervals is bounded by

    2 * sum_(s=s_min)^M binom(M,s) F^s.

Usually s_min=1. If all activations must have pair-closed support, a
single active RS row cannot activate a pair: one batch lies in one region
and contains at most one coordinate of that row. Then s_min>=2.

## Integer implementation and a useful limitation

For a region with u ordinary batches and activation weights whose integer
sum is v, write n=qu+v. Its mean weight is n/(16q). If R regions have
nonzero n, the one-row bound simplifies to

    F = [sum_(j=0)^(R-129) coefficient_of_z^j product_r(n_r+16z)]
        / (16^R q^128).

When R<=128 it is zero. Full regions contribute q+z after cancelling
their factors 16. This gives a small positive integer polynomial near
the active-row crossover, rather than a large-alphabet enumeration.

The helper maximizes polynomial coefficients over regional offsets and
over early lengths with the same defect budget. It keeps total length,
necessary internal defects, and region count. Replacing an overlap between
two regional fragments by two distinct regions only increases every
elementary symmetric polynomial at fixed sum. This permits a valid
split-region envelope when at most 256 regions result. Coefficientwise
maxima may still combine different geometries and can be loose.

The first six-defect replay using only this weighted enumerator reduced
137 residual coarse cells to 135. It did not close the family. In a
remaining case, many rows can avoid an activation batch while landing in
later batches of the same interval. The bound then pays no reset cost
for those rows. A real clean episode cannot ignore their later inputs.

## Clean supports constrain every input in the episode

Call a closed proper episode clean if its incoming supports only decrease
and every incoming pair is either empty or full. Its length l is at most
eight. Backward closure along the l distinct matching phases forces its
initial support U to satisfy

    |U| >= 2^l.

Every input coordinate used during that episode also lies in U. Indeed,
the current support is a union of current pairs, so A Q is supported
there. The next support is a subset, and X=Q_next-AQ. The activation
input itself equals the initial state. This is a deterministic containment
statement, not a claim that states are uniform.

If a region contributes s nonzero row coordinates across u batches of the
episode, those positions use only coordinates in U. A fixed state
coordinate occurs at most u times there. Therefore

    |U| >= ceil(s/u).

Together, these two inequalities retain input/state compatibility that
the activation-only weighted enumerator discards.

## A family that closes without a defect cutoff

Fix m in {5,...,64}. Consider histories with m closed proper episodes,
followed by one live interval extending through Q_(L-1). Require:

- Each closed episode is clean and has length l_i in {1,...,8}.
- Its l_i+1 input batches straddle exactly one regional boundary.
- The regions touched by different intervals are disjoint.
- The final input interval touches 129-2m regions.

Thus the allowed input support touches exactly 129 regions. No condition
is imposed on births or one-sided pairs inside the final live interval.
The argument below does not use an output-weight bound: it excludes every
nonzero message with this geometry except with the stated probability.

Fix s active outer rows. Each such row has a one-dimensional shortened RS
code and is nonzero in all 129 allowed regions. There are at most q^s
messages on those rows. In a closed episode split into u_i and v_i batches,
u_i+v_i=l_i+1, both regions therefore contain s nonzero inputs. Its initial
weight w_i satisfies

    w_i >= ceil(s/min(u_i,v_i)),    w_i >= 2^l_i.

The activation batch contains at most one position of each active row,
so w_i<=s. With d=max_i l_i, this implies s>=2^d.

For s specified rows, route containment across these two regions costs
at most (u_i*v_i/256)^s, using the regional falling-factorial inequality.
The terminal interval contributes a factor at most one. Set

    a_i = 32-ceil(log2 l_i),
    p_i = floor(a_i/min(u_i,v_i)).

The fixed-message reset bound satisfies

    (l_i/q)^w_i <= 2^(-p_i*s).

Balanced splits maximize the route factor and minimize this reset charge
simultaneously. For l in {1,...,8}, define u=floor((l+1)/2), v=l+1-u and

    g_l = (uv/256) * 2^(-floor((32-ceil(log2 l))/u)).

For fixed lengths, put y=Mq*product_i g_(l_i). The complete contribution
for fixed interval placements is at most

    2 * sum_(s=2^d)^M binom(M,s) (q*product_i g_(l_i))^s
      <= 2*y^(2^d)/(1-y),       provided y<1.

For each m and maximum length d, the verifier uses
g_d*(max_(l<=d) g_l)^(m-1). This bounds every length order without a
separate combinatorial factor: interval endpoints already specify their
ordered lengths. Every y is below one.

There are at most (L+1)^(2(m+1)) placements, bounded by 26(m+1) bits.
The exact rational check covers all 60*8=480 cases. The largest bound
after its placement union is -151 at (m,d)=(5,1). Adding nine bits for
the case union gives the entire family upper bound **2^-142**.

The probability includes the route and pair maps, with setup fixed for
all messages. It is not conditioned on a good-route event, and it covers
all active-row counts. The broad family above remains only one part of
the full code's possible histories.

## Relation to the six-defect frontier

For h<=6,a=0,c=6 with no full reset, all six defects are activations from
zero. The five early proper episodes are therefore clean. If T ranges
from 1893 through 1908, their regions plus the terminal interval can
reach at most 129 regions. The RS minimum forces equality, disjointness,
and all five early intervals to cross a boundary. The family theorem
therefore covers that entire count band, not just the previously checked
T=1893 example.

The replay also applies the same theorem to individual joint-length cells
whenever their terminal interval touches 119 regions. This is a proof
reuse gate, not a new construction or heuristic exception.

A completed replay with the broad count-band gate left 120 coarse residuals,
with upper 6491; it did not close six defects. After adding the joint-cell
gate, the next replay exceeded its 50-second cap. That run is incomplete;
no updated complete residual count is claimed from it.

Whole six-defect coverage still fails. In particular, the first remaining
families allow additional RS zeros or share regions between episodes.
An active row can then be absent from a region that the containment
argument wants to charge. The next step must retain those zero choices
and shared occupancies jointly; treating every active row as present in
every touched region would be invalid once more than 129 regions are used.

## Checks and preservation

Run serially:

    python -B diagnose_pair_reset_weights.py --test
    python -B diagnose_pair_reset_weights.py
    python -B diagnose_pair_reset_weights.py --clean-crossings
    python -B diagnose_pair_reset_weights.py --cover

Tests include 625 exact GF4 weighted RS enumerators, 2500 weighted
without-replacement inequalities, 128 fragment envelopes, 69900 coordinate
containment checks, and 1062 backward-closure checks. The family check
uses exact fractions, covers 480 cases, and takes under one second.
The six-defect diagnostic takes tens of seconds and rejects at its
50-second cap; it does not claim closure when residual cells remain.
Polynomial tables and caches have fixed size caps. Only the new source
and this compact note are retained. No raw files or Python caches are made.

The previous proofs, verifiers, and certified GF32 control are unchanged.
C: had over 22 GiB free before numerical work, but dropped to about 0.6 GiB
during the final replay. Local numerical runs stopped at that safety-floor
breach. Further runs require revalidated space or the configured remote host.
No unrelated files were deleted. Runs were serial; no
benchmark, production/paper edit, commit, push, or promotion occurred.

Continue the proof rather than redesign now. The construction already
provides the support-containment fact needed to repair this counting loss.
Next combine reset weights with the permitted RS zero patterns and shared
regions, and test an explicit remaining family before enlarging coverage.
If that loses essential information again, identify the missing interface
before deciding whether an additional inner operation is justified.
