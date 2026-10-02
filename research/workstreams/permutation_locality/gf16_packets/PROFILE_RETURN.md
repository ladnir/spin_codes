# Cancellation with a Constrained Expansion Profile

The earlier character bound permits every packet to take its worst
expansion pattern independently. That can effectively discard the
expansion's minimum distance. This bound enforces the minimum weight
before maximizing the character products. The encoder is unchanged.

## Weighted Feedback Coefficient

Use the iid comparison input from [the refresh bound](REFRESH_BOUND.md).
Fix a nonzero entering state a, and write y=Aa. There are W four-bit
packets. Each input packet is zero with probability 1-p and otherwise
uniform among the 15 nonzero values. Fix 0<z<=1.

For p<=15/16, define c=1-16p/15 and b=p/15. The packet distribution has
the positive decomposition

    Pr[X_i=x] = c 1_{x=0} + b       (x in F_2^4).

This is an identity for the comparison distribution, not new randomness
in the encoder. When p>15/16, the implementation uses the earlier bound.

For feedback character u, let r_i=wt(C_i^T u), where C_i contains the
four feedback columns at position i. Put v_i=wt(y_i) and
D_r=(1+z)^(4-r)(1-z)^r. The weighted character factor at that position is

    F_i(u) = E[z^wt(y_i+X_i) (-1)^(u dot C_i X_i)]
           = c z^v_i + b epsilon_i D_{r_i},
    epsilon_i = (-1)^((C_i^T u) dot y_i).

By Fourier inversion, the weighted lazy return to zero is

    R_a = E[z^wt(Aa+X) 1_{CX=a}]
        = 2^-s sum_u (-1)^(u dot a) product_i F_i(u).

In the product expansion, choosing c z^v_i at every position gives the
character-independent term c^W z^wt(y). Its contribution to R_a is zero
because a!=0. Apply the triangle inequality to the remaining terms. Since
c,b,D_r are nonnegative, their absolute sum is at most

    G(v,r) = product_i (c z^v_i + b D_{r_i}) - c^W z^(sum_i v_i).

Thus R_a<=2^-s sum_u G(v,r(u)). No independence between output weight
and feedback is used.

## Maximize Subject to Expansion Weight

Let d_A be the minimum nonzero expansion weight. Each v_i lies in [0,4]
and sum_i v_i>=d_A. The function G is coordinatewise nonincreasing and
convex in the real vector v. To see this, expand it over the nonempty
sets S of positions choosing the uniform term:

    G(v,r) = sum_{S nonempty} c^(W-|S|) b^|S|
               product_{i in S} D_{r_i} z^(sum_{i not in S} v_i).

Each summand is a nonnegative constant times an exponential of a linear
function. Therefore reducing the total weight to d_A can only increase
G. On the polytope 0<=v_i<=4, sum_i v_i=d_A, a maximum occurs at a
vertex. A vertex has floor(d_A/4) weights four, at most one weight
d_A mod 4, and all other weights zero.

It remains to place those weights. D_r decreases with r. If v_i>=v_j
and r_i<=r_j, pairing the larger weight with the smaller r changes the
two-factor product by

    c b (z^v_i-z^v_j)(D_{r_j}-D_{r_i}) >= 0.

The subtracted term is invariant under permutations. Repeated exchanges
therefore maximize G by pairing expansion weights in decreasing order
with character weights in increasing order.

Call that maximum G_d(r). It depends only on the character histogram,
whose multiplicity is already counted exactly. The common lazy-return
bound is

    R_a <= C_profile := 2^-s sum_u G_d(r(u)).

For the actual expansion, d_A=48 and W=32. The relaxed maximizing profile
contains twelve weight-four packets and twenty zero packets. It is a
worst-case upper bound, not an assumption about the entering state.

## Use and Verification

The `profile-return` kernel may replace the arbitrary-to-zero matrix
entry by alpha C_profile + beta H/L, if smaller. All other entries and
the three-coordinate domination invariant remain unchanged.

Profile pairings use exact integer counts. Character factors, their
products, and their final sum use outward Arb arithmetic. Floating
calculations only propose witnesses. Tests exhaust all weight assignments
for up to three packets and check the canonical relaxation directly.
They also enumerate every small input and entering state, check raw
lazy-return bounds, and compare eight consecutive bounded transitions
with exact transitions for one, two, and three updates.

This local bound does not alone certify the code. A whole-code result
requires a complete dense cover, matching sparse covers, and their
independently replayed sum.

## Actual Profiles and the Ten Exceptional States

The `shape-return` refinement retains the actual packet-weight histogram
instead of relaxing every state to the canonical profile. For each
histogram, rearrangement gives a bound shared by every state in that
class. The exact class multiplicity permits both a maximum over states
and an average over all nonzero states.

The canonical twelve-full-packet profile occurs in only ten states of
the actual expansion. Those states are handled separately. Fix one such
state a and partition the packet positions into S_0, where y_i=0000, and
S_4, where y_i=1111. The preceding character factors are then exactly

    F_0(r) = c + b D_r,
    F_4(r) = c z^4 + b (-1)^r D_r.

Define H_j(u)=product_{i in S_j} F_j(r_i(u)), for j in {0,4}. Its
weighted return is E_u[(-1)^(u dot a) H_0(u)H_4(u)], where u is uniform
over all feedback characters. Cauchy--Schwarz gives

    R_a <= sqrt(E_u[H_0(u)^2] E_u[H_4(u)^2]).

Both moments are evaluated from exact character histograms on the actual
position sets S_0 and S_4. Their correlation is not discarded by an
independence assumption: Cauchy--Schwarz explicitly allows that correlation.
These two marginal histograms are much cheaper to evaluate repeatedly
than the joint histogram. An optional signed joint census computes R_a
itself as a diagnostic. It is checked against exhaustive small-input
enumeration; the default kernel uses the cheaper Cauchy bound.

Let C_h bound a return for any state in ordinary histogram class h, and
let n_h count that class. Let C_a bound each exceptional state's return.
Every nonzero state belongs to exactly one of these disjoint sets. Hence
the two multipliers are

    c_shape = max({C_h}_h union {C_a}_a),
    cbar_shape = (sum_h n_h C_h + sum_a C_a)/L.

The kernel takes the minimum with the corresponding earlier zero-entry
bounds, after applying alpha and the refresh contribution. It does not
assume that the actual incoming state is uniform. The average applies
only to the component already dominated by a uniform-density measure.

## Condition on the Number of Active Packets

The `conditioned-return` refinement conditions on packet occupancy before
bounding a return. This can improve a bound taken after averaging over
occupancies. It changes neither the encoder nor the three-coordinate
domination invariant.

Let J count the nonzero packets in an epoch. Under the iid comparison
distribution, J has distribution Bin(W,p). Conditional on J=j, positions
form a uniform j-subset, and their values are independent uniform nonzero
packets. For each nonzero entering state a, define

    Y_a = wt(Aa+X),
    H_{j,k} = max_{a!=0} E[z^(k Y_a) | J=j],
    Hbar_{j,k} = (1/L) sum_{a!=0} E[z^(k Y_a) | J=j].

Here L=2^s-1. These moments come from the actual expansion histograms. For
packet weight v, its factor in the occupancy-generating polynomial is

    z^(k v) + x ((1+z^k)^4-z^(k v))/15.

Multiply the W factors, extract the coefficient of x^j, and divide by
binom(W,j). Taking the maximum or multiplicity-weighted average over
expansion histograms gives H_{j,k} or Hbar_{j,k}, respectively.

Let b_j bound max_{a!=0} Pr[CX=a | J=j], and let
nu_j=Pr[CX!=0 | J=j]. Both are supplied by the existing feedback census.
For any k>1, Holder's inequality gives

    E[z^Y_a 1_{CX=a} | J=j]
      <= H_{j,k}^(1/k) b_j^(1-1/k).

For the uniform-density coordinate, introduce an independent uniform
nonzero state U solely to evaluate its dominating measure. Then
Pr[CX=U | J=j]=nu_j/L, so

    E[z^Y_U 1_{CX=U} | J=j]
      <= Hbar_{j,k}^(1/k) (nu_j/L)^(1-1/k).

Neither inequality assumes independence between emitted weight and the
return event. The uniform auxiliary state describes the existing
dominating measure, not the actual conditioned state of the encoder.

For each j separately, take the minimum of these bounds for
k in {2,3,4,6,8,12,16}, the first-moment bound, and the existing
pointwise distance bound z^max(0,d_A-4j) times the event-probability bound.
Average the resulting bounds against Pr[J=j]. The two resulting values
may tighten the lazy contributions to M->Z and U->Z; refresh contributions
and all other entries remain unchanged. Taking these minima separately
for each j is valid because every candidate bounds that same conditional
expectation.

The implementation also takes these minima before maximizing over
expansion histograms. For histogram h of weight d_h, use its own moments
and the pointwise factor z^max(0,d_h-4j). This gives a bound c_{h,j} for
each state in that class. Thus max_h c_{h,j} bounds an arbitrary state,
and (sum_h n_h c_{h,j})/L bounds the uniform-density component. The latter
is combined with the direct uniform-state Holder bound above by taking
the smaller value. All class sizes n_h come from exact enumeration.

Outward evaluation uses exact rational packet probabilities, positive
polynomial coefficients, and Arb arithmetic. Tests check the moments and
raw return bounds against exhaustive small-input enumeration, then check
the resulting envelope over eight successive transitions.

## Exact Weight Tails at a Fixed Occupancy

The `trimmed-return` refinement uses the full output-weight distribution
instead of a few moments. Fix a nonzero entering state a and occupancy
J=j. There are D_j=binom(W,j)15^j equally likely inputs X: choose the j
active packet positions, then choose their nonzero values.

Let h describe the packet weights of Aa. Let N_{h,j,w} count those inputs
with wt(Aa+X)=w. These counts depend only on h. One packet of expansion
weight v contributes the integer polynomial

    y^v + x ((1+y)^4-y^v).

The coefficient of x^j y^w in the product over packets is N_{h,j,w}.
All coefficients are nonnegative integers, and their sum over w is D_j.
The implementation checks this identity for every j and every profile.

The event CX=a contains at most B_j=floor(D_j b_j) inputs, where b_j is
the existing bound on a nonzero feedback atom. Among all events of that
size, the largest sum of z^wt(Aa+X) selects the lightest outputs first.
Let T(N,B;z) denote that sum: consume the counts N_w in increasing w,
stopping after B outcomes. Then

    E[z^wt(Aa+X) 1_{CX=a} | J=j] <= T(N_{h,j},B_j;z)/D_j.

An exchange argument proves the bound. Replacing a selected heavier
outcome by an unselected lighter outcome cannot decrease the sum when
0<z<=1. Events with fewer than B_j outcomes are bounded by adding outcomes.
This is the sharp bound from these two marginals alone. It does not
assume that feedback and output weight are independent.

For the uniform-density coordinate, use the auxiliary independent uniform
nonzero state U. There are L D_j equally likely pairs (U,X). Their output
counts are sum_h n_h N_{h,j,w}. Exactly D_j nu_j pairs satisfy CX=U:
each X with nonzero feedback supplies one such state. Trimming the
aggregate output counts to this integer budget gives a second bound on
the uniform average. The implementation takes its minimum with the
average of the per-profile trimmed bounds. These bounds are then averaged
over J, or used directly in the fixed-occupancy operators.

## Use the Rank of the Feedback Restriction

The `rank-return` bound also constrains how many variable bits can remain
after fixing the feedback. Fix a and a j-element set S of active packet
positions. Write v_i=wt((Aa)_i), and let C_S contain the 4j feedback
columns at those positions. Put r_S=rank(C_S).

Outside S, the emitted output equals Aa. Inside S, relax the input domain
from nonzero packets to all 4j-bit words satisfying C_S X_S=a. This
solution set is empty or an affine space of dimension 4j-r_S. Choose free
coordinates for this affine space. They range over all bit patterns,
while the remaining output coordinates have weight factors at most one.
Consequently,

    sum_{X_S: C_S X_S=a} z^wt(Aa+X)
      <= z^(sum_{i notin S} v_i) (1+z)^(4j-r_S).

This includes relaxed inputs with zero packets, so it remains an upper
bound for the original nonzero-packet domain.

Let d_S=s-r_S. The annihilator of C_S has 2^d_S characters. Since
(1+z)^d_S <= 2^d_S, averaging over uniform j-subsets gives

    E[z^wt(Aa+X) 1_{CX=a} | J=j]
      <= (1+z)^(4j-s)/D_j
         sum_u sum_{S subset Z(u), |S|=j} z^(sum_{i notin S} v_i),

where Z(u)={i:C_i^T u=0}. The outer sum includes every feedback character,
including u=0. This identity for the annihilator count requires no
randomness assumption about the fixed feedback matrix.

For fixed profile h and |Z(u)|=c, the inner sum is maximized when Z(u)
contains the c largest v_i. Indeed, after factoring out z^(sum_i v_i),
the remaining sum is the j-th elementary symmetric polynomial in the
nonnegative numbers z^-v_i for i in Z(u). Replacing one available weight
by a larger weight cannot decrease that polynomial.

Let A_{h,j}(z) be the resulting upper polynomial, summed against the
exact number of characters with each value of c. Its coefficients count
choices of packet positions and are nonnegative integers. The resulting
profile bound is

    R_{h,j}(z) = (1+z)^(4j-s) A_{h,j}(z)/D_j.

For each h and j, take the minimum of R_{h,j} and the trimmed-output
bound before maximizing over h or averaging with class multiplicities.
This order matters: the rank relaxation can be loose at small j, whereas
the trimmed bound retains the zero event count there. All terms bound
the same conditional expectation, so these minima preserve domination.

The iid kernel averages the conditional bounds against Bin(W,p).
The `gf-rank` sparse kernel inserts the same bounds into each epoch's
fixed-occupancy operator, before the existing exact placement calculation.
It uses no iid approximation for sparse packet positions. Only the two
lazy return entries are tightened; the refresh entries, terminal weights,
and state-domination invariant remain unchanged.

Tests compare the integer weight counts with exhaustive small inputs,
verify trimming against every small event, and test the raw rank bound
even for degenerate feedback maps. Separate tests check the resulting
fixed-occupancy operators through eight transitions. These tests support
the local implementation; whole-code claims still require the full
occupancy cover and its independent outward replay.

## Exact Regional Placement for One Active Group

For one active four-row group, we can remove the Bernoulli comparison
across regions. Fix its union-support size u. The independent uniform
row permutations make that union a uniform u-subset of the 256 regions.
Conditional on this union, each occupied region receives one active
packet. Its position is uniform among that region's 2048 positions;
its GF16-randomized value is independently uniform nonzero.

Fix an output tilt lambda>0 and set z=exp(-lambda). Let R_0 and R_1
be the existing nonnegative region operators for zero and one active
packet. These operators already average the packet's position within
the region. They act on the same three state-domination coordinates.
For r regions and j active regions, let V_{r,j} be the row vector
obtained by averaging the ordered operator products over uniform
j-subsets of those regions. Initialize V_{0,0}=(1,0,0). Conditioning on
whether the last region is active gives the exact recurrence

    V_{r,j} = ((r-j)/r) V_{r-1,j} R_0
              + (j/r) V_{r-1,j-1} R_1.

Terms with indices outside 0..r-1 are omitted. Matrix order is retained;
the operators need not commute. The scalar V_{256,u}(1,1,1)^T bounds the
conditional output moment. The local operators remain envelopes, so
"exact placement" does not assert an exact output enumerator.

For bad-output cutoff D, define b_u as the minimum of 1 and
exp(lambda D) V_{256,u}(1,1,1)^T over the supplied positive tilts.
Each candidate bounds the same conditional failure probability, so
the tilt may depend on u. Let F(u) be the authenticated upper CDF on
the expected number of nonzero four-row messages with union size at
most u. Put bbar_u=max_{v>=u} b_v. Summation by parts then gives

    expected bad messages in one specified group
      <= sum_{u=1}^{256} (F(u)-F(u-1)) bbar_u.

The nonincreasing majorant is necessary: the differences of an upper CDF
are not shellwise count bounds. Multiplying this sum by 2048 accounts
for every choice of the sole active group. No independence between
different messages is used.

`single_group.py` evaluates this recurrence and sum with upward rounding.
At D=209715, its 256-bit run gives an upper-bound margin of
42.1225874867 bits. The older Bernoulli-comparison cover, even after
retuning and full subdivision, gave only about 35.85 bits. This removes
the occupancy-one obstruction at 10%; other occupancies remain to be
certified at that distance.

The optional `assemble.py --single-group-exact` path regenerates this
bound and the disjoint range q=2 through the dense handoff minus one.
It still requires the final sum over every occupancy to be below 2^-40.
Tests compare the recurrence against exhaustive noncommuting products,
check upward rounding against rational arithmetic, and check the CDF
majorant and the assembler's disjoint occupancy ranges.
