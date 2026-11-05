# Partial-packet shortening and mask symmetrization

Status: a proved route-moment bound and a mask-sum interface, not a
whole-code distance certificate. The proposed physical lane shuffle is
not implemented. A complete invariant mask sum can use the same bound
without changing the existing route.

The [one-diagonal inner](SINGLE_DIAGONAL.md) randomizes an arbitrary
coordinate set A. Unlike the sets in [Packet route](PACKET_ROUTE.md), A
can cut packets. The fixed lane order can then concentrate A on particular
outer rows. This note removes that concentration either by sampling lane
shuffles or by averaging a complete family of support masks.

## Route and exact shortening

Fix a field F and an even B. Put k=B/2. The outer is the direct sum of
pG fixed MDS [B,k] rows over F, grouped into G groups of p rows. Thus
N=pGB. Write c[g,l,j] for position j of row l in group g.

There are B regions, each containing G packets of p coordinates. For
each region j, sample an independent uniform permutation pi_j of G groups.
For each destination packet h in that region, independently sample a
uniform permutation sigma_(j,h) of its p lanes. Place

    X[jpG+ph+a] = c[pi_j(h),sigma_(j,h)(a),j].             (1)

Sample the entire route once and fix it for all messages. Its randomness
is independent of any inner setup. The p lanes are independent outer
rows, not the field lanes inside one payload record.

Fix A independently of this route. For each source row (g,l), let
T_(g,l)(A) be its number of coordinates routed into A. Write C_pi,sigma
for the routed outer and define

    D_(pi,sigma)(A) = dim_F(C_pi,sigma intersect F^A).

The MDS shortening formula gives, for every route and every scalar set A,

    D(A) = sum_(g,l) (T_(g,l)(A)-k)_+,
    D(A) = |A|-N/2+D(A^c).                                (2)

The second identity follows from (t-k)_+-(B-t-k)_+=t-k. Neither identity
requires whole packets. Without lane shuffling, the same exact identities
hold, but the moment bound below need not hold for a fixed A.

## A monotone comparison for assigning counts

Let a_1,...,a_G be fixed numbers in an ordered finite set. Let f_1,...,f_G
be nonnegative nondecreasing functions on that set. For a uniform
permutation pi,

    E_pi product_g f_g(a_pi(g))
      <= product_g (G^(-1) sum_h f_g(a_h)).                (3)

To prove (3), expand each f_g as a nonnegative constant plus nonnegative
multiples of upper-threshold indicators. It suffices to compare each
resulting product of threshold indicators. Suppose r factors remain.
Their eligible sets of labeled population positions are nested. Order
their sizes as n_1<=...<=n_r. If some n_i<i, the joint probability is zero.
Otherwise sampling the r labeled positions without replacement gives

    product_(i=1)^r (n_i-i+1)/(G-i+1)
      <= product_(i=1)^r n_i/G.                           (4)

The right side is the probability under independent draws. Expanding
back proves (3), including repeated population values. Monotonicity and
nonnegativity are essential; (3) is not asserted for arbitrary functions.

## Arbitrary-set moment

For tau>=1 and 0<=s<=N, define

    m = min(s,N-s),     z=m/N,
    H_B(z,t) = sum_(i=0)^B binom(B,i) z^i (1-z)^(B-i)
                              t^((i-k)_+),
    M_p(s,tau) = tau^((s-N/2)_+) H_B(z,tau^p)^G.           (5)

Then every fixed scalar set A satisfies

    E_(pi,sigma) tau^D(A) <= M_p(|A|,tau).                 (6)

First prove the uncomplemented bound, with z=|A|/N and no prefactor.
In region j, let a_(j,h) be the number of selected coordinates in
destination packet h. The group permutation assigns this multiset of
counts without replacement to the source groups. Conditional on a count a,
the lane shuffle selects a uniform a-subset of the p source rows.

Condition on every other region. The factor for a source group is
nonnegative and nondecreasing in each of its row-selection indicators.
Averaging over uniform a-subsets makes this factor nondecreasing in a:
use the first a lanes of one uniform lane ordering to couple all sizes.
Equation (3) therefore permits replacement of the current count assignment
by independent draws from its empirical count distribution. The selected
lane subsets are also independent between groups after this replacement.
Repeat the comparison for every region. The resulting groups are independent.

Within one group, write T_l for the occupancy of row l. The p rows need
not be independent. Holder's inequality nevertheless gives

    E tau^(sum_l (T_l-k)_+)
      <= product_l (E (tau^p)^((T_l-k)_+))^(1/p).          (7)

Each row has the same marginal distribution: a sum of B independent
Bernoulli variables with probabilities

    alpha_j = |A intersect region j|/(pG).

Their sum is |A|/(pG). The integer function t^((i-k)_+) is convex for
t>=1. Replacing two Bernoulli probabilities by their average increases
the expectation of this function: its mixed coefficient is the
nonnegative second difference, and their product increases. Repeated
pairwise averaging gives the Bin(B,|A|/N) upper bound in (7).
Multiply the G group bounds. Finally use (2) on A^c when |A|>N/2.
This proves (6).

The exponent p on tau and the number G=N/(pB) of groups both remain.
Equation (6) uses all B regions; it does not claim an interval endpoint
bound for arbitrary masks. It is valid at tau=1, with both sides one.

## Use without physically shuffling lanes

Let pi denote the existing packet route with fixed lane order. Let Gamma
be the group of all independent destination-packet lane permutations.
For any gamma in Gamma, pulling back coordinates gives

    D_(pi,gamma)(A) = D_pi(gamma A),                       (8)

where replacing gamma by its inverse only changes the action convention.

Consider a finite family C of pairs of masks (W,H), closed under the
joint action (W,H)->(gamma W,gamma H). Let a(W,H)>=0 be an invariant
weight. Then reindexing this finite sum gives the exact identity

    sum_C a(W,H) E_pi tau^D_pi(H)
      = sum_C a(W,H) E_(pi,gamma) tau^D_(pi,gamma)(H)
      <= sum_C a(W,H) M_p(|H|,tau).                       (9)

For inner batches containing whole packets, joint lane permutations
preserve each batch's w=|W|, h=|H|, and o=|W intersect H|. Thus a complete
family specified only by these statistics is invariant. A rank bound that
depends only on such statistics and dimensions can use (9), once that
bound is established. The matrices of the actual inner need not commute
with these permutations.

Equation (9) must be applied after replacing realizable masks by a
complete invariant family and using an invariant amplitude bound. Actual
state realizability need not be invariant. Nor does (9) assert (6) for
one fixed mask under the unshuffled route. For example, selecting lane zero
of every packet always selects all B coordinates of G complete outer rows.

For a class of M mask pairs with common |H|=s and an amplitude bound q^E,
the direct fractional-moment consequence is

    Pr[an eligible output in this class]
      <= min(1, M tau^(E-s) M_p(s,tau)),   1<tau<=q.        (10)

This uses the inverse law in SINGLE_DIAGONAL.md and (9). The mask
multiplicity M is outside the fractional power. Replacing M by
M^(log_q tau) in (10) is not justified: different masks have different
shortened spaces.

There is a useful sharper grouping when only W varies. Fix H first, and
suppose there are L eligible W masks, each with a uniform amplitude bound
q^E. All these outputs use the same shortened space D(H). Sum over W
before taking the fractional power. With nu=log_q(tau), the resulting
fixed-H bound is L^nu tau^(E-|H|+D(H)). If there are M_H eligible H masks
in a complete invariant family, with the same L, E, and |H|=s, (9) gives

    Pr[an eligible output in this class]
      <= min(1, M_H L^nu tau^(E-s) M_p(s,tau)).            (11)

For complete families with fixed per-batch h_t,w_t,o_t,

    M_H = product_t binom(b,h_t),
    L = product_t binom(h_t,o_t) binom(b-h_t,w_t-o_t).

Use only feasible overlap counts and an amplitude bound uniform over this
family. Equation (11) legitimately puts the W multiplicity inside the
fractional power, but leaves the H multiplicity outside. It does not
assume independence between output supports or between shortened spaces.
Closing the full mask sum remains a separate proof gate.

## A partial period class that the new moment distinguishes

Consider the fixed repeating schedule

    (h_t,w_t) = (b,b), (b,0), (1,0).

It passes the adjacent-weight and zero-output-run relaxations in
SINGLE_DIAGONAL.md. Its fresh-output rank bound R is zero. The fallback
D(H)<=|H|/2 therefore leaves only half a field-symbol charge per period.
This was an obstruction to that fallback, not a low-distance witness.

The partial-packet moment gives a different result for this schedule.
Temporarily count N=3rb coordinates in r complete periods. There are
only b choices of H per period; both full masks and all W masks are fixed.
Thus the mask multiplicity is b^r, not a global binomial estimate. The
output count is at most q^w, where

    w=N/3,       |H|=(2+1/b)N/3.

Choose tau^p=4, assuming tau<=q. The complemented density is
z=(1-1/b)/3. Since t^((i-k)_+)<=1+t^(i-k),

    H_B(z,4) <= 1+4^(-B/2)(1+3z)^B
               = 1+(1-1/(2b))^B < 2.

The class term from (10) is consequently at most

    b^(N/(3b)) tau^(-N/6) 2^(N/(pB)).

Its base-two exponent per coordinate is bounded above by

    -1/(3p) + 1/(pB) + log2(b)/(3b).                     (12)

For p=4,B=64,b=128, this coefficient is -47/768. For p=2,B=64,b=64,
it is -49/384. These are exact arithmetic evaluations of (12), not
numerical proof searches.

This schedule describes a live interior segment. The actual zero initial
state requires boundary batches; a constant number changes the exponent
by O(b), with parameters fixed. No finite bound for those boundaries or
union over all schedules is claimed here. The calculation establishes
only that the previous repeated-pattern loss is not intrinsic to the new
route bound. A useful next calculation should retain per-batch mask counts
and compare complete bounds for each class.

## Exact checks and scope

`theory/test_shuffled_packet_route.py` enumerates all 4096 labeled routes
for p=2, B=4, G=2. Symmetries reduce all 65536 scalar masks to 126 moment
representatives. The test checks exact rational moments, complements,
small generator-rank fixtures, the monotone comparison, and joint-mask
orbit sums. A negative control rejects (6) for an unshuffled fixed mask.
Root's bounded serial replay passes all 4096 labeled routes, 126 support
classes, 516096 complement identities, 252 rational moment bounds, 56 direct
rank fixtures, 4589 monotone comparisons, and six joint-mask orbit checks
covering 768 pairs. The negative control rejects the fixed-mask bound
without lane shuffles. All 65536 masks are covered by the checked support
symmetries; the inner's realizable masks are not asserted symmetric.

These checks support the stated finite interfaces. They are neither a
K20 distance certificate nor a performance measurement. The next step is
a sound amplitude-and-mask sum for the one-diagonal recurrence, retaining
the actual inverse support, all mask multiplicities, and the free terminal.
