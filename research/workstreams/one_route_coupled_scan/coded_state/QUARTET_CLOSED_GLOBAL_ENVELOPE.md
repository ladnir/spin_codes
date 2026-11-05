# All destinations, and the limit of a uniform one-step envelope

Checkpoint, 2026-10-30. The encoder, outer, routing, field multipliers, and
45%/40-bit objective are unchanged. This checkpoint covers the proof
interface for the remaining r=2100 class. It does not close that class or
the whole-code certificate.

The previous regional bounds retained only paths inside small population
sets. Here we first sum every destination. We then prove that replacing
the resulting transfer by one uniform product-weighted row bound is too
loose. This rejects a proof shortcut, not the construction or the exact
closed-set moment.

## Sum the destinations exactly

Use the seven source types (Z,V0,V1,B,L0,L1,L2). Let h_i>0 weight each
destination of type i. Write a=b^lambda, d=b^(-tau), y=b^(-6/5), and let
x mark original-inactive positions. Put T=1+xy and

    base = h_V0 + a*T*h_L0,
    N = sum_(j=1)^4 binom(4,j)*(x*h_Z)^j*base^(4-j).

For a labelled source quartet with type counts m, its destination-weighted
polynomial K_m has three cases:

- If it contains both a blocked and a selected coordinate, K_m=0.
- If it contains a selected coordinate but no blocked coordinate, K_m is
  T^4 times the product of the seven entries
  (h_V0+a*h_L0, h_V1, h_B, 0, a*h_L1, a*h_L2, d*h_B+a*h_L2),
  with multiplicities m.
- If it has no selected coordinate, K_m=N. When it also has no blocked
  coordinate, add the product of (base,h_V1,h_B) with the three ordinary
  source multiplicities.

The first case enforces the incompatible fullness requirements. In the
second case the whole quartet is touched. In the third case, N counts
non-full quartets; the additional term counts a full quartet without a
selected predecessor. No destination population is omitted.

Define the homogeneous source polynomial

    Q(z) = sum_(|m|=4) multinomial(4;m)*K_m*z^m/prod_i h_i^m_i.

For a source population n of size512, the exact product-weighted row sum is

    R(n) = [z^n] Q(z)^128 / multinomial(512;n).

This is a coefficient identity, not an independent-coordinate model.
The random quartet partition is accounted for by the source multinomial.
The factorization passes420 exact local comparisons. Eight width-eight
tests check the all-destination coefficient bound against direct expansion.

`quartet_closed_factored.py` also extracts the total selected and blocked
source counts by a two-dimensional positive recurrence. It keeps the exact
selected-age split and bounds only the remaining source counts by positive
evaluation. This is a reusable all-destination interface, but the present
full-width evaluations are FLOAT proposals.

## A bound covering all source populations

Evaluating at p=n/512 gives

    R(n) <= Q(p)^128 / Pr[Multinomial(512,p)=n].

The probability in the denominator is minimized by balanced category
counts. To see this, write its logarithm, apart from constants, as the sum
of f(n)=n*log(n)-log(n!). The increments are
f(n+1)-f(n)=n*log(1+1/n), which increase with n; the increment at zero is
zero. Balancing therefore minimizes the sum. The resulting probability
is an exact rational number. There are8550 exhaustive small checks.

A nonnegative quartic also admits Bernstein degree elevation. If its
normalized degree-four coefficients are k_m, the elevated coefficient at
|v|=D is

    B_v = sum_m k_m*prod_i binom(v_i,m_i)/binom(D,4).

The largest B_v bounds Q on the probability simplex. All691 tested
degree-elevation rows satisfy the exact Vandermonde normalization.

Balancing the three selected-age weight ratios reduces the source
polynomial to five categories. Its optimized FLOAT bound is still about
+987073 bits in the final r=2100 expression. Allowing all seven weights
independently gives about+1017571, including a larger multinomial prefactor.
Even deleting that prefactor diagnostically leaves about+914598.
The maxima in these proposals occur at homogeneous source populations;
increasing the tested Bernstein degree does not repair them.

These are failed upper bounds, not failure probabilities. A separate
directed argument below shows that optimizing this envelope further cannot
close the class in the stated parameter range.

## A directed obstruction for every product weighting

Consider any fixed positive product potential Phi(n)=prod_i h_i^n_i,
normalized by h_Z=1. Consider a common positive inactivity marker x and,
optionally, a payment-cap marker 0<eta<=1. The latter changes y to y*eta
and charges eta^(-255463) at the end.

Let Rmax be a uniform upper bound on every one-step weighted row sum.
For each source type i, the homogeneous population512*e_i is among these
rows. If Q_i is the weighted polynomial for four copies of that source,
then Rmax >= Q_i^128. This observation involves exact rows, without the
multinomial mean-atom relaxation.

Write the positive monomials of Q_i as

    c_j * b^(e_j/2560) * x^inactive_j * eta^paid_j
        * prod_k h_k^dest_(j,k) / h_i^4,

where e_j=707*activity_j-1905*deaths_j-3072*paid_j. Choose a rational
probability distribution p_(i,j) over the monomials of each row. Let
M_(i,k)=sum_j p_(i,j)*dest_(j,k)/4. Choose a stationary distribution alpha
for this stochastic matrix. Then weighted AM-GM, averaged using alpha,
cancels every log(h_k) term exactly.

The checker constructs two such rational witnesses from arbitrary rounded
numerical proposals. Each row probability is a normalized positive integer
vector. It solves the seven-state stationary equations over rationals and
checks every flow equality. A rational mixture of the two witnesses has
mean original inactivity exactly4*2100/8192. Their individual means are
approximately0.949375 and1.125951, bracketing the required1.025390625.
Thus the x term cancels exactly against the original regional quotas.

The mixed mean paid count is approximately0.260396 per quartet, below
255463/2^19, approximately0.487257. Consequently the eta contribution to
the final expression is nonnegative for every0<eta<=1. Retaining the
payment cap cannot evade this particular envelope obstruction.

There are132 positive monomials in the seven homogeneous rows. Rational
logarithm intervals bound their entropy contributions outward. The final
expression, even before its nonnegative terminal and multinomial-prefactor
costs, is bounded below by an affine function of ell=log2(b). Its directed
intercept is greater than964746.633048 and its exact rational slope is
approximately-22313.431086. For the entire interval1<=b<=8, a compact
trusted lower endpoint is

    3856048867539029 / 4294967296 > 897806.

At the currently selected b=2^(5/2), the trusted lower endpoint is

    3903966595927635 / 4294967296 > 908963.

The field multipliers lambda and tau remain fixed throughout this claim.
This rules out a *uniform, time-independent, one-step product-potential
envelope*, with the common markers just described. It does not rule out
population-dependent potentials, coupled multi-update bounds, or the
exact regional transfer. The stationary flow used in the AM-GM witness
balances types only on average; it is not a feasible codeword trajectory.

## Test a population-dependent repair

The natural missing interaction is between blocked and selected counts.
A selected coordinate requires a full next quartet; a blocked coordinate
forbids one. Ending runs creates blocked coordinates that must remain in
the next step's source population.

Test the potential

    Phi(n)=prod_i h_i^n_i * 2^(-gamma*n_B*n_S/512).

`quartet_closed_interaction.py` sums the joint outgoing(B,S) polynomial,
then applies this factor. It does not replace the product of counts by a
product of expectations. Eight small rational polynomial comparisons pass.
Full-width calculations are FLOAT only.

Optimizing against homogeneous source rows improves that restricted
diagnostic from about+778k at gamma1 to+402k at gamma6, still positive.
The gamma6 proposal is worse on the mixed population
(0,0,0,64,0,0,448): its repeated-row diagnostic is about+1.659 million.
Other sampled mixed rows also remain positive. These are location tests,
not an optimization over all populations or a new lower bound. Do not
claim that this quadratic potential closes the proof or exhausts the
population-dependent options.

## Consequence and next discriminating test

Do not continue tuning a uniform one-step product norm. Preserve the
exact all-destination factorization, but retain joint endpoint populations
through at least two updates. In particular, a death count in the first
update must become blocked source mass in the second, with the actual
quartet exclusion probability. Include transitions out of the previously
covered sets. Start with the homogeneous selected source and the mixed
blocked/selected source above; these expose different losses in the current
envelopes. Keep the original regional quota and the single terminal test.

This is proof-interface progress, not an additional certified distance or
margin. The full goal remains open, including other outer-row counts and
the final failure-budget audit.

## Reproduction and resources

The five new helpers are `quartet_closed_factored.py`,
`quartet_closed_global.py`, `quartet_closed_global_ages.py`,
`quartet_closed_global_obstruction.py`, and `quartet_closed_interaction.py`.
Run their main entry points serially with Python3 and NumPy on a host
meeting the20GiB free-space guard. The ages helper accepts `--paid-cap`.
The obstruction helper verifies the rational flows, quota mixture, payment
cap, and directed endpoints without requiring optimizer convergence.

All runs used one-thread, in-memory Peach sessions with180-second bounds;
selected computation times were under two seconds. They wrote no files.
An output lost during context rollover was rerun only after an authoritative
process check confirmed completion. No native benchmarks, construction
changes, production edits, paper edits, commits, or generated dumps.
