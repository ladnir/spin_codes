# A sharper bound on the actual BCH row-pair overlap

The convex reduction needs control of the whole row-overlap distribution,
not merely its variance. The first comparison law was unnecessarily broad.
This note constructs a tighter law that bounds every convex expectation.
The encoder, candidate messages, and setup distribution remain unchanged.

Let T80 be the set of weight-80 words in the fixed BCH [256,128] subcode.
Sample U,V independently from T80, and define

    a := |supp(U) intersect supp(V)|.

The set T80 is nonempty by the retained shell-size certificate.
Its translation symmetry gives E a=25. Minimum distance 38 implies
that a belongs to {0,...,61,80}. The atom at 80 has probability 1/|T80|.

The new comparison law, denoted nu_J below, has mean 25 and variance
approximately 42.39813614903227. For every convex g on [0,80],

    E g(a) <= E_{b~nu_J} g(b).

This is an upper-bound comparison, not equality in distribution or an
assertion that the actual overlap law is random-code-like. The full
second moment and the fixed setup-failure target remain unresolved.

## What the individual atom caps already imply

`weight80_fourth.json` contains a bound u_a on each probability Pr[a].
Those bounds cover all possible differences, not just weights through 80.
They use the outer shell caps and degree-14 two-block Christoffel bounds
on |T80 intersect (T80+d)| for each fixed difference word d.
`verify_weight80_dependence.py` replays every cap with rational arithmetic.

Under only these caps, total mass one, and mean 25, there is a common
maximizer for all convex expectations. Fill every atom outside [5,35]
to its cap, put zero mass strictly inside this interval, and determine
the masses at 5 and 35 from the mass and mean equations. Call this law nu_E.

The exact certificate verifies that both remaining masses lie within
their caps. They are approximately 0.1400234 and 0.3062666.
For any convex g, draw the affine line through (5,g(5)) and (35,g(35)).
The function g lies below that line inside [5,35] and above it outside.
Saturating the outside caps and assigning no interior mass therefore
maximizes the expectation. This proves the claimed convex comparison.

`certify_overlap_extremal_law.py` checks all 81 stop-loss functions exactly.
The variance of nu_E is approximately 225.58395272152455. It improves the
previous comparison variance of about 350.137, but is still too broad.

## Positivity constraints from the constant-weight space

Pair-distance distributions of constant-weight sets satisfy more than
nonnegativity and mass constraints. For 0<=j<=80 and 0<=i<=80, define
the normalized Hahn polynomial

    Q_j(i) := sum_{t=0}^{min(j,i)} (-1)^t
                binom(j,t) binom(257-j,t) binom(i,t)
                / (binom(80,t) binom(176,t)).

The Johnson distance between U and V is i=80-a, half their Hamming distance.
The associated projection kernels give

    E Q_j(80-a) >= 0.

This is the standard positivity condition for the Johnson scheme.
The normalization above divides the usual kernel by its positive multiplicity.
See [Chailloux and Debris-Alazard, Proposition 8 and Section 4.2](https://arxiv.org/html/2405.07666v2).

`johnson_overlap.py` evaluates the polynomials as exact rational numbers.
`test_johnson_overlap.py` verifies normalized orthogonality and scaled
idempotence of complete small Johnson kernels. Symmetry and positive scaled
idempotence prove that those test matrices are positive semidefinite.
The general positivity theorem, rather than the small tests, supplies
the premise at length 256.

## Exact certificates from numerical dual proposals

Fix a function g on the allowed overlap atoms. Choose rational lambda_0,
lambda_1, and nonnegative gamma_j. Let s_j be a fixed positive normalization
for the j-th Hahn constraint. Define

    L(a) := lambda_0 + lambda_1 a - sum_j gamma_j Q_j(80-a)/s_j,
    h_a := max(0, g(a)-L(a)).

Pointwise, g(a)<=L(a)+h_a. The mean equation, Hahn positivity, and atom
caps then give the exact upper bound

    E g(a) <= lambda_0 + 25 lambda_1 + sum_a u_a h_a.

The optimizer proposes the coefficients only. The checker computes every
hinge and the final upper bound with rational arithmetic. Thus the bound
does not depend on optimizer correctness or numerical feasibility tolerances.

This distinction mattered in this calculation. One unscaled numerical LP
reported variance about 47.75. Rechecking its dual against the exact matrix
gave an upper bound about 5375.74, which is unusable. Zeroing matrix entries
of magnitude at most 1e-9 reproduces the optimistic numerical dual value;
the large multipliers amplify the omitted entries. That failed proposal
is retained in `johnson_overlap_variance_attempt01.json`.

Scaling variables by their atom caps substantially improves the proposals.
The degree-20 proposal in `johnson_scaled_variance.json` gives an exact
variance bound of approximately 43.1966. Higher-degree numerical optima
were lower but their checked bounds were worse. The checker does not
promote those numerical optima to theorems.

## From stop-loss bounds to a single comparison law

For each integer t in [0,80], define the actual stop-loss expectation

    C(t) := E max(a-t,0).

The function C is convex, nonincreasing, and has slopes in [-1,0].
Its endpoints are C(0)=25 and C(80)=0. These properties follow directly
by averaging the piecewise-linear functions max(a-t,0).

For t=1,...,79, `certify_johnson_convex_law.py` proposes Hahn duals through
degrees 10,16,20,24. It checks them exactly and takes the smallest valid
upper bound, also allowing the bound from nu_E. Denote this bound by U_t.
Set U_0=25 and U_80=0.

The numbers U_t need not form a convex sequence. Let H be their lower
convex hull, interpolated linearly between its vertices. Despite the word
"lower," H remains an upper bound on the actual C. Between two hull
vertices, convexity puts C below the chord of its endpoint values, which
is below the chord of the corresponding U values. That chord is H.

Write Delta_t=H(t+1)-H(t). Define nu_J on {0,...,80} by

    nu_J(0)  := 1+Delta_0,
    nu_J(a)  := Delta_a-Delta_{a-1}       for 1<=a<=79,
    nu_J(80) := -Delta_79.

The exact certificate checks -1<=Delta_t<=0 and nondecreasing slopes.
It also checks nonnegative masses, total mass one, and mean 25.
The resulting stop-loss expectations are exactly H(t), hence dominate C(t).

Every convex function on the integer interval is an affine function plus
a nonnegative combination of these stop-loss functions. Equal means and
stop-loss domination therefore prove the comparison for every convex g.
The exact replay checks all 316 proposed duals, constructs the hull, and
checks all 81 stop-loss identities. No optimizer runs during replay.

Some diagnostic bounds illustrate the change:

| t | E_nu_E (a-t)+ | New certified upper before the hull |
|---|---:|---:|
| 10 | 16.992535 | 15.0000974 |
| 20 | 10.401500 | 5.7071573 |
| 30 | 3.810465 | 0.7974840 |
| 40 | 0.00080605 | 0.00043476 |

The stored rational values, not these rounded diagnostics, are authoritative.

## Consequence for the pending second moment

The function f(A) in `CONVEX_OVERLAP_REDUCTION.md` is convex in the total
core overlap A. Conditional on support intersection r, the actual row
overlaps are independent. We can therefore replace each by an independent
nu_J sample, one at a time, when bounding E f(A).

The remaining outer average may retain the exact hypergeometric law of r.
Its generating function is

    sum_r [binom(2610,r) binom(5579,2610-r) / binom(8189,2610)]
          * (sum_a nu_J(a) z^a)^r.

No occupancy tail or diagonal pair is removed by this replacement.
The current numerical global envelope is still too loose: selected bounds
on summands in the final sum exceed the required total threshold.
These are failed upper-bound proposals, not lower bounds on actual failure.
The next obligation is a sharper region envelope across every overlap.

There is also a valid way to choose a different region envelope for each
total overlap A. Conditional on the actual row overlaps a_1,...,a_r,
choose any positive log-convex envelope h_A that covers every retained
region type. The subset-to-multinomial coupling still gives an upper
bound f_{h_A}(A). However, the function A |-> f_{h_A}(A) need not be
convex, so it cannot be averaged by replacing the law of A directly.

Instead, bound each actual probability Pr[A=t] through its generating
function. For every z>0, the function a |-> z^a is convex. The established
comparison therefore bounds the actual generating function by the
displayed hypergeometric expression. Nonnegative coefficient bounds give

    Pr[A=t] <= z^(-t) sum_r p_r (sum_a nu_J(a) z^a)^r.

Multiplying this bound by f_{h_t}(t) and summing over every possible t is
valid. Both z and h_t may depend on t. This avoids requiring one envelope
optimized near the usual total overlap to remain sharp at every extreme.
No evaluation of this stronger all-total-overlap bound is yet certified.
