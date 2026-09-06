# A finite one-dimensional reduction for the remaining second moment

The corrected witness in `INVERSE_WEIGHT_AND_PARITY_CORRECTORS.md` removes
the global-parity penalty. Its unresolved second moment still involves
256 correlated region overlaps. This note reduces that expectation to
finite sums, provided a suitable bound on the region ratio is supplied.

The reduction is proved below. The required global bound on the region
ratio has not yet been certified, and the final finite sum is unevaluated.

## The bound that the reduction needs

Retain the notation R(x,y,k)=beta2(x,y,k)/(beta_x beta_y).
Let h be a positive sequence indexed by 0,...,208800, satisfying

    h(k+1)^2 <= h(k) h(k+2).

Thus h is log-convex; it need not be monotone. Require that

    R(x,y,k+d) <= h(k)

whenever x,y are even in [740,900], 0<=k<=897, 0<=d<=3,
and (x,y,k+d) is a feasible actual region type. Values of h beyond 897
are an extension used in the upper-bound argument, not extra region types.
The extension must remain positive and log-convex.

Here k denotes the overlap from core rows only. The correction rows
increase overlap by some d in {0,1,2,3}. The displayed bound is uniform
over their choices, including their dependence on core parity.

Let K_j be the core overlap in region j. Dropping retention indicators
from the normalized second moment gives

    E F^2/M^2 <= E product_{j=1}^{256} h(K_j).

This step does not require independence of the K_j.

## From uniform subsets to multinomial counts

Condition on the two core supports and their row words. Suppose their
supports intersect in r rows. For each common row i, let a_i be the
overlap of its two weight-80 words before permutation. Its contribution
to (K_1,...,K_256) is the indicator vector of a uniform a_i-subset.
These vectors are independent over common rows. Set A=sum_i a_i.

For a fixed integer a<=n, compare the indicator of a uniform a-subset
of {1,...,n} with the counts from a independent draws from that set.
The latter counts dominate the former indicator for every convex function
of the whole vector. The following coupling proves this claim.

First draw a times with replacement and record the distinct labels.
Add uniformly selected unused labels until the resulting set S has size a.
Permutation symmetry makes S a uniform a-subset. Conditional on S,
the draw counts vanish outside S and are exchangeable within S.
Their sum is a, so their conditional mean is exactly the indicator of S.
Conditional Jensen proves the claimed convex comparison.

Apply independent copies of this coupling to the common rows. Interpolate
log h linearly between integers. Log-convexity makes this interpolation
convex, and hence

    (z_1,...,z_n) |-> exp(sum_j log h(z_j))

is convex on the nonnegative box. Conditional Jensen now yields

    E[product_j h(K_j) | a_1,...,a_r] <= f(A),

where, for n=256,

    f(A) := E[product_j h(Z_j)],
    (Z_1,...,Z_n) ~ Multinomial(A;1/n,...,1/n).

The comparison uses convexity of a joint function. It does not assert
independence of the actual region counts or monotonicity of h.

## The multinomial expectation is convex in its total

The sequence f(A) is discretely convex. To prove this, fix a count vector z
with total A and put

    H(z) := product_i h(z_i),
    t_i := h(z_i+1)/h(z_i),
    u_i := h(z_i+2)/h(z_i+1).

Log-convexity gives u_i>=t_i. Add one or two independent uniform balls.
The conditional second difference of the expectations of H is

    H(z) * [((sum_i t_i)/n - 1)^2
              + (1/n^2) sum_i t_i (u_i-t_i)] >= 0.

Averaging over z proves f(A+2)-2f(A+1)+f(A)>=0. No sign condition on
the first difference of f is needed.

## An exact convex upper law for one row overlap

For independent uniform U,V in T80, write a=|supp(U) intersect supp(V)|.
Translation symmetry of the fixed BCH subcode makes each coordinate of
a uniform T80 word equal to one with probability 5/16. Thus E a=25.

For completeness, translations preserve this particular subcode, not just
the containing BCH code. Represent coordinates by F256. The containing
extended BCH code has vanishing power sums of degrees 0,...,36.
Expanding (x+b)^37 shows that translation leaves p37 unchanged.
It therefore preserves the restriction p37 in {0,...,31} and is transitive
on the 256 coordinates.

If U is not V, minimum distance 38 gives a<=61. Equality U=V gives a=80.
For a=40,...,61, let c_a be the retained upper bound on Pr[a], obtained
from difference weight d=160-2a:

    c_a := cap(A_d(C)) * cap(|T80 intersect (T80+dword)|) / lower(|T80|)^2.

The intersection cap is uniform over every C word dword of weight d.
Also put c_80=1/lower(|T80|). These are the exact terms used by
`certify_weight80_family.py` to bound Pr[a>=40]. Their sum is below
0.002564. Define a law nu by

    nu(a) = c_a                         for a=40,...,61 and a=80,
    nu(39) = (25 - sum_{a>39} a c_a)/39,
    nu(0) = 1 - nu(39) - sum_{a>39} c_a.

The exact certificate checks that every mass is nonnegative, their sum
is one, and the mean is 25. Its variance is approximately 350.137.

For every convex function g, E g(a)<=E_{b~nu} g(b). Indeed, replace each
value a<=39 by a mixture at 0 and 39 with the same mean. Convexity can
only increase E g. Then increase the mass at each a>39 to c_a, adjusting
the masses at 0 and 39 to preserve total mass and mean. The increase
per unit added mass is

    g(a) - (a/39)g(39) + (a/39-1)g(0) >= 0.

The last inequality follows from the increasing slopes of a convex function.
All intermediate endpoint masses are nonnegative: mass at 39 decreases
to its certified nonnegative final value, and mass at 0 increases.
This proves convex order, not stochastic order.

## The remaining finite sum

Conditional on support intersection r, the actual a_i are independent
copies of a. Discrete convexity of f permits replacement of these a_i,
one at a time, by independent samples from nu. If nu^{*r}(A) denotes
their sum distribution, then

    E F^2/M^2
      <= sum_{r=0}^{2610} p_r sum_{A=0}^{80r} nu^{*r}(A) f(A),

    p_r := binom(2610,r) binom(5579,2610-r) / binom(8189,2610).

This is the exact support-intersection distribution for two independently
sampled 2610-subsets of 8189 positions. It includes coincident supports
and all diagonal messages through the overlap law's atom at 80.

For any A and positive real z, nonnegative coefficient extraction gives

    f(A) = A!/256^A [z^A] (sum_{k=0}^A h(k) z^k/k!)^256
         <= A!/(256 z)^A (sum_{k=0}^A h(k) z^k/k!)^256.

The finite truncation at A is exact for the coefficient. No infinite
series or unproved Poisson tail estimate is needed. Numerical optimization
may propose z; a final certificate must verify the chosen bound outward.

The reduction would suffice if its finite sum were below (9/16)*2^40.
The central atlas supplies only part of h's required domain. Constructing
a global log-convex h, and evaluating the entire sum without losing too
much, are the next two obligations.

An optional further bound removes the hypergeometric sum. Let p=2610/8189
and let W equal zero with probability 1-p, or an independent nu sample
with probability p. Then the preceding bound is at most E f(sum_i W_i)
for 2610 independent copies W_i. To justify this, first observe that
r |-> E f(sum_{i=1}^r b_i), with b_i independent from nu, is convex.
For nonnegative integer a,b, convexity of f gives

    f(s+a+b)-f(s+a)-f(s+b)+f(s) >= 0.

Averaging over a,b proves the assertion. The same sampling coupling
used above compares the hypergeometric intersection size with
Bin(2610,p) for every convex function, by drawing from a population
of 2610 ones and 5579 zeros. This proves the optional bound.

In particular, put P(z)=1-p+p sum_a nu(a) z^a. The optional bound equals

    sum_{A=0}^{208800} [z^A] P(z)^2610 * f(A).

For any positive t, each nonnegative coefficient is at most
P(t)^2610/t^A. A separate tilt may be used for each A. Combining this
coefficient bound with the bound for f(A) gives a fully finite route
using only one-variable evaluations and positive summation. The numerical
loss from both coefficient bounds must be assessed, not assumed small.

`certify_overlap_convex_law.py --verify` replays the exact law certificate.
`test_overlap_convex_reduction.py` checks subset comparison, convexity,
endpoint replacement, and coefficient identities in exhaustive toy cases.
One test uses U-shaped h specifically to exclude a monotonicity shortcut.
