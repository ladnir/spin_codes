# Scalable bit-route SPIN with the Turbo inner: proof draft

This document assembles the distance argument before changing the
manuscript. It uses the paper's one-bit regional route and shared growing
Golay--BAA outer. The inner is the actual Turbo `t=64,s=22` expansion,
transpose feedback, and independent uniform `GL(22,2)` update at every
step. All setup is sampled once and reused across messages. Probability
below is over that one setup.

## Claim and imported outer facts

For `L=128m`, choose the least multiple of 24 satisfying
`b >= (39/4) log2(Lb)` and put `N=Lb`. Let `E_m` be the resulting
rate-half encoder. The proposed theorem is

```
Pr[d_min(E_m) <= floor(0.11 N)] = o(1)  as m -> infinity.
```

The outer theorem already used by Scalable SPIN supplies an event `G_b`
depending only on the two BAA permutations. It has probability `1-o(1)`.
On `G_b`, every nonzero outer word has relative weight in
`W=[0.104,0.896]`, and its realized weight multiplicities obey the
certified concave majorant `a_hat_BAA`, with uniform `O(log b)` remainder.
In particular, the ratio of its weight-`w` multiplicity to `C(b,w)` is
at most

```
exp[b(-ln(2)/2 + ell + epsilon_b)],
ell=0.01281,  epsilon_b=O(log b/b),
```

uniformly on `G_b`. The one outer is reused at all `L` positions, so
this event is paid once. Condition on an arbitrary realized outer in
`G_b` for the route and inner bounds below.

## Exact state update and a weighted envelope

Write `M=2^22-1`. For a fixed nonzero entering state `q`, a uniform
`GL(22,2)` map makes `M_i q` uniform over the `M` nonzero states. It is
independent of the step input `X_i`. The step emits
`Y_i=X_i+Aq` and leaves state `M_i q+A^T X_i`. The expansion `A` is
injective; every column of `A^T` is nonzero. Its image spectrum and the
weight enumerator `K_j` of `ker(A^T)` are reconstructed by
`verify_pair_return_map.py` and `shared_inner_dense_screen.py`.

Take iid `Bernoulli(beta)` input bits for one step and weight an output
by `z^wt(Y_i)`, where `0<z<=1`. Define `g0,g1,f_w,F,U,H,u` and the
nonnegative matrix `T(beta,z)` exactly as in
[SHARED_INNER_FEASIBILITY.md](SHARED_INNER_FEASIBILITY.md).
Here is the domination argument with the measure being bounded made
explicit. An unnormalized state measure `mu` has envelope `(z0,f,u0)`
when `mu(0)<=z0` and there exists a nonnegative measure `nu` on the
nonzero states with `mu(q)<=nu(q)+u0/M` and
`sum_{q!=0}nu(q)<=f`. The zero-state part of one weighted step has
mass `H` at zero and at most `g0^64-H` arbitrary live mass. For any
nonzero entering `q` and any target live state `q'`, the weighted
transition to `q'` is at most `f_{wt(Aq)}/M<=F/M`; its weighted
transition to zero is at most `min(f_{wt(Aq)},u)/M`. Summing these
inequalities against `nu` gives the second row of `T`.

For the uniform envelope, each entering live state has mass at most
`u0/M`. Summing `f_{wt(Aq)}` over those states gives `u0 U`; each target
live state receives at most `u0 U/M`, and zero receives at most
`u0 min(U,u)/M`. This gives the third row. The inequalities are
entrywise and remain true after arbitrary earlier steps; induction gives

```
E[z^wt(I(X))] <= e_0 T(beta,z)^(N/64) 1
```

for iid Bernoulli input `X`. If a fixed positive column `v` satisfies
`T v<=rho v`, the right side is at most
`(v_0/min_j v_j) rho^(N/64)`. The vector and prefactor may depend on
the chosen finite witness, but not on `N`. No assumption that
`A^T X_i` is uniform enters this argument.
`test_shared_inner_transfer.py` independently enumerates a small
`GL(2,2)` instance and checks each envelope row with exact rational
arithmetic; it is a cross-check of the transition orientation, not a
proof of the general statement.

## Positive occupation

Let `Q` outer positions be active, `alpha=Q/L`, and let `x` be the mean
relative weight of their realized outer words. The existing route
domination lemma is independent of the inner law. Conditional on the
row weights, it bounds every nonnegative function of the routed word
by `(b+1)^Q(L+1)^b` times its iid Bernoulli-`alpha*x` expectation.
Changing the Bernoulli density to `beta=p*y` costs at most
`exp[N(D(alpha||p)+alpha D(x||y))]`. Concavity of `a_hat_BAA` bounds the
outer multiplicities at mean `x`. For `alpha>=10^-4`, all polynomial,
position-count, and route-comparison factors are `exp(o(N))`
uniformly. Therefore the exponent per output bit is at most

```
alpha*a_hat_BAA(x) + D(alpha||p) + alpha*D(x||y)
       + log(rho(T(p*y,exp(-lambda))))/64 + 0.11*lambda.
```

For fixed `p,y,lambda` and one affine outer segment, this function is
convex in `(alpha,alpha*x)`: its outer term is affine, `D(alpha||p)` is
convex, and `alpha*D((alpha*x)/alpha||y)` is a convex perspective.
The rectangle in `(alpha,x)` maps to a quadrilateral in these variables,
so its four vertices suffice. The exact-rational partition and 256-bit
Arb checks in `shared_inner_dense_certificate.py` cover
`[10^-4,1] x W` in 46,317 accepted boxes, all with outward upper
exponent below `-1.00147e-7`. A 512-bit replay gives the same number of
leaves and passes. The certificate computes a positive rational
Collatz vector for each witness's three-state matrix; floating-point
optimization only proposes witnesses. The finite witness set has a
finite maximum moment prefactor. Thus the positive-occupation sum is
`L exp(-eta N+o(N))=o(1)` for any fixed
`0<eta<1.00147e-7`.

## Growing sparse occupation

For `4096<=Q<=10^-4 L`, use the paper's fair-row comparison instead of
the `(L+1)^b` route factor. Each of the `Q` fair active rows contributes
one uniformly located mark per region, carrying an independent fair bit.
Generate independent marks with probability `p=(8/5)alpha`, condition
each region's mark count to equal `Q`, and then give each mark a fair bit.
Before conditioning, every input bit is iid Bernoulli
`beta=(4/5)alpha`. The inverse mark-count probability is at most
`8 sqrt(Q) exp[L D(alpha||p)]` per region, as in the existing proof.

Set `z=1-(8/5)alpha`. The exact-rational check in
`shared_inner_dense_screen.py --sparse-gate` proves that for
`0<alpha<=10^-4`,

```
T(beta,z) v(alpha) <= (1-48 alpha) v(alpha),
v(alpha)=(1,(1+24 alpha)/1024,1/1024).
```

The check uses only the kernel's distance-eight property, the
expansion's minimum live weight 16, and its mean live weight at least
32. Here is its analytic enclosure. With `c=4/5`, `d=8/5`, and
`beta=c alpha`,

```
T00 <= (1-beta)^64 + C(64,8) beta^8,
T0F <= 64 beta,
TF0,TU0 <= 64 beta/M,
TFU <= (1-d alpha+cd alpha^2)^16,
TUU <= 1-32d alpha + (64cd+C(64,2)d^2)alpha^2.
```

The first line follows from the kernel minimum weight and the eighth
factorial moment of a binomial variable. The last two follow from
`g0<=1`, `g1=1-d alpha+cd alpha^2`, the minimum/mean expansion
weights, and the second Bonferroni bound for `(1-u)^w`. Substituting
the displayed `v(alpha)` and dividing each residual by `alpha` yields
the exact rational upper bounds `-3.020856`, `-1.439532`, and
`-2.663212` at `alpha=10^-4`; each bound is monotone in `alpha` on
the interval. Hence the moment prefactor is at most 1024.

The outward logarithmic check in
`shared_inner_remaining_certificate.py` combines that contraction
with the selected outer and the schedule. Its leading coefficient per
`Qb` is below `-0.0134042`. The remaining terms are no larger than
in the current IMT proof: `log(8 sqrt(Q))/Q<0.002` for `Q>=4096`,
the other uniform remainders are eventually below `0.004`, and the
new moment prefactor 1024 is smaller than the former 2048. Thus for
all sufficiently large native lengths,

```
E[Z_{floor(0.11N),Q} | G_b] <= exp(-0.006 Qb)
```

throughout this range, and the geometric sum tends to zero.

## Every fixed occupation

Fix `Q<4096` independently of `L`. Each region receives at most `Q`
one-bit impulses. Conditional on their originating rows and on their
number `a`, the independent regional shuffle places them at a uniform
`a`-subset of the `L` bit slots. Put `n=L/64` steps per region and
`z_L=exp(-theta/L)` for fixed `theta>0`.

For an empty step starting from a live state `q`, the successor is
exactly uniform nonzero, independent of its emitted weight. Let
`phi_L=(1/M)sum_{q!=0} z_L^wt(Aq)` and
`p0=2^21/M`. An empty gap of `g>=1` steps has exact live-to-live
kernel

```
z_L^wt(Aq) phi_L^(g-1) / M
```

for every entering live `q` and exiting live `q'`. Its zero-to-zero
entry is one and both cross-class entries are zero. Since
`log(phi_L)=-theta*64*p0/L+O(theta^2/L^2)`, uniformly for
`1<=g<=n` this kernel differs multiplicatively by `1+O(theta/L)`
from the live component of
`D_gamma(64g/L)`, where `gamma=theta*p0` and
`D_gamma(u)=diag(1,exp(-gamma*u))`.

For a fixed `a<=Q`, the probability that two impulses fall in one or
adjacent steps, or that the first or last impulse touches a region
boundary, is `O_Q(1/L)`. On the complementary event, every gap has at
least one empty step, so it refreshes a live state before the next
impulse or the region boundary. An impulse's emitted-weight factor
lies in `[exp(-64 theta/L),1]`. For the unweighted impulse, a zero
state becomes live because every feedback column is nonzero; a
uniform live state returns to zero with probability `r=1/M`. Thus
the limiting zero/live impulse matrix is

```
P0 = [[0,1],[r,1-r]].
```

Ordered uniform impulse positions converge by a Riemann sum to the
simplex of spacings. For each fixed `a<=Q`, the full-state regional
weighted kernel therefore converges entrywise, uniformly over all
entering and exiting states, to the lift of

```
K_a(theta) = a! integral_{u_i>=0, sum_i u_i=1}
  D_gamma(u_0) P0 D_gamma(u_1) ... P0 D_gamma(u_a) du.
```

The error is `O_{Q,theta,M}(1/L)`: bad placements have that
probability and weighted kernels have row sums at most one; on good
placements, each finite gap and impulse differs by `1+O(1/L)`, and
the integrand has bounded derivatives. This is a *full-state* bound,
so it handles the state carried across regions. Entries that vanish
in the limit vanish identically in the finite kernel: for `a=0` the
class cannot change, and for `a=1` zero cannot return to zero. Every
other coarse entry is positive at fixed `Q,theta`. Dividing the
additive error by the positive limiting entries gives an entrywise
comparison `F_{a,L} <= (1+C_{Q,theta,M}/L) J K_a R` for all `a<=Q`,
where `J` classifies the entering state and `R` lifts a live exiting
class uniformly over its `M` states. Across `b=O(log L)` regions the
loss is `exp(O_Q(b/L))=1+o(1)`.

Replacing the live-to-live entry of `P0` by one gives the same
conservative `P_plus` used in the current fixed-occupation proof.
`shared_inner_remaining_certificate.py` recomputes the `Q=1,2`
continuum exponents with Turbo's `r,p0`, obtaining outward upper
bounds `-0.109827` and `-0.109812` per active-row bit. For every
fixed `Q>=3`, it directly recomputes all 62 fugacity endpoints with
the new parameters; the largest upper bound is below `-0.00086852`
per active-row bit. The old coefficient extraction and weighted
product-norm argument depend only on `P_plus,p0` and the outer
majorant, so they apply to these updated values. Consequently each
fixed `Q<4096` has a bound `exp(-c_Q Qb+o_Q(b))` with `c_Q>0`.

## First-moment assembly and status

Let `Z_{d,Q}` count nonzero messages with `Q` active outer positions
and output weight at most `d=floor(0.11N)`. For an arbitrary outer in
`G_b`, the preceding route and inner bounds are uniform. The finite
sum over `1<=Q<4096` tends to zero. The range
`4096<=Q<=10^-4L` has a uniform geometric bound. The remaining
range contributes at most `L exp(-eta N+o(N))`. Hence

```
sum_{Q=1}^L E[Z_{d,Q} | G_b] = o(1).
```

Markov's inequality and `Pr[G_b^c]=o(1)` give the proposed distance
claim. The outer is injective, the route is a permutation, and the
inner recurrence is bijective for every fixed schedule; thus every
realized encoder has rate one half. The shared outer uses `O(N)` work
at this schedule, and the fixed-width route and inner use `O(N)`
ordinary and transposed work.

This is a **proof draft**, not a replacement for the current theorem.
The new analytical steps are the three-state domination and the
fixed-region comparison. They have been derived above and their
numerical gates are outward checked, but they still need independent
mathematical review and a formal replay or Lean integration before the
manuscript should present this as its proved scalable instantiation.
