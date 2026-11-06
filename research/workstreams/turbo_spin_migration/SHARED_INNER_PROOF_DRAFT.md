# Scalable bit-route SPIN with the Turbo inner: conventional proof

This document assembles the distance argument before changing the
manuscript. A LaTeX version of the new full-state lemma and theorem
assembly is in [shared_inner_paper_proof.tex](shared_inner_paper_proof.tex).
The argument uses the paper's one-bit regional route and shared growing
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

Fix `Q<4096` independently of `L`. Write `t=64`, `n=L/t`, and
`z_L=exp(-theta/L)` for a fixed `theta>0`. Conditional on which `a<=Q`
active rows send a bit to a region, its independent shuffle puts their
distinguishable bits at uniformly chosen distinct slots among the `L`
slots. Let `F_{a,L}(z_L)` be the resulting `(M+1)`-state weighted region
kernel, with the averaging over those slots and all step mixers included.
It does not average the entering state: that state is whatever the
previous region left behind.

**Lemma (uniform full-state region comparison).** Define
`p0=2^21/M`, `gamma=theta*p0`,
`D_gamma(u)=diag(1,exp(-gamma*u))`, and
`P0=[[0,1],[1/M,1-1/M]]`. Let `J(q,c)` indicate whether state `q` is
zero or live, and let `R(c,q)` put the zero class at zero and distribute
the live class uniformly over the `M` nonzero states. Thus `RJ=I_2`.
Set

```
K_a(theta) = a! integral_{u_i>=0, sum_i u_i=1}
  D_gamma(u_0) P0 D_gamma(u_1) ... P0 D_gamma(u_a) du.
```

For every fixed `Q` and `theta`, there is a finite constant
`C_{Q,theta,M}` such that, for all sufficiently large multiples `L` of
64, every `0<=a<=Q`, and all entering and exiting states,

```
F_{a,L}(z_L) <= (1+C_{Q,theta,M}/L) J K_a(theta) R
```

entrywise. The constant is independent of the active-row identities,
the region index, and the entering state.

*Proof.* Let `B_g=W_0(z_L)^g` be the empty-step kernel. An empty step
leaves zero at zero. From a live state `q`, the uniform `GL(22,2)` mixer
sends the next state uniformly to a nonzero state, independently of the
current output `Aq`. Consequently, for every `g>=1` and live `q,q'`,

```
B_g(q,q') = z_L^wt(Aq) phi_L^(g-1)/M,
phi_L = (1/M) sum_{v!=0} z_L^wt(Av).
```

The other entries of `B_g` are zero except `B_g(0,0)=1`. The uniform
live mean weight is `t*p0`. Hoeffding's lemma, applied to a random
variable in `[0,t]`, gives

```
0 <= log(phi_L) + theta*t*p0/L <= theta^2*t^2/(8L^2).
```

Since `g<=n=L/t` and `0<=wt(Aq)<=t`, put
`c_theta=theta*t+theta^2*t/8`. Then `B_g` is bounded above and below,
entrywise on its nonzero entries, by `exp(+-c_theta/L)` times
`J D_gamma(t*g/L) R`. This also covers `g=1` and arbitrary entering
live states; no stationary-state premise is needed.

Let `V_j(z_L)` be the full-state kernel of a step with one bit in local
coordinate `j`. Its output has weight at most `t`, so, entrywise,

```
exp(-theta*t/L) V_j(1) <= V_j(z_L) <= V_j(1).
```

Every column `A^T e_j` is nonzero. A zero entering state therefore
becomes live under `V_j(1)`. A uniform live entering state becomes zero
with probability exactly `1/M`: the uniform mixer maps each nonzero
state uniformly over the `M` nonzero states, and exactly one such state
cancels `A^T e_j`. Therefore `R V_j(1) J=P0` for **every** `j`.

Call a placement good when its impulse steps are pairwise separated by
at least one empty step and neither the first nor the last step is an
impulse step. A union bound gives, for `L>1`,

```
Pr[bad] <= 2*a*t/L + 3*t*binom(a,2)/(L-1).
```

Indeed, an impulse hits either boundary step with probability `2t/L`;
conditional on one slot, at most `3t` slots are in its step or an
adjacent step. On a good placement there are `a+1` empty gaps, each
of length at least one. Replacing each gap by its full-state lifted
matrix, each impulse by `V_j(1)`, and using nonnegative matrix
multiplication gives a componentwise two-sided factor
`exp(+-[(a+1)c_theta+a*theta*t]/L)` around

```
J D_gamma(t*g_0/L) P0 D_gamma(t*g_1/L)
    ... P0 D_gamma(t*g_a/L) R.
```

This factorization is valid for arbitrary entering and exiting states.
In particular it retains the state carried between regions. The slot
coordinates within an impulse step have disappeared because
`R V_j(1) J=P0` holds for every coordinate.

To compare the discrete gaps to the integral, take `a` independent
uniform points on `[0,1]`, sort them, and quantize them to the `L`
slots. The probability of a repeated slot is at most
`binom(a,2)/L`. Conditional on no repeat, the ordered slots have
exactly the required uniform-without-replacement law. Sorting and
quantization change every ordered normalized position by at most
`1/L`; converting positions to spacings changes their total absolute
error by at most `2(a+1)/L`. Replacing each impulse step by its slot
position changes the spacing vector by at most `2t(a+1)/L` more.
The integrand above is Lipschitz in spacings with entrywise difference
at most `gamma` times their total absolute difference: each `D_gamma`
has derivative at most `gamma`, while `P0` is stochastic. The bad
placements, including repeated quantized slots, contribute at most
their `O_Q(1/L)` probability to any entry, because every weighted
full-state kernel is substochastic. It follows that

```
max_{a<=Q,q,q'} |F_{a,L}(q,q')-(J K_a(theta) R)(q,q')|
    <= C'_{Q,theta,t}/L.
```

For `a=0`, both cross-class entries are identically zero in both
kernels. For `a=1`, the zero-to-zero entry is identically zero:
one impulse activates the state and all subsequent empty mixers are
invertible. The other coarse entries of `K_1` are positive. For
`a>=2`, all four coarse entries are positive since `1/M>0`,
`1-1/M>0`, and all diagonal gap factors are positive. Thus the
minimum positive entry of the finite collection
`{J K_a(theta)R : 0<=a<=Q}` is positive. Dividing the additive bound
by this minimum proves the claimed componentwise multiplicative bound;
the exactly zero entries require no division. This proves the lemma.

The regional shuffles and mixers are independent between regions.
Conditional on each region's contributing rows, apply the lemma to
each full-state kernel and multiply the entrywise inequalities. The
loss is `(1+C_{Q,theta,M}/L)^b=exp(O_Q(b/L))=1+o(1)` because
`b=O(log L)`. It remains valid when coefficient extraction later
fixes each row's total weight: the inequality holds for each
contribution pattern before those nonnegative terms are summed.

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

The conventional mathematical argument closes the 11% asymptotic
claim for this setup law, using the already proved outer-selection and
bit-route lemmas. The new analytical steps are the three-state
domination and the full-state fixed-region comparison above. The
numerical gates pass with outward arithmetic at both 256 and 512 bits.
This proof has not been translated into Lean. The revised manuscript
states the uniform-mixer Scalable theorem and gives its full argument in
`research/paper/structured_proof.tex` and
`research/paper/structured_uniform_appendix.tex`. The earlier IMT proof
remains archived outside the compiled paper.
