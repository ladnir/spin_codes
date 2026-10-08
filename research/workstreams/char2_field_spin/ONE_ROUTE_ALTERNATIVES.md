# One-route alternatives: concrete proof screens

The target is one permutation pass, characteristic-two field operations,
and a jointly credible proof and implementation. The certified two-scan
control with two independent permutations is not a candidate winner here.
No encoding benchmark was run for this note.

## 1. A uniform global route needs a substantially larger outer

This restriction is independent of the inner arithmetic. Let a causal
inner start in zero state. If an outer row has at least b-k+1 of its
b coordinates in the last e positions, the [b,k] MDS shortening to
those coordinates is nonzero. A nonzero source word supported there
produces output only in that suffix, so the sampled code fails distance
e+1 for every choice of the causal inner's coefficients.

For one fixed row, the exact routing probability is

    p = sum_(j=b-k+1)^b binom(e,j) binom(N-e,b-j) / binom(N,b).

At K=65536, N=131072, e=13107, the following are **expected bad-row
counts Lp**, not probabilities obtained by assuming row independence:

| Outer | Lp from the causal suffix event |
|---|---:|
| GRS[16,8] | 0.048426 |
| GRS[32,16] | 5.211e-6 |
| GRS[64,32] | 1.499e-13, about 2^-42.60 |
| GRS[128,64] | 3.153e-28, about 2^-91.36 |

Thus replacing the regional route by a uniform global route does not
repair the existing small-outer construction. GRS[64,32] leaves little
room relative to a 40-bit target even before other failures are counted.
GRS[128,64] avoids this particular boundary obstacle, but its extra
outer arithmetic must earn its cost.

## 2. One weighted accumulator: simple, but the current bound does not close

After one uniform permutation, compute y_i=c_i x_i+y_(i-1), where
c_i is independently uniform in a nonzero coefficient set of size Q.
This is one coefficient product and one XOR per symbol, with no
second permutation. Its transpose is a suffix sum followed by local
coefficient products. The companion's generic return-rank argument
applies with g=1/Q+(1-1/Q)/tau.

For a single scan, let H_g(s,e) count support paths with s active
inputs, output weight at most e, and one factor g per return. Its
two-state generating matrix is

    [[1, xz], [gx, (1+x)z]].

For a repeated [b,k] MDS outer define the shortened-dimension polynomial

    Q_outer(u,tau) = 1 + sum_(w=b-k+1)^b binom(b,w) u^w tau^(w-b+k).

The same quotient-rank proof gives the one-permutation bound

    (1/tau) sum_(s=1)^(2e)
      [u^s] Q_outer(u,tau)^L / binom(N,s) * H_g(s,e).

There is only one binomial normalization. The inverse scan limits the
source support to 2e. This is a valid whole-code formula, not a numerical
certificate.

The local screen uses a cheap upper bound on H_g. A path with r closed
runs contributes at most

    binom(s-r-1,r-1) binom(N,r) binom(e,s-r) g^r.

A path ending in an open run contributes at most

    binom(s-r,r-1) binom(N,r-1) binom(e,s-r+1) g^(r-1).

These follow by summing run lengths and replacing binom(N-h,r) by
binom(N,r), or the analogous open-run term. At Q=255 and 10% distance,
a floating-point search over 801 tau values gives one-active-row upper
bounds greater than one for each rate-half outer of length 16, 32, 64,
and 128. This screen does not disprove those codes or exclude sharper
analysis. It does mean that a larger outer plus one ordinary weighted
accumulator is not an immediately certifiable replacement using this
simple bound. Deprioritize it until a tighter run calculation earns a win.

## 3. Two same-order weighted scans without a second permutation

A more substantial alternative is

    u_i = c_i x_i + u_(i-1),
    y_i = d_i u_i + y_(i-1),

with independent coefficient sets and one route applied only before
the first scan. In matrix notation the inner is A D_d A D_c, not the
two-permutation control. Its transpose is D_c A^T D_d A^T. Both scans
are contiguous passes; the products can be placed outside each prefix
or suffix dependency chain. This design has two local coefficient
products and two XOR scans per symbol, but only one permutation pass.

There is a direct four-state support operator. A state records whether
the previous u and y are nonzero. For the first scan, 0 plus input 0
stays 0; 0 plus input 1 becomes 1; 1 plus input 0 stays 1; and 1 plus
input 1 either returns to 0 with factor g or stays 1. Apply the same
rule to the second scan using the newly computed u as its input.
Mark original active inputs by x and final nonzero outputs by z.
The resulting T_g(x,z), with states 00,01,10,11, is constructed exactly
by `cascade_operator` in the accompanying script.

The coefficient of x^s z^h in [1,0,0,0] T_g^N [1,1,1,1]^T bounds
joint support paths. The first scan's return coefficient is exposed
before the second's; required-nonzero forms must be kept in the alive
branch condition. The fresh coefficients justify separate return
charges even though the scans use the same order. For a uniform global
route, substitute this joint coefficient into the one-permutation
formula, with source cutoff 4e. Do not insert an intermediate binomial
denominator: there is no second permutation.

For a regional route, use the four-state regional coefficient matrix
normalized by binom(L,h) in each region and retain state across regions.
The source's joint column support and shortened dimension must remain
in the count. A scalar normalization by binom(N,total weight) is invalid.

There is an immediate warning in this otherwise simple support bound.
At x=0, its active-first-scan submatrix has Perron root

    lambda = (z + sqrt(z^2+4gz))/2.

This exceeds one when z>1/(1+g). Thus the crude rank-return count
accumulates slack during long stretches with no original input. A
successful implementation may need a rank-aware refinement, just as
the current batched inner does. The formula is a concrete fallback to
screen, not a proof or performance win already obtained.

## 4. Circular scans repair a boundary, not the whole small-outer problem

A tail-biting weighted accumulator can enforce a cyclic state instead
of starting at zero. For example y_i=alpha*y_(i-1)+c_i*x_i with
y_(-1)=y_(N-1) is invertible when alpha^N!=1. At power-of-two N,
any alpha other than 0 or 1 in GF256 satisfies that condition, since
GF256* has odd order. One global linear aggregate determines the
initial state; encoding then needs the scan and a contiguous correction
pass, not a second permutation.

This removes the automatic zero initial state for a suffix input.
However, imposing zero initial state is only one extra homogeneous
constraint. A shortened outer subspace of dimension at least two can
still satisfy it. More generally, closing an s-symbol boundary adds
at most s constraints, so a suffix containing b-k+s+1 coordinates of
one MDS row still guarantees a bad word. Small-state circular variants
therefore do not erase the small-outer concentration problem.

## 5. A cheap high-dimensional refresh with an exact projective law

For a larger state in F_q^s, independently sample u uniformly from the
nonzero vectors and then v uniformly subject to v^T u != 1. In
characteristic two,

    T = I + u v^T

is invertible. This costs a dot product and a scaled vector addition,
or 2s scalar products per step, rather than the s^2 products of a dense
matrix. Setup rejection has probability 1/q for each fixed nonzero u.

For a fixed nonzero state x, the exact image probabilities are

    Pr[T x = x] = q^(s-1)/(q^s-1),
    Pr[T x = c x] = 1/(q^s-1),  c != 0,1,
    Pr[T x = y] = (1-1/q)/(q^s-1),  y outside span(x).

Multiplying T by an independent uniform scalar lambda in F_q* yields
the particularly simple distribution

    lambda*T*x ~ (1-1/q)*Uniform(F_q^s minus zero)
                    + (1/q)*Uniform(span(x) minus zero).

This is a distributional identity, not merely an atom bound. A separate
scalar multiplication adds s products per step, for about 3s in total.
At q=256 its projective lazy part has weight 1/256. It is a concrete
candidate for a larger state such as s=8 with t=32, while retaining
GF256-coefficient arithmetic and one regional routing stage.

The useful fixed-vector law does **not** establish a full quotient-rank
return bound. For a matrix H of state forms and a prescribed target H,
the event lambda=1 leaves

    (lambda*T-I) H = u (v^T H),

which has rank at most one, even if H has rank eight. This event has
probability 1/(q-1). It therefore cannot be analyzed as uniform GL(s,q)
refresh merely because its action on one vector is nearly uniform.
The proof must retain the scalar-line branch and shared-form dependence.
No complete distance certificate or performance result is claimed for
this larger-state alternative. The independent GF4 tests verify the
image law, scalar symmetrization, and low-rank return exception exactly.

## Recommendation

Keep the regional route and prioritize the sharper current rank-return
operator. Its new state-direction constraint removes a demonstrated
source of proof loss without adding encoder work. The same-order
two-scan operator is the most concrete fallback with one permutation,
but its extra arithmetic and long-run rank slack need an explicit
screen before implementation. Uniform-global routing is not a safe
shortcut for the existing GRS[16,8] outer.

## Reproduction

`theory/one_route_alternatives.py --causal-suffix` evaluates exact
hypergeometric probabilities before displaying decimal summaries.
`--screen` runs the explicitly non-rigorous one-active-row search.
The script also contains exact rational run bounds and the four-state
cascade operator. No full certificate is claimed for these alternatives.
