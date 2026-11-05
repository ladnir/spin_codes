# An exact positive action on arbitrary population vectors

2026-11-01. The adaptive population operator admits an exact Horner
recurrence on homogeneous polynomials. It applies the operator to an
arbitrary positive population vector using two sparse quartic
multiplications per degree step. It requires neither a dense
population matrix nor a general polynomial convolution or FFT. This
removes the two function-comparison losses of the minimum-exponential
method, while retaining the adaptive Hall majorant and the separate
regional quota comparison.

## The two local polynomials

Fix a marker `y>1`, a positive label evaluation `x`, and one of the
networks from
[QUARTET_GENERAL_FLOW_COVER.md](QUARTET_GENERAL_FLOW_COVER.md).
The old coordinate types are `W`, `E=A outside W`, and `Z=A^c`.
Let `T(X_W,X_E,X_Z)` be the positive homogeneous quartic for choosing
a hole, and `S(X_W,X_E,X_Z)` the corresponding quartic for selecting
the current quartet. Their coefficients include the old-type
multinomial, all current label multiplicities, the marker `x`, the
local cover cost, and respectively the probabilities `g,p_select`.

For example, in the original `(1,1,2)` network, a type with
`w+e+z=4`, `u=w+e` contributes to `T` the coefficient

```
multinomial(4;w,e,z)
 *sum_(j,r) binom(w,r)binom(4-w,j-r) x^j y^(w-r) g(u,j,r),
```

and `S` replaces `g` by `p_select`. For partial capacities use the
local gain and cost in
[QUARTET_PARTIAL_TEMPORAL_FLOW.md](QUARTET_PARTIAL_TEMPORAL_FLOW.md).
The probabilities depend on the old local type and current labels,
but not on the population vector to which the operator will be applied.

Write `n=4m`, with production values `n=512,m=128`. If exactly `k`
quartets are holes, the new hole population is `4k`. A monomial
`X_W^r X_E^s X_Z^(4k-r-s)` within the hole factor `T^k` gives new
intersection population `r+s`: both old `W` and old `E` belong to
old `A`. The selected factor contributes nothing to either new
population.

## Exact arbitrary-vector formula

Let `f(4k,v)` be any nonnegative function on the safe physical states
`0<=v<=4k`, `0<=k<=m`. For a homogeneous degree-`4k` polynomial,
define `D_k` to multiply its monomial coefficient at old exponents
`(r,s,4k-r-s)` by `f(4k,r+s)`. Then

```
G_f=sum_(k=0)^m binom(m,k) D_k(T^k) S^(m-k).         (1)
```

For an old physical state `(ell,omega)`, its exact action is

```
(K_x f)(ell,omega)
 =[X_W^omega X_E^(ell-omega) X_Z^(n-ell)]G_f
   /multinomial(n;omega,ell-omega,n-ell).            (2)
```

The binomial coefficient chooses which of the `m` quartet factors
are holes. The final denominator averages the old labels over the
fresh uniform partition. There is no new-population denominator and
no intermediate division by an old multinomial.

Crucially, `D_k` acts before multiplication by `S^(m-k)`. Weighting
the final total old degree would confuse the new intersection with
the entire old population.

## Horner recurrence

Set

```
P_0=1,       H_0=f(0,0).
For k=1,...,m:
    P_k=T*P_(k-1),
    B_k=binom(m,k) D_k(P_k),
    H_k=S*H_(k-1)+B_k.                              (3)
```

An induction gives

```
H_k=sum_(j=0)^k binom(m,j) D_j(T^j) S^(k-j).
```

Thus `H_m=G_f`. Every polynomial at stage `k` is homogeneous of
degree `4k`. Both multiplications in (3) are by a fixed quartic with
15 monomials, so they use the same positive shifted-array operations
as the existing one-exponential engine. The binomial coefficient is
`binom(m,k)` with the final `m` fixed throughout the recurrence.

After setting `X_Z=1`, the implementation needs triangular arrays of
old exponents `(r,s)` and degree at most `n`. The arithmetic cost is
`O(15*m*n^2)` and memory is `O(n^2)`; there is no dense 33,153-square
transition matrix. Relative to one quartic-power action, the principal
work is a second sparse multiplication and a coefficientwise vector
weight/addition at each stage. Runtime must be measured by the parent;
no timing is asserted here.

The powers `P_k` and accumulators `H_k` need independent common
scales. When forming `S H_(k-1)+B_k`, align those scales before the
positive addition. An implementation must detect loss of a positive
term through underflow instead of clipping or silently dropping it.
The logarithm of `binom(m,k)` and any input-vector scale are included
once in `B_k`. Floating range guards alone do not provide directed
rounding certification.

Two exact identities provide independent checks. For `f=1`, (1)
reduces to `(S+T)^m`, hence `K_x 1` is the known label/cover-cost
product. For `f(ell,omega)=alpha^ell beta^omega`,

```
G_f=[S+alpha^4*T(beta X_W,beta X_E,X_Z)]^m,
```

which must agree with the already checked exponential-action engine.
Small-width physical partition enumeration can additionally test
non-exponential vectors.

## Regional quotas and an optional endpoint shortcut

For the half-row class, let `h=4096`, `C0=binom(8192,4096)`, and
`B=[x^h]K_x^16/C0`. Given a fixed positive physical vector `phi`,
compute all sixteen exact actions at one fixed evaluation `x` before
applying the Cauchy coefficient bound. Then

```
B phi <= F(phi),
F(phi)(s)=min_x x^-h (K_x^16 phi)(s)/C0.             (4)
```

The minimum is pointwise in the common physical state space. There is
no bin maximum, moved function minimum, or exponential re-envelope.
The quota comparison in (4) remains an upper bound, not exact
coefficient extraction.

The exact first-region endpoint is
`b=[x^h]K_x^15 q_x/C0`, with the appropriate original or partial
capacity endpoint `q_x`. Optionally, for the same vector `phi`, set

```
mu_x=max_s q_x(s)/phi(s).
```

Positivity gives the valid shortcut

```
b(s)<=b_upper(s)
 :=min_x x^-h mu_x (K_x^15 phi)(s)/C0.               (5)
```

The fifteenth vector is already available from the chain used in
(4), so (5) needs no second fifteen-step chain. It may be weaker than
applying the exact operator directly to `q_x`; its convenience does
not assert sharpness. For the proposed original-network initial
vector `phi(ell,omega)=(sqrt(2))^(ell+omega)` and `y=2`,
`mu_x=(1+x)^512` exactly, attained at the empty state.
For a later arbitrary witness, recompute `mu_x` on all physical
states; this special endpoint identity does not automatically persist.

Define physical ratios

```
rho=max_s F(phi)(s)/phi(s),
nu=max_s b_upper(s)/phi(s).
```

The resulting probability comparison is

```
C0*y^-d_star *phi(full)*rho^255*nu.                 (6)
```

It pays the row-set union once and contains no extra terminal
`y^512` factor. The same interface applies to other eligible row
counts with their non-`H` quota and binomial normalizer.

## Geometrically damped physical witnesses

The map `F` in (4) is positive, monotone, and homogeneous. If
`rho=max_i F(phi)_i/phi_i`, then for every `0<=theta<=1` the witness

```
psi_i=phi_i^(1-theta)*F(phi)_i^theta
```

has bulk Collatz ratio at most `rho`. Indeed, with
`r_i=F(phi)_i/phi_i`, one has `psi<=rho^theta phi`, hence

```
F(psi)_i/psi_i
 <=rho^theta*r_i^(1-theta)
 <=rho.                                             (7)
```

This uses only positivity, monotonicity, and homogeneity. It assumes
neither linearity nor logarithmic convexity or concavity. Normalizing
`psi` by a common scalar leaves the ratio unchanged. The half-step
update can be formed directly in logs as
`log2(psi)=(log2(phi)+log2(F(phi)))/2`, followed by a common shift.

Half damping also removes the elementary two-cycle of the matrix
`[[0,a],[b,0]]`: from any positive vector it produces a vector
proportional to `(sqrt(a),sqrt(b))`, exactly its Perron vector.
This example is not a general convergence theorem. Nor does (7)
guarantee that the complete probability bound (6) decreases, because
its full-state and endpoint factors also change. Recompute endpoint
ratios for each new witness, retain previous valid endpoint bounds,
and keep the best complete bound.

When a regional comparison has not been excluded by a lower
certificate, half damping is a principled next witness update, not a
new marker or capacity search. The lemma has passed independent
mathematical review. The sourcewise Cauchy obstruction recorded below means
that further witness iteration of that comparison is now
parked, despite this valid contraction result.

## Checkpoint and next action

The factorization, Horner induction, destination indexing, final
normalizer, and endpoint shortcut received independent mathematical
audits. The frozen implementation in
[quartet_adaptive_vector.py](quartet_adaptive_vector.py) also passed
independent static review of its direct `T/S` sums, fixed-`m` binomial,
`W+E` destination degree, independent scales and positive addition,
single old-multinomial division, and input/output scale bookkeeping.
The parent ran the small checks serially: 240 split coefficients,
420 physical actions, 168 non-exponential input rows, and 84 each
for linearity, unit vectors, exponential vectors, and external scale
handling passed.

The width-512 exponential control then passed all 33,153 rows against
the existing engine, with maximum logarithmic difference
`1.11e-15`. The Horner action took 3.714 seconds; the reference
exponential action took 1.569 seconds. A separate unit-vector control
also passed every row, with maximum logarithmic difference
`2.72e-15`, in 3.745 seconds. These validate one-update evaluations
in floating arithmetic, not a regional tail certificate.

The first authorized regional control used one positive witness with
three fixed-quota-evaluation chains of sixteen exact actions, using
the saved fifteenth iterates for (5). Its caller was uncached with a
300-second limit, reflecting roughly 180 seconds for 48 actions.
An optional cache of the powers `T^k`, which are independent of the
input vector, is a separate performance change and does not alter
the comparison.

The frozen caller
[quartet_horner_region.py](quartet_horner_region.py) passed independent
static review of normalized core calls, external scale accumulation,
the saved fifteenth iterate, fixed-`x` ordering, physical quota minima,
endpoint multipliers and reuse, start/union factors, and witness scale
invariance. The parent ran its tiny checks serially: 261 two-step
coefficients; 42 each for quota, endpoint, and scale rows; 126
initial-endpoint comparisons; 84 floating comparisons; four full
comparisons; four fixed-evaluation traces; and the deadline guard
passed.

The parent completed the full uncached regional control in 168.38
seconds. It returned `log2(rho)=737.582633`,
`log2(nu)=737.582464`, and logarithmic tail upper bound
`+137000.328`. The full physical row was worst, and every physical
row chose `x=2`. The bound did not close a codeword class. It retained
the exact physical operator at each fixed `x`, so this result is no
longer attributable to a bin or exponential-envelope projection.

A subsequent exact three-state argument proves a limitation of this
specific grid. For `x in {1,2,4}`, even the regional full-to-full
entry has Cauchy lower bound greater than `2^518`. It follows that
every positive potential for the unsplit comparison, and even every
destination-split comparison on this grid, gives a complete tail
bound greater than `2^79750` at threshold 60518. See
[QUARTET_REGIONAL_CAUCHY_OBSTRUCTION.md](QUARTET_REGIONAL_CAUCHY_OBSTRUCTION.md).
This rules out further potential iteration for that fixed comparison,
not the exact quota-conditioned operator or the construction.

The parent's tiny floating continuous-`x` diagnostic found a
three-state full-to-full minimum near `x=2.756067424`, with logarithm
`177.10420`, below the roughly `205.2268` bulk budget. That diagonal
alone therefore did not establish a continuum obstruction.

A subsequent exact certificate does exclude the original network's
entire sourcewise Cauchy comparison, even when the minimum in (4)
ranges over every positive `x`. On the three states it uses
`v=(1,2^314,2^677)`, extended by zero elsewhere, and proves
`F(v)>2^400 v`. Opposite derivative signs at rational brackets and
outward rational logarithm bounds give certified row-minimum lower
logs `410.755936`, `410.358917`, and `411.295089`; the proof asserts
the weaker integer threshold 400 exactly. Positive comparability and
the endpoint lower bound then force every positive-potential complete
comparison at target 60518 to exceed `2^49660`.

The parent ran this certificate serially, and independent source
review passed its integer matrix/derivative powers, logarithm
enclosures, convex tangent bound, and probability normalization. Its
source is
[quartet_cauchy_phase_obstruction.py](quartet_cauchy_phase_obstruction.py),
with theorem and scope in the regional-obstruction note linked above.
Thus further potential iteration or finer quota grids cannot close
this sourcewise comparison. The theorem does not refute exact quota
conditioning, destinationwise quota minima, partial-network variants,
or the construction.

The continuous entrywise minimum of the restricted three-state
regional matrix had floating Perron logarithm `177.1042`. This is a
diagnostic distinction, not a proof that destination splitting will
close the full operator. That distinction should be resolved before
another full-size control.

## Destination groups before quota extraction

Let the physical destination states be partitioned into fixed groups
`J_1,...,J_r`. For a fixed `0<epsilon<1`, put

```
w_j(s)=1_{s in J_j}+epsilon*1_{s not in J_j},
kappa=1+(r-1)*epsilon,
f_j=phi*w_j.
```

The weights do not depend on the quota variable or on the random
setup, and `sum_j f_j=kappa*phi`. Consequently the exact regional
operator satisfies

```
B phi <= U(phi)
       :=kappa^-1 sum_j min_{x in X_j}
           x^-h K_x^16 f_j/C0.                         (8)
```

Each minimum is pointwise in the physical source state. Different
groups may use different positive evaluation sets `X_j`. Equation
(8) follows by decomposing the input vector before coefficient
extraction, applying the positive Cauchy bound separately to each
piece, and then summing. It retains every physical source row and
does not aggregate source states. The sourcewise obstruction above
does not apply to this different order of operations.

The endpoint can use the same saved fifteenth iterates without a
factor of `epsilon^-1`. Define the positive-polynomial pieces
`g_j(x)=q_x*w_j`, so `sum_j g_j(x)=kappa*q_x`. With the common
multiplier `mu_x=max_s q_x(s)/phi(s)`, one has
`g_j(x)<=mu_x*f_j` for every group. Thus

```
b <= kappa^-1 sum_j min_{x in X_j}
       x^-h mu_x K_x^15 f_j/C0.                        (9)
```

Here the decomposition precedes coefficient extraction. The
pointwise multiplier `mu_x` is used only after evaluating at a
positive `x`; it need not be polynomial. In particular, the required
multiplier is `max q_x/phi`, not `max q_x/f_j`.

For any fixed positive witness, take `rho=max U(phi)/phi` and let
`nu` be the corresponding endpoint ratio from (9). Positivity of the
physical operator gives `B^255 b<=rho^255*nu*phi`, hence (6) remains
valid. This argument requires no convergence or iteration claim.

The next bounded comparison uses the fixed groups `w<256` and
`w>=256`, where `w` is persistent overlap, and the initial witness
`log2(phi(a,w))=(315*a+372*w)/512`. Its respective evaluation sets
are `{1.6,1.86,2.04}` and `{2.04,2.45,2.76}`. The restricted
three-state floating diagnostics distinguish this split from the
split by zero versus positive hole population; they do not certify
the full comparison.

The positive regularization `epsilon=2^-8192` has a direct bound for
this witness. Here `1<=phi<=2^687`, `C0>2^8178`, and
`K_x 1<=(2+x)^512 1`. Throughout the proposed evaluation interval,
`(2+x)^2/x<9` and `9^4<2^13`. Therefore each regularized off-group
contribution to a quota-scaled bulk action is less than

```
2^-8192 *2^687 *9^4096 /2^8178 <2^-2371.
```

The endpoint has the same bound using
`mu_x<=(2+x)^512`. Since `phi>=1`, these are also bounds relative to
the witness; summing the two groups costs at most a factor of two.
This addresses approximation sharpness only. An implementation must
still preserve positive terms and verify its numerical range.

The frozen caller
[quartet_destination_quota.py](quartet_destination_quota.py) passed
independent static review. The parent ran its tiny checks serially:
261 physical coefficients; 168 each for grouped quota, direct
endpoint, and shortcut endpoint rows; 84 each for regularization
identities and conservative omission of `kappa`; four composition
and fixed-evaluation traces; 84 floating rows; and 42 scale rows
passed. The implementation omits division by `kappa`, which weakens
(8)--(9) conservatively.

The full bounded control completed in 197.291 seconds. It returned
`log2(rho)=642.63667` at the empty physical state,
`log2(nu)=642.233`, and logarithmic tail upper bound
`+112868.7587`. The minimum of the comparison ratios over all physical
states was `518.937`, at `(a,w)=(460,1)`. These are floating
diagnostics, not certified lower or upper bounds. They give no reason
to continue potential iteration of this comparison. The next
discriminator is a lower bound on an exact quota-conditioned regional
return coefficient, which can distinguish a quota-comparison loss
from a limitation of the adaptive majorant itself. No codeword class
was closed by this control.
