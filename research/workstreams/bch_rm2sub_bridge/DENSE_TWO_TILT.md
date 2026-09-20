# A second tilt for the selected (128,19) dense bound

The mean-only Poisson-binomial comparison admits input tails that are too
expensive to use in parts of the (128,19) proof. A second tilt changes the
bounding distribution before that comparison. It does not change the encoder.
The resulting bound keeps a two-dimensional cover: occupancy and mean label
probability. It does not enumerate every all-one count or every band type.

The construction, random setup, and bad-distance event are those in
`CURRENT_UNDERSTANDING.md`. Write L=K/128, N=256L, and H=floor(N/10).
Let T(r,z) be the selected four-state epoch upper transfer for independent
Bernoulli(r) inputs. Define its output moment bound by

\[
M(r,\lambda)=e^{\lambda H}e_Z^{\mathsf T}
 T(r,e^{-\lambda})^{N/128}\mathbf1.
\]

All bounds below concern the expected number of bad messages under the
shared setup. They require a union over every occupancy before they imply
a full distance/setup certificate.

## Fix the label coordinates before choosing the tilts

Partition the ordinary nonzero outer weights into bands g. Keep the all-one
word in a separate label *. Let A_w^cap be the certified outer shell caps.
Fix one coordinate q_g in (0,1) for each ordinary band, and set q_*=1.
These coordinates stay fixed across a complete density cover.

Fix Q occupied positions and an ordered assignment a of labels to those
positions. Define nu_a=Q^-1 sum_i q_(a_i) and theta_a=Q nu_a/L.
Inactive rows have coordinate zero. Different assignments with the same
nu_a are counted together; no distribution over assignments is assumed.

Choose an input tilt x>0. The actual row reference probabilities are

\[
p_g=\frac{q_g}{x+(1-x)q_g}.
\]

They may depend on the witness x because the assignment's cover coordinate
nu_a depends on the fixed q_g, not on p_g. Define

\[
C_g(x)=\max_{w\in g}
\frac{A_w^{\rm cap}x^w}
 {\binom{256}{w}q_g^w(1-q_g)^{256-w}},
\qquad C_*(x)=x^{256}.
\]

The all-one term is exact. Its original density factor is one; x^256 is
only the normalizer introduced by changing the input law.

## Change the input law before comparing with independent bits

For one region, let P_p be the product Bernoulli law with the row reference
probabilities p_i. Let P_q use the corresponding q_i. For a bit vector b
of weight j, direct multiplication gives

\[
P_p(b)=\left(\prod_i(1-p_i+p_ix)\right)x^{-j}P_q(b).
\]

The identity includes deterministic zero and all-one coordinates. A uniform
shuffle preserves j, so it also preserves this identity.

The shuffled P_q density is at most L+1 times the density of independent
Bernoulli(theta_a) bits, by the comparison proved in
`K30_CLOSURE_PROGRESS.md`. Multiplying the latter density by x^-j gives

\[
x^{-j}\theta_a^j(1-\theta_a)^{L-j}
=B(\theta_a,x)^L r^j(1-r)^{L-j},
\quad
B(\theta,x)=1-\theta+\theta/x,
\quad
r=\frac{\theta}{x+(1-x)\theta}.
\]

For an ordinary row of weight w, its density cost times the 256 region
normalizers equals

\[
\frac{A_w^{\rm cap}(1-p_g+p_gx)^{256}}
 {\binom{256}{w}p_g^w(1-p_g)^{256-w}}
=\frac{A_w^{\rm cap}x^w}
 {\binom{256}{w}q_g^w(1-q_g)^{256-w}}.
\]

The product row reference law makes the regions independent. Applying the
comparison to all 256 regions therefore bounds the fixed assignment by

\[
(L+1)^{256}\left(\prod_i C_{a_i}(x)\right)
B(\theta_a,x)^N M(r(\theta_a,x),\lambda).
\]

The transfer keeps the state across regions. The argument does not assume
independence of those states or independence between emitted weight and
the update syndrome.

## Sum the labels and cover the remaining two coordinates

For an auxiliary slope eta, define

\[
S(x,\eta)=x^{256}e^{-\eta}
 +\sum_g C_g(x)e^{-\eta q_g}.
\]

Expanding S^Q sums every ordered label assignment with its tilted cost.
For any interval I of nu values, the preceding assignment bound implies

\[
U_{Q,I}\le\binom LQ(L+1)^{256}S(x,\eta)^Q
\sup_{\nu\in I}
\left\{e^{\eta Q\nu}
B(Q\nu/L,x)^N M(r(Q\nu/L,x),\lambda)\right\}.
\]

Indeed, each term inside I is bounded by the supremum times
exp(-eta Q nu_a), and enlarging the remaining nonnegative sum gives S^Q.
Every interval may choose its own x, lambda, and eta. The coordinates q_g
remain fixed, so the intervals still cover every assignment. Overlapping
closed endpoints can only overcount. Bounds for all intervals are summed.

At x=1 this is the earlier label-sum bound. At the deterministic all-one
endpoint, C_*^Q B(1,x)^N=1 when Q=L. Thus the new input tilt adds no
spurious all-one multiplicity or penalty.

## Numerical implementation and audit boundary

`ladder_dense_tilt.py` constructs exact row costs with Arb and bounds every
entry of T over a density interval. It reuses the exact monotone decomposition
of the Bernstein coefficients, then powers the nonnegative upper matrix.
The choice of witness is a binary64 search; its result is not an accepted bound.

The scalar exponent is Q log S + eta Q nu + N log B(Q nu/L,x).
The checker uses a tangent upper bound for the concave function log B.
After that substitution the exponent is bilinear in Q and nu, so its maximum
over a rectangle is attained at a corner. This retains cancellation between
the label term and the input normalizer. Separately maximizing both terms
would discard that cancellation.

`search_ladder_dense_tilt.py` subdivides integer Q intervals and rational nu
intervals. A retained binary partition tree proves complete coverage, including
unresolved leaves. Each accepted leaf bounds the sum over its integer Q values.
The final dense union is an exact sum of retained dyadic upper bounds.
A separate 512-bit run checks every leaf and reconstructs that union.
An unresolved checkpoint is not a dense certificate, even if all sampled
points have positive margins. A dense certificate still excludes the sparse
occupancies below its stated minimum.

Tests check the change of measure by exact enumeration, the row-cost identity,
all-one normalization, the x=1 reduction, the scalar tangent, and rejection of
incomplete or malformed partition trees. They supplement the derivation;
they do not replace its universal probability comparisons.

## Hold the matrix reference fixed within a box

The first interval checker still loses precision by maximizing each matrix
entry over an input-probability interval. `ladder_dense_fixed_reference.py`
avoids that loss, at the cost of a polynomial factor in Q.

Fix r in (0,1) for the entire box. For each theta=Q nu/L, choose

\[
x(\theta)=\frac{\theta(1-r)}{r(1-\theta)}.
\]

This makes the reference probability in M equal to r, independent of theta.
The preceding label sum cannot be applied unchanged: its costs now vary
with nu. Instead, apply it separately to each distinct nu value. There are
at most binomial(Q+12,12) such values for the 13 labels, because a label
composition determines nu. Hence a valid interval bound is the supremum
of the previous expression with x=x(theta), multiplied by this composition
count. No enumeration is needed. This factor is included in every retained
box and is not inferred from sampled types.

It remains to bound the scalar part while keeping M(r,lambda) fixed. Write
u=log x. The function log S(exp(u),eta) is convex: each band's log cost is
a maximum of affine functions of u, and log-sum-exp preserves convexity.
Over a finite u interval, choose a line alpha+beta u that dominates both
endpoints, with beta in [0,256]. Convexity makes it an upper bound throughout.
For an interval ending at theta=1, use beta=256 and the line through an upper
bound at the lower endpoint. It remains valid to infinity because all
exponents in S are at most 256.

After the substitution, the scalar exponent is at most

\[
Q\alpha+\eta Q\nu-Q\beta\log\frac r{1-r}
-N\log(1-r)
+Q\beta\log\theta+(N-Q\beta)\log(1-\theta).
\]

Both log coefficients are nonnegative because Q<=L and beta<=256.
Replace log theta and log(1-theta) by their tangent upper bounds at an
interior theta_0. The resulting expression is affine in nu for each Q.
For each nu it is convex in Q: its Q-squared coefficient is

\[
\frac{\beta\nu}{L}
\left(\frac1{\theta_0}+\frac1{1-\theta_0}\right)\ge0.
\]

Its maximum over the Q,nu rectangle is therefore at a corner. At theta=1,
the feasible endpoint has Q=L and nu=1. The coefficient of log(1-theta)
then vanishes when beta=256; the bound extends by continuity.

The checker reconstructs alpha with outward endpoint evaluations; a rounded
slope only proposes beta. It evaluates the fixed-r matrix directly with Arb.
Tests compare the scalar bound with interior evaluations and check that a
singleton box with x=1 reproduces the earlier bound plus the composition cost.
`search_ladder_dense_fixed.py` uses this checker with the same exact partition
tree audit and a balanced subdivision rule. The earlier search implementations
and their unresolved receipts remain unchanged.
