# A high-distance gate for one-permutation field codes

**The full-uniform-permutation causal-scan prototype below is withdrawn.**
[HIGH_DISTANCE_GLOBAL.md](HIGH_DISTANCE_GLOBAL.md) proves that a suffix
preserves every shortened outer subspace. With repeated GRS[16,8], the
30% target fails with probability at least 1-2^-282 already at K=2^16,
for every scan width and every coefficient choice. The local polynomial
moment remains valid, but it cannot fix this structural obstruction.

The new target is rate one half and at least 30% symbol distance, with a
longer-term target approaching 50% over large characteristic-two fields.
A further 10% scalar-extension certificate would not address this target.
This note gives construction obstructions, a concrete proof prototype, and
a new local coefficient bound. It gives no new whole-code certificate or
encoding-time claim.

## Changes that are necessary

Let N be output length and R=K/N. Suppose an inner processes t symbols per
step with s state symbols. Its zero-state path emits the input unchanged.
Suppose the outer repeats a fixed [b,k] linear code. With T=N/t steps, take
a=floor(s(T-1)/k)+1 outer blocks, when this many blocks are available.
Their message space has dimension ak. Requiring zero feedback in the first
T-1 steps imposes at most s(T-1) equations. A nonzero message survives and
has output weight at most ba. Consequently,

    d/N <= s/(R t) + O(b/N).

The final state is discarded, so no last-step feedback equation is needed.
Input diagonals, arbitrary refresh coefficients, and larger coefficient
alphabets do not remove this witness. This is the dimension argument in
[BATCHED_PROOF.md](BATCHED_PROOF.md), applied to the new target.

At R=1/2 and t=16, s=2 has ceiling 25%, and s=3 has ceiling 37.5%.
Thus the two-state inner cannot reach 30%. The three-state inner cannot
approach 50%. At least s=4 is necessary at t=16 for the latter objective;
this removes the stated obstruction but proves no lower bound.

The coefficient field is a separate restriction. A generator over GF16 or
GF256 has exactly the same minimum symbol distance after scalar extension.
Wider records alone therefore cannot give the large-field improvement.
Native polynomial coefficients can remove this restriction while preserving
linearity over the consumer's field. Arbitrary binary mixing inside records
does not preserve that interface.

Finally, a union over all bad messages is the wrong global count.
[SUBFIELD_PROOF.md](SUBFIELD_PROOF.md) exhibits one rare routing event with a
large subspace of bad words. Its probability is tiny while its contribution
to the expected number of bad projective words exceeds 2^83. Increasing
the state dimension only shifts this example. The proof must count support
constraints or deficient cuts, retaining the dimension of each surviving
subspace.

## The withdrawn one-permutation prototype

Fix F=GF(2^m) in its declared polynomial basis, with m in {8,32,64,128}.
Write theta for the residue class of X modulo the field polynomial.
Use a genuine characteristic-two GRS[16,8] outer B independently on each
eight-symbol message block. Evaluation points may be the polynomial bit
patterns 0,...,15. This is the existing small-coefficient outer from
[SMALL_COEFFICIENTS.md](SMALL_COEFFICIENTS.md), not an integer Newton matrix
reduced modulo two.

Sample one uniform permutation P of all N coordinates. Independently, for
each scan position i and tap j=1,...,min(w,i), sample a coefficient a_(i,j)
from U_d, the set of all native polynomials of degree less than d. Zero is
included. Put Q=2^d and compute

    u = P B x,
    y_i = u_i + sum_(j=1)^min(w,i) a_(i,j) y_(i-j).

Output all N coordinates of y. The triangular scan has diagonal one and
is invertible even when a tap is zero. Thus every setup has rate one half.
The outer is local, P is the only global permutation, and the scan reads and
writes a sequential stream. A second block shuffle is not part of this
proposal. Unlike regional routing, the full uniform P permits the single
normalization 1/binom(N,s) for source weight s.

The original proof target was w=4 and d=m/2, with distance 30%. The suffix
obstruction now rules out this prototype. Increasing d changes cancellation costs;
increasing w changes both arithmetic cost and a structural distance ceiling.
The set U_d contains theta and generates F over GF2 when d>=2. For m>8 it
is not confined to the embedded GF16 or GF256. At m=8 the consumer's field
itself is GF256; the same distinction cannot be requested within that field.

There is another hard ceiling for this architecture. Let H_0 be the outer
block-diagonal check matrix. The final check matrix is H_0 P^-1 D, where
D is the inverse width-w scan. Each column of D has at most w+1 nonzeros,
so each final check column touches at most w+1 outer check groups. The
companion's exact group-count lemma therefore gives the asymptotic ceiling

    d/N <= (1/2) 2^(-1/w) + o(1).

The ceilings are 25% at w=1, 35.36% at w=2, 42.05% at w=4, and 45.86%
at w=8, rounded upward. No coefficient choice improves these ceilings.
Reaching 49% requires at least w=35 by this necessary condition. Thus w=4
passed this necessary width screen; the stronger suffix obstruction rules
out the proposed 30% target. It is not a path to 50% at fixed cost.
The exact finite group-count lemma is in the companion's
[width/coefficient analysis](../../../../gen-BAA/field-universality/paper/width-coefficient-tradeoff.tex).

## A new exact moment for polynomial coefficients

Put q=|F|=2^m and assume m>=2d-1. For independent A,B sampled from U_d,
and a projective normal [u:v] in P1(F), define

    alpha([u:v]) = max_(z in F) Pr[u A + v B = z],
    T_r = sum_([u:v] in P1(F)) alpha([u:v])^r.

The two axis normals have alpha=1/Q. For a nonzero ratio c, the remaining
atom is |U_d intersect c U_d|/Q^2. Every nonempty fiber is an affine coset
of the kernel because U_d is an additive GF2-subspace. Thus the maximum
is attained at zero and depends only on this intersection size.

Every ratio with a nonzero intersection has a unique reduced polynomial
representation a/b with a,b nonzero and degrees below d. Indeed, equality
of two ratios gives a polynomial cross-product of degree at most 2d-2.
This degree is below m, so reduction modulo the field polynomial cannot
create a new equality.

If h=max(deg(a),deg(b)), the intersection has size 2^(d-h): every solution
is a common polynomial multiple of (a,b), including the zero multiple.
There is one reduced pair at h=0 and 3*2^(2h-1) at each h>=1. To count
these pairs, the number of coprime ordered nonzero pairs of degree at most
h is 2^(2h+1)-1. Gcd decomposition proves this by subtracting twice the
total pair count at degree h-1 from the total at degree h. Subtracting the
coprime counts at h and h-1 gives the displayed height count.

There are Q^2/2-1 attained nonzero ratios. All other nonzero ratios have
intersection size one. Consequently the exact moment, for integer r>=2,
is

    T_r = 2/Q^r + (q-Q^2/2)/Q^(2r)
          + Q^(-2r) [Q^r
              + sum_(h=1)^(d-1) 3*2^(2h-1) * 2^(r(d-h))].

In particular,

    T_2 = (1+3d/2)/Q^2 + q/Q^4,
    T_3 = (9/2)/Q^3 - (7/2)/Q^4 + q/Q^6.

The q/Q^(2r) term matters when the coefficient alphabet is small relative
to the ambient field. It must not be dropped by reasoning only about
ratios attained by nonzero pairs.

## Charging a complete batch of fresh zero equations

The following transfer follows the companion's
[affine-plane batch argument](../../../../gen-BAA/field-universality/expander-lift/CONVOLUTION_AFFINE_PLANE_BATCHES.md),
with the exact characteristic-two moments above.

Fix a finite-dimensional F-space V and r>=2 equation families in V*:

    f_i = ell_i + sum_j a_(i,j) v_(i,j),  i=1,...,r.

All offsets ell_i and feedback forms v_(i,j) must be jointly fixed before
any coefficient in this batch is revealed. Conditional on that common
history, all coefficients are independent samples from U_d. Assume each
row's feedback forms span at least two dimensions before the whole batch,
and q>2r. These hypotheses exclude adaptive choices of a later row's basis
based on an earlier sampled equation.

Choose two independent feedback forms in each row. There exists an F-linear
projection V* -> F^2 injective on every chosen pair. A random projection
fails on one pair with probability less than 2/q, so a union bound proves
existence under q>2r. Fix a successful projection before exposing coefficients.
Expose all coefficients outside the selected pairs. The projected rows are
then independent affine images of U_d^2 with invertible two-by-two bases.

For a fixed line through zero, each row's membership probability is bounded
by alpha at its transported normal. Each invertible basis permutes P1(F).
AM-GM, followed by summation over every line, bounds the probability that
the projected rows have rank at most one by T_r. Their joint probability
of rank zero is at most Q^(-2r). Projection cannot increase rank. Hence,
if R is the rank of the original sampled equations and tau>=1,

    E[tau^(2-R)] <= K_r(tau),
    K_r(tau) := 1 + (tau-1) T_r + (tau^2-tau) Q^(-2r).

This is a bound on two rank gains from the entire batch. It does not assert
one gain per equation. For dimension potential Phi(D)=tau^(D-1) when D>=1
and Phi(0)=0, imposing the batch costs at most tau^-2 K_r(tau) times Phi(D).
Discarding branches that violate earlier nonzero requirements only helps.

For d=m/2 and the illustrative choice tau=q, these expressions imply

    K_2(q) < 4+3d/2,       K_3(q)-1 < (9/2) 2^(-d).

Thus at m=128,d=64, pairs cost less than 100 and triples have excess moment
below 2^(-61.83007). At m=32 and 64, the pair bounds are 28 and 52.
At m=8 the pair bound is 10. These are local moment bounds, not setup-failure
exponents for a complete code.

The global proof must keep tau tunable or group deficient subspaces
explicitly. Fixing tau=q globally can recreate the rare-window overcount:
a surviving subspace of dimension D contributes q^(D-1). The local batch
lemma does not by itself remove that global multiplicity.

Smaller alphabets trade arithmetic representation against rank cost.
At m=128,d=48, the exact triple bound has excess below
(9/2)*2^-16 + 2*2^-32. At d=40 that bound is already large. In particular,
the pair bound at d=m/4 has a term of order 2^(m+1), so half-width results
cannot be attributed to byte coefficients. The exact T_3 bound is stronger
than substituting the coarser T_2/Q into K_3.

## The global proof still needed

The scan permits the required nonadaptive exposure order. Fix a proposed
output support. In the recurrence, formally set every prescribed zero
output to zero. Expose the permutation and all coefficients on nonzero
output rows first. All input forms and remaining nonzero output forms are
then fixed. The equations on zero-output rows have independent, still
unexposed tap coefficients. Their offsets and feedback families are fixed
together, so eligible pairs or triples may be selected before those taps
are exposed. Required nonzero forms must remain alive after quotienting.

This local improvement does not complete the global argument. A useful
certificate must do all of the following:

- Count every source support using its shortened dimension. For the GRS
  outer the source polynomial is
  1+sum_(j=9)^16 binom(16,j) z^j tau^(j-8).
- Control zero-equation families whose feedback rank is zero or one, and
  count structural duplicate rows. Charging a triple for two gains leaves
  one structural deficit; the extra row is not free globally.
- Preserve subspace dimensions after previous constraints. Scalar minimum
  distance alone does not control the shortened subspaces needed here.
- Cover sparse supports, endpoint effects, dense profiles, and finite
  coefficient specializations with one finite failure bound.

The companion's
[bottleneck analysis](../../../../gen-BAA/field-universality/expander-lift/CONVOLUTION_BOTTLENECK_THEOREM.md)
explains why large backward neighborhoods do not suffice: internal cuts
may retain too much source dimension. Its generic-rank characterization
does not itself protect minors after small-coefficient specialization.
The present moment supplies a concrete specialization ingredient. It does
not license importing a two-permutation certificate into this one-permutation
construction.

## Operation cost and the next decision

The width-w prototype uses w native field-by-coefficient products and w
field additions per output, apart from the outer. At m=128,d=64, one such
product uses two 64-by-64 carryless data products. For modulus
X^128+X^7+X^2+X+1, the high part has degree at most 62, so multiplying it
by X^7+X^2+X+1 needs only one fold, with no subsequent field reduction.
One implementation uses a third PCLMUL for that fold; another uses fixed
shifts and XORs. VPCLMUL can batch independent lane products, but the scan's
output dependency still has to be scheduled explicitly.

At m=32,d=16 and m=64,d=32, the data product fits one 64-by-64 CLMUL,
followed by a fold. At m=8 the native AES-field product can use GFNI.
These counts exclude loads, shuffles, reductions between SIMD lanes, and
the outer. For 128-bit records, three CLMUL operations are not automatically
cheaper than an optimized general 128-by-128 Karatsuba product plus its
reduction. No speed advantage is claimed from the phrase small coefficient.

Storing w=4 independent 64-bit taps takes 32 bytes per output, plus the
permutation indices. This substantial traffic makes width four a proof
prototype, not an approved performance choice. The batched t=16 alternative
would need at least four state symbols before approaching 50%; a three-stage
two-by-two block-shear refresh uses 12 native coefficient products per step.
Its own joint output/return analysis is still required, and GFNI for fixed
subfield maps requires an explicit compatible field representation.

The subsequent [global analysis](HIGH_DISTANCE_GLOBAL.md) supplies an exact
obstruction to this prototype. The next proof must change the route or the
outer's shortened subspaces before considering a native kernel. A complete
30% bound for a replacement ensemble would be new: the binary rate-half lift
cannot inherit that guarantee, whereas another 10% field lift does not
advance the stated distance objective.

## Exact checks and independent review

[test_polynomial_box_rank_moment.py](theory/test_polynomial_box_rank_moment.py)
has five passing standalone tests. They enumerate ratio intersections over
GF4, GF16, and GF256; check the closed second and third moments; enumerate
affine pairs and triples with offsets and different feedback planes over
GF16; and verify the native-size rational bounds. They reject use of the
histogram when m<2d-1. These are exact local checks, not a numerical
distance search or an encoding benchmark.

An independent review confirmed the ratio-count proof and the affine-batch
transfer with the common-history condition. It specifically rejected the
weaker condition that each later basis is fixed only before its own random
coefficients: that permits adaptive orientations and loses independence.

    python -B research/workstreams/char2_field_spin/theory/test_polynomial_box_rank_moment.py -v
