# Mathematics of a `g=2` triangular profile cover

This note isolates the finite convexity argument needed for a theorem-level
cover of the `g=2` packet-profile simplex.  It does not certify any numerical
witness.  Its purpose is to state exactly what a future outward-checked atlas
must prove and to rule out several tempting but invalid shortcuts.

## Profiles, multiplicities, and one frozen witness

Put

```text
M = N/2 = 1048576,
a = (a0,a1,a2),
a0+a1+a2 = M.
```

Here `aj` is the number of two-bit atoms of Hamming weight `j`.  The class
multiplicities are

```text
c0 = 1, c1 = 2, c2 = 1.
```

Consequently the number of ordered binary atom strings with profile `a` is

```text
Q(a) = M!/(a0! a1! a2!) * 2^a1.                 (1)
```

This is the normalization used by `normalization_log2(2, a)`.  The factor
`2^a1` must be present: the middle class contains the two concrete atoms
`01` and `10`.

Fix every choice made in one valid Cauchy witness: inner output pole, the
three per-concrete-atom fugacities, Perron domination vector, and one globally
valid frozen outer branch.  Its per-profile log bound has the form

```text
F_W(a) = C_W - <q_W,a> - log2 Q(a).              (2)
```

The signs in (2) agree with the current scripts: if `f_j` is the fugacity of
one concrete weight-`j` atom, then `log2(f_j)` occurs in `q_W`, while the
factor `c_j^a_j` remains in `Q(a)`.  Equivalently, one may use aggregate
class fugacity `t_j=c_j f_j`; the two conventions must not be mixed.

The outer branch in (2) must remain frozen throughout the cell.  In
particular, either a fixed linear-BL witness or the globally affine output of
`fixed_total_weight_outer` is admissible.  Reoptimizing the pole, taking the
pointwise minimum of outer branches, or using the clamped optimized value
from `total_weight_outer` at each vertex does not produce one affine witness.

## Convexity lemma

Extend (1) to real `a_j >= 0` by replacing factorials with Gamma functions:

```text
log2 Q(a)
  = log2 Gamma(M+1)
    - sum_j log2 Gamma(a_j+1)
    + a1.
```

**Lemma 1 (strict convexity).**  For a fixed witness `W`, the extension of
`F_W` is continuous on the closed nonnegative simplex and strictly convex on
its affine plane `sum_j a_j=M`.

**Proof.**  Apart from constants and affine terms,

```text
F_W(a) = sum_j log2 Gamma(a_j+1).
```

Its Hessian is diagonal with entries

```text
psi_1(a_j+1) / ln(2),
```

where the trigamma function `psi_1` is positive on `(0,infinity)`.  Hence
the Hessian is positive definite, including on every nonzero tangent vector
whose coordinates sum to zero.  Continuity at `a_j=0` follows from
`Gamma(1)=1`.  QED.

**Corollary 2 (one-witness triangle test).**  Let `T` be any real triangle in
the profile plane and suppose the same finite witness `W` is valid throughout
`T`.  If

```text
F_W(v) <= tau
```

at all three vertices `v` of `T`, then `F_W(a) <= tau` for every real point
of `T`, and therefore for every integer profile in `T`.

Indeed, for `a=sum_i lambda_i v_i`, convexity gives

```text
F_W(a) <= sum_i lambda_i F_W(v_i) <= tau.
```

The triangle need not be unimodular and its vertices need not be integer.
For a clean exact certificate, integer vertices are preferable because (1)
then uses factorials only.  If fractional vertices are used, their Gamma
values need rigorous outward enclosures.

The phrase **the same witness** is essential.  It is invalid to check vertex
`v_i` with a separately optimized witness `W_i`: the pointwise minimum of
convex functions need not be convex.  A fixed convex mixture is allowed, but
it has additional domain and arithmetic obligations stated next.

## Fixed convex mixtures

Let `Z(a)>=0` be the same per-profile quantity for every component, and
suppose each witness `W_i` proves

```text
log2 Z(a) <= F_i(a)
```

throughout a common domain `D`.  Use `log2(0)=-infinity`.  If
`lambda_i>=0` and `sum_i lambda_i=1`, then pointwise on `D`

```text
log2 Z(a)
  <= min_i F_i(a)
  <= sum_i lambda_i F_i(a).                       (3)
```

The second inequality is just “minimum at most weighted average.”  It does
not use Jensen's inequality and it fails in general if a weight is negative
or if the weights do not sum to one.

If every component uses the same profile normalization `Q(a)` and

```text
F_i(a) = C_i - <q_i,a> - log2 Q(a),
```

then its fixed mixture is

```text
F_lambda(a)
  = C_lambda - <q_lambda,a> - log2 Q(a),          (4)
C_lambda = sum_i lambda_i C_i,
q_lambda = sum_i lambda_i q_i.
```

Thus Lemma 1 applies without change.  The coefficient of `-log2 Q(a)` is
exactly one only because the weights sum to one.  Components expressed using
different class-fugacity conventions must first be rewritten against the
same `Q(a)`.  More importantly, all components must bound the same `Z(a)`:
one cannot average bounds for different conditional events or different
outer ensembles merely because their displayed formulas look alike.

**Theorem 3 (fixed-mixture triangle criterion).**  Let `T` be a triangle and
let `I` be a finite set of witnesses such that every `W_i`, `i in I`, is
finite and valid on all of `T`, bounds the same `Z(a)`, and has the common
normalization `Q(a)`.  Let one vector `lambda` satisfy
`lambda_i>=0` and `sum_i lambda_i=1`.  If

```text
sum_i lambda_i F_i(v) <= tau
```

at each vertex `v` of `T`, then `log2 Z(a)<=tau` at every integer profile in
`T`.

**Proof.**  Equation (3) bounds `log2 Z` by (4).  Equation (4) is convex by
Lemma 1, so its vertex values bound it throughout `T`.  QED.

The vector `lambda` may be chosen separately for different triangles, but it
must be fixed across all points and all vertices of any one triangle.  A
profile-dependent LP choice recreates a pointwise minimum and loses the
convexity argument.

### Mixture sparsity

For triangle vertices `v_1,v_2,v_3`, map witness `i` to

```text
b_i = (F_i(v_1),F_i(v_2),F_i(v_3)) in R^3.
```

Caratheodory's theorem says that an arbitrary point of `conv{b_i}` can be
represented using at most four witnesses.  The actual certificate LP is
slightly sharper:

```text
minimize    t
subject to  sum_i lambda_i F_i(v_k) <= t,  k=1,2,3,
            sum_i lambda_i = 1,
            lambda_i >= 0.                         (5)
```

An optimal extreme solution of (5) uses at most three positive `lambda_i`.
Indeed, if `r<=3` vertex inequalities are active and `p` weights are
positive, an epigraph extreme point has `p+1` free variables (the weights and
`t`) constrained by the normalization and those `r` active equalities;
extremality requires `p+1<=r+1`, hence `p<=3`.  Equivalently, this is the
standard support bound for a mixed strategy against three columns.  For a
degenerate edge or point the corresponding bounds are two and one.

Therefore a successful triangle never needs more than three mixture
components if (5) itself is the only selection LP.  The safe generic
Caratheodory statement remains “at most four” for a preselected feasible
mixture that is not replaced by a minimax optimum.

### From binary64 LP weights to an exact certificate

Binary64 weights returned by `linprog` are discovery data, not certificate
coefficients.  A robust hardening procedure is:

1. Use the floating LP only to choose a support of at most three witnesses.
2. Compute outward upper endpoints `U_ik >= F_i(v_k)` for every selected
   witness and vertex, and a downward endpoint `tau_lo <= tau`.  Store these
   endpoints as exact dyadic rationals.
3. Find exact rational weights `lambda_i>=0`, `sum_i lambda_i=1`, satisfying

   ```text
   sum_i lambda_i U_ik <= tau_lo,  k=1,2,3.       (6)
   ```

   With at most three weights this can be done by exact rational linear
   algebra or a tiny rational LP.
4. Re-evaluate (6) exactly; do not infer it from the floating LP residual.

One may instead round the first `p-1` weights to dyadics and define the last
as `1-sum`, provided all weights remain nonnegative and (6) is recomputed.
Rounding each weight independently and then approximately renormalizing is
not sufficient.  A useful diagnostic error estimate is

```text
|sum_i (lambda_i-lambda_i') F_i(v_k)|
  <= ||lambda-lambda'||_1 max_i |F_i(v_k)|,
```

but the final proof should use exact evaluation of (6).

If the floating solution has strict margin beyond all interval and rounding
errors, density of the rationals guarantees a nearby rational solution.  At
zero margin this need not be true for an arbitrarily prescribed rounding;
the exact rational feasibility problem must decide the issue.  Because the
coefficients in (6) are rational interval endpoints, every nonempty hardened
LP has a rational basic feasible solution.

## The exact weight-21 clipping

The physical Hamming weight represented by a profile is

```text
w(a) = a1 + 2 a2.
```

The relevant integer set is

```text
Lambda = {a in Z_{>=0}^3 : sum a_j=M and w(a)>=21}.   (7)
```

Using all `binom(M+2,2)` profiles in the union allocation is safe.  If the
exact count is wanted, precisely

```text
sum_{a2=0}^{10} (21-2 a2) = 121
```

profiles have physical weight below 21, so

```text
|Lambda| = binom(M+2,2) - 121.                    (8)
```

The continuous cut `w>=21` has the fractional edge point
`(a1,a2)=(0,10.5)`.  It is safe but not exact to cover that larger continuous
polytope.  The exact integer hull has an additional valid inequality

```text
a1+a2 >= 11,
```

because 21 physical one-bits require at least 11 nonzero two-bit atoms.  In
`(a1,a2)` coordinates its five vertices are

```text
(M,0), (0,M), (21,0), (1,10), (0,11).            (9)
```

Lifted back to `(a0,a1,a2)`, they are

```text
(0,M,0),
(0,0,M),
(M-21,21,0),
(M-11,1,10),
(M-11,0,11).
```

Thus the clipped integer domain is a pentagon, not a triangle.  A triangular
atlas must cover or triangulate this pentagon.  Exact geometric coverage can
be checked with integer orientation determinants or rational barycentric
coordinates; binary64 point-in-triangle tests are not certificate evidence.

## Exact edges and zero fugacities

The three support-two faces reduce to the following integer segments:

```text
a2=0:  a1=21,...,M,
a1=0:  a2=11,...,M,
a0=0:  a1+a2=M.
```

A witness with fugacity `f_j=0` is valid only on the face `a_j=0`.  Adopt the
extended-value convention

```text
-a_j log2(f_j) = 0        if a_j=0 and f_j=0,
                   +infinity if a_j>0 and f_j=0.
```

Therefore one zero-fugacity edge witness may cover an edge interval by
checking the two interval endpoints, but it cannot cover a triangle that has
even one point with positive omitted-class count.  In particular, continuity
from an edge does not justify the adjacent `a_j=1` lattice layer.

The same restriction applies to mixtures.  If `lambda_i>0` and component
`i` has zero fugacity in a class that is positive somewhere in the triangle,
then `F_i=+infinity` there and the weighted average in (3) is infinite.  A
finite full-support component does not “rescue” that average.  A component
with `lambda_i=0` may simply be deleted.  Hence every positive-weight
component of a triangle mixture must be valid on the intersection of all
profiles in that triangle.  Face-restricted witnesses may be mixed only for
triangles lying wholly in the same face (or a smaller common face).

This does not prevent using the pointwise minimum in the initial inequality
of (3): a finite full-support witness may still bound a point where another
witness is inapplicable.  What fails is replacing that minimum by one finite
convex average containing an inapplicable component.

The current `support` eligibility rule in the cell-cover code captures the
right mathematical condition only when it requires the omitted coordinate's
*upper* bound to equal zero.  Merely having one or two vertices on that face
is insufficient.

## Full support and thin boundary layers

Every integer profile is either on a zero-coordinate face or has
`a0,a1,a2>=1`; there is no limiting region between those cases.  The integer
hull of the full-support part of (7) has four vertices:

```text
(M-20,19,1),
(1,M-2,1),
(1,1,M-2),
(M-11,1,10).                                  (10)
```

It is a quadrilateral.  A positive-fugacity witness that succeeds at the
vertices of either triangle in a triangulation of (10) covers that entire
triangle, including the thin `a_j=1` layers.  More generally those layers
may be partitioned into smaller triangles, but every witness used there must
have strictly positive fugacity in all three classes.

Replacing zero by a representable positive floor, such as `2^-16`, does give
a mathematically full-support witness.  It is then valid on the thin layer,
but its actual affine charge must be used and all three triangle vertices
must pass.  There is no theorem that a sufficiently thin layer inherits the
edge margin: the term `-a_j log2 f_j` can be large already at `a_j=1`.

## Hybrid convex-triangle and singleton certificate

A finite cover need not close every leaf by convex interpolation.  It may use
fixed-witness triangles where possible and explicitly certify every lattice
profile in the remaining terminal triangles.

Let `P=conv(Lambda)` be the integer hull described above.  A hybrid
certificate has two kinds of leaves:

* a **convex leaf** `(T,W)`, where `T` is a rational triangle and `W` is one
  fixed witness or fixed rational mixture satisfying Theorem 3 on all of
  `T`; and
* a **terminal leaf** `(T,{W_a})`, where every lattice point
  `a in T intersect Lambda` is explicitly listed or generated exactly, and
  `W_a` is any valid pointwise witness proving `F_{W_a}(a)<=tau`.

The witnesses `W_a` in one terminal leaf may differ from point to point.
This is valid because no interpolation is claimed there.  Checking only the
vertices of a terminal triangle with different witnesses is not enough; all
of its lattice profiles must be certified individually.

**Theorem 4 (hybrid completeness criterion).**  Suppose there is a finite
family of leaves and an exact ownership map

```text
owner : Lambda -> leaves
```

such that:

1. `a` belongs to the triangle of `owner(a)` for every `a in Lambda`;
2. if `owner(a)` is a convex leaf, that leaf's one fixed witness or mixture
   satisfies the outward vertex inequalities of Theorem 3;
3. if `owner(a)` is a terminal leaf, its exact singleton table contains `a`
   and proves a pointwise outward inequality at `a`; and
4. every witness used at `a` is finite and valid at `a`, including all
   support restrictions.

Then `log2 Z(a)<=tau` for every `a in Lambda`.

**Proof.**  Fix `a`.  Its owner exists by condition 1.  For a convex owner,
Theorem 3 applies.  For a terminal owner, condition 3 is exactly the desired
pointwise inequality.  QED.

### Exact partition and boundary ownership

An exact real triangulation of `P` is a convenient sufficient construction
of `owner`, but the theorem only requires an exact partition of its lattice
points.  A verifier may establish this in either of two ways:

* prove with rational orientation determinants that the closed triangles
  cover `P` and have pairwise-disjoint relative interiors, then assign every
  shared edge and vertex to one incident leaf by a deterministic rule; or
* verify recursively that each parent lattice set is the disjoint union of
  its children.  For an integer affine split this is naturally written as
  `ell(a)<=c` versus `ell(a)>=c+1`, which has neither a gap nor an overlap.

The outer pentagon boundary, the physical-weight boundary, all internal
diagonals, and zero-coordinate faces must participate in this check.  A sum
of reported leaf cardinalities proves completeness only after disjointness
or canonical ownership has been established.

Closed triangles are allowed to overlap on edges and vertices.  Such overlap
does not invalidate the mathematics: define `owner(a)` to be, for example,
the least-indexed valid leaf containing `a`.  Duplicate singleton rows may
likewise be deduplicated by the exact profile tuple.  What is invalid is to
count duplicate rows toward a claimed coverage cardinality without first
showing that every distinct member of `Lambda` occurs.

For a compact machine certificate, each terminal leaf should record its
rational triangle, an exact lattice enumeration rule, the number of distinct
owned profiles, and a digest of the sorted owned profile/witness records.
The verifier should recompute the ownership counts and check that their sum
is exactly `|Lambda|` from (8).

### Union-bound accounting

The probabilistic branches are the profiles `a`, not the proof leaves,
witnesses, mixture components, or duplicate geometric memberships.  With a
uniform allocation

```text
tau = -40 - log2 L,
```

where `L>=|Lambda|`, Theorem 4 gives

```text
sum_{a in Lambda} Z(a)
  <= |Lambda| 2^tau
  <= 2^-40.                                      (11)
```

Each profile appears exactly once in (11), through its canonical owner.  A
profile covered by several triangles or several witnesses receives the
minimum available certified bound, but it is not a new event each time and
must not consume multiple profile allocations.  Similarly, the components
of a convex mixture are geometric-mean proof devices in (3), not separate
union-bound events.

One may instead use nonuniform per-profile thresholds `tau(a)` provided the
certificate proves

```text
sum_{a in Lambda} 2^tau(a) <= 2^-40.
```

Again the sum is over distinct profiles after ownership/deduplication.  If an
implementation deliberately sums per-leaf contributions without canonical
ownership, duplicates are an additional safe overcount, but then its target
allocation must pay for those duplicate terms and it is not using the sharper
profile-count accounting above.

### Cell-local weighted aggregation

The uniform allocation is convenient for constructing a cover, but it need
not be used in the final union ledger.  For every convex leaf `T_r`, let
`U_r` be an outward upper bound on the maximum of its one fixed witness or
fixed mixture at the three vertices, and let

```text
L_r = |T_r intersect Z^2|
```

be its exact closed-triangle lattice count.  Convexity gives the same bound
`U_r` at every integer profile in the cell.  If `S` is the globally
deduplicated set of singleton profiles in the terminal leaves, with outward
point bounds `u_a`, then

```text
sum_{a in Lambda} Z(a)
  <= sum_r L_r 2^U_r + sum_{a in S} 2^u_a.       (12)
```

Equation (12) remains valid if adjacent closed triangles count shared edge or
vertex profiles more than once: every summand is nonnegative, so retaining
those duplicates is a conservative overcount.  It is therefore unnecessary
to build a delicate canonical boundary owner merely to improve the union
ledger.  The verifier computes every `L_r` independently from doubled area
and boundary gcds using Pick's theorem, checks any producer-reported count,
and evaluates the final log-sum-exp with outward intervals.

This aggregation can be much sharper than
`|Lambda| 2^(max_a tau(a))`, while consuming exactly the same geometric cover
and witness inequalities.  In particular, it distinguishes a nearly tight
witness on a small terminal region from one that is nearly tight throughout
the entire profile simplex.

## Certificate criterion and plan audit

Let `tau=-40-log2 L`, where `L` is the chosen number of profile branches
(either the safe all-profile count or the exact value (8)).  A theorem-ready
triangular certificate consists of a finite list `(T_r,W_r)` in which `W_r`
is either one witness or one exact rational convex mixture satisfying:

1. every integer point of `Lambda` belongs to at least one `T_r`;
2. every positive-weight component of `W_r` bounds the same per-profile
   quantity, uses the same `Q(a)`, and is valid on all of `T_r`;
3. the mixture weights are nonnegative exact rationals, sum exactly to one,
   and remain fixed throughout `T_r`;
4. every zero fugacity in a positive-weight component corresponds to a
   coordinate identically zero throughout `T_r`;
5. `F_{W_r}(v) <= tau` is proved with outward arithmetic at every vertex of
   `T_r`;
6. the profile normalization is exactly (1), with no duplicate class factor;
7. the finite triangle-union coverage is checked with exact integer/rational
   geometry.

A hybrid certificate may replace any unresolved `(T_r,W_r)` by a terminal
leaf satisfying Theorem 4.  Its completeness check is the exact ownership
partition, not the number of optimization calls or the number of rows before
deduplication.

Under these conditions Corollary 2 and a union bound prove the desired
profile sum.

The triangular-cover idea itself is mathematically sound.  The plan would be
invalidated by any of the following shortcuts:

* selecting a different optimized witness at each vertex;
* using floating LP weights, independently rounded weights, or weights that
  vary with the profile inside one triangle;
* treating the pointwise minimum of the outer, inner, or bijection branches
  as convex;
* mixing witnesses that do not bound the same per-profile quantity or use the
  same orbit normalization;
* extending a zero-fugacity edge witness into a positive-count boundary
  layer, including by averaging it with a finite full-support witness;
* using the four vertices of the continuously clipped simplex while claiming
  they are the exact integer hull, or omitting the `(1,10)` integer-hull
  vertex;
* dropping `2^a1` from `Q(a)`, or including it both in `Q(a)` and in an
  aggregate fugacity charge;
* relying on sampled or binary64 triangle coverage rather than exact finite
  geometry and outward vertex inequalities;
* treating differently witnessed terminal-triangle vertices as a convex
  certificate for unenumerated interior lattice points;
* summing leaf cardinalities across shared boundaries without exact ownership
  or deduplication, or charging duplicate proof coverage as distinct profile
  events in the union bound.

No flaw was found in the core convexity reduction.  The substantive remaining
obligation is computational: construct enough fixed-witness or fixed-mixture
triangles, outward-certify every singleton left in terminal triangles, and
verify one exact ownership partition of all profiles.
