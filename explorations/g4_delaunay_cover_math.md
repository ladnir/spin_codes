# Theorem conditions for a `g=4` anchor-mesh cover

This note states the minimum exact checks needed to turn the diagnostic
Qhull anchor mesh into a finite first-moment certificate.  Qhull, Delaunay
selection, and binary64 LP output remain discovery data.

## Profile polytope and convex witnesses

For `g=4`, put `M=N/4=524288` and

```text
Lambda = {a in Z_{>=0}^5 : sum_j a_j=M,
                              a1+2a2+3a3+4a4 >= 21}.
```

The concrete-atom orbit size is

```text
Q(a) = M!/prod_j a_j! * 4^(a1+a3) * 6^a2.       (1)
```

There are

```text
L = binom(M+4,4) - 717                           (2)
```

feasible profiles; `717` is the number of nonnegative `(a1,...,a4)`
with weighted sum at most `20`.

Use coordinates `x=(a1,a2,a3,a4)` and recover
`a0=M-sum x_j`.  The continuous clipped profile polytope is

```text
P = {x in R_{>=0}^4 : sum x_j <= M,
                         x1+2x2+3x3+4x4 >= 21}.
```

Its eight rational vertices are the four nonzero pure profiles and the four
intersections of the weight-21 plane with the edges from the all-zero class:
`a_j=21/j`, `a0=M-21/j`, for `j=1,...,4`.

For one frozen valid witness,

```text
F(a)=C-<q,a>-log2 Q(a).                           (3)
```

After the Gamma extension of (1), the nonaffine part of `F` is
`sum_j log2 Gamma(a_j+1)`.  Its Hessian is diagonal with positive entries
`psi_1(a_j+1)/ln(2)`.  Hence `F` is convex on `P`.

It follows that one fixed witness passing the five vertices of a rational
4-simplex passes throughout that simplex.  The same holds for one fixed
convex mixture of witnesses that all bound the same per-profile quantity,
use the same `Q(a)`, and are finite throughout the simplex.  The mixture
weights must be nonnegative exact rationals, sum exactly to one, and be fixed
across the simplex.  A minimax LP against five vertices has an optimal basic
solution with at most five positive witnesses (the generic Caratheodory bound
in `R^5` is six).  Binary64 weights may select the support, but the final
weights and all outward vertex inequalities must be recomputed exactly.

Fractional boundary vertices require outward Gamma evaluations.  Applying
integer factorial formulas to rounded versions of these vertices is invalid.

## Exact four-dimensional geometry

Let the stored rational anchors be `V`, including all eight vertices of `P`.
A candidate cell is a five-element subset
`sigma={v0,...,v4}`.  Its exact coordinate matrix is

```text
B_sigma = [v1-v0 | ... | v4-v0].
```

The cell is nondegenerate exactly when `det(B_sigma) != 0`.  Exact
barycentric membership is given by

```text
lambda_1..lambda_4 = B_sigma^-1 (x-v0),
lambda_0 = 1-sum_{i=1}^4 lambda_i.
```

All determinants and barycentric values must be integer/rational, not
binary64 tolerances.

The following is a sufficient exact mesh certificate:

1. every anchor lies in `P`, and the convex hull of the anchors is exactly
   `P` (including all eight extreme vertices);
2. every cell is nondegenerate and consistently oriented;
3. the cells form a simplicial complex: the intersection of two cells is
   exactly the convex hull of their shared vertices;
4. every interior tetrahedral facet belongs to exactly two cells with
   opposite induced orientations, and every boundary facet belongs to one
   cell and lies in the correct facet of `P`; and
5. the resulting complex is connected and has no missing component.

An alternative, often easier to verify, is an exact **regular triangulation**
certificate.  Assign exact rational lifting weights to the anchors and give,
for every cell, the affine lower-hull hyperplane through its five lifted
vertices.  Every other lifted anchor must lie strictly above that hyperplane.
The verifier must also establish that the ledger lists **every** lower-hull
facet, for example by exact lower-hull enumeration or by the exact facet and
boundary-incidence audit above.  Together with `conv(V)=P`, the complete set
of simplicial lower facets certifies a triangulation of `P`.
A deterministic rational symbolic perturbation can break cospherical ties.

The `QJ` option joggles coordinates.  Therefore the present Qhull simplex
list is not automatically a triangulation of the unjoggled rational anchors.
It may be used as a candidate only after one of the exact checks above.
Matching total 4-volume alone is insufficient: holes and overlaps can have
equal volume.  Floating `covered_volume_fraction=1` is not a proof.

## Support faces

A witness with zero fugacity in class `j` is finite only where `a_j=0`.
Thus it is eligible for a simplex precisely when every vertex of that simplex
has coordinate `a_j=0`.  Mixing it with a full-support witness does not make
the average finite off that face.

A full-dimensional cell therefore needs components with positive fugacity in
all five classes.  Stronger zero-fugacity witnesses may be used in separate
lower-dimensional face meshes.  Boundary profiles can instead be assigned
to adjacent full-dimensional cells, but then those cells' full-support
witnesses must certify them.  Every support stratum used separately needs an
exact ownership rule so that no boundary profile is omitted.

### Exact support partition

For a profile `a`, let `supp(a)={j:a_j>0}`.  The sets
`Lambda_S={a in Lambda:supp(a)=S}` are pairwise disjoint and their union is
`Lambda`.  Every nonempty `S subset {0,1,2,3,4}` is feasible except `{0}`:
if `S` contains a positive class, the remaining mass can be put in its
largest class.  Thus exactly 30 supports are feasible: 15 not containing
zero and 15 equal to `{0} union P` for nonempty
`P subset {1,2,3,4}`.

Fix `S`, put `s=|S|`, and shift

```text
a_j=1+b_j  (j in S),   a_j=0  (j notin S).
```

Then `b_j>=0`, `sum b_j=R=M-s`, and

```text
sum_{j in S} j b_j >= h_S,
h_S=21-sum_{j in S}j.                            (15)
```

If `0 notin S`, the cut is redundant because every profile has physical
weight at least `M`.  The exact integer hull is the shifted simplex with
vertices `R e_j`, and

```text
|Lambda_S|=binom(M-1,s-1).                       (16)
```

If `S={0} union P`, write `y_i=b_i` for `i in P` and recover
`b_0=R-sum_i y_i`.  Put `beta_P=sum_{i in P}i` and `h_P=21-beta_P`.
Here `1<=h_P<=20`, and the continuous shifted hull is

```text
H_P={y>=0:sum_i y_i<=R, sum_i i y_i>=h_P}.       (17)
```

Its lattice points are **exactly** `Lambda_S`: (17) is the original mass and
physical-weight condition after an integral change of variables.  Hence the
continuous hull is a sound, complete cover of all integer profiles in the
support.  It is generally larger than their integer hull as a real polytope,
so it can be conservative for convex interpolation.

The exact support count is

```text
|Lambda_{0 union P}|=binom(M-1,|P|)-E(P,h_P),    (18)

E(P,h)=sum_{t=0}^{h-1}[z^t] product_{i in P}(1-z^i)^(-1).
```

The subtraction counts nonnegative `y` of weighted sum below `h`; `b_0` is
automatically nonnegative because `h<=20<<R`.  Exact exclusion counts are:

```text
P       h_P   E       P       h_P   E
1        20   20      12       18   90
2        19   10      13       17   57
3        18    6      14       16   40
4        17    5      23       16   27
24       15   20      34       14   12
123      15  147      124      14  100
134      13   61      234      12   27
1234     11   94
```

These exclusions sum to 716; adding the infeasible `{0}` profile gives the
717 subtraction in (2).

### Exact integer-hull kinks for supports containing zero

The continuous hull (17) has far vertices `R e_i` and edge points
`(h_P/i)e_i`, which may be fractional.  Its integer hull often has additional
low boundary vertices.  Below, coordinates are ordered by the displayed
`P`; the complete vertex set is `B_P union {R e_i:i in P}`.  In addition to
`y_i>=0` and `sum y_i<=R`, the displayed inequalities are the nontrivial
covering facets.

```text
P      B_P                                      covering facets
1      (20)                                     y1>=20
2      (10)                                     y2>=10
3      (6)                                      y3>=6
4      (5)                                      y4>=5
12     (0,9),(18,0)                             y1+2y2>=18
13     (0,6),(2,5),(17,0)                       y1+3y3>=17; y1+2y3>=12
14     (0,4),(16,0)                             y1+4y4>=16
23     (0,6),(2,4),(8,0)                        2y2+3y3>=16; y2+y3>=6
24     (0,4),(8,0)                              y2+2y4>=8
34     (0,4),(2,2),(5,0)                        3y3+4y4>=14; 2y3+3y4>=10
123    (0,0,5),(0,6,1),(0,8,0),(1,7,0),
       (15,0,0)                                 y1+2y2+3y3>=15;
                                                  y1+y2+2y3>=8
124    (0,0,4),(0,1,3),(2,0,3),(0,7,0),
       (14,0,0)                                 y1+2y2+4y4>=14;
                                                  y1+2y2+2y4>=8
134    (0,0,4),(0,3,1),(1,0,3),(0,5,0),
       (1,4,0),(13,0,0)                         y1+3y3+4y4>=13;
                                                  y1+y3+2y4>=5;
                                                  y1+y3+y4>=4
234    (0,0,3),(0,4,0),(6,0,0)                 2y2+3y3+4y4>=12
1234   (0,0,0,3),(0,0,1,2),(0,0,4,0),
       (0,1,3,0),(1,1,0,2),(0,4,1,0),
       (2,0,3,0),(3,0,0,2),(0,6,0,0),
       (1,5,0,0),(11,0,0,0)                    y1+2y2+3y3+4y4>=11;
                                                  y1+2y2+3y3+3y4>=9;
                                                  y1+2y2+2y3+3y4>=8;
                                                  y1+y2+2y3+2y4>=6
```

The `(vertices,facets)` counts are:

```text
P=1,2,3,4: (2,2);  12,14,24: (4,4);  13,23,34: (5,5);
123,124: (8,6);  134: (9,7);  234: (6,5);  1234: (15,9).
```

There is a short exact verification.  Each covering inequality is checked
against the finite set where its positive left side is below its constant;
no point there satisfies the original weighted cover.  Conversely, exact
vertex enumeration of the displayed rational polytope gives precisely
`B_P union {R e_i}`, all feasible integers.  The polytope therefore both
contains every feasible integer point and lies in their convex hull.  This
also exposes kinks such as `(2,5)` for `P=13`, absent from the continuous
one-halfspace hull.

Use (17) first: it has few vertices and exactly the desired lattice points,
so it cannot omit a profile.  Use the tabulated integer hull only when
fractional continuous vertices make a witness fail or add material slack;
integer vertices use factorials rather than fractional Gamma evaluations.

### Support-local union ledger

Since supports partition `Lambda` exactly, one may aggregate without overlap:

```text
sum_{a in Lambda} Z(a)
 <= sum_{S feasible} |Lambda_S| 2^U_S,            (19)
```

where counts come from (16) and (18), and `U_S` is the outward convex maximum
of one fixed witness or mixture on a hull containing `Lambda_S`.  If a
support is further meshed into closed cells, replace its term by the
overlap-safe cell sum (13).  Duplicate cell boundaries are harmless within a
support; distinct exact supports never overlap.

## Lattice ownership and union accounting

An exact continuous triangulation of `P` already proves geometric coverage
of every profile in `Lambda`; it is not necessary to enumerate all `L`
profiles.  Closed cells overlap on lower-dimensional faces.  Either:

* define a canonical owner using exact barycentric membership and the least
  cell identifier; or
* deliberately retain all closed-cell memberships as a nonnegative
  overcount.

Let `U_sigma` be an outward upper bound for (3) throughout cell `sigma`, and
let `n_sigma` be an exact or rigorous upper bound on
`|Lambda intersect sigma|`.  Then the overlap-safe cell-local ledger is

```text
sum_{a in Lambda} Z(a)
  <= sum_sigma n_sigma 2^U_sigma.                 (4)
```

No ownership is needed for (4); profiles on shared faces are merely counted
more than once.  Euclidean 4-volume or the diagnostic determinant proxy is
not a lattice count.  Exact cell counts require a fixed-dimensional lattice
point algorithm (for example a verified Barvinok/generating-function
calculation), including rational boundary vertices.

The lowest-complexity safe choice is

```text
n_sigma = L
```

for every cell and an outward evaluation of

```text
log2 L + log2(sum_sigma 2^U_sigma) <= -40.        (5)
```

This loses at most the log of the number of cells relative to ideal ownership
but avoids four-dimensional lattice counting.  The still coarser condition
`max U_sigma + log2 L + log2 C <= -40`, with an exact upper bound `C` on the
number of cells, is also valid.  A reserved `24` cell bits is sound only if
the final, adaptively refined cell count is proved at most `2^24`.

For sharper accounting, provide exact canonical ownership counts
`n_sigma^own`, prove they sum to (2), and use them in (4).  Duplicate proof
coverage is never a new probabilistic event.  Convex-mixture components are
also proof devices, not union-bound branches.

If a failed terminal cell is handled by singleton witnesses, every owned
lattice point in that cell must be exactly enumerated and certified.  One
rounded centroid per failed 4-simplex is only an anchor-generation heuristic;
it provides no completeness statement.

## Cheap rigorous cell counts

The coarse choice `n_sigma=L` is always available, but three inexpensive
exact bounds can be materially sharper.  Take their minimum.

### Coordinate box plus stars and bars

For every one of the five profile coordinates, compute from the rational cell
vertices

```text
l_j = ceil(min_sigma a_j),
u_j = floor(max_sigma a_j).
```

If some `l_j>u_j`, the cell contains no integer profile.  Otherwise define

```text
B(l,u) = #{a in Z^5 : l_j<=a_j<=u_j and sum_j a_j=M}.
```

This is an exact upper bound on the number of cell profiles.  It is much
sharper than a four-coordinate box because it also enforces the omitted
coordinate's interval and the exact total mass.

To evaluate it, omit any coordinate `r`, put `d_j=u_j-l_j` for `j!=r`, and
set

```text
L_r = M-u_r-sum_{j!=r} l_j,
S_r = M-l_r-sum_{j!=r} l_j.
```

For four caps `d=(d_1,...,d_4)`, let

```text
Phi(t;d)
 = sum_{T subset {1,2,3,4}} (-1)^|T|
     binom(t-sum_{i in T}(d_i+1)+4,4),            (6)
```

where `binom(n,4)=0` for `n<4`.  Inclusion-exclusion gives

```text
B(l,u) = Phi(S_r;d)-Phi(L_r-1;d).                 (7)
```

The value is independent of which coordinate is omitted; checking that fact
for all five choices is a useful implementation invariant.  Formula (7) uses
only exact big integers and 32 binomial terms.

The physical-weight cut can be inserted almost for free.  Exactly `717`
global profiles have weight at most `20`.  Enumerate those small quadruples
once and count how many also satisfy the five intervals `[l_j,u_j]`.  Then

```text
n_box = B(l,u) - #{box profiles of physical weight <=20}              (8)
```

is still an exact upper bound for `|Lambda intersect sigma|`.

### Unimodular strip boxes

Use the lattice chart `x=(a1,a2,a3,a4) in Z^4`.  For any stored integer
matrix `U` with `|det U|=1`, put `z=Ux`.  For each row, its exact minimum and
maximum over a simplex occur at vertices.  Therefore

```text
n_U = product_i max(0,
        floor(max_sigma (Ux)_i)-ceil(min_sigma (Ux)_i)+1)              (9)
```

is a rigorous upper bound.  The transformation is a lattice bijection, so no
Jacobian or determinant division is needed.  A small fixed bank is cheap:
identity; elementary differences `x_i-x_j`; the total nonzero count
`x1+x2+x3+x4`; and physical weight `x1+2x2+3x3+4x4`, each completed to a
unimodular matrix.  Every stored matrix should be verified to have determinant
`+1` or `-1` exactly.  Discovery may propose LLL/shear matrices, but only the
integer matrix and determinant check enter the certificate.

The minimum of (8) and all values (9) is valid because it is the minimum of
upper bounds.  Intersecting their numerical counts by multiplying or dividing
them is not valid.

### Cubically padded determinant bound

Plain simplex 4-volume, even divided by a lattice determinant, is **not** an
upper bound on its lattice points: boundary and thin-cell effects can
dominate.  A safe determinant-sensitive replacement is

```text
#(sigma intersect Z^4)
 <= vol_4(sigma + [-1/2,1/2]^4)
 =  sum_{S subset {1,2,3,4}} vol_|S|(proj_S sigma).                    (10)
```

For each lattice point of `sigma`, its half-open unit cube is disjoint from
the cubes of the other points and lies in the Minkowski sum in (10), proving
the inequality.  The equality follows by adding the four unit coordinate
segments one at a time; extrusion adds the corresponding coordinate
projection volume.  The empty projection has volume one.

All vertices are rational, so every projection volume is rational.  Exact
projected convex-hull computation is small (at most five projected points).
An even simpler safe implementation upper-bounds a `k`-dimensional projected
hull by the sum of the absolute volumes of all its `binom(5,k+1)` projected
`k`-simplices.  Caratheodory shows these simplices cover the projected hull,
although they may overlap.  Let the resulting exact rational upper endpoint
be `R_sigma`; then

```text
n_pad = floor(R_sigma)                                                (11)
```

is an integer upper bound.  Applying (10) after any verified unimodular `U`
is also sound and can improve alignment.

The recommended per-cell count is therefore

```text
n_sigma = min(L, n_box, min_U n_U, min_U n_pad,U).                    (12)
```

Every term is independently checkable with exact integer/rational arithmetic.

### Overlap-safe weighted ledger

For a cell `sigma`, let one fixed exact witness or rational mixture have
outward vertex upper endpoints `V_sigma,k`.  Convexity gives the closed-cell
maximum

```text
U_sigma = max_k V_sigma,k.
```

If the exact regular mesh covers `P`, then nonnegativity and (12) give

```text
sum_{a in Lambda} Z(a)
 <= sum_sigma sum_{a in Lambda intersect sigma} Z(a)
 <= sum_sigma n_sigma 2^U_sigma.                                    (13)
```

Closed cells may share facets, edges, and vertices: (13) deliberately counts
those profiles repeatedly, which is safe.  No ownership or deduplication is
needed for this ledger.  The final computation must be an outward-rounded
log-sum-exp of all nonzero terms

```text
log2sumexp_sigma(log2(n_sigma)+U_sigma) <= -40.                       (14)
```

Cells with `n_sigma=0` contribute nothing.  Mixture components are proof
devices and do not create extra union terms.  Support-restricted components
remain eligible only when finite on the entire closed cell.

For implementation, start with (8) and the fixed unimodular strip bank; both
are tiny compared with witness verification.  Add (10)--(11) only for cells
whose contribution materially affects (14).  Exact per-cell Barvinok counts
or canonical ownership should be a last resort.

## Continuous versus discrete cover

The theorem concerns only `Lambda`, so a discrete cover is logically enough.
Nevertheless, for the first `g=4` certificate an exact rational continuous
triangulation is preferable:

* it proves completeness without enumerating roughly `2^71.4` profiles;
* convex interpolation is native to its five vertices; and
* the coarse ledger (5) avoids exact simplex lattice counts.

A discrete/hybrid cover becomes preferable only if many continuous cells
cannot be closed but their lattice residue is demonstrably small.  In four
dimensions that residue must be counted and enumerated by an exact algorithm;
sampling centroids or grids cannot replace it.

## Fatal shortcuts and recommended certificate

The following invalidate an end-to-end claim:

* treating Qhull `QJ`, binary64 determinants, or volume fractions as exact;
* optimizing a different witness at each vertex, or using a pointwise minimum
  as though it were convex;
* retaining binary64 mixture weights in the final proof;
* using a zero-fugacity face witness in a cell that leaves that face;
* rounding fractional hull vertices before evaluating `Q`;
* using simplex volume as its number of integer profiles;
* reporting one residual centroid per failed simplex instead of all residual
  lattice profiles; or
* applying a fixed cell-budget reserve after adaptive refinement without
  checking the final cell count.

The lowest-complexity sound certificate is therefore:

1. reconstruct one exact rational regular triangulation of the continuous
   clipped polytope `P` from the diagnostic anchors;
2. choose one fixed witness or at-most-five-component exact rational mixture
   per 4-simplex, respecting support domains;
3. outward-check its five vertex values, including fractional Gamma values;
4. avoid terminal singleton cells in the first pass; and
5. use the outward cell-local ledger (5), or the max-plus-cell-count version
   if its margin suffices.

Only if (5) misses should the implementation add verified per-cell lattice
counts or a discrete terminal-cell layer.
