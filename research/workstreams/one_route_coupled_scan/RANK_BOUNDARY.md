# Keep residual rank between regions

The online route operator retains rank and support inside a region. This
note gives a five-state linear envelope that retains the residual state
rank between regions as well. Its states are quiet and ranks 1, 2, 3, 4.
The construction is unchanged. The envelope is a proof interface, not a
numerical closure or a whole-code distance certificate.

## 1. Boundary functions and the linear envelope

Use the common zero-form space W and retained nonzero requirements from
[PROOF.md](PROOF.md). At a boundary let H be the state-support mask and

    r = dim((span(Q)+W)/W).

The valid states are (0,empty), called quiet, and (r,H) with
1<=r<=|H|<=4. This r is the residual rank of the current state forms, not
the remaining dimension of the whole message space. In particular, do not
reset the global factor tau^(-dim W) at a boundary.

Fix any positive function v(r,H) on the live states. Define five disjoint
basis functions:

    f_0(quiet)=1,                  f_0(live)=0,
    f_b(r,H)=v(r,H) 1[r=b],       f_b(quiet)=0,    1<=b<=4.

For a nonnegative coefficient column c, write L(c)=sum_b c_b f_b. Let T
be a sound regional continuation operator, uniform over every incoming
quotient, retained history, and remaining labeled source forms in its
scope. Assume T is monotone, positively homogeneous, and subadditive.
The [online route recurrence](ONLINE_ROUTE.md) has these properties:
each local rule is a maximum or sum of nonnegative linear continuations,
and the regional recurrence uses positive sums and compositions of those
rules.

Construct a nonnegative matrix M, with incoming state as its row index
and terminal state as its column index:

    M[0,b] = T(f_b)(quiet),
    M[a,b] = max_(H: |H|>=a) T(f_b)(a,H)/v(a,H),    1<=a<=4.    (1)

These definitions imply T(f_b)<=L(M e_b). For c>=0, subadditivity and
positive homogeneity therefore give

    T(L(c)) <= sum_b c_b T(f_b) <= L(M c).                     (2)

Monotonicity extends (2) to any terminal function bounded by L(c).
Thus chronological regional matrices M_1,...,M_s compose as

    T_1(T_2(...T_s(L(c))...)) <= L(M_1 M_2 ... M_s c).         (3)

The support mask is forgotten only at these boundaries; its effect within
each region remains inside T. The matrix entries are accounting bounds,
not transition probabilities. They need not have row sum at most one.

This refinement is not guaranteed to dominate the two-state envelope.
Splitting terminal ranks can let the maxima in T(f_b) choose incompatible
geometries for different b. Hence sum_b T(f_b) can exceed T(sum_b f_b).
Retaining rank may compensate for this loss, but that is a numerical
question, not a monotonicity theorem for the refinement.

## 2. Empty columns must retain rank too

An empty batch preserves the span of the state forms before new zero
constraints: Q'=(R+D)Q and R+D is invertible. Imposing zero equations can
lower its residual rank. A retained live branch cannot become quiet in
an empty stretch, and its residual rank cannot increase.

For L consecutive empty batches, apply (1) to the actual L-fold empty
operator. Its matrix E_L obeys

    E_L[0,0]=1,
    E_L[0,b]=E_L[a,0]=0                 for a,b>0,
    E_L[a,b]=0                         unless 1<=b<=a<=4.     (4)

Alternatively, the L-th power of a sound one-batch rank envelope is safe,
though it can lose more information. One physical outer column contains
four regions, so an empty-column matrix uses L=4m, with m=2048.

The scalar inequality T_empty(v)<=lambda v does not justify the matrix
diag(1,lambda,...,lambda) in this basis. A rank-two state can lose one
independent rank while remaining live; its continuation into f_1 need not
be zero. The scalar inequality controls the combined live terminal v,
not each rank-restricted terminal f_b. Likewise, lambda^(4m) cannot be
placed blindly on the four live diagonal entries of an empty column.

## 3. Source-pivot partitions before a late column

The [lazy outer argument](LAZY_OUTER_COEFFICIENTS.md) partitions each
group's pre-mixer vector by its first nonzero coordinate. Surviving pivot
branches occur at distinct prefix-rank gains i in {0,1,2,3}; there is at
most one such branch for each i. Its relative global potential is tau^-i.

A branch with gain i leaves a quiet state quiet. On a live incoming
state (r,H), it preserves H and leaves residual rank s satisfying

    max(1,r-i)<=s<=r.

Killing a required nonzero state coordinate rejects the branch. Thus the
following operator bounds this source partition for every f>=0:

    B(f)(quiet) = g_tau f(quiet),
    B(f)(r,H) = sum_(i=0)^3 tau^-i
                    max_(max(1,r-i)<=s<=r) f(s,H),
    g_tau = sum_(i=0)^3 tau^-i.                              (5)

This is a uniform upper envelope. It does not claim that all four pivot
branches survive, or that every allowed s is attainable. It is again
monotone, positively homogeneous, and subadditive.

Applying (1) to B gives an explicit rank matrix P:

    P[0,0]=g_tau,     P[0,b]=P[a,0]=0 for a,b>0,
    P[a,b]=(sum_(i=a-b)^3 tau^-i)
                max_(H: |H|>=a) v(b,H)/v(a,H),
                       for 1<=b<=a<=4,
    P[a,b]=0           for b>a.                              (6)

In particular P[a,a]=g_tau, not one. The sum in (6) includes rank gains
that do not lower the state rank: an equation can constrain source
directions outside the current state span.

Apply P once per group's raw occupied column, before its first component
region. For the two-group class this gives P^2 before each late occupied
column. No independence of the two partitions is needed: (5) is uniform
in the common quotient after the preceding partition. Do not apply P
per region, per active batch, or per coefficient charge. It has no output
weight factor and introduces no second initial-dimension factor.

## 4. Regional and column products for the adjacent-pair class

Fix two cyclically adjacent groups and the same h complete occupied
columns in both, 33<=h<=64. Put b=h-32. The shortened message dimension
is D=8b over F. The first b occupied columns are the fresh information-set
prefix; the remaining 32 are the late raw-Toeplitz columns. Empty physical
columns do not advance this prefix.

Let M_P^F and M_S^F be the rank envelopes of the online regional operators
for, respectively, one full source pair and two singleton source pairs,
using the fresh-source rule. Let M_P^L and M_S^L be their lazy-outer
counterparts. The singleton-pair operator must keep sampling without
replacement and the same-batch collision; see [ONLINE_ROUTE.md](ONLINE_ROUTE.md).
The late operators also require that note's enlarged filtration, delayed
source nonzero checks, and proved Toeplitz diagonal exposure order.

For an adjacent edge whose even regions contain the full pair, define

    C_F^even = (M_P^F M_S^F)^2,
    C_L^even = P^2 (M_P^L M_S^L)^2.                         (7)

For the other edge parity, interchange P and S in the four regional
factors, but keep the source-pivot factor P^2 at the left. Here the letter
P on M_P denotes the full-pair route, whereas the standalone matrix P is
the source partition (6). Let C_0 be the empty-column matrix from (4).

For each fixed h and parity alpha, define a matrix recurrence that counts
the shared occupied-column sets once. Start A[0,0]=identity and all other
entries zero. Processing one more physical column gives

    A[j+1,l] = A[j,l] C_0
                 + A[j,l-1] C_sel(l-1),

    C_sel(u) = C_F^alpha if u<b, and C_L^alpha otherwise.     (8)

Omit invalid indices. This recurrence uses chronological matrix order;
neither the regional nor the column matrices may be commuted. Its entry
A[64,h] already sums every h-subset of the 64 columns. Do not multiply
it by another binomial coefficient, and do not square the column-set
count for the two groups: the current class requires the same set.

## 5. Free terminal and initial source accounting

The final state is unrestricted. A coefficient vector bounding the
constant-one terminal is

    t_0=1,
    t_r=max_(H: |H|>=r) 1/v(r,H),       1<=r<=4.             (9)

The encoder starts quiet, so its initial row is e_0^T. With
e=floor(0.35 N), a bound for this fixed h class, summed over all cyclic
adjacent pairs, is

    (256/255)^64 * 4096 * tau^(8(h-32)-1) * z^-e
       * sum_(alpha in {even,odd}) e_0^T A_alpha[64,h] t.    (10)

There are 4096 edges of each parity. The factor tau^(D-1) is applied
once for the entire shortened message space; regional rank states do
not replace it. The 64 raw matrices are conditioned back to invertible
once through (256/255)^64. Their two pivot partitions per late column
are already accounted for in (7), and are not another conditioning cost.

Formula (10) is conditional on the local and exposure lemmas used to
construct every envelope. It covers only adjacent groups with identical
complete-column supports. Partial columns, unequal column sets,
nonadjacent pairs, and larger group unions require separate accounting.
One may sum complete h-class bounds, or select the minimum of several
independently valid bounds for each h before summing. This is not a
whole-code certificate until all remaining source classes are covered.

## Verification obligations

Before interpreting a numerical value as a bound, check the basis images,
normalizations, all maxima, and chronological products with rigorous upper
arithmetic. Use one consistent positive basis v throughout, including
its terminal vector. Inverse factors, especially rank prices and basis
normalizations, require their own rounding direction; an upward-rounded
denominator is not by itself an upper bound for a quotient.

Check empty-column rank changes and source-pivot rank changes explicitly.
Keep the online route and lazy-outer filtration hypotheses; retaining more
boundary states does not justify lookahead or independence of the four
components of a shared Toeplitz matrix. Compare the result with the older
two-state bound before attributing an improvement to rank retention.

No numerical run or new closure is claimed here. The next discriminating
test is the same large-h adjacent-pair class at the existing markers,
changing only this boundary representation and its empty/pivot operators.
