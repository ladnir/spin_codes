# Every support in one outer group

The [33-column argument](PARTIAL_COLUMN_MOMENTS.md) uses a pair whose rank
is constant across columns. That property fails on larger column sets.
This note gives an all-column-count replacement: exact joint support
excludes zero pair columns, and the degree of a nonzero minor bounds the
number of rank-one columns in a rational-rank-two pair.

The selected finite calculation bounds setup failure below 2^-147 at 35%
distance for every message using exactly one outer group. It does not
cover two or more active groups, and it makes no implementation-time claim.

## 1. Exact pairs bound the shortened dimension

Let E=GF(q), q=256, and let F be the degree-16 extension in the
[construction](CONSTRUCTION.md). Retain the global rank marker tau=q^2.
Fix an outer realization and a nonempty exact source support S within one
256-coordinate group. Let V_S be the E-linear space of outer words
supported within S, and write D_S=dim_E V_S. Scalar extension gives the
same shortened dimension over F.

Suppose S is realizable as an exact support over F. Then every coordinate
in S is a nonzero linear functional on V_S. Otherwise that coordinate
would also vanish on its scalar extension. For an independent uniform pair
(v,w) in V_S times V_S, each specified coordinate vanishes in both words
with probability q^(-2). Thus the number A_pair(S) of ordered pairs whose
joint coordinate support is exactly S satisfies

    A_pair(S) >= q^(2 D_S) (1-|S|/q^2).

Consequently, for D_S>=1,

    1[exact S is realizable] tau^(D_S-1)
        <= A_pair(S)/(q^2-|S|)
        <= A_pair(S)/(q^2-256).                            (1)

Unrealizable supports and D_S=0 contribute zero. This is a pointwise
statement for each outer realization, before averaging its randomness.
The members of the pair need not satisfy any inner zero equation or
represent the eventual F-valued bad message. They only bound its initial
shortened-dimension factor. The inner continuation remains uniform over
the realized shortened space.

Unlike a containment pair, a pair counted in A_pair(S) is nonzero in every
occupied outer column. This removes the rank-zero column case below.

## 2. Polynomial pairs after shortening

Fix an occupied-column set H of size h, where 33<=h<=64, and put b=h-32.
Before the column mixers, each shortened RS row is

    f_H(alpha_j) p(alpha_j),       deg p<b,       j in H,

where f_H vanishes on the 64-h excluded evaluation points and is nonzero
on H. A pair of four-row words is therefore a 4-by-2 polynomial matrix

    P(X)=(v(X),w(X)),       deg P_(i,j)<b.

The nonzero factors f_H(alpha_j) do not affect column rank or coordinate
support. Every counted pair has evaluation rank one or two on each j in H.
The rational rank means rank over E(X), not rank over E at one column.
The nonzero matrix P has rational rank either one or two.

### Rational rank one

Every rational-rank-one P has a factorization

    (v,w)=u(X)(a(X),c(X)),

where a and c are coprime polynomials and u is a nonzero four-component
polynomial vector. Fix a representative of the common nonzero constant
scaling of (a,c); this makes the factorization unique. Coprimality ensures
that u has polynomial, rather than rational, entries.

Write d=max(deg a,deg c), using the nonzero polynomial if one is zero.
Then 0<=d<b and deg u_i<b-d. There are at most
(q+1)q^(2d) normalized scalar pairs of degree d: choose their projective
leading coefficient pair in q+1 ways and their lower coefficients in
q^(2d) ways. Ignoring coprimality only increases this count.
There are at most q^(4(b-d))-1 nonzero choices for u. Hence a uniform
upper bound for the number of rational-rank-one polynomial pairs is

    R_1(b)=sum_(d=0)^(b-1) (q+1) q^(2d)
                              (q^(4(b-d))-1).              (2)

This count also includes pairs that vanish at some column of H. It remains
an upper bound for the required class, which does not. Every pair in the
required class has rank one at all h occupied columns.

### Rational rank two

For a rational-rank-two P, some 2-by-2 minor is a nonzero polynomial of
degree at most 2b-2. Every rank-one evaluation makes this minor vanish.
If R is the exact set of rank-one positions and r=|R|, then

    r<=2b-2.                                               (3)

In particular, at least h-(2b-2)=66-h occupied columns have rank two.
This conclusion uses distinct evaluation points; it does not assume
independent column ranks.

Fix H and R. Let s=min(r,b), and choose a specified s-position subset of R,
for example its first s positions. Evaluation at these s points is a
surjective linear map on degree-less-than-b polynomial coefficients.
A rank-one 4-by-2 matrix over E has exactly

    N_1=(q^4-1)(q+1)

possible values. Therefore the number of polynomial pairs with this exact
rank profile is at most

    R_interp(b,r)=q^(8(b-s)) N_1^s,       s=min(r,b).

This bound imposes only the selected rank-one evaluations. Dropping the
other rank requirements and the rational-rank-two condition increases
the count. The selected subset is fixed from R, so there is no factor for
choosing an information set. However, this bound discards every rank-one
constraint after the first b positions. The next argument retains them.

### Every rank-one evaluation can be counted

At each j in R, fix the line in E^2 containing the four rows of P(alpha_j).
Choose a nonzero annihilator (a_j,c_j) for that line. Let V be the space
of polynomial row pairs (v,w) of degree less than b that satisfy

    a_j v(alpha_j)+c_j w(alpha_j)=0,       j in R.

If V contains two linearly independent rows over E(X), these r equations
are independent over E. In particular, dim_E V=2b-r.

For r<=b, evaluation at the r distinct points is surjective, proving
independence. Suppose r>b and the equations have a nonzero dependency
(t_j). For each 0<=k<b, the dependency implies

    sum_j t_j a_j alpha_j^k = sum_j t_j c_j alpha_j^k = 0.

Put mu_j=1/product_(l in R,l!=j)(alpha_j-alpha_l). The kernel of this
b-row Vandermonde matrix consists of vectors
mu_j p(alpha_j) with deg p<r-b. Indeed, Lagrange interpolation gives
sum_j mu_j alpha_j^m=0 for m<r-1. These evaluations form a space of
dimension r-b, equal to the kernel dimension.

Thus there are polynomials A,C of degree less than r-b with

    t_j a_j=mu_j A(alpha_j),
    t_j c_j=mu_j C(alpha_j).

The pair (A,C) is nonzero because (t_j) is nonzero and each annihilator is
nonzero. Every (v,w) in V satisfies A v+C w=0 at all r evaluation points.
This also holds where t_j=0: both A and C vanish there. The polynomial
A v+C w has degree at most r-2, so it vanishes identically. Hence all of
V lies in one line over E(X), a contradiction.

There are (q+1)^r possible line patterns. For any pattern admitting a
rational-rank-two P, the preceding lemma gives exactly q^(4(2b-r)) choices
for four rows in V before imposing their remaining rank requirements.
Patterns whose constraint space has smaller rational rank contribute no
such P. Therefore

    R_line(b,r)=(q+1)^r q^(4(2b-r))

is also a uniform upper bound. An actual rank-one evaluation is nonzero,
so its line is uniquely determined. The sum over line patterns loses no
objects through an ambiguous line choice. Use the stronger of the bounds:

    R_2(b,r)=min(R_interp(b,r),R_line(b,r)).                (4)

Both are exact integer expressions. For r>b, the line bound retains the
extra constraints that the information-set bound discarded. The exponent
2b-r is at least two in the range (3).

## 3. Column probabilities under the same sampled mixer

For a fixed nonempty coordinate mask M in an occupied column, write k=|M|.
The actual mixer is uniform among invertible 4-by-4 Toeplitz matrices over
E. Retain its atom cap

    rho_T=1/((q-1)q^3).

At a rank-one column, the joint support of the pair is the support of the
image of one fixed nonzero vector. Exactly (q-1)^k target vectors have
support M. Thus its exact-support probability is at most

    xi_k=(q-1)^k rho_T.                                    (5)

At a rank-two column, it suffices to bound containment in M. The
[Toeplitz product-span argument](PARTIAL_COLUMN_MOMENTS.md#3-the-actual-toeplitz-probability-for-a-column-subspace)
gives

    theta_1=0,
    theta_2=1/((q-1)q^2),
    theta_3=1/((q-1)q),
    theta_4=1.                                            (6)

For each fixed polynomial pair, multiply these bounds across columns.
The column mixers are independent. The two words in the pair always use
the same mixer; no independent-copy approximation is made.

## 4. A positive ordered matrix sum

Use the common two-state region envelopes B_0 and B_1 from
[the empty-span argument](EMPTY_SPAN_RANK.md), verified at the same tau,z,
and lambda. B_0 treats a zero source coordinate; B_1 treats one nonzero
source coordinate uniformly routed within its region. For one active
group, the actual pair route gives precisely that uniform location,
independently between regions. The fixed-route linear envelope is proved
before taking this average.

For M contained in {0,1,2,3}, define its ordered column matrix

    C_M=B_(1[0 in M]) B_(1[1 in M])
        B_(1[2 in M]) B_(1[3 in M]).

Define

    A_0=B_0^4,
    A_1=sum_(M nonempty) xi_(|M|) C_M,
    A_2=sum_(M nonempty) theta_(|M|) C_M.                   (7)

The order of the four scalar regions is retained. These matrices need not
commute. Scalar indeterminates u and v commute with the matrix entries.
Set

    P_(h,r)=[u^h v^r] (A_0+u(A_2+v A_1))^64,
    P_(h,all1)=[u^h] (A_0+u A_1)^64.                      (8)

The first expression sums all occupied-column choices H of size h and all
rank-one subsets R of size r in chronological column order. The second
has rank one in every occupied column. The sums over M in (7) also count
each exact scalar source mask. No additional binomial or mask-choice
factor belongs outside (8).

Let ell=(1,0) select the initial quiet state. With the verified live
potential v_min, let t_free=(1,1/v_min)^T include every terminal state.
For N=2^21 and e=floor(7N/20), define

    U_one = 8192 z^(-e)/(q^2-256) ell
              sum_(h=33)^64 {
                  R_1(h-32) P_(h,all1)
                  + sum_(r=0)^(2h-66) R_2(h-32,r) P_(h,r)
              } t_free.                                 (9)

Then U_one bounds the ideal-setup probability that a nonzero message using
exactly one outer group has output weight at most e. All 8192 group choices
are included. The terminal state is free.

To justify the quantifier order, first fix the exact source support and
outer realization. Apply the uniform inner continuation with its actual
shortened dimension. Replace that factor using (1), then average the
outer by counting polynomial pairs and using (5)--(6). Finally sum source
masks, occupied columns, rank profiles, and the active group. A message
may still depend on every setup coin.

All coefficients in (2), (4), and (7)--(9) are nonnegative. A direct
dynamic program initializes a row vector p_(0,0)=ell and all other
p_(d,r)=0. At each of 64 column positions, update

    p'_(d,r)=p_(d,r) A_0
             +1[d>0] p_(d-1,r) A_2
             +1[d>0 and r>0] p_(d-1,r-1) A_1.

The resulting p_(h,r) equals ell P_(h,r). In particular p_(h,h) also
evaluates ell P_(h,all1). This is a bounded positive computation with
0<=r<=d<=64, suitable for directed upward arithmetic.

## 5. Selected finite calculation and scope

Root evaluated [the one-group checker](theory/check_one_group.py) with
tau=65536 and four pairs of exact rational markers:

| z | lambda |
|---|---|
| .999 | .9962 |
| .998875 | .9956 |
| .99875 | .9951 |
| .998625 | .9946 |

Each marker pair has its own positive envelope whose local inequalities
are verified with exact rational arithmetic. For each occupied-column
count h, first evaluate its complete bound U_h for each marker pair.
Then take the smallest of those four bounds and sum over 33<=h<=64.
Both upward-rounded 80- and 120-digit calculations give approximately

    sum_(h=33)^64 min_markers U_h = 5.5103561928e-45
                                 < 2^(-147).

The largest h-class contribution is h=44. The 8192 group choices, every
scalar support mask, and the free terminal state are included. Both replays
explicitly verify the strict integer-scaled comparison with 2^-147.

The bounded exact tests cover 2343 GF4 projective-line constraint systems,
including 910 rank-deficient kernels whose row-pair determinants all vanish;
36 exact-support pair-bridge fixtures; and 20 ordered profile-DP coefficients
against direct enumeration. These are regression checks of the analytic
lemmas, not an exhaustive proof over GF256. An independent read-only audit
also checked the polynomial counts, matrix order, source factors, per-h
minimum, free terminal, and directed arithmetic. It found no discrepancy;
the numerical executions were root's serial runs.

The per-h minimum is valid because every complete U_h bounds the same
fixed source-message class. Different h values may use different proof
markers before the final union bound. This does not take entrywise minima
of transfer matrices, mix incompatible live potentials during a product,
or choose markers according to a sampled setup.

The initial R_interp-only calculation failed badly at large h: its h=64
term was about 7.30e147 at z=.999, lambda=.9962. It discarded rank-one
evaluation constraints after the first b positions. The R_line argument
recovers those constraints by separating rational-rank-one matrices from
full-rank line-constraint spaces. That argument received an independent
algebraic review. No encoder change was needed.

This closes the complete one-active-group support class, not the whole
code. The [two-group route](TWO_GROUP_ROUTE.md) now retains physical pairs
and batch collisions, but it still leaves most multigroup classes open.
Those inputs cannot be replaced by independent singleton insertions.
