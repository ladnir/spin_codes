# Exact supports for the GF256 outer

The rank continuation counts shortened message spaces over the large field.
The established outer proof instead counts GF256 codewords. This note connects
the two using **exact** coordinate supports. It avoids introducing all
supersets of a GF256 word's support.

The bridge is a proved counting inequality. It is not a distance certificate,
and it does not transfer the old inner certificate to the coupled scan.

## 1. One fixed outer group

Let E=GF(q), let F be an extension of E, and fix an injective E-linear group
code C of length n and dimension k. Write C_F for its scalar extension to F.
For a nonempty coordinate set S, define

    C_S = {c in C : supp(c) is contained in S},
    d_S = dim_E C_S,
    A_S = #{c in C : supp(c)=S}.

The same zero-coordinate equations define the shortened space over F.
Their rank is unchanged by field extension, so its F-dimension is d_S.

Call S realizable over F if some word of C_F has exact support S. In that
case d_S>=1, and every coordinate in S is a nonzero linear form on C_S.
Conversely, a coordinate form identically zero on C_S stays zero after
extension. Such an S is not realizable over F.

All lower bounds on A_S below assume that S is realizable over F. For
nonrealizable supports the rank enumeration must assign weight zero, not
tau^(d_S-1).

### A hyperplane bound when the group length is at most q

Put s=|S| and d=d_S. If d=1, every nonzero vector of C_S has exact support S,
and therefore

    A_S=q-1.

Suppose d>=2. Each coordinate in S vanishes on a proper hyperplane of C_S.
Fix one such hyperplane H_1. Its complement contains (q-1)q^(d-1) vectors.
Any other distinct coordinate hyperplane intersects this complement in
exactly (q-1)q^(d-2) vectors. A repeated hyperplane removes no new vectors.
A union bound inside the complement of H_1 consequently gives

    A_S >= (q-s+1)(q-1)q^(d-2).                    (1)

For s<=q, this lower bound is positive. Thus, if **each group** has length
n<=q, exact supports over E and F coincide. The claim does not require the
complete encoder length N to satisfy N<=q.

For 1<=tau<=q, (1) implies

    tau^(d-1) <= f(s) A_S/(q-1),
    f(s) = q/(q-s+1).                             (2)

For d=1, the factor can instead be exactly one. The factor in (2) ranges
from near one to q. At n=q=256, its worst case alone would cost eight bits
per active group. That worst case is avoidable for the actual rate-half
outer.

This argument cannot be applied to an arbitrary long binary code. For
example, the binary [3,2] even-parity code has no full-support binary word,
although its scalar extension can have one. Its group length three exceeds
q=2, so it does not contradict (1).

## 2. The actual redundancy gives a stronger uniform factor

The following elementary bound uses the group redundancy n-k rather than
only its length. It also applies when n>q, provided the displayed factor
remains positive.

Because the coordinate map is injective on C_S, select d coordinate forms
that form a basis of its dual. There are (q-1)^d vectors on which these d
forms are all nonzero. Parameterize those vectors by their d nonzero
coordinate values.

Each of the remaining s-d coordinate forms is nonzero. Fix one variable
whose coefficient in that form is nonzero. For each assignment to the
other d-1 variables, at most one value of the selected variable makes the
form zero. Thus each remaining form excludes at most (q-1)^(d-1) vectors.
The union bound gives

    A_S >= (q-1)^(d-1)(q-1-s+d).                 (3)

Shortening a length-n, dimension-k code gives d>=s-(n-k). For the actual
group, n=256, k=128, and q=256. Consequently

    q-1-s+d >= 127,
    A_S >= 127 * 255^(d-1).                      (4)

Also d<=128. Combining (4) with tau<=256 yields the uniform bridge

    tau^(d-1) <= C_* A_S/255,
    C_* = (255/127)(256/255)^127 < 3.31.          (5)

The explicit rational expression, not a rounded decimal, is the proof
constant. Its loss is less than 1.73 bits per active group. No property of
the random mixers was used in (3)--(5); these hold for every realized
injective [256,128] outer whose support S is realizable over F.

For a known support size, use min(f(s),C_*) rather than C_* alone. If a
calculation also retains d, then for d>=2 it can use the smaller of

    q/(q-s+1) * (tau/q)^(d-1),
    (q-1)/(q-1-s+d) * (tau/(q-1))^(d-1).          (6)

Both are factors multiplying A_S/(q-1). For d=1 use one. A dimension-
dependent factor cannot be pulled outside the outer expectation without
retaining that dimension or taking a justified uniform upper bound.

## 3. Several active groups

Fix r active groups and exact supports S_1,...,S_r that are realizable for
their respective outer setups. Write d_i for their shortened dimensions
and D=sum_i d_i. The direct sum of the group spaces gives

    tau^(D-1)
      = tau^(r-1) product_i tau^(d_i-1)
      <= tau^(r-1) / (q-1)^r
         * product_i [c_i A_{S_i}],              (7)

where c_i may be C_* or another justified factor above. The factor
tau^(r-1) is required. Replacing the product denominator by a single
projective factor while omitting it would be incorrect.

Group setups are independent. Hence an outer expectation of the right
side factors for fixed supports when the c_i are deterministic bounds.
This remains valid after multiplying by any nonnegative continuation
weight depending on the exact labeled supports and the independent route,
provided that weight is uniform over the remaining outer coefficients.
This is the interface supplied by a sound support-only rank transfer.

The total support need not be realized by one GF256 scalar multiple of
an F-valued witness. The argument instead uses the shortened dimension
and independent exact-support counts in each group.

## 4. The established Toeplitz measure applies directly

The actual outer uses four RS[64,32] constituents over E. Each of their
64 four-symbol columns receives an independent uniform invertible
Toeplitz4 matrix. For a fixed nonzero source column and target column,
the established atom bound is

    rho_T = 1 / ((q-1)q^3).

See [the Toeplitz argument](../char2_field_spin/HD_GROUPED_OUTER.md#invertible-toeplitz-column-mixers)
and [the sharper group density](../char2_field_spin/HD_SHIFTED_PAIR_SWAPS.md#use-the-actual-toeplitz-outer-bound).
The latter is implemented symbolically in
[hd_shifted_outer_density.py](../char2_field_spin/theory/hd_shifted_outer_density.py).
These references concern this GF256 outer distribution, not a uniformly
random F-valued matrix.

For completeness, fix a nonzero target word x in E^256, occupying h
four-symbol columns. Invertibility of the mixers preserves zero columns.
If h<33, no nonzero source can produce x. Otherwise, each constituent
shortened to these h columns has dimension h-32. Fix any h-32 occupied
columns and expose their mixers. Their inverse target values determine
at most one source group. The remaining 32 occupied-column mixers are
still independent. Each must hit its specified target, at cost at most
rho_T. Therefore

    Pr_outer[x belongs to C] <= rho_T^32.        (8)

Because the outer is injective, membership counts at most one message.
For an exact support S, let h(S) be the number of occupied four-coordinate
columns. Summing (8) over its (q-1)^|S| nonzero value assignments gives

    E_outer A_S = 0                       if h(S)<33,
    E_outer A_S <= rho_T^32 (q-1)^|S|    otherwise.       (9)

Combining (5), (7), and (9) bounds the expected source factor for r fixed
active groups and fixed labeled supports by

    tau^(r-1) [C_* rho_T^32/(q-1)]^r
      * (q-1)^(sum_i |S_i|),                       (10)

with zero weight if any group occupies fewer than 33 columns. The left
side here includes the indicator that every exact support is realizable.
Equation (10) can be multiplied by the corresponding nonnegative rank
continuation weight and summed over supports. It is not yet that sum.

## 5. What can be compressed without losing the route

For one group, introduce a marker z_r for coordinate r and u for an
occupied column. The polynomial

    P(u,z) = product_{j=0}^{63}
      [1 + u (product_{a=0}^{3}(1+(q-1)z_(4j+a)) - 1)]

has coefficient (q-1)^|S| u^h(S) on the exact labeled support S. Retaining
only u-degrees 33 through 64 supplies precisely the measure in (9), apart
from the scalar rho_T^32. This preserves the four-region column structure.

For multiple groups, keep their coordinate labels until the route has
been accounted for. Whether two active coordinates share a source pair
depends on both region parity and neighboring group indices. A scalar
weight enumerator, or a count of active groups alone, loses this information.
Similarly, a uniform distribution over occupied width-four batches does
not follow from the uniform permutation of physical pairs.

A useful first sparse gate is therefore one active group. Each occupied
region places exactly one input coordinate uniformly among its 8192
positions; different regions use independent choices. The source side
needs only its occupied-column counter and four labeled support bits per
column, rather than a list of 2^256 supports. The inner continuation must
still carry its state and any lookahead across region boundaries.

The next discriminating case is two adjacent active groups. Alternate
regions pair them together or place them in distinct physical pairs.
Retain both group column counters and the four possibilities for the two
input-support bits in each region. Passing the one-group test does not
replace this neighboring-group gate or the later larger-cluster coverage.

## 6. The exact-support bridge is not automatically tighter

The containment bridge in [rank continuation](RANK_CONTINUATION.md) remains
valid. For a fixed code its scalar generating identity is

    sum_S u^|S| (q^(dim C_S)-1)
      = sum_{c!=0} u^wt(c) (1+u)^(n-wt(c)).

Supersets can add mass, but that does not by itself make containment an
inferior bound. At tau=q its dimension inequality is nearly tight for a
fixed S. The new bridge trades superset bookkeeping for the factor c_i.

For example, drop the occupied-column gate and apply the same point cap
rho_T^32. With one uniform marker u>=0, containment is bounded by

    rho_T^32/(q-1) * [(1+qu)^256 - (1+u)^256],

whereas the uniform exact-support bridge gives

    C_* rho_T^32/(q-1) * [(1+(q-1)u)^256 - 1].

For these two uniformly relaxed expressions, containment is coefficientwise
tighter: (256/255)^256<C_*, while each coefficient uses an exponent at most
256. Do not replace containment by the uniform C_* envelope for a numerical
improvement. Retain the exact-support formulation when its dimension-dependent
factor or labeled geometry helps, and preserve the occupied-column gate.
The exact-support bridge's immediate benefit is a clean, finite-loss
interface to labeled regional supports. Its numerical benefit depends on
the continuation and on which refinements survive aggregation.

The next acceptance gate is a sound width-four sparse calculation using
(9), the actual pair route, and the common-quotient continuation. No native
code, benchmark, or whole-code numerical certificate was produced for
this note. The existing certified baseline remains unchanged.
