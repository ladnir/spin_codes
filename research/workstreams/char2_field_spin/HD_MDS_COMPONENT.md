# Exact MDS weights in the component transfer

Exact MDS feasibility can be retained by a positive component transfer.
It removes substantial mass in bounded examples, without tracking field
values individually. It does not repair the full raw first moment: the
independent window/rank obstruction described below survives the actual
MDS checks, not merely their earlier envelope. The useful scope is sparse
classes or a future proof that bounds setup events rather than all witnesses.

## Partition by exact active rows

Fix the placement, with one coordinate of each source row in each of b
regions of length R. In the current construction, b=64, R=32768, and each
row belongs to GRS[64,32] over GF65536. Prefix weights are independent
uniform nonzero scalars. Prefix state and both post states carry globally.
Regional permutations are mutually independent and uniform.

For a fixed final output, invert the post-inner, regional permutations,
and unweighted prefix. Its source differences have a fixed support before
exposing the independent inverse prefix weights. Let

    gamma(w) = M_w / (q-1)^w,

where M_w counts MDS row codewords with one prescribed exact support.
Thus gamma(0)=1, gamma(w)=0 for 1<=w<=32, and gamma(w)<=z=(q-1)^(-32)
for w>=33. The exact source-membership probability is the product of
these row factors.

Partition nonzero inverse words by their EXACT active row set J. This is
different from the permitted-row Bernoulli expansion. At every coordinate
of J, the prefix difference is free; all other differences are zero.
The corresponding free component values form t in F^(b|J|). Let w_r(t)
count the nonzero differences in row r, and let C_J t be the output after
the actual regional route and post-inner.

For v in (0,1], define

    E_J(v) = E_(route,post) sum_t
                 [every r in J has w_r(t)>0]
                 product_(r in J) gamma(w_r(t)) * v^wt(C_J t).

The zero assignment is absent. Projective division gives the valid target

    Pr[distance<=e] <= 1/(q-1) sum_(a=1)^R
                           v_a^(-e) sum_(|J|=a) E_J(v_a).

Marks may be chosen separately for exact-active classes. There is no
(1-z)^(R-a) coefficient in this decomposition.

The already computed sparse receipts can be reused conservatively.
If W_J(v) is the unfiltered component enumerator, then

    E_J(v) <= z^a (W_J(v)-1).

The published singleton and pair numerical uppers already dropped their
(1-z)^(R-a) coefficients. The same numbers therefore bound these exact-active
classes by the displayed inequality. This is a separate justification,
not an identification of permitted rows with active rows.

## Nonzero differences need only local equality information

At one free prefix cut, let u be the preceding component value and u' the
new value. The source difference is u'-u. It is not determined by whether
the two component values are individually nonzero: two live values may be
equal.

Nevertheless, full field-value tracking is unnecessary. Conditional on
any fixed preceding value with the indicated live flag, the possible new
values have these multiplicities. Put Q=q-1.

| Old component | New component | Nonzero difference | Multiplicity |
|---|---|---|---:|
| zero | zero | no | 1 |
| zero | live | yes | Q |
| live | zero | yes | 1 |
| live | live | no | 1 |
| live | live | yes | Q-1 |

In the fourth row, the new value equals the old value. In the fifth, it is
one of the other nonzero values. These counts are uniform in the particular
old nonzero value. They can therefore be composed along the global prefix.

Each cut has a source-row label. Multiply its contribution by a row variable
x_r when the difference is nonzero. The regional transfer becomes a positive
matrix polynomial in x_1,...,x_a. Its other state records the carried
component live flag and the post-state support class.

For each sequence of component live flags, the interval lengths determine
the prefix occupancy h in that region. Apply the normalized regional post
kernel F_h, preserving post state across boundaries. The local post support
law is uniform over actual nonzero input values. Thus repeated component
values introduce no extra post-state assumption. Regional routing still
uses normalization C(R,h), not a full-length interleaver.

After multiplying the b regional polynomials, the coefficient of
x_1^(w_1)...x_a^(w_a) is the component contribution with those exact row
degrees. Multiply it by product gamma(w_r), excluding zero row degrees.
This gives E_J(v) exactly when F_h is exact. A positive upper kernel gives
an upper bound instead.

The helper uses the conditioned four-product post law in its exact tests.
The shared two-state majorant can be substituted to cover either that law
or the three-product post ensemble. Their distributions are not equated.

## Two positive relaxations

For feasibility alone, replace each feasible active row factor by z and
require every counter to reach d=b-k+1. Counters may saturate at d, since
larger degrees receive the same factor. This needs (d+1)^a degree states.
For the current outer and a=3, that is 34^3=39304, rather than 65^3.
The region's cut order still uses the actual source-row labels.

A cheaper relaxation uses a single degree mark lambda>=1. Every feasible
active assignment has total degree at least da. Hence

    E_J(v) <= z^a lambda^(-da)
                 sum_(t!=0) lambda^(sum_r w_r(t))
                            E[v^wt(C_J t)].

The component frontier remains four states with the shared two-state post
control. Its per-cut weights, in old/new live-flag order, become

    00: 1;    01: Q*lambda;
    10: lambda;    11: 1+(Q-1)*lambda.

The 11 weight must retain its equality term; replacing it by Q*lambda
loses the distinction. Subtract the all-zero assignment before applying
lambda^(-da). This bound retains only total degree, not every row's minimum.
It can be weaker than the saturated-counter bound.

## Bounded exact checks

Six tests pass. The positive transfer agrees with two independent oracles:
enumeration of all GF4 component values, and forward enumeration of all
active GRS messages and all nonzero prefix diagonals. These checks include
changing cut order between regions and global prefix carry.

At v=1, the source-weighted transfer counts exactly (q^k-1)^a active
outer words. The tests also verify that saturated counters preserve the
feasibility sum and that the degree-mark bound dominates it.
Directed-Arb kernels and rational marks are checked against the same exact
moments; scalar conversion does not turn rational factors into floats.

At v=1/4, the following ratios compare new bounds with z^a(W_J(v)-1).
These are small GF4 examples, not evidence of a K20 certificate.

| Fixture | Feasibility-only ratio | Exact MDS ratio |
|---|---:|---:|
| One GRS[4,2] row; common cuts | 0.0208991 | 0.0197031 |
| One GRS[4,2] row; moving cuts | 0.00389509 | 0.00372109 |
| Two GRS[3,1] rows; changing order | 0.0806542 | 0.0806542 |

Positive degree marks 2 and 4 also improve the unfiltered bounds in these
fixtures. No K20 mark optimization or broad parameter search was run.

```text
python -B research/workstreams/char2_field_spin/theory/hd_mds_component_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_mds_component.py
```

## Why this does not close the whole first moment

The dense component investigation found a stronger obstruction than
infeasible MDS row supports. Its dimension argument applies to actual
codewords and every choice of nonzero prefix weights.

Fix a row set J. Where consecutive physical prefix cuts both belong to J,
the single coordinate between them is a unit free component. Select m such
coordinates, with h_c in region c. Their prefix differences are supported
only on J. Inverse prefix diagonals preserve that support.

Consider the routing event that these coordinates all land in the first L
positions of their respective regions. Its probability is exactly

    product_c C(L,h_c)/C(R,h_c).

Inside this m-dimensional input space, impose the two post-state exit
equations after each regional window. They cost at most 2b dimensions.
Impose the ACTUAL MDS parity checks on all rows of J. They cost at most
(b-k)|J| further dimensions, independently of the prefix coefficients.
The remaining dimension is at least

    D = m - 2b - (b-k)|J|.

The global post map is injective. Every nonzero surviving source therefore
produces a nonzero output supported in at most bL positions. If bL<=e and
D>0, the raw expected number of low-weight projective codewords is at least

    product_c C(L,h_c)/C(R,h_c) * (q^D-1)/(q-1).

The [independent fixed-profile calculation](HD_COMPONENT_DENSE.md) gives a
huge value for the two current placements, despite an astronomically small
routing-event probability. Its directed helpers are `hd_component_dense.py`
and `hd_component_dense_tests.py` in the theory directory.
Exact MDS membership cannot remove these witnesses: the dimension argument
already imposed it. A bound that counts every codeword can therefore fail
even when the corresponding bad setup event is negligible.

The next all-support proof needs to control rank defects or condition on
typical routes before counting witnesses. The new degree transfer remains
useful for sparse classes and for such conditional bounds. It should not
trigger a larger raw-first-moment sweep. Only the new files in this checkpoint
were changed; no native build, remote call, benchmark, or commit was made.
