# Compressing persistent group types with Hölder's inequality

Row conditioning permits useful packet laws, but different groups can have
different laws. Their types persist across all 256 regions. Resampling a
type independently in each region is not an equivalent experiment.

This note gives a product-measure upper bound that avoids retaining every
original type in the placement calculation. Several such bounds can serve
as the coarse types in [typed placement](typed_placement.py). The argument
provides a coverage mechanism, not a completed numerical certificate.

## The weighted measure of one active group

Write m = 256 for the number of regions and L = 2048 for the number of
four-row groups. Partition the admissible positive BCH row weights into
disjoint intervals. Add a separate label for a row fixed to zero. For an
active interval I and a fixed Bernoulli witness p, define

    beta(w,p) = binomial(m,w) p^w (1-p)^(m-w),
    Gamma_I(p) = max_{w in I, Abar_w > 0} Abar_w / beta(w,p).

Here Abar is the authenticated coefficientwise BCH spectrum cap. An
interval with no positive capped mass can be omitted. A witness must give
positive mass to every permitted weight. Thus p = 1 is allowed for a row
type supported only at weight 256. A fixed zero row instead uses p = 0
and scalar factor 1.

A group type t is an unordered multiset of four row labels, excluding the
all-zero multiset. Fix its row witnesses once. Let n_I count its copies of
active interval I and let n_0 count its zero rows. Define

    h_t = 4! / (n_0! product_I n_I!)
            * product_{active rows i} Gamma_{I_i}(p_i),
    theta_t(b) = [x^b] product_{rows i} (1-p_i+p_i x),  b = 0,...,4.

The factorial factor counts distinct row-label assignments. It appears
once, inside h_t. Each theta_t is a probability distribution on packet
weights. Conditional on weight b, the independent lane permutation gives
a uniform mask among the binomial(4,b) possible masks.

Let mu_t denote the setup-averaged counting measure of all messages in
type t, restricted to this group. Pointwise row domination gives

    mu_t <= h_t theta_t^(tensor m).

Both sides include the same conditional uniform lane masks. Equivalently,
the displayed inequality can be read on packet-weight trajectories after
averaging any nonnegative downstream function over those masks. The
reference allows weight-zero packets, including an entirely zero
reference group. This is harmless domination of a nonzero-message class;
theta_t(0) must not be discarded.

The actual all-zero group type is excluded. An inactive group is kept
identically zero outside this active-group measure. This distinction
preserves the meaning of the occupancy q.

## One common product measure

For a finite collection T of active group types, define

    nu_b = (sum_{t in T} h_t theta_t(b)^m)^(1/m),  b = 0,...,4.

For every packet-weight trajectory (b_1,...,b_m), generalized Hölder gives

    sum_t h_t product_{c=1}^m theta_t(b_c)
      <= product_{c=1}^m (sum_t h_t theta_t(b_c)^m)^(1/m)
       = product_{c=1}^m nu_{b_c}.

To apply Hölder, use the m sequences
`h_t^(1/m) theta_t(b_c)`, indexed by t, with all exponents equal to m.
The type is fixed throughout every term on the left. The product measure
on the right dominates this persistent-type mixture; it is not obtained
by asserting that types are resampled.

Set z = sum_b nu_b and phi_b = nu_b/z, omitting an empty family with z = 0.
Then the full averaged measure of one active group is at most

    z^m phi^(tensor m).

Apply this domination separately to each of q selected groups. Let M_q(phi)
be the region transfer operator for q selected groups with packet law phi.
Selected groups whose reference packet is zero still occupy their selected
positions; their zero packets can subsequently be marginalized as in
[the homogeneous driver](row_verify.py). The resulting occupancy bound is

    binomial(L,q) z^(mq) exp(lambda H)
       e_zero M_q(phi)^m terminal,

Here lambda > 0, H = 209715, e_zero initializes the zero inner state, and
terminal sums the mass coordinates without double-counting density or
tail coordinates. The transfer operators bound the tilted output moment.
This includes all original type assignments within T. No extra
factor |T|^q, row-position factor, or Gamma factor belongs outside z.

The product envelope is coordinatewise minimal for the relaxed measure
`sum_t h_t theta_t^(tensor m)`. Indeed, the constant trajectory
`(b,...,b)` requires any dominating product measure eta to satisfy

    eta_b^m >= sum_t h_t theta_t(b)^m.

This minimality concerns the relaxed reference mixture, not the true BCH
counting measure. Better bounds can still come from different witnesses,
tighter row classes, several families, or retained inter-region structure.

## Several weighted families

Choose nonnegative allocations alpha_(c,t) to a finite set of families c.
Require sum_c alpha_(c,t) >= 1 for every active type t being covered. A
disjoint partition uses indicator allocations. Define

    nu_(c,b) = (sum_t alpha_(c,t) h_t theta_t(b)^m)^(1/m),
    z_c = sum_b nu_(c,b),
    phi_(c,b) = nu_(c,b) / z_c.

Discard families of zero mass. Apply the preceding argument separately
to each weighted family. The complete active-group measure then obeys

    sum_t mu_t <= sum_c z_c^m phi_c^(tensor m).

Expand the product of these sums for q groups. For family counts n_c with
sum_c n_c = q, let M_(n_c) be the region operator from typed placement.
Its selected groups retain their family laws across all m regions. Its
local operators average zero packets as part of those laws. The
contribution for these counts is at most

    binomial(L,q) * q! / product_c n_c!
      * product_c z_c^(m n_c)
      * exp(lambda H) e_zero M_(n_c)^m terminal.

Sum this expression over every family-count vector. Group locations are
counted by binomial(L,q); the multinomial assigns families to those
locations. Original row labels and original types were already absorbed
into the family measures.

Two families need only a two-type placement calculation. Full occupancy
coverage additionally requires that the row-label partition cover every
possible nonzero group, that every family-count vector be included, and
that the sum be bounded outward. Selected homogeneous examples do not
establish these conditions or the desired numerical bound.

For outward implementation, first round each unnormalized nu_(c,b)
upward to an exact dyadic value. Sum and normalize those dyadics exactly
to obtain a rational law phi_c. Increasing the unnormalized coordinates
preserves pointwise domination. The existing exact-rational packet-law
interface can then be used without treating an interval midpoint as a
probability certificate.

`holder_exact.py` implements this normalization. Its root proposals are
checked by exact rational powers, including the preceding dyadic grid
point. Six tests cover every short trajectory for several toy families,
degenerate laws, and the constant-trajectory constraints. It accepts exact
coefficients, not binary64 estimates of h_t.

## Scalar partition diagnostic

The binary64 pilot uses row categories `{0}`, [38,62], [64,192], [194,218],
and `{256}`. These cover every weight with positive authenticated BCH cap.
There are 69 nonzero unordered four-row types. The all-middle type has
product-envelope cost `256 log2(z) = 516.448` bits per group.

Putting all exceptional types in one family raises that scalar cost to
about 619.28 bits. Even a two-family low/high split costs about 535 bits
per family. These coarse partitions discard too much information.

A finer partition first fixes the number of exceptional rows. Within each
count, classify the sum of row scores as negative, zero, or positive. The
five row categories receive scores -2,-1,0,1,2. Omitting empty classes gives
11 exceptional families, in addition to the all-middle family. Optimizing
their reference probabilities gives at least 47.49 bits per group of
scalar advantage over the all-middle cost; several families give more.

These are binary64 costs of product envelopes, not distance margins or
outward certificates. Their packet laws differ. Recombining these families
into one would lose the benefit again; the inner analysis must preserve
the family information or pay a justified additional bound.

```sh
python -B research/workstreams/permutation_locality/independent_rows/holder_pilot.py --optimize
python -B research/workstreams/permutation_locality/independent_rows/test_holder_exact.py
```

`holder_verify.py` converts the 69 types to exact coefficients and replays
a selected exceptional family against an all-middle background. At witness
p = 1/4, its outward scalar costs are 516.448286505 for the middle type
and 469.925117095 for `exactly_1_low_strict` (values rounded upward).
This exceptional family includes both one-zero-row and one-light-row types.
The envelope-only check is:

```sh
python -B research/workstreams/permutation_locality/independent_rows/holder_verify.py --envelope-only --family exactly_1_low_strict --p 1/4
```

The full replay uses exact rational packet laws in `typed_inner.py`, with
the scalar and multinomial factors derived above. At q = 64, tilt .052,
and p = 1/4, selected mixtures now have outward bounds:

| All-middle groups | `exactly_1_low_strict` groups | Margin, rounded down |
|---:|---:|---:|
| 63 | 1 | 77.901375 |
| 62 | 2 | 83.497137 |

These bounds include every original type choice within the exceptional
family and all group/row locations. For example, the second row includes
two light-row groups, two zero-row groups, and one of each. It does not
cover high-row families or other exceptional-group counts.

```sh
python -B research/workstreams/permutation_locality/independent_rows/holder_verify.py --groups 64 --exception-counts 1 2 --family exactly_1_low_strict --p 1/4 --tilt .052 --joint-cancellation
```

Independent review found no factor error. Three further tests check the
exact sum of type masses, the disjoint twelve-family partition, and
deterministic endpoint types.

## Two shortcuts that do not follow

Persistent types cannot simply be resampled. For one group and two
regions, take two packet laws `(1/4,3/4)` and `(3/4,1/4)`, each with
unnormalized weight 1. The persistent-type measure of trajectory `(0,0)`
is `1/16 + 9/16 = 5/8`. Resampling the average law in each region, while
preserving total mass 2, gives only `2(1/2)^2 = 1/2`. Hölder gives `5/8`
and therefore handles this example correctly.

Nonnegative transfer operators alone do not show that homogeneous types
are worst. For two selected groups, consider scalar region operators

    B_AA = B_BB = 1/4,    B_AB = B_BA = 1.

Over two regions, homogeneous assignments have moment 1/16, whereas a
heterogeneous assignment has moment 1. These are nonnegative operators
bounded by 1. Thus the proposed extremal shortcut does not follow from
the abstract transfer interface. This example is not a counterexample
to an additional IMT-specific theorem; such a theorem would require
structure beyond nonnegativity and exchangeability.
