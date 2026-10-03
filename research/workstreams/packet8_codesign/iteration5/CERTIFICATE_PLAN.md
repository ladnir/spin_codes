# Outward certificate plan for the wider actual24 construction

This retained plan has now been implemented. The [completed checkpoint](README.md)
records the 68.8937-bit whole-code certificate and independent global replay.
The text below preserves the design and arithmetic audit that preceded them.

The positive g4 pipeline can certify the existing construction without a new
local census model. Rebase the chosen upper operators using exact rational
mixture weights. Then evaluate chronological block products, fractional
powers, and placement with outward positive arithmetic.

This note specifies a verifier design, not a completed certificate. The
construction and conditioning argument are those audited in
[`AUDIT_WIDER.md`](../iteration4/AUDIT_WIDER.md). They use 256 outer groups,
64 byte regions, 32 physical steps per region, and one continuous 24-bit
state. The cutoff is 13107. No arithmetic below resets state between regions.

## Local interface and arithmetic audit

Fix an exact dyadic weight parameter z in (0,1). For active occupancy k,
let K_k be the actual positive kernel on inner states, weighted by
z to the emitted Hamming weight. Let E be the matrix whose rows are the
normalized measures delta0, U, P1,...,P8. Here U is uniform on nonzero states;
P_i is the genuine weighted nonzero birth measure for active occupancy i,
divided by its genuine total mass. The required upper matrix satisfies

```
E K_k <= U_k E                         (pointwise on output states).
```

The local producer supplies finite nonnegative binary64 entries U_k,
interpreted as exact dyadic numbers. Only row zero may have birth columns.
Its U column is zero. For k=0, the zero row is delta0 and every nonzero row
has an exact zero in the zero-target column. These structural zeros must
survive conversion and thinning.

The read-only audit of `local_outward24.py` found no mathematical blocker.
Its profile polynomials are exact integers over explicit denominators.
Uniform means and single-packet births are exact before endpoint conversion.
For larger occupancies, the normalized Walsh transform has error at most

```
delta_input + gamma_b * max_abs_rounded_input,
gamma_b = b / (2^53-b),                 b=24.
```

Each butterfly contributes one rounded sum or difference per stage. Negation
is exact. The common dyadic lattice survives cancellation, and the checked
lattice exponent excludes nonzero subnormal intermediates after scaling.
The magnitude check excludes overflow, including positive profile sums.
Projection onto nonnegative values cannot increase error from the true masses.
The subsequent gamma_(S+chunks) allowance conservatively covers bincount and
chunk accumulation for S=2^24 states.

Birth normalization uses the exact total input moment minus an enclosed
zero-syndrome mass. The full-rank restrictions through three packets justify
setting that mass to zero for k<=3. Dividing enclosed numerators by the
opposite positive denominator endpoints therefore bounds means in the
genuine normalized P_i families. The returned lower matrix encloses the
filled comparison matrix; it is not a pointwise lower physical kernel.

Five exact GF4-cubed tests passed independently. This audit accepts the
producer's stated IEEE binary64 and NumPy operation contract. The final
receipt must record that contract and authenticate the producer and census
sources; the tests alone do not establish a platform-independent arithmetic
theorem.

## Exact rebase of selected upper operators

Potential occupancy j first gives the operator

```
V_j = sum_{k=0}^j C(j,k) 255^k / 256^j * U_k.
```

Compute each coefficient outward and regard the resulting dyadic matrix V_j
as fixed. Its coefficients may exceed the true ones. That does not invalidate
a subsequent change of normalized comparison measures.

For j=1,...,8, define exact rational numbers from these selected dyadics:

```
c_jk = V_j[0,P_k],     b_j = sum_k c_jk,
Q[B_j,P_k] = c_jk/b_j.
```

Set Q's zero and uniform rows to their identity rows; all other entries are
zero. Require b_j>0. Then Q is nonnegative, Q1=1, and e0 Q=e0 exactly.
The new measures E'=Q E are normalized mixtures of genuine P_k. They need
not equal the genuine potential-j birth distributions.

Define W_j before rounding as follows:

* Its zero row has V_j[0,0] at zero and b_j at B_j when j>0.
* Every nonzero row has only zero and U targets, with coefficients from Q V_j.
* For j=0, the zero row remains exactly delta0.

The structural form gives the exact identity W_j Q=Q V_j. Therefore

```
E' K_j^potential <= Q V_j E = W_j E'.
```

Rounding W_j upward preserves this inequality. Equality need not survive
that final rounding. In particular, there is no requirement that the selected
upper birth masses approximate true masses from below.

Do not round the rows of Q independently and then assert Q1=1. Instead,
evaluate every new birth-row coefficient directly as

```
W_j[B_i,t] = (sum_k c_ik V_j[P_k,t]) / b_i,       t in {0,U},
```

using exact dyadic sums/products and exact rational division followed by an
upward dyadic conversion. The denominator b_i is the exact sum defining Q,
not a separately inflated bound. This is inexpensive for nine 10-by-10
matrices. No inverse of Q is used; linearly dependent mixture rows are valid.

This route is preferable to rebuilding genuine potential birth families.
The latter is also sound with outward numerator/denominator intervals, but
it introduces unnecessary dependence on additional true normalization data.

## Extract a known occupancy factor

Let a=((1+z)/2)^8, defined as an exact rational. Replace each selected
rebased operator by W_j/a^j, evaluating the division upward. Every fixed-q
regional occupancy sequence has total potential count q. Across 64 regions,
the extracted factor is consequently a^(64q), independent of hidden states.
After the fractional exponent alpha, the scalar prefactor is

```
[ beta^q * z^(-13107) * a^(64q) ]^alpha,
beta = 2^512/(2^32-1)^8.
```

Exact per-occupancy scaling commutes with the exact rebase, because it cancels
from each normalized mixture row. Upward rounding can alter the selected
mixtures, but the preceding domination argument still applies. This scaling
reduces dynamic range; it does not remove the need for outward arithmetic.

## Fractional powers without trusting a library root

Write alpha=p/d in lowest terms, with 0<p<=d. The selected denominators are
small. For an exact positive rational x, a proposed dyadic y is a certified
upper endpoint for x^alpha precisely when

```
y >= 0 and y^d >= x^p.
```

Check this condition with integer powers and cross multiplication. A library
root may propose y, but cannot certify it. Raise a failed proposal until the
integer condition holds. Handle x=0 exactly. This procedure also certifies
the scalar prefactor, which is rational before the fractional power.

For the millions of block entries, use a positive tangent table instead.
At an exact t>0, certify u>=t^alpha once. Concavity gives the globally valid
upper function

```
x^alpha <= (1-alpha)u + (alpha*u/t)x.
```

Both coefficients are nonnegative and can be converted upward. Any table
entry is valid for any x>=0; selecting a nearby entry affects tightness only.
Preserve x=0 explicitly to avoid adding a spurious constant to structural
zeros.

For scale independence, write x=m*2^e with 1<=m<2. Let ep=kd+r with
0<=r<d, including when e is negative. Evaluate an upper tangent at m, multiply
by an integer-checked upper endpoint for 2^(r/d), then scale by 2^k outward.
Use 256 dyadic mantissa bins or a denser table if a small replay needs it.
The table can record exact root checks and coefficient inequalities once.
There is no need to certify an approximate logarithm or exponential.

## Positive products, sums, and underflow

A verifier must bound the actual operation sequence, not add one nextafter
to an arbitrary library dot product. For a nonnegative n-term dot evaluated
with ordinary binary64 products and additions, a conservative gamma_(2n)
allowance covers either sequential or tree reduction. If all intermediates
are normal, the exact dot is at most the computed dot divided by
1-gamma_(2n).

With gradual underflow, let tau=2^-1074. A sufficient bound is

```
exact_dot <= (computed_dot + n*tau)/(1-gamma_(2n)).
```

The n*tau term covers at most 2n-1 rounded operations, each with absolute
underflow error at most tau/2. Evaluate the correction and inflation upward.
Reject nonfinite intermediates. This formula assumes round-to-nearest,
gradual underflow, and no unaccounted reassociation or extra arithmetic.
FTZ/DAZ must be disabled and checked; otherwise use an independently justified
fallback. Normal-input tests do not establish absence of FTZ/DAZ.

Determine structural zeros from Boolean operand support. A computed zero
with a positive product term may be underflow and must receive the correction.
Conversely, an exactly empty support stays zero. This preserves the local
operator structure required for rebasing.

Use exact powers of two to rescale matrix or row-vector intermediates, carrying
an integer exponent. If scaling itself loses a positive entry, account for
that loss outward or fall back to exact/MPFR/Arb arithmetic. Positive summation
over fine tuples needs its own operation-count allowance. Exact combinatorial
ratios must be rounded upward; floating factorials and log-gamma values are
not certificate inputs.

## Global verifier and acceptance condition

For each selected pair (z,alpha), the verifier performs these steps:

1. Authenticate maps, local receipts, arithmetic contract, and dyadic endpoints.
2. Thin active operators and rebase the chosen upper family exactly as above.
3. Optionally extract a^j, retaining its exact scalar contribution.
4. Multiply each ordered four-step tuple before taking entrywise alpha powers.
5. Average the powered tuples with weights prod C(8,j_i)/C(32,sum j_i).
6. Compose eight macros with exact hypergeometric placement for 256 regional slots.
7. Compose all 64 regional matrices, starting at zero and summing every terminal state.
8. Multiply by C(256,q) and the certified scalar prefactor for each q>=2.

The fine-count weights stay outside alpha. The outer subset count C(256,q)
also stays outside alpha. The beta and cutoff factors are inside alpha.
Only the final row vector is collapsed to a scalar; collapsing each regional
matrix would incorrectly reset state.

For q=1, use unpowered active operators and exact outer shell coefficients.
Compose the degree-one regional placement and the 64-region support
polynomial with positive outward arithmetic. Divide each support coefficient
by C(64,v), apply z^(-13107), and cap at one. A different certified z may be
chosen for each support size v. Multiply by the exact shell counts and 256.

Every q=1,...,256 must have an authenticated upper endpoint. Taking the
minimum among independently valid witnesses for a fixed q is valid. Sum the
selected dyadic endpoints exactly, then compare that sum with 2^-40 using
integers. The final Boolean result must follow from this comparison, not a
rounded display of the margin in bits.

The next implementation step is a small positive-arithmetic verifier plus one
complete witness replay. Check alpha=1 against a separately evaluated ordinary
placement, exercise exact-zero and underflow fixtures, and compare a few
entries with a high-precision outward backend. Only then replay all selected
witnesses and the exact occupancy union. No new construction search or full
local census model is needed for this plan.

## Independent tiny checks

`test_certificate_plan.py` tests exact rebasing and all short products,
including dependent mixture rows. It also tests occupancy scaling, integer-
checked roots, positive tangents, and an underflow case that defeats a purely
relative error allowance. These checks support the algebraic design; they
do not constitute a whole-code certificate.
