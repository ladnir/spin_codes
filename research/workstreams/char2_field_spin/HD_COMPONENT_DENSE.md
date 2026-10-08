# Dense classes and an obstruction to the codeword first moment

The entire dense tail of the component expansion, comprising permitted-row
classes 29943 through 32768, is below 2^-322.58291595 at K20. An information
set gives this bound without assumptions on the fixed placement geometry.

The middle range has a different problem. For both individual-row placement
and the native four-row-group placement, rare routes contain large subspaces
of actual low-weight codewords. Their expected number of scalar lines exceeds
2^420595 and 2^441457, respectively. This rules out an uncapped codeword
first-moment proof, even with exact outer-code membership and class-specific
Chernoff marks. It does not rule out a 40-bit bound on setup failure.

## Construction and probability space

Use R=32768 outer rows, a GRS[64,32] code in each row, and the field F of
order q=65536. Write b=64, rho=32, Q=q-1, K=32R, and N=64R. The forbidden
weight cutoff is e=floor(2N/5)=838860.

A fixed placement orders the rows within each of the b coordinate regions.
One input-weighted prefix scan follows this placement. The route then applies
an independent uniform permutation within each region. Finally, the post-inner
has two additive states, initialized to zero, with global carry across regions.

For the four-product post-inner, each position uses nonzero c,d1,d2,f with
c+d1+f*d2 nonzero, and performs

    y = s1 + f*s2 + c*x,
    s1' = s1 + d1*x,
    s2' = s2 + d2*x.

The three-product variant sets c=d1, samples d1,d2,f independently from F*,
and therefore computes y=s1'+f*s2 using the old second state. Both maps are
triangular with a nonzero input-to-output diagonal. Both are injective.

The obstruction below is pointwise in every admissible prefix and post
coefficient setup. Its only randomness is the independent uniform regional
route. Averaging over either coefficient ensemble preserves the lower bound.
The uniform route is a mathematical ensemble; no claim about a finite PRNG
seed distribution is made here.

The dense-tail upper additionally uses the established source law from
[HD_AFFINE_PROOF.md](HD_AFFINE_PROOF.md): independent nonzero prefix weights
give MDS membership probability at most z=Q^(-rho) for each nonzero source
row. Its probability space therefore averages over the prefix diagonal.

## A rare route with a large actual-code kernel

We first give a finite linear-algebra argument. Fix a permitted source-row
set J of size a. Let pi_c(j) be the source row at physical prefix position j
of region c. Select the positions

    H_c = {j : 0 <= j < R-1, pi_c(j) in J, pi_c(j+1) in J},
    h_c = |H_c|,          m = sum_c h_c.

Every selected coordinate gives a unit vector in the prefix-output space.
Its inverse prefix difference is supported only at j and j+1. Both source
rows belong to J. Multiplication by the inverse nonzero prefix weights does
not change that support. These m unit vectors are independent. Excluding
j=R-1 avoids cross-region boundary terms; the prefix itself is not reset.

Fix a window of length L at the start of every output region, with h_c<=L.
Let E be the event that every selected position in H_c routes into that
region's window. Independent uniform regional permutations give exactly

    Pr_route[E] = product_c C(L,h_c)/C(R,h_c).

Fix any route in E and any admissible coefficient setup. On the m-dimensional
selected input space, impose two post-state exit equations per region. There
are at most 2b equations. Their kernel has zero output outside the windows:
the first region starts at zero, and each zero exit state persists through
the remaining zero-input positions. This argument also supplies zero initial
states for the following region.

Now impose the actual GRS parity checks. Source rows outside J already vanish.
The rows in J contribute at most rho*a further linear equations. Consequently
the intersection has dimension at least

    d = m - 2b - rho*a.

All maps after the source placement are invertible. Thus this intersection
gives at least (q^d-1)/Q distinct scalar lines of actual codewords, each of
weight at most bL. The argument uses the actual MDS equations, not the
relaxation that assigns probability z to every nonzero row.

Let Z_e be the number of actual output-code lines of weight at most e. When
d>=2 and bL<=e, the preceding argument proves

    E[Z_e] >= Pr_route[E] * (q^d-1)/Q
            > Pr_route[E] * q^(d-2).

This inequality holds for either post ensemble and for every prefix diagonal.
It remains valid if all coefficient randomness is averaged afterward.

## The two frozen K20 witnesses

The helper fixes a=19660 and L=13107, so bL=838848<e. It constructs J with
Python's Random(0).sample(range(32768),19660). The resulting row-mask digest is

    eab8a631796d80d096c309dff05dca4985e3947b732733e264d89126749cd436

The helper checks this digest before using the witness. It also freezes and
checks every one of the 64 unit-component counts, rather than relying on
an average degree calculation.

Both placements use the 64 fixed odd strides in
`theory/hd_two_state_regional_strides.py`. Their physical row orders are

    individual: pi_c(j) = stride_c*j mod R,
    group four: pi_c(j) = 4*((stride_c*floor(j/4)) mod (R/4)) + (j mod 4).

These orders invert the outer store addresses in
`spin/experiments/char2_field_spin/hd_affine_two_state.cpp`. In particular, the
second formula is the native four-row-group construction, not the original
individual-row construction under a different name.

| Quantity | Individual rows | Four-row groups |
| --- | ---: | ---: |
| Minimum h_c | 11693 | 11771 |
| Maximum h_c | 11940 | 11867 |
| m | 754354 | 756089 |
| Actual-code dimension lower d | 125106 | 126841 |
| log2 Pr_route[E] | -1581068.941247536344 | -1587966.187611903450 |
| log2 E[Z_e] lower | 420595.058752463655 | 441457.812388096549 |

The logarithms are directed Arb evaluations at 160-bit precision. Integer
counts and the dimension inequality are exact. The script verifies the
conservative statements E[Z_e]>2^420595 and E[Z_e]>2^441457.

These rare events prove actual distance failure on their occurrence. Their
probabilities, about 2^-1581069 and 2^-1587966, are far below 2^-40. Therefore
the result does not demonstrate failure of the desired distance guarantee.
It demonstrates that many codewords on rare bad setups make the untruncated
expectation unsuitable for proving that guarantee.

Exact source feasibility cannot remove this obstruction: the counted
subspaces already satisfy all source parity checks. An exact-active-row
partition may distribute their words among subsets of J. We do not identify
the witness with the exact-active class a=19660.

The positive permitted-row expansion has an additional direct obstruction
in its a=19660 class. Before imposing the MDS equations, the same event gives
at least q^(m-2b)-1 nonzero assignments in the permitted component space.
For every mark 0<v<=1, each low-weight assignment contributes at least one
after multiplying by v^(-e). Moreover,

    Q^(-rho*a) >= q^(-rho*a),
    (1-Q^(-rho))^(R-a) > 1/2.

Thus this one permitted-J term exceeds half the conservative first-moment
lower bound in the table, for every mark. Choosing a separate mark for each
class, or even for each J, cannot repair that permitted-class inequality.
This is stronger than the common-mark obstruction in
[HD_COMPONENT_RANK.md](HD_COMPONENT_RANK.md).

## A complete dense-tail upper from an information set

The dense tail can still be bounded without summing component profiles.
Fix a permitted set J of size a and fix the route and post setup. The
component output space has dimension m_a=ba, because it is the invertible
image of the source space supported on J.

Choose an information set of m_a output coordinates. Projection onto that
set is injective and cannot increase weight. The number of nonzero vectors
of output weight at most e is therefore at most

    sum_{w=1}^e C(m_a,w) Q^w <= e*C(m_a,e)*Q^e,

where the second inequality assumes e<=m_a/2. Summing over J and multiplying
by z^a/Q bounds the class contribution by

    T_a = e*C(R,a)*C(ba,e)*Q^(e-1-rho*a).

For the positive permitted-row expansion, we have dropped the factor
(1-z)^(R-a)<=1. For the exact-active partition, each word in that class has
source membership probability at most z^a, and the permitted component space
contains the exact-active space. These are distinct justifications for the
same upper bound; permitted and exact-active classes are not equated.

Let A satisfy e<=bA/2. For every a>=A, the ratio of consecutive bounds is at
most

    r_A = (R-A)/(A+1) * Q^(-rho)
          * ((bA+1)/(bA-e+1))^b.

Indeed,

    C(m+b,e)/C(m,e) = product_{i=1}^b (m+i)/(m-e+i)
                    <= ((m+1)/(m-e+1))^b.

Both ratio factors decrease with a. If r_A<1, the entire tail is at most
T_A/(1-r_A). At A=29943, a directed 160-bit computation gives

    log2(sum_{a=29943}^{32768} class-upper_a) < -322.5829159528.

This covers all 2826 classes, beginning at 91.37878418% of the outer rows.
It is uniform over every fixed regionwise placement and every injective post
setup. It is not a complete distance certificate: the middle classes cannot
be filled by the uncapped first moment just analyzed.

## Reproduction and next proof step

The new files are `theory/hd_component_dense.py` and
`theory/hd_component_dense_tests.py`. The following commands reproduce the
receipts and eight tests:

```text
python -B research/workstreams/char2_field_spin/theory/hd_component_dense.py
python -B research/workstreams/char2_field_spin/theory/hd_component_dense_tests.py
```

The tests check the native address inverse, unit-vector support, and exact
regional containment probabilities. Small GF4 instances also construct actual
codewords with nontrivial prefix weights, both placements, and both post laws.
They check every source parity equation and every zero-state window exit.
The remaining tests cover the information-set bound, geometric tail ratio,
and both directed K20 receipts. All eight pass. No native build, remote call,
performance benchmark, or existing-file edit was performed for this branch.

The next proof should bound the probability of bad routing or the tail of
a rank defect, rather than charge every line in a large bad kernel. A
conditioning event must cover the relevant family of source sets J and
output windows, including their union cost. Removing only this fixed witness
does not address that family. The selected unit components and one window
per region are only one mechanism; nonunit components and general low-weight
output sets remain outside this argument. Sparse exact-source filtering is
still useful, but it cannot close the full uncapped moment by itself.
