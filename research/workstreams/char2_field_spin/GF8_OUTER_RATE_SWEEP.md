# Outer-rate sweep with the eight-state inner fixed

Increasing the outer dimension from 87 to 88 gives a new complete K18
certificate at 10% symbol distance with 72 bits of setup-failure margin.
The requested rate rises from 2/3 to 32/47, approximately 68.085%.
Output length falls by 2.083%. This is a proof and size improvement;
no new encoding time was measured.

The same change remains unclosed at K20. Dimension 89 also has substantial
middle-class gaps under the current bounds. Thus the sweep does not justify
a large rate increase with the existing inner. These are limitations of
the tested bounds, not proofs that the corresponding codes have low distance.

## Construction and changed parameter

Keep GF256, the 128-coordinate outer length, one regional routing stage,
and the width-16/eight-state symmetric Cauchy inner. Their maps and setup
distribution are specified in the [existing construction](FIELD128_K20_CERTIFICATE.md#exact-construction-and-setup-distribution).
Only the outer dimension d and its natural padding change.

For requested dimension K, put

\[
L=64\left\lceil\frac{K}{64d}\right\rceil,\qquad
K_{\rm pad}=dL,\qquad N=128L.
\]

Each outer row is the evaluation of a degree-below-d polynomial at the
128 byte labels. The fixed maps are unchanged as d varies. Setup samples
independent uniform regional permutations and independent uniform nonzero
GF256 input and state diagonals, then fixes them for every message.
The bad event is a nonzero padded source with output weight at most
`e=ceil(N/10)-1`.

Using additive-polynomial coefficients instead of systematic source values
preserves the full outer code image for every d. Zeroing unused source
coordinates gives a subcode, so both source interfaces inherit the full-space
distance guarantee. The 128-bit implementation applies the scalar code to
sixteen byte lanes; a nonzero lane witnesses record support, without an
additional union over lanes.

## New whole-code K18 certificate

For d=88 and K=262144, the geometry is L=3008, K_pad=264704, and N=385024.
The nominal rate is 11/16; after padding the requested rate is 32/47.
The forbidden threshold is 38502, so the certified minimum distance is
38503, at least 10% of N.

The [aggregate replay](theory/gf8_outer_rate_certificate.py) recomputes the
following components. No saved numerical receipt is trusted for acceptance.

| Component | Complete range | Margin, bits |
|---|---|---:|
| Sparse active rows | A=1 through 375 | 72 |
| Middle inverse-input weights | W=15401 through 92405 | 129 |
| Dense inverse-input weights | W=92406 through 385024 | 15031 |

Every nonzero GRS[128,88] row has weight at least 41. Outside the sparse
range, W>=41*376=15416, which lies inside the middle range. The middle
and dense intervals are adjacent. These facts cover every nonzero source.

The aggregator matches the construction and geometry of all receipts,
checks the range joins, and sums their actual probability bounds as exact
fractions. The upward-rounded total is

`1.2528267213880597358333991877892946739763396629303542176730365764542430696110719E-22`,

which is below 2^-72. Probability is over the ideal setup distribution
above. This proves the whole-code statement, not just selected weights.
The proof uses the existing two-state transfer; the new rank-history
control is not needed for this result.

The [sparse](GF8_OUTER_RATE_SPARSE.md), [middle](GF8_OUTER_RATE_MIDDLE.md),
and [dense](GF8_OUTER_RATE_DENSE.md) notes specify their arithmetic, tests,
and complete coverage. The aggregate has six tests, including a full replay
and rejection of mismatched geometry, missing coverage, and insufficient
probability sums.

```text
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_certificate.py
python -B research/workstreams/char2_field_spin/theory/gf8_outer_rate_certificate_tests.py
```

## Rate frontier and finite-size effects

| Size | Existing d=87 output length | d=88 output length | Status of d=88 |
|---|---:|---:|---|
| K16 | 98304 | 98304 | No output-size reduction; not separately certified in this sweep |
| K18 | 393216 | 385024 | Complete 10% certificate, 72 bits |
| K20 | 1548288 | 1531904 | Middle/dense closed; sparse trough remains unclosed |

At K20, the d=88 two-state sparse screen bottoms near A=256 with about
4.59 optimized bits. The independently audited rank-history control raises
the directed bound at the frozen marks to 23 bits, and to 30 bits at the
optimized two-state marks. Neither reaches 40. These are pointwise class
bounds, not a 30-bit whole-code certificate.

At d=89, refined K20 middle screens give approximately -7295 and -7391
bits at densities .13 and .14. Raising the sparse cutoff does not remove
these middle-weight gaps. Dimensions 90 and 92 have larger sampled gaps;
the sweep does not claim an exhaustive impossibility result for higher rates.

At K16, d=87 through 93 all retain L=768 and N=98304. There is no reduction
in routed data or inner work, although outer scheduling can still differ.
The first smaller output occurs at d=94, with N=90112 and requested rate
8/11. That profile fails the tested middle bounds by thousands of bits.

## Engineering decision

Preserve the new d=88 K18 certificate as an option. The 2.083% reduction in
routed/output data is useful but small. It is not an end-to-end speedup
measurement. In particular, the previously measured .244/.956 ms timings
belong to the d=87 coefficient-basis implementation, not this new profile.
No production parameter or native kernel changed during the sweep.

The next substantial gain needs a cheaper or stronger inner, or a materially
tighter middle-class analysis. Merely changing d does not currently support
the 3/4-rate target. Do not spend a large implementation campaign on the
small d=88 size gain before choosing the next construction direction.
Existing certificates and benchmarks remain intact; no commit, push, or
native benchmark was performed for this sweep.
