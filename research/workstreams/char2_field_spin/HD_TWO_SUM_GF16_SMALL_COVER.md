# Small active-group counts under the Group16 outer

This screen broadens the 3.398-ms Group16 candidate beyond one matching
profile at r=17. It covers every source location for 1 through 17 active
groups, with either zero or one dense group. A dense group has 3 through
16 active RS constituents. The other groups have one or two.

The first FLOAT run passes all 320 cohorts at r=2 through 17, with union
log2 bound -563.3281943682605. Both sparse cohorts at r=1 also pass. The
remaining one-group dense cohort fails this bound at +2867.498534517523.
A positive log bound is uninformative, not a demonstrated low-weight word.
Fresh 256/384-bit replay now verifies all 322 passing cells at
-563.328194368258 bits, including both r=1 sparse cells. The
[directed record](HD_TWO_SUM_GF16_SMALL_CERTIFICATE.md) supplies exact
ceilings, arithmetic checks, and pins. The r=1 dense cell remains excluded
from that receipt. No whole-code guarantee follows from this result.

## One envelope for every fixed source support

Let r be the number of active outer groups. There are binomial(8192,r)
possible source supports. For a fixed regional support of size m and f
complete source pairs, the paired-route analysis supplies a nonnegative
two-state matrix R[m,f]. Define E[m] as its entrywise maximum over
0 <= f <= floor(m/2).

E[m] bounds every possible regional subset, including subsets of any
fixed r-element source support. Consequently, it applies under both
alternating matchings. This comparison loses matching information, but
does not resample the source support or assert independent regional masks.
All sixteen regional factors in an outer column use E. A bank through
occupancy 17 supplies the smaller-r calculations without changing their
source counts.

## Sparse and dense groups

For n1 singleton groups and n2 two-row groups, fixed row choices contribute
16^n1 120^n2 once. The sparse support polynomial counts each active symbol
with a factor 255. Its sixteenth power retains total support at least
9(n1+2n2), then divides by 255^[8(n1+2n2)]. These are the RS[16,8]
support floor and prescribed-support normalizer.

For a fixed dense row set A, let N_A(T) count ambient input vectors whose
mixed output support is T. Coordinates of these ambient inputs may vanish.
The shortened-RS comparison and the exact closed kernel give

    W(w)^16 = sum_{a=3}^{16} sum_{|A|=a} N_A(T)^16 / 255^(8a),
    where |T|=w.

The right side is independent of the positions in T. The implementation
forms this rational exactly, using common denominator 255^128. Hölder
bounds the sum over all 65,399 fixed row sets across all sixteen columns
by the product of W(|T_j|). No additional arity multiplicity is charged.
These weights are a comparison measure, not a sampled per-column row set.

The ordered weight-and-hole DP mixes the dense masks and the exact sparse
one-hole distribution before taking the common-potential quotient. An
additional variable marks nonempty dense columns. The sixteenth power
retains at least nine such columns, because every nonzero dense group
contains a nonzero RS row. This condition does not replace the separate
sparse support floor.

For h in {0,1}, the source and coloring prefix is

    binomial(8192,r) binomial(r,h) binomial(r-h,n1)
      * 16^n1 120^n2 / (255 p z^734003),

where r=n1+n2+h. The sparse row normalizer is applied separately as above;
dense row normalizers occur inside W. There is one projective factor and
one endpoint-potential factor for the complete word.

## Frozen witnesses and results

The screen tries three existing rational (z,log2 p) pairs. For each cohort,
it selects the smallest of the three bounds before summing cohorts.
Each witness bounds the same event, so this selection needs no union
factor for the number of witnesses. No witness was optimized in this run.

The source envelope and merged-arity weights are intentionally conservative.
The separate [two-dense-group control](HD_TWO_SUM_GF16_DENSE_TWO.md) shows
that merging arities can lose hundreds of bits: empty columns favor small
arity, while dense columns favor large arity. A row set cannot actually
change arity between columns. The separate
[r=1 certificate](HD_TWO_SUM_GF16_SINGLE_DENSE_CERTIFICATE.md) now retains
both arity and the full two-state matrix. It closes that excluded class
at 53.223053203418 bits, without changing the encoder. An independent
1024-bit sum with the 322 rounded-cell certificate gives the same rounded
upper bound for their disjoint union. That union covers all r=1..17
source supports with at most one dense group, not the whole code.

The screen took 4.894 seconds of reported computation. Five focused tests
check the exact combined moments, pointwise source-envelope domination,
the pair-only polynomial control, sparse support multiplicities, and the
once-only prefix. Two independent read-only reviews found no algebra or
counting issue before the run.

Reproduce from `theory` with an unused output filename:

```text
python -B -m unittest hd_two_sum_gf16_small_cover_tests
python -B hd_two_sum_gf16_small_cover.py --max-seconds 300 --output NEW_OUTPUT.json
```

The frozen screen is `hd_two_sum_gf16_small_cover_r1_17_allsource_float.json`.
Its source pin is
`3f6ade5acce2e005dbf65f8ecfbbb0fe7fe527c09eddfcca2581eb38595c52cd`.
The implementation and artifact remain research-only. Cohorts with two
or more dense groups and active-group counts above 17 are not covered here.
