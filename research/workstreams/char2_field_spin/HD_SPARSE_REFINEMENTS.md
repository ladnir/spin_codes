# Bounded sparse-class refinements

This checkpoint records unsuccessful three-row bounds and two reusable
positive identities. It supplies no new complete distance certificate.
The one-row and two-row bounds in `HD_COMPONENT_RANK.md` and
`HD_COMPONENT_PAIRS.md` remain valid.

## More live components are not enough by themselves

For a permitted set of a rows, each region has a-1 interior components
and one new carried component. The interior lengths and final length are
positive; the initial carried length can be zero. Selecting k live
interior components has multiplicity C(a-1,k)*(q-1)^k.

The helper `hd_component_sparse.py` extends the pair matrix using the
minimum possible occupancy of each branch. It also tightens the
one-return zero-to-zero bound using the range of a uniform m-subset,
where m is the minimum positive occupancy. A larger uniform subset has
a stochastically larger range, so this replacement is an upper bound.

At K20 and a=3, marks with denominators 600, 800, 1000, 1200, 1400, 1450,
1500, 1550, 1600, 1800, 2000, 2500, and 3000 all fail. The smallest
tested log2 upper is about +339.01, at v=1399/1400. This is a loose
upper bound, not evidence that the construction fails with that probability.
The bound lets incompatible short-gap choices coexist across branches.

## Independent cyclic offsets: a changed ensemble

A possible construction change adds an independent uniform cyclic offset
of row groups in each region before the prefix. It retains within-group
order. The offset can be folded into the existing placement address, but
no native implementation or performance measurement was made.

Independent offsets allow regional component matrices to be averaged
before multiplying them. The carried component and post state remain in
the matrix frontier; neither state is reset. The test helper shifts and
reorders all cuts, including wraparound.

Write Q=q-1, p=1/Q, beta=v+(1-v)p. Let B_h be the no-return birth moment
for a uniform h-subset of R positions. Summing over occupancies gives

    sum_(h=1)^R B_h
      = v * sum_(l=0)^(R-1) beta^(R-1-l)*(R+1)/((l+1)*(l+2)).

The identity follows by summing the probability of first activation at l:

    sum_(h=1)^(R-l) C(R-l-1,h-1)/C(R,h)
      = (R+1)/((l+1)*(l+2)).

The sum of exactly-one-return birth moments is

    v^3*p^2 * sum_(l=0)^(R-3)
      beta^(R-3-l)*6*(R+1)*(R-l-2)/((l+2)*(l+3)*(l+4)).

Higher-return terms have the positive factorial tail used in the earlier
run-kernel proof. For unequal incoming/outgoing component flags, occupancy
changes injectively as one fixed circular gap crosses the region boundary.
For group size L, each occupancy is reached at most a times among R/L
offsets. Thus L*a/R times the summed positive kernel entries bounds their
mean. The helper also takes the minimum with the existing uniform bound.

For L=4, a=3, and denominators 800, 1000, 1200, and 1400, the tested
log2 uppers are approximately +494.79, +231.06, +107.64, and +277.50.
The changed ensemble therefore does not pass this bounded gate. These
results do not justify a native implementation or promotion.

Five tests verify both positive identities, the minimum-range bounds,
and all three-cut configurations at R=4,8,12 for two marks. They compare
the shifted matrix with explicit averages over all group-aligned offsets.

```text
python -B research/workstreams/char2_field_spin/theory/hd_component_sparse_tests.py
python -B research/workstreams/char2_field_spin/theory/hd_component_sparse.py
python -B research/workstreams/char2_field_spin/theory/hd_component_shift.py
```

The [geometry refinement](HD_COMPONENT_GEOMETRY.md) separately preserves
the frozen affine orbit for within-group triples. Its failed bound is
dominated by the unbinned fixed profile J=(0,1,2). At v=1149/1150, the
directed upper for that profile is 2^18.09317428. Applying the total-degree
MDS filter from `HD_MDS_COMPONENT.md` with lambda=9/8,5/4,3/2,2,3,4,8
increases the log2 upper to 33.85,47.96,72.39,110.96,165.35,203.94,296.94.
This bounded check does not justify a larger degree-mark sweep. It does
not evaluate the stronger per-row degree filter.

The [actual-code first-moment obstruction](HD_COMPONENT_DENSE.md) now
rules out closing the existing ensemble by source-count refinements alone.
The next priority is a capped routing/rank argument, as started in
[the conditioning gate](HD_ROUTING_CONDITION.md), not a larger mark sweep.
