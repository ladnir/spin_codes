# Packet-8 drive-stratified inner transfer

## Scope

This note keeps the frozen packet-8 Riffle construction unchanged.  A local
step has a uniformly permuted incoming 64-bit state support, eight binary
packets, the fixed 64-bit accumulator, and the systematic `[128,64,22]` EBCH
parity map.  The purpose is to replace the old global point-mass cap by a
rigorous cap conditional on the accumulator-drive weight.

The formulas below are finite combinatorics.  The current probe evaluates
them in binary64; the proof certificate must evaluate the same formulas with
outward error control.

## Exact local masses and point caps

Let `R` be the incoming state support, uniformly distributed among the
`q`-subsets of 64 coordinates.  Let `U` be the current 64-bit input, grouped
as eight ordered bytes, and give it the nonnegative weight

```
f(U) = product_i f_[wt(U_i)].
```

Write

```
W = U xor R,       V = Acc(W),
d = wt(W),         y = wt(V),       r = wt(PV),
```

where `Acc` is the fixed accumulator bijection and `P` is the systematic EBCH
parity map.  Define the exact drive/emission mass

```
H[q,d,y] = (1 / C(64,q))
           sum_(R,U) f(U) 1[wt(R)=q, wt(W)=d, wt(Acc(W))=y].
```

For a fixed drive word `W`, define

```
a_q(W) = (1 / C(64,q)) sum_[R:wt(R)=q] f(W xor R).
```

For each drive weight, use the point cap

```
c[q,d] = max_[W:wt(W)=d] a_q(W).
```

Both tables are exact finite objects.  `H` is computed by an eight-byte
transfer that tracks incoming accumulator parity, outgoing parity, selected
state weight, drive weight, and emitted weight.  To compute `c`, observe that
the one-byte polynomial for a drive byte depends only on the byte's Hamming
weight.  Consequently a 64-bit drive is represented, for this maximization,
by a multiset of eight numbers in `{0,...,8}`.  Enumerating the
`C(16,8)=12870` multisets and multiplying their one-byte polynomials finds the
maximum exactly; no sampling or independence assumption is used.

## Shared EBCH-column inequality

Fix `q` and `y`, and let

```
x[d,r] = sum_[W:wt(W)=d, wt(Acc(W))=y, wt(P Acc(W))=r] a_q(W).
```

Then

```
sum_r x[d,r] = H[q,d,y].                         (1)
```

Let `B[y,r]` be an upper bound on the number of systematic EBCH words whose
data half has weight `y` and parity half has weight `r`.  Because `Acc` is a
bijection, the words indexed by different drive weights are disjoint.  Since
`a_q(W) <= c[q,d]`,

```
sum_d x[d,r] / c[q,d] <= B[y,r].                 (2)
```

Equation (2) is the key improvement over the earlier drive-stratified probe:
one EBCH split column is shared by every drive stratum, rather than reused
independently for each `d`.

## Rank-one transportation bound

For a nonnegative next-state test vector `v`, maximize

```
sum_(d,r) x[d,r] v[r]
```

subject to (1), (2), and nonnegativity.  Put

```
n[d,r] = x[d,r] / c[q,d].
```

The source supply is `H[q,d,y]/c[q,d]`, the sink capacity is `B[y,r]`, and
the profit of one transported unit is the rank-one product `c[q,d] v[r]`.
The rearrangement inequality (or the elementary two-by-two exchange
argument) therefore maximizes the objective by pairing source caps in
nonincreasing order with test-vector entries in nonincreasing order.  A
two-pointer greedy fill evaluates this maximum, denoted `Phi[q,y](v)`.

Define the positive homogeneous, monotone upper operator

```
F(v)[q] = sum_y z^y Phi[q,y](v),       0 < z < 1.
```

If `A` is the exact averaged one-step transfer, the preceding inequalities
give `A v <= F(v)` componentwise for every nonnegative `v`.  Monotonicity then
gives `A^B 1 <= F^B 1` for all block counts `B`.

For any strictly positive vector `v`, if an outward-certified number
`lambda` satisfies

```
F(v) <= lambda v,
```

and `gamma >= max_q 1/v[q]`, then the inner MGF from initial state zero is at
most

```
gamma lambda^B v[0].
```

This is the nonlinear Collatz witness used by the probe.  It is a direct
finite inequality; convergence of the power iteration is irrelevant once
`v`, `lambda`, and `gamma` are frozen and checked outward.

## Hard-face diagnostic

At the existing sharp hard-face witness

```
profile = (0,12483,12483,12483,0,34328,34328,34329,121710)
z       = 0.3,
```

the shared-drive operator reduces the fixed inner MGF bound by
`81113.8020` bits, or `2.47540` bits per inner group.  Since every outer term,
profile charge, and orbit normalization is held fixed, the former combined
value `+23165.9309` becomes `-57947.8711` in binary64.  The per-profile target
is `-168.7010`, leaving `57779.1701` bits of diagnostic room.

This large margin is why the proof-hardening step should certify one frozen
Collatz vector with deliberately coarse outward guards instead of attempting
to certify an optimized eigenvalue.

## Remaining-face audit

An exhaustive audit of all 198580 nested-prefix triangle barycenters found
5322 residuals under the previous atlas.  Reusing the hard-face witness closes
49 of them and leaves 5273.  This is expected support dependence: the hard
witness assigns near-zero fugacity to its inactive packet classes, so it is a
poor certificate when those classes become active.

Retuning without changing the transfer lemma closes the worst residual tested:
the face `100<1c0<1f9` moves from the old-atlas value `+90563.70` to a robust
retuned value `-4038.25`, then to `-67516.36` after the shared-drive
improvement.  A structurally different sparse residual `100<101<1a5` moves to
`-164678.16`.  These are binary64 atlas diagnostics, not outward face-cover
certificates.  They identify the next finite task precisely: retune shared-
drive witnesses over the residual faces, construct a convex finite cover, and
outward-check only the selected witnesses and cover inequalities.
