# Structured geometry for the full-support `g=8` shard

The full-support root has dimension eight and dominates the profile census.
An effective refinement must keep three properties simultaneously:

1. every integer profile has one owner;
2. every proof cell has few exact rational vertices;
3. the verifier can derive a count or an independently valid count upper bound.

This note compares four structured decompositions.  A fixed-denominator
stick-breaking chain gives the best contract.  Its proof cells have at most
256 vertices at every depth.

## Common root

Fix a permutation `pi=(pi_0,...,pi_8)` of the classes.  For a full-support
profile, define

```text
b_i = a_i-1,                 b_i >= 0,
T = M-9,                     sum_i b_i = T.
```

The physical-weight constraint is redundant on this root.  All decompositions
below therefore act on the shifted simplex

```text
Delta_T={b in R_{>=0}^9 : sum_i b_i=T}.
```

An ownership cell is a set of integer profiles.  A proof cell is a closed
rational polytope that contains those profiles.  Proof cells may overlap.
Only ownership cells enter count aggregation.

## Candidate comparison

| decomposition | vertices before arbitrary cuts | exact ownership | count method | decisive limitation |
| --- | ---: | --- | --- | --- |
| nested prefix bands | at most 21,728 | integer cumulative bounds | `O(8T)` prefix-sum DP | substantially more vertices |
| ordered chambers | 9 per chamber | stable sorting with integer strictness | restricted-partition DP | 362,880 root chambers |
| fixed-denominator stick chain | at most 256 | primitive integer ratio cuts | `O(8T)` range-sum DP | requires separate ownership and proof constraints |
| balanced recursive joins | at most 256 | primitive integer ratio cuts | tree convolution | exact counting can cost quadratic time |

Arbitrary affine cuts destroy the stated vertex bounds.  A producer must
remain inside the selected structured family until it explicitly changes to
the generic BSP contract.

## Nested prefix bands

Define cumulative variables

```text
y_k=sum_(h=0)^k b_(pi_h),       0<=k<=7.
```

A prefix-band cell fixes integer intervals `L_k<=y_k<=U_k` and also has

```text
0<=y_0<=... <=y_7<=T.
```

A split `y_k<=t` versus `y_k>=t+1` partitions the integer points exactly.
Refinement replaces one endpoint, so it introduces no new facet direction.

Every vertex partitions `y_0,...,y_7` into maximal constant blocks.  Each
block must contain a coordinate fixed at one interval endpoint.  Otherwise,
that block admits a local translation.  For a composition
`8=s_1+...+s_m`, there are at most `prod_r 2s_r` anchored block choices.
Consequently the number of vertices is at most

```text
V_0=1,        V_d=sum_(s=1)^d 2s V_(d-s),        V_8=21728.
```

The integer count is the number of nondecreasing integer sequences within
the eight intervals.  Prefix sums compute this count in `O(8T)` integer
additions.  Coordinate-box inclusion--exclusion also gives an independent
upper bound.

This family is sound and countable.  Its vertex bound is too large for the
first full-support design.

## Ordered chambers

For a permutation `pi`, the closed chamber is

```text
b_(pi_0)<=...<=b_(pi_8).
```

Define gaps `z_0=b_(pi_0)` and
`z_k=b_(pi_k)-b_(pi_(k-1))` for `k>0`.  Then

```text
sum_(k=0)^8 (9-k) z_k=T,        z_k>=0.
```

Thus every closed chamber is an eight-simplex with nine vertices.  To obtain
unique integer ownership, sort the pairs `(b_i,i)` lexicographically.  An
adjacent inversion of class indices requires a gap of at least one.  After
subtracting these fixed gaps, the chamber count is a coefficient of

```text
product_(w=1)^9 (1-x^w)^(-1).
```

The coefficient is exact and admits an `O(9T)` dynamic program.  However,
the root has `9!=362880` owned chambers.  Witness-dependent affine refinement
also turns each simplex into a general polytope whose vertex count grows with
the number of cuts.  Chamber symmetry cannot identify cells because witness
charges depend on the class labels.

Ordered chambers are useful for diagnostics, not for the primary shard.

## Fixed-denominator stick-breaking chain

Fix one dyadic denominator `Q=2^k` in the manifest, with `1<=k<=52`.  Define

```text
R_i=sum_(h=i)^8 b_(pi_h),       0<=i<=8,
R_0=T,                          b_(pi_i)=R_i-R_(i+1).
```

A split chooses a stage `i<8` and a numerator `p` with `0<p<Q`.  Define the
integer affine form

```text
A_(i,p)(b)=Q*b_(pi_i)-p*R_i.
```

The ownership children are

```text
left:  A_(i,p)(b)<=0,
right: A_(i,p)(b)>=1.
```

These sets are disjoint and contain every integer parent profile.  Reusing a
single denominator is essential.  It lets a path retain only the strongest
lower and upper numerator at each stage.

The proof domain relaxes every right constraint by one:

```text
right proof constraint: A_(i,p)(b)>=0.
```

This relaxation contains the owned right child.  It can overlap the left
proof domain only on `A_(i,p)=0`.  The verifier must never use proof-domain
overlap for counting or ownership.

For each stage, the resulting proof constraints have the form

```text
alpha_i R_i <= R_(i+1) <= beta_i R_i.
```

At a vertex with `R_i>0`, one endpoint constraint is active.  If `R_i=0`, all
later residual masses are zero.  Choosing one endpoint at each nonzero stage
therefore generates every vertex.  The proof cell has at most

```text
2^8=256
```

exact rational vertices, independently of tree depth.

The exact owned count also has a linear-time dynamic program.  Let `c_i(r)`
count prefixes that reach `R_i=r`.  The path bounds determine one integer
interval

```text
ell_i(r) <= R_(i+1) <= u_i(r).
```

For an upper ratio numerator `p_hi`,

```text
ell_i(r)=ceil((Q-p_hi)r/Q).
```

For a strict lower ratio numerator `p_lo`,

```text
u_i(r)=floor(((Q-p_lo)r-1)/Q).
```

Missing bounds use `0<=R_(i+1)<=R_i`.  Range additions followed by one prefix
sum compute `c_(i+1)` from `c_i` in `O(T)` operations.  Summing `c_8(r)` gives
the leaf count.  The recurrence uses only exact integers.

Exact leaf counts are optional for collapsed aggregation.  A verifier may
instead recompute the root census or a coordinate-box count as an upper
bound.  Expanded aggregation requires verifier-recomputed stick-chain
counts.

## Support-recursive joins

Let a rooted binary tree have the nine classes as leaves.  For each internal
node `v`, define `m_v` as the sum of `b_i` over its leaves.  Constrain the
left-child mass by a fixed-denominator ratio interval relative to `m_v`.

After the same one-unit proof relaxation, each of the eight internal choices
contributes two homogeneous endpoint constraints.  Hence the proof cell has
at most 256 vertices.  Exact ownership again uses the unrelaxed integer gap.

For exact counts, define

```text
C_v(m)=sum_n C_left(n) C_right(m-n),
```

where `n` ranges over the owned interval at node `v`.  A balanced tree needs
general restricted convolutions.  Their direct cost can be quadratic in
`T`.  The stick chain is the special case in which one child is a leaf, so
each convolution becomes a range sum.

Balanced joins offer no vertex advantage over the chain and complicate
independent counting.  They should remain a fallback for a demonstrated
witness-alignment benefit.

## Recommended producer and verifier contract

Use `stick-chain-dyadic-v1` for the first structured full-support producer.
The immutable root record contains `pi`, `Q`, and `T`.  Each split contains
`stage`, `numerator`, its canonical nine-coordinate integer affine form, and
the canonicalization child-swap bit.

The producer must satisfy these conditions.

1. Use the manifest's single `pi` and `Q` at every node.
2. Emit only numerators strictly inside the current owned interval.
3. Canonicalize the affine form under the full-support gauge.
4. Reject children without an exact integer progress witness or positive
   verifier-recomputable count.
5. Mark binary64 witness comparisons as discovery diagnostics only.

The verifier independently performs these checks.

1. Reconstruct `A_(i,p)` from `pi`, `Q`, `stage`, and `numerator`.
2. Canonicalize it and compare every serialized integer field.
3. Route ownership with `A<=0` and `A>=1`.
4. Build proof geometry with `A<=0` and `A>=0`.
5. Retain only the strongest lower and upper numerator per stage.
6. Enumerate at most 256 exact rational endpoint vertices.
7. Recompute the declared count method or replace it by an independent upper
   bound.
8. Outward-evaluate one fixed selector at every proof vertex.

The current generic BSP verifier uses the exact right constraint `A>=1` for
both ownership and proof geometry.  It must not claim the 256-vertex bound.
Add structured verifier support only after the manifest versions the new
geometry rule and declares that proof cells are ownership relaxations.

## Fatal traps

- Changing `Q` below a node prevents one-bound-per-stage normalization.
- Using the exact `A>=1` constraint for proof geometry loses the 256-vertex
  guarantee.
- Counting overlapping proof cells double-counts threshold profiles.
- Treating a binary64 ratio as conservative makes ownership unsound.
- Adding a witness-difference cut silently changes the cell family to a
  generic BSP.
- Accepting a nonempty relaxed proof cell does not prove that its owned
  integer cell is nonempty.
- Using producer counts without replay invalidates expanded aggregation.
- Reordering `pi` below the root breaks the chain count recurrence and the
  immutable ownership contract.
