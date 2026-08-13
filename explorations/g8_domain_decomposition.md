# Domain decomposition for the `g=8` certificate

The `g=8` profile is

```text
a=(a0,...,a8),       sum_j a_j=M=262144,
sum_j j*a_j >= 21.
```

The certificate must cover these integer profiles without constructing a
global eight-dimensional triangulation.  The proposed decomposition uses
exact support masks, integer slabs, and local lattice BSP trees.  A fixed
witness or fixed rational mixture certifies each terminal proof polytope.

## Exact support census

For an integer profile `a`, define `supp(a)={j:a_j>0}`.  Every nonempty
support is feasible except `{0}`.  Hence there are

```text
2^9 - 2 = 510
```

feasible exact-support strata.  A support of size `s` has dimension `s-1`.

| dimension | support size | feasible strata | feasible integer profiles |
| ---: | ---: | ---: | ---: |
| 0 | 1 | 8 | 8 |
| 1 | 2 | 36 | 9,437,096 |
| 2 | 3 | 84 | 2,886,184,992,415 |
| 3 | 4 | 126 | 378,293,710,105,607,109 |
| 4 | 5 | 126 | 24,791,478,291,771,024,604,728 |
| 5 | 6 | 84 | 866,508,443,723,541,949,201,514,342 |
| 6 | 7 | 36 | 16,224,627,887,200,131,391,413,327,495,204 |
| 7 | 8 | 9 | 151,895,545,730,963,601,519,675,551,057,510,391 |
| 8 | 9 | 1 | 553,017,927,440,720,481,221,689,864,611,271,442,433 |

The first three rows contain 128 strata and can use low-dimensional methods.
Only one stratum has dimension eight.  Nevertheless, that full-support
stratum contains about `99.972538%` of all feasible profiles.  Profile count
alone cannot dismiss sparse strata because their per-profile bounds can be
larger.

Here is the exact count within one support.  Put `s=|S|` and write

```text
a_j=1+b_j for j in S,       sum_j b_j=M-s.
```

If `0` is not in `S`, then

```text
|Lambda_S| = C(M-1,s-1).
```

If `S={0} union P`, define

```text
h_P = 21-sum_(j in P) j,
E(P,h) = sum_(t=0)^(h-1) [z^t] product_(j in P) (1-z^j)^(-1),
```

with `E(P,h)=0` for `h<=0`.  Then

```text
|Lambda_S| = C(M-1,|P|)-E(P,h_P).
```

Across support sizes `s=1,...,9`, the physical-weight cut excludes

```text
1, 52, 437, 957, 582, 70, 0, 0, 0
```

profiles.  The first exclusion is the infeasible profile supported on `{0}`.
Thus the exact feasible count is

```text
L8 = C(M+8,8)-2099
   = 553169839211945865258921061892182603726,
log2(L8) = 128.70099010348517...
```

The generic `profile_count(8,N)` helper returns `C(M+8,8)`.  That value is a
safe upper bound, but a final support ledger can subtract the 2,099 excluded
profiles exactly.  The uniform per-profile target based on `L8` is about
`-168.7009901035` bits.  The final weighted ledger need not impose that target
on every profile.

## Certificate tree

Each support mask is an independent root.  A root domain is the shifted
simplex, optionally intersected with the physical-weight halfspace.  The
certificate processes a root through three layers.

1. **Root witness.** Try one support-eligible witness or one fixed rational
   mixture on the complete root polytope.
2. **Laminar slabs.** Split difficult roots by integer band masses, such as
   `a0+...+a3` and `a4+...+a8`.  Refine a band only after its parent fails.
3. **Local lattice BSP.** Split a residual node by a primitive integer affine
   form.  Use coordinate splits first.  Use witness-dominance cuts only where
   they materially reduce the certified contribution.

For an affine form `A(a)` with integer coefficients and an integer threshold
`t`, the children are

```text
left:  A(a)<=t,             right: A(a)>=t+1.
```

These closed proof domains contain a disjoint partition of the parent's
integer profiles.  The open real slab between them contains no integer
profile.  Rational discovery cuts must be scaled to primitive integer
coefficients before they enter an artifact.

A terminal node stores its exact rational constraints.  The verifier
reconstructs the vertices from those constraints.  It then evaluates one
fixed witness or one fixed rational mixture at every vertex.  Convexity
bounds the witness throughout the terminal polytope.  Floating-point vertex
sets or Qhull incidence claims are not certificate inputs.

Zero-fugacity eligibility is checked at the support root.  A witness whose
class-`j` fugacity is zero may occur only when every profile in the node has
`a_j=0`.  Refinement within a positive-support root cannot make such a
witness eligible.

## Ownership and counts

Ownership is procedural and exact.

1. `supp(a)` selects one support root.
2. At each tree node, the integer comparison selects one child.
3. The resulting root-to-leaf bit string is the unique leaf identifier.

No least-cell rule or floating barycentric test is needed.  A profile on an
affine threshold goes left.  Worker artifacts may overlap as real polytopes,
but their owned integer sets cannot overlap.

Coordinate and laminar-band trees should carry exact counts.  For a
coordinate box on `s` active classes, let `l_j<=a_j<=u_j`, put
`T=M-sum_j l_j`, and put `d_j=u_j-l_j`.  Inclusion-exclusion gives

```text
B(l,u) = sum_(R subset S) (-1)^|R|
         C(T-sum_(j in R)(d_j+1)+s-1,s-1),
```

where the binomial is zero when its upper argument is too small.  Subtract
the globally enumerated weight-at-most-20 profiles that lie in the box.
Because the split tree is disjoint, exact child counts must sum to the exact
parent count.

Laminar band-mass counts are also exact.  For a band with `q` active classes
and positive mass `r`, its unconstrained count is `C(r-1,q-1)`.  Products and
finite convolutions of these terms count nested band intervals.  The same
2,099-profile table handles the physical-weight exclusion.

For an arbitrary affine BSP leaf, the verifier may instead use the minimum
of independently valid upper bounds:

```text
parent count, coordinate-box count, and verified unimodular-strip counts.
```

Two aggregation modes are valid for a refined node `v`.

```text
expanded:   sum_(leaf ell below v) n_ell * 2^U_ell,
collapsed:  n_v * 2^(max_(leaf ell below v) U_ell).
```

Here `U_ell` is an outward upper bound for the leaf's per-profile term.
Expanded mode uses owned-leaf counts.  Collapsed mode needs only complete
leaf coverage and a valid count for the parent.  The `g=4` BSP verifier uses
the collapsed rule for each replaced root.  A `g=8` verifier can compute both
valid bounds and retain the smaller one.

The final ledger is an outward log-sum-exp over the chosen root or leaf
terms.  Witness-mixture components and BSP branches are proof devices; they
are not additional probabilistic events.

## Divide-and-conquer artifacts

Use one immutable manifest for parameters, support-count formulas, witness
source hashes, and arithmetic conventions.  Shard geometry by support mask.
Each shard contains a complete tree for its mask and no references to another
mask's local identifiers.

A practical work split is:

- dimensions 0--2: 128 direct low-dimensional strata;
- dimensions 3--5: 336 adaptive support shards;
- dimensions 6--7: 45 dense boundary shards;
- dimension 8: one dedicated full-support shard.

Discovery jobs may propose cuts and mixtures independently.  The verifier
must rebuild every count, vertex, ownership decision, and outward inequality.
Aggregation begins only after every one of the 510 feasible masks has one
complete shard.  Long performance measurements remain serialized; support
sharding does not require concurrent benchmarks.

## Staged stopping criteria

Every node has exactly one state: `EMPTY`, `CERTIFIED_LEAF`, `SPLIT`, or
`UNRESOLVED`.  A stage may freeze a shard only when it has no unresolved
node.

1. **Census gate.** Verify all 511 nonempty masks, reject `{0}`, recover the
   dimension table, and check that support counts sum to `L8`.
2. **Root gate.** Attempt one outward root certificate per feasible support.
   Freeze successful roots.  Rank failures by `log2(n)+U`, not by their worst
   per-profile value alone.
3. **Slab gate.** Apply exact laminar or coordinate splits to the dominant
   failures.  Freeze a node only after one eligible fixed witness or mixture
   passes all exact vertices.
4. **BSP gate.** Use local affine cuts only for the remaining contribution.
   A terminal singleton fallback must enumerate and certify every owned
   profile; one residual centroid is never sufficient.
5. **Integration gate.** Require zero unresolved nodes, exact support
   coverage, valid count conservation where claimed, and an outward global
   upper endpoint at most `-40`.

For planning, a complete diagnostic near `-50` bits supplies ten bits for
outward hardening and ledger overhead.  It is not a theorem result.  The only
final stopping condition is a complete independently replayed ledger at or
below `-40` bits.

This architecture confines eight-dimensional work to the full-support shard
and its local residual nodes.  It never constructs a global eight-dimensional
mesh, and it preserves the exact ownership and cell-local accounting used by
the successful `g=4` certificate.
