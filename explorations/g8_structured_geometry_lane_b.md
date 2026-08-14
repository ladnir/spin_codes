# Structured geometry for the `g=8` profile simplex

The certificate needs cells with exact ownership, few vertices, and tractable
counts.  The witness bound is convex on every fixed cell.  Thus an outward
check at all cell vertices bounds every profile in that cell.

This note compares four coordinate systems.  The cumulative-mass grid gives
the best scaling while preserving the verifier's integer split interface.
All artifacts in this lane bind `G8_SUPPORT_MANIFEST.json` at SHA-256
`cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616`.

## Exact cumulative coordinates

Fix an exact support

```text
S={s_0<...<s_d}.
```

For each active class, write `a_{s_i}=1+x_i`.  Define

```text
T=M-(d+1),
y_k=sum_{i=0}^k x_i       for 0<=k<d.
```

The map from `a` to `y` is a unimodular bijection between integer profiles
on `S` and integer chains

```text
0<=y_0<=...<=y_{d-1}<=T.
```

Choose `B=2^h` balanced integer intervals that partition `{0,...,T}`.  Each
`y_k` has one interval owner.  Because the chain is nondecreasing, the tuple
of interval indices is also nondecreasing.  Hence the grid has exactly

```text
C(B+d-1,d)
```

nonempty cells.  Independent coordinate boxes would have `B^d` cells.

The ownership test uses canonical prefix splits.  If interval `m` ends at
`U_m`, then

```text
y_k<=U_m
```

is equivalent to

```text
sum_{i=0}^k a_{s_i}<=U_m+k+1.
```

The split has zero inactive coefficients.  Its pivot coefficient is zero,
and its active coefficient gcd is one.  Therefore, the split already uses
the verifier's support-canonical form.

## Vertices and counts

Consider one nondecreasing interval tuple.  Equal consecutive interval
indices form runs.  A run of length `r` in interval `[L,U]` is the chain
polytope

```text
L<=z_1<=...<=z_r<=U.
```

The run has `r+1` vertices:

```text
(L,...,L,U,...,U).
```

Different runs use disjoint ordered intervals.  Therefore, one cell is a
product of these run polytopes.  Its vertex count is

```text
prod_run (r+1)<=2^d.
```

All vertices are integral before the optional physical-weight intersection.
For full support, the base weight is `0+...+8=36`.  The physical-weight cut at
21 is redundant, so every full-support cell has the explicit product form.

The same run decomposition gives the exact integer count.  An interval of
integer width `w` contributes

```text
C(w+r-1,r)
```

for a run of length `r`.  The cell count is the product over runs.  If class
zero is active and the base weight is below 21, subtract the weight-below-21
profiles owned by the cell.  The global excluded set has only 2,099 profiles.
The prototype enumerates this set exactly.  It uses the verifier's rational
intersection code for clipped supports of dimensions one and two.

Collapsed aggregation needs only the exact root count and complete cell
coverage.  Expanded aggregation may use the product count for every cell.

## Comparison

| Geometry | Full-support cells | Vertices per cell | Counting | Main limitation |
| --- | ---: | ---: | --- | --- |
| Braid chambers | `9! = 362880` | 9 | Partition counts with tie ownership | Large fixed cell count |
| Raw dyadic boxes | `B^8` | at most 256 | Inclusion-exclusion | Sum constraint couples boxes |
| Dyadic stick boxes | `B^8` | at most 256 | Nested `O(dM)` dynamic program per cell | Poor uniform scaling |
| Cumulative boxes | `C(B+7,8)` | at most 256 | Product over equal-bin runs | Natural class order is fixed |

Braid chambers have the fewest vertices per cell.  However, stable tie
ownership complicates counts, and class-dependent witnesses prevent a simple
quotient by permutations.

Stick-breaking splits have exact integer ownership.  Their ratio constraints
are triangular, but a uniform depth `h` creates `2^{8h}` cells.  Cumulative
coordinates exploit the chain order and remove most of these cells.

For full support, cumulative-grid scaling is:

| Levels `h` | Bins `B` | Nonempty cells | Maximum vertices | Total vertex incidences |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 2 | 9 | 25 | 165 |
| 2 | 4 | 165 | 81 | 6,435 |
| 3 | 8 | 6,435 | 256 | 490,314 |
| 4 | 16 | 490,314 | 256 | 61,523,748 |

Levels one and two are small enough for discovery.  Level three remains
plausible only if most cells reuse a small fixed witness bank.  Level four is
not a sensible global starting grid.  Level two is the practical first
full-support grid: it needs 165 ownership cells and 6,435 vertex inequalities.

## Proof conditions

A cumulative-grid certificate is sound if all following conditions hold.

1. The balanced intervals partition every integer in `{0,...,T}` exactly.
2. Every profile maps to one nondecreasing interval tuple.
3. Every emitted prefix split uses the support-canonical integer form.
4. The verifier reconstructs every cell vertex from exact constraints.
5. One fixed witness or rational mixture passes every vertex of its cell.
6. Exact cell counts sum to the manifest's exact root count when expanded
   aggregation is used.
7. Collapsed aggregation uses the manifest root count and complete cell
   coverage.

The structural prototype proves no witness inequality.  It emits no shard
completion claim.

## First bounded run

Start with the clipped dimension-two support `{0,1,2}` and four bins.  The
command creates ten cells and checks exact count conservation.

```powershell
python scripts/prototype_packet_group_g8_cumulative_boxes.py `
  --manifest G8_SUPPORT_MANIFEST.json --support-mask 0x007 --levels 2 `
  --output out/g8_cumulative_boxes_007_l2.json
```

If the ten-cell artifact passes review, run full support at level one.  That
run has nine cells and at most 25 vertices per cell.

## Full-support h=2 diagnostic

The bounded h=2 run used all 510 frozen base witnesses and 32 supplementary
witnesses.  It evaluated one fixed singleton and one rational minimax mixture
on every exact vertex of each cell.  The run found no binary64 pass among the
165 cells.  Every mixture improved its cell's best singleton.

The best mixture cell has bin tuple `(0,0,0,0,3,3,3,3)`.  Its upper score is
`768190.7452450581`, and its selector has seven components.  Its maximizing
vertex is

```text
(1,65534,1,1,131070,1,65534,1,1).
```

The worst mixture cell has bin tuple `(0,0,0,0,0,0,1,1)`.  Its upper score is
`1454007.7514306216`, and its selector has eight components.  Its maximizing
vertex is

```text
(1,1,1,1,65534,1,2,65534,131069).
```

This cell contributes `1454125.2595775214` to the log2 union diagnostic.  The
full mixture aggregate has the same displayed value.  Thus, this cell is the
first exact refinement target.

Do not refine the entire grid to h=3 first.  That grid needs 6,435 minimax
problems and 490,314 vertex incidences.  Instead, intersect the worst cell
with its local h=3 grid.  Split `y_0,...,y_5` at `32766` and `y_6,y_7` at
`98300`.  In profile coordinates, the canonical thresholds are

```text
32767, 32768, 32769, 32770, 32771, 32772, 98307, 98308.
```

Chain monotonicity leaves 21 nonempty exact children: seven refinements of
the length-six run and three refinements of the length-two run.  The existing
run-product formula counts each child exactly.  Reoptimize one fixed rational
mixture on each child before refining another h=2 cell.

The diagnostic artifact is
`out/g8_full_support_cumulative_h2_diagnostic.json`.  Its SHA256 digest is
`47308693e73546e249f21c1a2edef4a820c4017820b6d9f8f670e95d720151b3`.
Its binary64 scores still require independent outward replay.  The artifact
is not a certificate or a probability bound.
