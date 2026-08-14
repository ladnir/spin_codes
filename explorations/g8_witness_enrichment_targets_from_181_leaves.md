# Witness-enrichment targets from the 181-leaf checkpoint

The two-wave checkpoint has 181 active leaves.  Geometry refinement reduced
the worst diagnostic contribution, but the frontier remains far above the
target.  The next witness batch should target active leaf maxima while
retaining diversity across the original h2 cells.

## Acquisition score

For each active leaf `ell`, independently reconstruct its exact vertex hull.
Evaluate its current fixed rational selector at every vertex.  Define

```text
v_ell = canonical maximizing vertex,
E_ell = log2(n_ell)+U_ell.
```

Resolve a maximum tie by the canonical profile vector.  Rank leaves by
decreasing `E_ell`, then by leaf identifier.  This ranks the contribution
that the new witness bank must reduce, rather than only the per-profile score.

Deduplicate exact profile vectors.  If one profile maximizes several leaves,
assign it the log-sum-exp of those leaf contributions for acquisition ranking.

Exclude every exact profile already present in a frozen tuning-target source.
Measure proximity by transported mass

```text
d(a,b)=(1/2) sum_j |a_j-b_j|.
```

Treat candidates with `d<2048` as one neighborhood.  Retain the candidate
with larger acquisition contribution.  A lower-ranked near duplicate may
survive only if it covers another top-16 leaf.

Apply a cap of two active maxima per original h2 cell.  Within that cap, use
farthest-first selection under `d` after choosing the highest contribution.
This prevents one narrow descendant family from consuming the batch.

## Recommended core batch

The following 12 profiles are the best core batch.  Each is the active worst
vertex of the named leaf.  None duplicates a previous supplementary target.

| leaf | profile `(a0,...,a8)` |
| --- | --- |
| `h2:162/L/L` | `(131069,1,80565,1,1,1,1,50504,1)` |
| `h2:162/L/R` | `(131069,65534,15032,1,25073,1,1,25432,1)` |
| `h2:161/R/L` | `(196602,1,1,15032,1,21595,1,1,28910)` |
| `h2:158` | `(131069,1,1,1,1,65534,2,1,65534)` |
| `h2:160/R` | `(131069,1,65534,1,2,65534,1,1,1)` |
| `h2:004/R` | `(65534,1,1,1,1,1,2,65534,131069)` |
| `h2:163/R` | `(196602,2,1,1,65534,1,1,1,1)` |
| `h2:163/L` | `(177353,19251,1,1,65534,1,1,1,1)` |
| `h2:010/R` | `(28957,1,1,1,1,102112,1,1,131069)` |
| `h2:010/L` | `(28956,1,1,1,36579,2,1,1,196602)` |
| `h2:159/L` | `(131069,1,1,1,1,117424,1,13645,1)` |
| `h2:159/R` | `(131069,1,1,1,65534,51892,1,1,13644)` |

The closest core profile to an old target has transported distance 30,005.
Thus the core batch is not a local replay of the previous 32 targets.

If only eight tuning slots are available, use the first eight rows.  This
keeps the two-per-h2-cell cap and covers six original h2 cells.

## Optional transition-neighborhood batch

With 16 slots, add both endpoints from two transition neighborhoods that
created the current dominant descendants.

| transition | profile `(a0,...,a8)` |
| --- | --- |
| `h2:162/L`, left regime | `(196602,1,2,15031,1,1,50504,1,1)` |
| `h2:162/L`, right regime | `(196602,1,2,15031,1,50504,1,1,1)` |
| `h2:161/R`, left regime | `(196602,1,1,2,65534,1,1,1,1)` |
| `h2:161/R`, right regime | `(196602,1,1,65535,1,1,1,1,1)` |

The first pair straddles the `y_5` transition at 236,703.  The second pair
straddles the `y_3` transition at 211,632.  Tuning both sides tests whether a
new witness family can resolve the regime change instead of merely improving
one endpoint.  These four profiles are also exact-new relative to the prior
supplementary targets.

Do not replace a core profile by the transition profile nearest to it.  A
core profile targets the active fixed-mixture maximum.  A transition endpoint
tests whether the bank lacks a witness between two leader regimes.

## Tuning artifact

Each target row should bind:

- the 181-leaf checkpoint digest;
- leaf identifier or transition identifier;
- exact profile and its canonical digest;
- acquisition contribution and rank;
- nearest prior target and transported distance;
- target role: `active-worst` or `transition-endpoint`.

Tune one content-bound checkpoint per target.  Deduplicate any identical
tuned witness row before constructing the candidate atlas.  The candidate
atlas remains diagnostic until a manifest binds it.

## Replay after tuning

Freeze the 181-leaf geometry.  Reconstruct every exact hull and count.  Add
the candidate witnesses to the old bank, then recompute one rational minimax
selector per leaf under the same denominator and tolerances.

The enlarged bank should not worsen an exact binary64 LP optimum.  Allow only
the declared rationalization tolerance in the emitted selector score.  Report:

```text
old and new U_ell,
old and new E_ell,
new-witness selector weights,
global expanded log-sum-exp,
global collapsed maximum,
worst leaf identifier.
```

Declare the enrichment batch a discovery success only if all conditions hold.

1. At least 8 of the 12 core leaves improve by 32,768 bits.
2. At least 6 distinct new witnesses receive positive weight on a core leaf.
3. The worst frontier contribution decreases by at least 65,536 bits.
4. The global expanded diagnostic decreases by at least 65,536 bits.
5. No leaf worsens beyond the recorded rationalization tolerance.
6. The lower-quartile positive core improvement projects closure within eight
   further enrichment or refinement stages.

For condition 6, let `Delta_25` be that lower quartile and let
`G=max_ell(E_ell+40)` after replay.  Require

```text
Delta_25>0,                 ceil(G/Delta_25)<=8.
```

If conditions 1--4 fail, stop geometry expansion and revisit the witness
family or outer model.  If they pass, run one further finite checkpoint wave
with the enriched bank before selecting another target batch.

These gates are discovery criteria.  A theorem result still requires a new
manifest, fixed rational selectors, and independent outward evaluation at
every reconstructed leaf vertex.
