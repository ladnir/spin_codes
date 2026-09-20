# Collision-surcharge attempt

## Question

The exact adjacent-pair calculation refutes entrywise domination by the iid
lane model.  The remaining question is whether each multi-active group can be
charged a small scalar surcharge that composes across groups.

For a profile with (a) active blocks in (m) active groups, define the defect
count

\[
 D:=a-m.
\]

Each defect is an active block beyond the first block in its group.

## Universal likelihood bound

Fix all active-lane patterns.  In one transposed region, the exact profile has
(q^m) possible shift tuples.  The iid model has (q^a) labeled lane tuples.
Therefore, every nonnegative entry satisfies the likelihood-ratio bound

\[
 R_{\mathrm{profile}}(z)
 \le q^D R_a^{\mathrm{iid}}(z).
\]

Fresh shifts are independent across 256 regions.  The complete moment can
therefore lose (256D\log_2q) bits.  For (q=32), this charge is (1280D)
bits.  Even one defect consumes more than the low-occupation margin.  The
universal likelihood bound is rigorous but unusable.

## Exact one-defect calculation

Fix one shared group with pattern (\{0,d\}).  Let (s:=a-2) singleton groups
choose independent lanes.  For a fixed shift of the shared pair, an
exponential generating function sums the singleton occupations exactly.  The
calculation then averages all 32 shifts and raises the region transfer to the
256th power.

The following table reports the largest scalar surcharge over all pair
separations.  Each row uses the Chernoff tilt selected by the iid calculation
for the same total occupation.

| Total occupation (a) | Worst separation | Scalar surcharge |
|---:|---:|---:|
| 2 | 1 | 0.054578 bits |
| 4 | 1 | -0.037505 bits |
| 8 | 1 | -0.034907 bits |
| 16 | 1 | -0.027177 bits |
| 32 | 1 | 0.104204 bits |
| 48 | 1 | -0.011678 bits |
| 64 | 1 | -0.007661 bits |

Negative values mean the shared pair improves the bound.  The adjacent pair
is worst in every tested row.

## Two defects and denser groups

The calculation also places two independent adjacent-pair groups in the same
word.  At occupation 32, one pair costs 0.104204 bits and two pairs cost
0.209581 bits.  The other tested occupations improve under both one and two
pairs.  This behavior is consistent with an additive scalar surcharge.

Selected denser patterns give the same qualitative result.  At occupation 32,
the scalar surcharges are:

| Shared-group pattern | Group weight | Surcharge | Surcharge per defect |
|---|---:|---:|---:|
| contiguous | 3 | 0.28 bits | 0.14 bits |
| contiguous | 4 | 0.50 bits | 0.17 bits |
| contiguous | 8 | 1.49 bits | 0.21 bits |
| contiguous | 16 | 2.16 bits | 0.14 bits |
| full group | 32 | -2.26 bits | -0.07 bits |

Evenly spaced patterns improve the iid bound.  At occupations 8, 16, and 64,
the tested dense patterns also improve it.

These values are floating-point diagnostics.  The scripts compute each stated
profile exactly within the floating-point transfer model, but they do not use
outward-rounded arithmetic.  The selected dense patterns are not exhaustive.

## Available margin

The iid calculation has large pointwise margins beyond the first few active
blocks:

| Occupation (a) | Iid pointwise margin |
|---:|---:|
| 16 | 981.53 bits |
| 32 | 2034.17 bits |
| 48 | 3163.07 bits |
| 64 | 4294.29 bits |

The observed scalar surcharges are negligible relative to these margins.
They provide positive evidence for the construction.

## Obstruction

The calculation does not prove a composable surcharge.  A scalar comparison
after 256 regions does not permit insertion of another defect group.  The
entrywise comparison composes, but its charge is much larger.  Proving an
arbitrary-background insertion inequality remains the missing step.

Consequently, the current cyclic-shift construction remains unproved at 9%.
The evidence suggests that the construction has sufficient distance, but the
present proof route lacks a rigorous compression across an arbitrary number
of defect groups.

The receipts are `receipts/one_pair_collision_surcharge_small_delta09.json`,
`receipts/one_pair_collision_surcharge_large_delta09.json`,
`receipts/two_pair_collision_surcharge_delta09.json`, and
`receipts/collision_surcharge_selected_profiles_delta09.json`.

