# Exact EBCH-128 block-outer certificate

Updated: 2026-08-08

This is a separate exact certificate lane for the finite full-split
construction.  It replaces the direct sum of 4096 `RM(4,9) [512,256,32]`
outer blocks by 16384 copies of the committed extended BCH `[128,64,22]`
weight enumerator.  The length, dimension, inner construction, and target
distance are unchanged:

- `N=2^21`, `K=2^20`;
- outer: 16384 extended BCH `[128,64,22]` blocks;
- inner: full-split extended BCH `[128,64,22]`, with `b=64`;
- target `d=floor(.09*N)=188743`.

## Result

The regenerated rational/outward artifact reports:

| Family | Diagnostic log2 upper bound | Theorem-safe gate |
|---|---:|---:|
| `1 <= h <= 500` | `-20.805420115939` | `588/2^30` |
| `501 <= h <= N` | `-50.102265480007` | `2^-50` |
| complete | `-20.805420113752` | `616562689/2^50 <= 2^-20.80` |
| random codimension-20 subcode | `-40.805420113752` | `616562689/2^70 <= 2^-40.80` |

Thus the same first-moment argument proves that some realization is a binary
`[2^21,2^20,d_min >= 188744]` code when this outer is used.  The smaller
block size has 20.80 theorem-safe bits of full-dimension margin.  Independently
choosing a uniformly random codimension-20 outer subcode multiplies every
nonzero-word first-moment row by less than `2^-20`.  It therefore proves the
stronger-margin existence statement
`[2^21,1048556,d_min >= 188744]` with 40.80 theorem-safe bits.

The initially positive bounds for the critical and post-prefix windows were
pole-choice slack.  After rational pole retuning, the limiting post-prefix
component is the exact `h=501..2000` cap at `-53.424193575894`; the limiting
overall component is the exact/outward `h<=500` family.

The low-weight improvement subdivides the first 4000-block gap bucket into
250-block rational subbuckets for weights through 150.  The remaining weights
retain the previous coarser upper bound.  This changes the first bucket from
`-20.7464` to `-33.1587`; the exact ultra-late sum at `-20.8057` then becomes
the full-dimension bottleneck.

## Why full dimension cannot reach 40 bits in this ensemble

The artifact also contains a lower-bound audit.  The weight-22 outer words
whose permuted support lies in the final 4000 inner blocks have placement
expectation `2^-34.858313...`.  Conditional on this event, the expected output
weight is at most 4000 times the exact per-live-block mean.  Markov's
inequality gives conditional probability at least `0.321832747...` that the
output weight is at most 188743.  Consequently this family alone contributes
more than `2^-36.493931` to the actual first moment.

Therefore no upper-bound cleanup can give the unchanged full-dimension,
uniform-interleaver ensemble 40 bits.  The codimension-20 result is a genuine
expurgation step, not merely more favorable reporting.  Its rate loss is only
20 dimensions out of `2^20`.

## Reproduction

Fast artifact and input-fingerprint check:

```powershell
python scripts\certify_ebch128_outer_fullsplit.py
```

Full regeneration and byte-for-byte artifact comparison:

```powershell
python scripts\certify_ebch128_outer_fullsplit.py --recompute-artifact
```

Canonical files:

- `scripts/certify_ebch128_outer_fullsplit.py`
- `scripts/ebch128_outer_fullsplit_rational.json`
- `scripts/ebch128_64_spectrum.csv`
- `scripts/EBCH128_64.wd`
- `scripts/fullsplit_h500_gap_sums_exact.json`
- `scripts/fullsplit_h500_inner_bounds_dyadic.json`

The verifier structurally checks the spectrum length, total multiplicity,
minimum positive weight, and complement symmetry, and authenticates all four
input artifacts by SHA-256.  As in the existing RM/EBCH theorem lane, the
mathematical provenance of the committed EBCH weight enumerator is a declared
input: this certificate validates and consumes that table but does not derive
the enumerator from a generator polynomial.

The manuscript theorem remains on the frozen RM-512 lane until this alternate
certificate receives a proof review and is deliberately promoted.
