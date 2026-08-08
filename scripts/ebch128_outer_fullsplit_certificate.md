# Exact EBCH-128 block-outer certificate

Updated: 2026-08-07

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
| `1 <= h <= 500` | `-19.775758215947` | `19/2^24` |
| `501 <= h <= N` | `-50.102265480007` | `2^-50` |
| complete | `-19.775758214875` | `1275068417/2^50 <= 2^-19.75` |

Thus the same first-moment argument proves that some realization is a binary
`[2^21,2^20,d_min >= 188744]` code when this outer is used.  The smaller
block size loses about 17.5 bits of certificate margin relative to the frozen
RM-512 row, but it remains far below one.

The initially positive bounds for the critical and post-prefix windows were
pole-choice slack.  After rational pole retuning, the limiting post-prefix
component is the exact `h=501..2000` cap at `-53.424193575894`; the limiting
overall component is the exact/outward `h<=500` family.

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
