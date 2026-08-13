# First `g=8` witness-reuse reconnaissance

The first `g=8` profile probes found large margins when each tested profile
received local tuning. The next question was whether a small frozen witness
bank could cover distant profiles without repeated optimization.

These experiments are binary64 diagnostics. They do not certify a profile
region or contribute to an end-to-end union bound.

## Low-support representatives

The low-support survey selects one deterministic representative from each
feasible support stratum of dimension zero, one, or two. This gives 128
representatives. The baseline uses three inner seeds and five conditioned-row
outer seeds. They form 15 frozen affine witnesses.

```powershell
python scripts\probe_packet_group_g8_low_support_witness_bank.py `
  --seed-bank baseline --outer-bank-mode frozen-regression `
  --screen-iterations 1 --final-iterations 4 `
  --output out\g8_low_support_recon_smoke.json
```

Only 6 of 128 representatives meet the uniform per-profile target. The worst
representative is the pure class-7 profile. Its diagnostic deficit is about
1.64 million bits. The failures concentrate on supports containing classes 7
and 8, which the small bank does not seed directly.

Artifact SHA-256:

```text
4af47c5680497672e884543689e96bacb5c58b302488c5d0d1365ee6836ce1ba
```

Two enrichment experiments separate local strength from witness transfer.
Adding six high-class inner seeds, without matching outer seeds, still covers
only 6 representatives. It reduces the worst deficit by about 541,000 bits.

The paired-local experiment tunes both components at nine selected supports.
Every selected representative then passes by at least 289,000 bits. However,
the 81 cross-combinations cover only 12 of 128 representatives. The other
support strata remain far from the selected anchors.

```powershell
python scripts\probe_packet_group_g8_low_support_witness_bank.py `
  --seed-bank enriched-high-class --outer-bank-mode paired-local `
  --screen-iterations 1 --final-iterations 3 `
  --outer-max-iterations 60 --outer-max-evaluations 800 `
  --output out\g8_low_support_recon_paired_local_smoke.json
```

The paired-local artifact has SHA-256
`6fefd7d715959b27abb4ca635cf3361af0c11a29e446443bff5c34ffbb86efc1`.
Its large off-anchor values are transfer diagnostics. They are not local
profile bounds.

## Full-support profiles

The full-support survey optimizes five density-family seeds. It evaluates the
frozen affine witnesses on 240 deterministic full-support profiles.

```powershell
python scripts\survey_packet_group_g8_full_support_witness_reuse.py `
  --output out\g8_full_support_witness_reuse_smoke.json `
  --seed-limit 5 --random-per-concentration 12 `
  --ridge-ratios 10,50,90 --inner-screen-iterations 1 `
  --inner-final-iterations 4 --outer-max-iterations 80 `
  --outer-max-evaluations 1200
```

The bank covers 62 of 240 profiles. It covers 46 of 63 profiles on the
binomial curve but only 6 of 108 near-boundary ridge profiles. None of the
sampled near-vertex or seeded-skew profiles pass. All five witnesses lead at
least one profile; the largest leader share is 36.7 percent. Thus witness
reuse is real but local.

The worst sampled profile lies near the equal `0/8` ridge. Its diagnostic
deficit is about 1.78 million bits.

Artifact SHA-256:

```text
03bb0e4952dcf23010bac4237868b7752ed835c963955d19f39a88b12c54ad72
```

An enriched bank adds seven ridge, near-vertex, and skew seeds. It evaluates
224 profiles and covers 83, compared with 66 for the density subset on the
same catalogue. The added witnesses newly cover 17 profiles. They improve 54
profiles and reduce one exponent by about 937,435 bits. The worst remaining
ridge still misses the uniform target by about 615,609 bits.

The enriched artifact has SHA-256
`18c26b526ec62efefe13ab53b32198bd3254f2130ad3b4f3211d6fdec8ba6ccc`.

## Consequence for the proof search

The data supports two conclusions. First, the current inequality has large
local margin at every tuned sparse anchor tested so far. Second, a small
global witness bank does not span the profile domain.

The next stage therefore tuned one representative in every exact support
stratum. These 510 witnesses form a support-local atlas, not a certificate.
The proof must now test affine reuse with an adaptive exact slab or BSP
decomposition. Refinement should target cells whose cell-local union
contribution remains largest.

## Complete support-local atlas

The support-local atlas stage is complete. It tunes paired inner and outer
witnesses at one deterministic representative of each feasible exact support.
The 16-worker Peach run completed in 202.47 seconds.

All 510 representatives pass the uniform per-profile target. The worst local
margin is 58,470.196 bits. The worst anchor has full support. This establishes
large pointwise diagnostic slack at every support anchor; it does not cover
the remaining profiles in any support.

```powershell
python scripts\build_packet_group_g8_support_seed_atlas.py `
  --workers 16 --checkpoint-dir out\g8_support_seed_checkpoints `
  --output out\g8_support_seed_atlas.json `
  --inner-screen-iterations 1 --inner-final-iterations 3 `
  --outer-max-iterations 60 --outer-max-evaluations 800
```

The atlas SHA-256 is
`c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e`.
The frozen manifest SHA-256 is
`cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616`.
