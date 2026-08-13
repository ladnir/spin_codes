# Riffle `g=8` starting point

The `g=4` certificate is frozen at tag `riffle-g4-certified-v1`. Start `g=8`
work from that tag. Do not modify the `g=4` artifacts to make a `g=8`
diagnostic pass.

## Baseline to preserve

The local code remains the binary systematic extended BCH code
`[128,64,22]` in both layers. Only the packet width changes from four to eight
bits. The `g=8` proof must retain the declared lane bijections, packet
permutation, and recursive state permutations unless the construction itself
is revised explicitly.

The `g=4` theorem proves

```text
Pr[d_min <= 188743] <= 2^-61.78815155534398186100.
```

The proof has `21.7881515553` bits of margin beyond the 40-bit target. Run
`python scripts\verify_g4_savepoint.py` before using the checkpoint.

## Machinery that transfers

The following components are independent of the four-bit profile dimension.

- The packet-orbit composition and the fixed-witness convexity argument.
- Exact support stratification and support eligibility for zero fugacities.
- Exact rational mixtures with independently hardened vertex evaluations.
- Dual column generation over a finite witness atlas.
- Cell-local profile-count bounds and outward log-sum-exp aggregation.
- The one-conditioned-row outer bound, including the graph-hole treatment.
- Exact BSP ledgers and their independent outward verification.

The one-conditioned-row inequality is the decisive `g=4` improvement. Test it
before developing a new outer bound.

## What does not transfer directly

A global anchor triangulation does not scale from four to eight profile
dimensions. Its cell count can grow too quickly, and a floating-point mesh is
not a certificate. The `g=8` effort therefore needs an adaptive decomposition.

The initial proof search should proceed in this order.

1. Generalize the witness evaluator and conditioned-row outer to nine profile
   classes. Verify them on `g=4` regression cases.
2. Measure the bound on pure supports, low-dimensional support strata, and a
   representative sample of full-support profiles.
3. Identify whether the leading gap comes from the outer bound, the inner
   witness, or profile-domain aggregation.
4. Choose geometry only after locating the residual region. Prefer support
   strata, analytic slabs, and exact BSP cells over a global 8D mesh.
5. Freeze each successful discovery witness. Certify it with outward arithmetic
   before adding it to a complete ledger.

The first stopping point is a quantitative reconnaissance report. It should
state the worst profiles, their certified or diagnostic margins, the active
witness branches, and the predicted cost of completing the residual domain.

## First reconnaissance checkpoint

The packet-width interface now supports `g=8`. The regression gate reproduces
the decisive frozen `g=4` outer value bit for bit. It also checks the `g=8`
total-mass identity and five frozen `g=8` discovery vectors:

```powershell
python scripts\check_packet_group_conditioned_row.py
```

The exact domain census has 510 feasible supports and

```text
binom(262144+8,8)-2099
=553169839211945865258921061892182603726
```

feasible profiles. The census gate is:

```powershell
python scripts\check_packet_group_g8_census.py
```

A bounded five-profile run found large diagnostic margins. Its worst tested
profile was the `Binomial(8,1/2)` profile, with a combined exponent near
`-58640.97`. This result is evidence that the local inequality remains strong.
It says nothing about complete coverage of the profile domain.

The proposed decomposition is documented in
`explorations/g8_domain_decomposition.md`. It uses 510 support shards, exact
integer slabs, and local BSP trees. Only the full-support shard has dimension
eight.

The conditioned-row outward evaluator now supports `g=8`. Its regression gate
checks five frozen vectors, the exact mass identity, the canonical `g=4`
leader, and an independent 180-digit reconstruction:

```powershell
python scripts\check_packet_group_g8_conditioned_row_outward.py
```

The next proof-critical gap is reusable witness coverage across the support
shards. Profile-local optimization can seed an atlas, but it cannot serve as
the final coverage argument.

The first low-support reuse experiment used three inner seeds and five outer
seeds across one representative from each of the 128 support strata of
dimension at most two. Only six representatives passed the uniform target.
The worst failures concentrate on classes 7 and 8. This is a transfer failure
of the small witness bank, not evidence against profile-local bounds. The next
low-support bank should include support-local seeds for the high classes before
any geometric subdivision.

The enriched surveys confirm that local adaptation is effective. Nine paired
inner/outer sparse anchors all pass locally by at least 289,000 bits, but they
cover only 12 of 128 low-support representatives. On 224 full-support samples,
an enriched 14-witness bank covers 83 profiles. The density subset covers 66.
One targeted witness improves its profile by about 937,435 bits. These results
favor a 510-support local atlas followed by exact adaptive slab or BSP coverage.
They do not favor a larger fixed global bank.

The commands, artifact hashes, and detailed results are in
`explorations/g8_witness_reuse_recon.md`. The proposed shard artifact contract
is in `explorations/g8_support_shard_schema.md`.

The complete support-local reconnaissance tunes one paired inner/outer witness
at a deterministic representative of every feasible support. All 510 anchors
pass the uniform target. The worst local margin is about 58,470 bits and occurs
at the full-support anchor. Peach produced the atlas in about 202 seconds with
16 workers.

The frozen atlas is `out/g8_support_seed_atlas.json`, with SHA-256
`c48644fa62ad74619a4db30b07950b14aaf1419d040760f7d659f43288fca45e`.
`G8_SUPPORT_MANIFEST.json` binds the atlas rows, exact support counts, spectra,
and evaluator sources. The manifest's SHA-256 is
`dafff5c978d51515355960293832223a8a4736ca4b2c0405496c89d395ac3890`.
The manifest includes the graph spectrum, EBCH weight distribution, systematic
split slices, the pure-Python point-cap policy, and the frozen Git tree for
transitive proof code. These are discovery inputs. They do not establish
coverage within a support.

## Proof dependencies

The EBCH and graph spectrum tables remain authenticated mathematical inputs.
The graph-hole analysis assumes the construction's graph-syndrome uniformity.
These assumptions must remain visible in any `g=8` theorem statement.

The complete `g=4` proof and reproduction command are in `PROOF_STATUS.md`.
The mathematical development is in
`explorations/riffle_group_chain_proof.md`. The machine-readable checkpoint is
`G4_CERTIFICATE_MANIFEST.json`.
