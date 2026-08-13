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
shards. Profile-local optimization is too expensive to serve as the coverage
algorithm.

## Proof dependencies

The EBCH and graph spectrum tables remain authenticated mathematical inputs.
The graph-hole analysis assumes the construction's graph-syndrome uniformity.
These assumptions must remain visible in any `g=8` theorem statement.

The complete `g=4` proof and reproduction command are in `PROOF_STATUS.md`.
The mathematical development is in
`explorations/riffle_group_chain_proof.md`. The machine-readable checkpoint is
`G4_CERTIFICATE_MANIFEST.json`.
