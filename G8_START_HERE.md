# Riffle `g=8` starting point

The current `g=4` certificate uses the globally sampled puncture lane and is
frozen by `G4_GLOBAL_LANE_CERTIFICATE_MANIFEST.json`. Start `g=8` from that
savepoint. Do not modify the `g=4` artifacts to make a `g=8` diagnostic pass.

## Baseline to preserve

The local code remains the binary systematic extended BCH code
`[128,64,22]` in both layers. Only the packet width changes from four to eight
bits. The `g=8` proof must retain the declared lane bijections, packet
permutation, and recursive state permutations unless the construction itself
is revised explicitly.

For the global-lane construction, the `g=4` theorem proves

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
- The one-conditioned-row outer bound, including the global-lane graph-hole treatment.
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
`cc446d0c6a0c986a8712ef78c4aa7c9d6073687665b7e0f315164339a5b81616`.
The manifest includes the graph spectrum, EBCH weight distribution, systematic
split slices, the pure-Python point-cap policy, and the frozen Git tree for
transitive proof code. These are discovery inputs. They do not establish
coverage within a support.

## First adaptive coverage results

Direct root witnesses or fixed minimax mixtures close 17 of the 128 supports
of dimension at most two. Exact adaptive splits are therefore necessary even
in low dimension.

For full support, coordinate BSPs with 128 and 1,024 terminal cells close no
cell. The best bound remains about 1.18 million bits above zero. Thirty-two
witnesses tuned at the largest residual vertices all pass locally. They reduce
the best terminal bound to about 183,000 bits, but they still close no existing
coordinate cell.

Dominance-informed affine cuts are sound discovery tools, but the first exact
polyhedral implementation does not scale. A seven-node tree finishes quickly;
a 15-node tree takes about 17 seconds; a 63-node tree does not finish within
the bounded run. The number of exact vertices grows after each arbitrary cut,
and no tested dominance leaf closes.

The next geometry should retain explicit low-complexity vertices. Use laminar
or ordered chambers and let witness dominance choose their structured
boundaries. Do not continue the current arbitrary-halfspace BSP by increasing
its node budget.

## Structured-geometry results

Three bounded experiments tested refinements with exact ownership and small
vertex sets.

An ordered decomposition by one band mass produced seven cells with at most
40 vertices. It closed no cell. Even the best band leaves incompatible pure
corners together, so adding more thresholds along the same mass cannot close
the root.

A balanced level-two cumulative grid produced 165 cells, at most 81 vertices
per cell, and 6,435 total vertex incidences. Rational minimax mixtures improved
every singleton bound. However, the grid closed no cell. Its worst diagnostic
contribution was about `+1,454,125` bits.

Reflection blocks `{0,8}`, `{1,6}`, `{2,7}`, `{3,5}`, and `{4}` jointly control
their pure corners. Majority-mass caps do not preserve that property. The caps
introduce cross-block boundary vertices, and their ten leaves also closed no
cell. Retain only the exact maximum-mass faces as possible boundary lemmas.

The first adaptive cumulative wave refined the eight largest level-two
contributions. Each parent received one exact prefix split at a discovered
witness transition. All eight splits improved their parent by between
22,533 and 66,713 bits. The global worst contribution fell by about 38,959
bits, to `+1,415,166`. The 16 children used 540 exact vertex incidences and
the run took about 22 seconds.

This wave passes the continuation gate, but it is not close to a certificate.
A linear extrapolation suggests roughly 31--63 comparable waves. Continue
with a resumable cumulative-prefix engine. Stop if improvement plateaus or
the exact vertex sets grow beyond the declared cap. The diagnostic selectors
also use 32 supplementary witnesses that are not yet in the frozen manifest;
regenerate the manifest before outward certification.

The corrected engine ran two waves and stopped at its finite gate. Both waves
accepted eight splits, but the lower-quartile rate projected 53 more levels
and 4,269 leaves. Geometry-only refinement is therefore a no-go.

The inner and outer witness components can be recombined independently. A
factorized minimax replay improved every active leaf and reduced the frontier
by about 257,500 bits, to `+1,134,921`. One further adaptive wave under this
stronger family reduced it by about 49,550 bits, to `+1,085,371`. The proof
selector should store separate rational inner and outer marginals; explicitly
materializing their Cartesian product would be unnecessarily quadratic.

This factorization needs one theorem statement made explicit. The inner bound
must hold uniformly after conditioning on the complete outer setup used by
the conditioned-row lemma. The construction retains independent inner lane,
packet, and state randomization, but the final proof and manifest must bind
that conditional statement.

Sixteen additional witnesses tuned at active worst vertices did not change
the factorized worst aggregate. Pointwise tuning is therefore exhausted for
the controlling cells. The next witness jobs should use the minimax LP dual
barycenter of a failed cell and tune inner and outer components separately.
Accept a component only when it has a verified negative reduced cost against
the current component hull.

The first cell-aware round tuned eight inner and eight outer components at
exact dual barycenters. The fixed 189-leaf replay reduced the frontier by
about 400,367 bits, from `+1,085,371` to `+685,003`. Every new component
passed its reduced-cost test, and all 189 leaves improved.

An audited second round tuned only eight outer components. It reduced the
frontier by another 51,671 bits, to `+633,333`. Every component transferred
outside its source cell, and 159 non-source leaves improved. The slowdown is
material: the lower-quartile gain projects about six comparable rounds, which
exceeds the two-round continuation gate. Stop pure pricing and alternate with
one exact adaptive geometry wave under the enlarged factorized bank.

Independent checks found no affine-sign or normalization error. They
reconstruct both component identities, subtract the orbit normalization once,
and retain the old-selector fallback. The result remains diagnostic until the
conditional inner lemma and outward factorized replay are complete.

The factorized proof state is now independently outward-replayed. The v2
manifest binds 92 inner components, 34 outer components, and exact rational
selectors for 197 leaves. Pure-Python interval replay checks 8,687 vertex
inequalities. Its expanded interval is approximately

```text
[580922.1202845245, 580922.1203803732].
```

Thus the ledger is rigorous but misses the `-40` target by about 580,962 bits.
The interval width is below `0.0001` bit, so outward arithmetic is not the
source of the gap.

Two complementary outer families were then added with the audited
128-replacement graph correction. Eight three-band BL2 rows plus eight exact
total-spectrum rows reduce the diagnostic frontier to about `+514,758`. A
fresh round at the shifted worst leaves reduces it to about `+486,266`, but
improves only 28,492 bits and fails the continuation gate. Stop outer-column
generation at this point.

The remaining high-leverage omission is combinatorial. The current outer
columns discard the certified three-band nonclumping properties: pair capacity
one, absence of Pasch/intercalate patterns, exact band balance, and dense
four-block pattern counts. The next proof family should retain enough pattern
state to use these constraints. Changing `N`, adding more mixtures, or running
more prefix waves has much lower expected leverage.

The nonclumping refinement must act on an augmented outer fiber. The nine
packet-weight counts do not determine the three collision partitions. A sound
first lemma conditions on the active outer blocks and retains the 195 band
occupancy counts. Finner's inequality then yields an exact coefficient bound
that uses band balance and pair capacity. Pasch exclusion requires a separate
motif-sensitive majorant; inserting unconditional motif probabilities into a
packet-profile fiber would be unsound.

The full active-set scan has now evaluated the pair-capacity coefficient bound
at `h2:073`. It sums all 16,385 possible support sizes and uses eight fixed
collision tilts. The optimized sum equals the zero-tilt sum:

```text
pair-capacity saving = 0 bits.
```

Only the isolated support-size-one term improves, by about `0.0022` bits.
Dense support terms dominate the outer moment and select zero collision tilt.
Do not extend the pair-capacity histogram scan to the other seven leaves.

The motif-sensitive majorant also fails its theorem gate. The fitted
interactions are near 0.0182 and 0.0186. The uncentered polymer criterion is
about 9.6, where convergence requires at most one. A conservative
centered-Ising/Kotecky--Preiss criterion is about 19.3. Dobrushin contraction
holds but does not certify the required truncated log-partition remainder.
Stop the current nonclumping family. A continuation would need a new
block-resummed loop theorem, with little evidence that its gain could approach
the remaining hundreds of thousands of bits.

## Conditioned-row layout audit

The exact-graph conditioned-row branch has a construction-binding defect. Its
63-free-row decomposition requires one perfect matching of retained data
blocks across all three band-tile classes. The frozen sampler chooses an
independent puncture lane in each selected band-zero tile. An allowed pair of
punctures can collide in band one, so no matching can contain every punctured
block. Exact spectra and outward arithmetic do not repair this global defect.

The proof-only repair conditions a fixed diagonal matching on the unpunctured
data word and pays a pointwise 128-replacement packet-fugacity tax. The first
diagnostic costs about 558 bits at the canonical `g=4` leader and 515 bits at
`g=8` leaf `h2:073` before the trivial outer cap.

The zero-hot-path-cost construction repair is to sample one global puncture
lane and use it in all 128 selected tiles. A fixed lane class is a perfect
matching in every band. This rule preserves every block's puncture marginal
and makes the existing exact-graph one-row formula applicable. It also permits
a separated rank-two probe with a second lane class. The construction must
choose this rule explicitly before the old `g=4` or current `g=8` receipts are
called end to end.

The full audit is in
`explorations/conditioned_row_sloped_layout_audit.md`.

The structured-geometry proofs and commands are in:

- `explorations/g8_full_support_structured_geometry.md`;
- `explorations/g8_structured_geometry_lane_b.md`;
- `explorations/g8_full_support_ordered_band_geometry.md`;
- `explorations/g8_adaptive_cumulative_prefix_refinement.md`;
- `explorations/g8_independent_inner_outer_recombination_audit.md`;
- `explorations/g8_cell_aware_factorized_column_generation.md`;
- `explorations/g8_cell_aware_second_wave_gate.md`;
- `explorations/g8_three_band_nonclumping_outer_lemma.md`;
- `explorations/three_band_nonclumping_profile_fiber_audit.md`.

## Proof dependencies

The EBCH and graph spectrum tables remain authenticated mathematical inputs.
The graph-hole analysis assumes the construction's graph-syndrome uniformity.
These assumptions must remain visible in any `g=8` theorem statement.

The complete `g=4` proof and reproduction command are in `PROOF_STATUS.md`.
The mathematical development is in
`explorations/riffle_group_chain_proof.md`. The machine-readable checkpoint is
`G4_GLOBAL_LANE_CERTIFICATE_MANIFEST.json`.
