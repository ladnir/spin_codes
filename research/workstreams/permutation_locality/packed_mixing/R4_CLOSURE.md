# Packed GL32 with Four Inner Updates

2026-09-30. **Closed: greater than 10% relative distance with more than
52.05 bits of whole-code setup-failure margin**, at K=2^20 and rate
one-half. Fresh 384-bit replay covers every occupancy and passes the exact
aggregate and unchanged-source checks. This is the separate R4 operating
point, not a replacement for the [retained R2 certificate](FIRST_CLOSURE.md).

The latest precomputed transposed implementation takes **5.334629 ms** for
128-bit XOR elements; the matched confirmation spans **5.32--5.35 ms**.
[Persistent packed GFNI state](../gfni_inner_REPORT.md) reduces time by
9.41% against a matched 5.889024 ms scalar-fused control, preserving the
same 19-state map and its distribution. Setup and allocation are excluded.

The preceding [wide BCH stores and compact GL32 coefficients](../packed_bch_tune_REPORT.md)
measured 5.854159 ms, reduced time by 4.91% against that campaign's
matched control, and halved coefficient storage to 8 MiB.

The earlier [exact-map inner fusion](../packed_fused_r4_REPORT.md) measured
6.212096 ms, versus a matched sequential R4 control of 6.753111 ms. Both
optimization stages, and the subsequent GFNI implementation, preserve the
sampled encoder exactly. Each report
keeps its own controls, checks, and source hashes; the historical timings
are not mixed into the new matched percentage. The earlier
[update-count comparison](../packed_updates_REPORT.md) also remains available.

## Construction

The stable identifier is `canonical-gl32-width8-shared4-r4`. The forward
encoder maps K=2^20 input bits to N=2^21 output bits:

1. Encode 8,192 rows with the fixed BCH[256,128] constituent, then collect
   each four adjacent rows into one of 2,048 groups.
2. Divide each group into 32 canonical blocks of eight adjacent columns.
   Independently sample a uniform GL(32,2) map for each four-row/eight-column
   block and apply it to those 32 bits.
3. Independently shuffle the 256 columns in each group, sharing that
   permutation across its four rows. Each column is one four-bit packet.
4. Route column j to region j. Independently shuffle the 2,048 packet
   positions in each region, then concatenate the 256 regions.
5. Apply IMT with t=128, s=19, zero initial state, no terminal flush, and
   four independent sampled transvections per step.

All setup choices are independent and are fixed for every message after
setup. There is no additional GF16 packet randomizer. The expansion and
feedback maps A and C are the fixed `Map128S19` maps in
[the selected maps](../../../../spin/src/kernels/generated/SelectedMaps.h).
For input block X_i and entering state Q_i, the recurrence is

    Y_i = X_i + A Q_i,
    Q_(i+1) = T_(i,4) T_(i,3) T_(i,2) T_(i,1) Q_i + C X_i,
    Q_0 = 0.

Each transvection is I+u v^T, with u uniform nonzero and v uniform in its
orthogonal complement; v=0 gives an identity update. The raw-input feedback C X_i is
added after the four transvections. The output uses the entering state.
The measured transposed encoder applies the adjoint of this forward map.

Every setup realization has rate exactly one-half: the BCH outer is
injective, the local maps and routing are invertible, and the inner can
be inverted successively from its outputs and the known initial state.

## Claim and Coverage

The bad-weight cutoff is 209715 = floor(N/10). Let Z count nonzero messages
whose output has weight at most this cutoff. A whole-code first-moment
upper bound U gives Pr[d_min <= 209715] <= E[Z] <= U over the stated ideal
setup distribution. Outside that event, minimum weight is at least 209716,
so relative distance is strictly greater than 10%.

The sparse portion covers q=1..32 and the dense portion q=33..2048, where
q is the number of nonzero four-row groups. Sparse support partitions cover
their complete domains. The dense witnesses form a disjoint exhaustive
337-interval partition of the auxiliary comparison domain
[2079/7829504,1]. These are not merely selected comparison points.

Only fresh outward bounds may enter the final sum. The replay authenticates
the BCH premises, reconstructs the actual R4 state operators, evaluates
every saved witness at 384-bit precision, and sums the resulting dyadic
endpoints exactly before one final upward rounding. Search scores and saved
numerical endpoints are not proof inputs. The four dense workers each
construct their own ordinary model; no search cache or pickled model is
imported. The coordinator checks their exact scope, assignment, precision,
and witness hashes before aggregation.

Both fresh child replays and the whole-code coordinator have passed.
An independent audit of every one of their 32 sparse and 337 dense
endpoints reproduces the coordinator's exact sum:

| Covered messages | Margin from the outward upper bound (bits) |
|---|---:|
| Sparse, q=1..32 | 58.024133122014 |
| Dense, q=33..2048 | 52.076065452225 |
| Whole code | 52.052884355478 |

The whole-code upper bound is approximately 2.14052570818e-16. Its exact,
upward-rounded dyadic representation is

```text
U = 163139234043523590988548289792810425783243537599046889881897963564631206129321543877737667526548377970320243027327188234527379 * 2^-468.
```

Integer arithmetic verifies U^20 < 2^-1041, so the conservative margin
**exceeds 52.05 bits** without relying on a rounded logarithm. The requested
40-bit target is therefore exceeded by more than 12.05 bits. This bounds
setup-failure probability under the stated distribution, not the actual
distance of every sampled code or a cryptographic attack work factor.

## What Changed the Bound

The [hill-climb ledger](HILL_CLIMB.md) records the experiments in order.
Three changes matter:

- Exact production-coordinate enumeration excludes BCH words confined to
  five canonical octets. Incidence bounds and an intersecting-family
  refinement further tighten the number of four-row tuples supported on
  few octets. These are code-specific counting improvements, not an assumed
  BCH weight spectrum.
- Joint normalized regional-count bounds improve the R2 dense analysis,
  but the tested improvements are much smaller than its remaining gaps.
  They do not provide the present R4 closure route.
- Increasing the actual number of inner updates to four reduces those
  dense gaps. R3 improves them substantially but still fails selected
  points at 10%. Each update count is analyzed separately; more mixing is
  not assumed to improve every low-weight tail monotonically.

The earlier R2 result remains >9.5% distance with >36.68 bits of margin.
Its stronger 10% sparse result is not a whole-code certificate. Neither
that R2 claim nor its timing is transferred to R4.

## Retained Witnesses and Replay

Numerical receipts remain local under ignored `tmp/`; no raw data is added
to version control. The complete search inputs are:

| Input | SHA256 |
|---|---|
| `tmp/packed-hill/dense-r4-ekr-d10-complete-p256.json` | `a44142155d23553b589983cfb87d3ed3f8a0514217e2e46a0efb4477065bf01c` |
| `tmp/packed-hill/sparse-r4-ekr-d10-p256.json` | `1bfede612f7d5534848adfcd83b4998b62b474a24ac5274e5458bf72ab323331` |

The final receipt and its fresh child receipts are under
`tmp/packed-hill/r4-d10-whole-parallel-p384/`:

| Receipt | SHA256 |
|---|---|
| `sparse.json` | `92a51555bc1865709652050865c67c53f74a749e9b52735235284cd6a3317e77` |
| `dense.json` | `fd33f911e084ba7c868a87954c644a41ee3fd0a806320859ca297faeee7db9e5` |
| `whole.json` | `637a1ac049a7aa81441705e5a4fbeebf79e78fa80e63fe54f93cfe6dad05a8f0` |

Run from the repository root, choosing a new output directory:

```text
python -B research/workstreams/permutation_locality/packed_mixing/whole_hill_replay.py --dense tmp/packed-hill/dense-r4-ekr-d10-complete-p256.json --sparse tmp/packed-hill/sparse-r4-ekr-d10-p256.json --output-dir tmp/packed-hill/r4-d10-whole-rerun-p384 --precision 384 --target-bits 40 --dense-workers 4
```

The sparse and dense stages run sequentially. Only the independent dense
checks run in parallel. Source-byte checks reject changes during replay;
the final receipt records source hashes and both fresh child receipts.
The stopped serial directory `tmp/packed-hill/r4-d10-whole-p384/` is a
partial verification run, not a whole-code certificate.

All 260 packed-mixing tests pass, including the parallel replay checks.
Independent review also checked the complete search partition and a real
cold-spawn import/serialization path. The final read-only audit reproduces
all three aggregate endpoints, checks all 479 proof-source manifest files
against the current source set, and verifies 545 reference-hash links across
504 unique files. The old R2 whole-code and search receipt hashes remain
unchanged. Production defaults and paper claims remain unchanged.

## Exact-Map Optimization and Next Step

The optimized kernel fuses four transvections without changing their
product. For T_i=I+u_i v_i^T,

    T_4 T_3 T_2 T_1 = I + sum_(i=1)^4 a_i v_i^T,
    a_i = T_4 ... T_(i+1) u_i.

The suffix is the identity for i=4. Precomputing the a_i permits four dot
products against the original state followed by one combined update. The
transpose uses the adjoint rank-at-most-four correction. This preserves each
sampled map exactly; it is not a weaker setup distribution. The scalar
subset-table implementation saves 0.541016 ms (8.01%) in the matched run.
The AVX512 grouped-table candidate performs similarly but needs more setup
storage. Both pass every-epoch complete basis-map tests, inner adjoint
checks, and full-encoder equality against the original sequential masks
for two seeds at K=2^14 and K=2^20.

Preserve this operating point before changing the construction again.
The subsequent [packed-stage optimization](../packed_bch_tune_REPORT.md)
reduces the GL32/BCH phase from roughly 3.4--3.5 ms to roughly 3.2 ms.
Its next step is integration into the reusable precomputed path. The old
R2 and sequential/fused R4 implementations remain intact as controls;
no production default or paper claim has changed.
