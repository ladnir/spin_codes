# IMT inner

**Independent-Map Transvection (IMT)** names the inner only. Describe the
complete construction as **SPIN with a BCH outer, randomized bit-transpose
permutation, and IMT inner**. Specify the BCH instance, message size, and inner
parameters separately; do not introduce a combined construction name.

IMT retains the recursive structure and RM-derived expansion of the earlier
design. It separates expansion from feedback and replaces field multiplication
with a sampled transvection. RM2Sub remains the name of the earlier inner,
not an alternative name for IMT.

## Interface

Over F_2, fix an expansion map A from s state bits to t output bits and a
feedback map B from t input bits to s state bits. Here, independent maps means
separately specified maps: it does not mean freshly sampled A and B.
The certified instances use fixed maps and do not require B=A^T or BA=0.

For input block X_i and state q_i, the inner computes

```text
Y_i     = X_i + A q_i
q_(i+1) = M_i q_i + B X_i,    q_0 = 0.
```

The shared setup independently samples each epoch's mixer M_i=I+u_i v_i^T.
It samples u_i uniformly from nonzero s-bit vectors, then v_i uniformly from
the orthogonal complement of u_i, including zero. Thus the sampler permits
the identity update. State persists across regions; the final state is discarded.

For transpose input U_i and reverse state r_(i+1), the adjoint computes

```text
V_i = U_i + B^T r_(i+1)
r_i = A^T U_i + M_i^T r_(i+1).
```

The terminal reverse state is zero. A^T consumes raw U_i, not emitted V_i.
The [transfer argument](asymmetric/TRANSFER_ARGUMENT.md) derives the bounds
for these separate maps.

## Certified instances

All instances below have t=128 and s=19. The rate belongs to the
complete code; IMT itself preserves length. The margins bound the probability
of a bad shared setup, not the failure probability of an individual message.

| BCH outer | Rate | log2(K) | Fixed IMT feedback | Relative distance threshold | Certified margin (bits, approximate) |
|---|---|---:|---|---|---:|
| [128,32,32] | 1/4 | 20 | Weight three, `greedy3_2` | 16.5% | 41.0481676058 |
| [128,32,32] | 1/4 | 20 | Same map | 19% | 30.0334910634 |
| [256,128,d>=38] | 1/2 | 16 | Weight five, `weight5_seed0` | 10% | 41.8183267960 |
| [256,128,d>=38] | 1/2 | 18 | Weight five, `weight5_seed0` | 10% | 50.1890763816 |
| [256,128,d>=38] | 1/2 | 20 | Same map | 10% | 50.0620882264 |
| [256,128,d>=38] | 1/2 | 22 | Same map | 10% | 48.3904916829 |
| [256,128,d>=38] | 1/2 | 24 | Same map | 10% | 46.4572549297 |

See the [quarter-rate certificate](asymmetric/CERTIFICATE_RESULT.md),
[half-rate certificate](asymmetric/bch256/weight5/README.md),
and the [finite migration records](finite_migration/README.md) for exact scope
and replay commands. The supported default remains RM2Sub while migration
is incomplete. All five selected half-rate lengths now have complete IMT
certificates. The K16 dense cover passed 512-bit replay, and its exact union
is bound to the retained 0.523367 ms transpose measurement.
The asymptotic IMT argument has passed an
in-session analytic review and is now in the manuscript at 11%, with passing
numerical evidence. Existing finite RM2Sub claims do not transfer by renaming.

## Code and record names

Existing source paths, symbols, target names, and receipt identifiers retain
their historical spelling. Verification binds their exact bytes, so this
terminology change does not rewrite them or relabel old benchmark records.

| Existing identifier or path under this directory | Meaning in current prose |
|---|---|
| `asymmetric/AsymmetricInner.h`, `asymmetricReverse` | IMT transpose kernel used by the quarter-rate implementation |
| `generated/asymmetric_greedy3_2_sparse/` | SPIN implementation with the selected weight-three IMT feedback map |
| `asymmetric/bch256/weight5/implementation/Weight5Inner.h` | IMT transpose kernel with weight-five feedback and compile-time emission policies |
| `asymmetric/bch256/weight5/implementation/Weight5Spin.cpp` | Complete half-rate SPIN implementation using that IMT kernel |
| `sparse_pages` in the half-rate implementation | Selected implementation policy, not a new inner or construction name |
| `ASYMMETRIC_*` receipts | Historical identifiers for IMT evidence; each receipt retains its original scope |

Use IMT in new prose and new public inner identifiers. Keep outer choice,
permutation choice, feedback-map identity, and implementation policies separate.
Historical symmetric-map experiments are not renamed wholesale.

Next: integrate the two certified instances as explicit opt-ins, retaining
RM2Sub as the baseline. The [migration plan](asymmetric/MIGRATION_PLAN.md)
tracks the remaining proof and implementation work.

The current [asymptotic investigation](imt_asymptotic/README.md) targets an
IMT-only paper. All occupancy regimes now have bounds at 11% distance,
with the same 39/4 outer-growth constant. The reviewed proof now replaces the
paper's asymptotic instantiation. Remaining finite cells still gate the global
default switch; the review is not independent external verification.
