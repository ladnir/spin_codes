# Fixed affine placement: forward cost checkpoint

Status: parked by the user's 2026-10-06
[blockwise restart](HD_BLOCKWISE_RESTART.md). The
[preservation index](ACCUMULATOR_DETOUR_CHECKPOINT.md) records the verified
source archive and final proof status. The measurements below are unchanged.

The four-row placement control takes **8.327--8.360 ms** at K=2^20.
It is approximately twice as fast as the retained 16.622 ms forward RS
baseline, but it has no complete distance certificate. The sub-5-ms,
40%-distance, 40-bit goal remains open. Production and earlier certified
winners are unchanged.

## Construction being measured

Each 128-bit application record contains eight GF65536 lanes. The same
scalar code acts independently on those lanes. There are R=K/32 source
rows, each encoded by GRS[64,32]. All 2K output records are retained.

The forward pipeline is:

1. Encode the outer rows and store their coordinates in 64 planes.
   For coordinate c, affinely permute groups of four consecutive rows,
   preserving order inside each group. All offsets are zero.
2. Apply one input-weighted prefix. Every position has an independent
   nonzero GF65536 weight. The carry continues between planes.
3. Independently permute each plane's R records, then apply the two-state
   inner. Both states continue between planes.

The post-inner at position i is

    y = s1 + f_i*s2 + c_i*x;
    s1' = s1 + d1_i*x;
    s2' = s2 + d2_i*x.

The four coefficients are nonzero, with the joint distribution conditioned
on c_i+d1_i+f_i*d2_i being nonzero. Setup prepares
cp_i=c_i+d1_i+f_i*d2_i, so encoding can use the updated states. This
changes the schedule, not the map.

There is one sampled regional routing stage. Fixed affine placement is
fused into outer stores; it adds no separate full-record permutation pass.
Changing from individual-row placement to four-row groups nevertheless
changes the construction. The individual-row spectral guarantee does not
transfer. [The cut audit](HD_AFFINE_CUTS.md) records bounded tests for both,
not a certificate for either.

## Measurements

All entries are full forward calls at K20 on one Peach CPU2 core. Setup,
allocation, input generation, and checksums are outside timing. Every
timed call receives freshly regenerated input. There are three warmups.
Except for the selected checkpoint, controls use 31 timed calls.

| Placement / schedule | Median ms | Interpretation |
|---|---:|---|
| Individual rows, initial schedule | 13.083 | First full cost control |
| Individual rows, 16-record scan batches | 12.160--12.172 | Same map, improved scheduling |
| Individual rows, separate route pass | 17.463 | Reject extra staging |
| Four-row groups, ordinary stores | 9.273 | Changed placement, full cache-line stores |
| Four-row groups, streaming stores, tile64 | 8.469 | Same grouped map |
| Four-row groups, streaming stores, tile256 | 8.665 | Reject larger tile |
| **Four-row groups, streaming stores, tile16** | **8.327--8.360** | Selected cost checkpoint |

The selected schedule was rerun in two process launches, with 101 calls
each. Their p10--p90 intervals were 8.338--8.409 ms and 8.300--8.368 ms.
A separate 101-call instrumented launch measured:

| Phase | Median ms |
|---|---:|
| Outer encoding and fixed placement | 2.980 |
| Weighted prefix | 1.711 |
| Regional reads and two-state post-inner | 3.658 |
| Full instrumented call | 8.347 |

Phase medians need not sum to the median full time. Costs are distributed
across all three phases. The map still performs three state scans and
five full-field coefficient products per retained record.

GCC 15.2, Release, AVX512/GFNI, and a 4.5-GHz userspace frequency request
were used. Boost was disabled. Native builds and timings were serial under
shared benchmark locks. Large buffers requested transparent huge pages;
physical huge-page backing was not verified. The RS comparison uses the
retained forward baseline, not a new matched RS rerun.

## Validation and reproduction

The implementation is `spin/experiments/char2_field_spin/hd_affine_two_state.cpp`.
The default control still uses individual-row placement. The selected
research variant adds these compiler definitions:

    AF_GROUP=4 AF_NT=1 AF_TILE=16 AF_SPLIT_ROUTE=0

The streaming-store path dispatches once on alignment and fences before
prefix reads. Unaligned calls use ordinary stores. Sixteen-record scan
batches expose independent products before the prefix chains. The timed
encoder performs no allocation.

Release and AddressSanitizer/UndefinedBehaviorSanitizer checks pass at
K=128,256,1024,4096. The independent scalar reference directly evaluates
the additive polynomial basis and uses the original old-state recurrence.
Tests cover affine-versus-contiguous outer equality, exact input/output
aliasing, four record offsets, dirty destinations, and buffer guards.
These are implementation checks, not distance evidence.

The full K20 checksum is `18284413787472232435` for every grouped schedule
above. The individual-row checksum is `1848557650131273058`.

Checkpoint source SHA256:

    bc71ce8fe05b6948a69456c6f279f58e58ae1d713a6fe68ddeb05c01ae3338e4

Selected Release binary SHA256:

    9ad57a9a166c4fd4f224ac3c8a8216c987ac837db537699608957f88b8eb2bdc

```sh
cmake -S spin/experiments/char2_field_spin -B /tmp/spin-affine \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=g++-15 \
  '-DCMAKE_CXX_FLAGS=-DAF_GROUP=4 -DAF_NT=1 -DAF_TILE=16'
cmake --build /tmp/spin-affine --target spin_char2_hd_affine_two_state
/tmp/spin-affine/spin_char2_hd_affine_two_state --check
taskset -c 2 /tmp/spin-affine/spin_char2_hd_affine_two_state --bench 101 full
```

Use the shared benchmark lock on a shared machine. Never run full and
phase measurements concurrently. This target requires AVX512/GFNI;
it is not a portable production API.

## Proof status and next step

The [source law](HD_AFFINE_PROOF.md) and component contraction are valid,
but subsequent work proves that the [uncapped first moment cannot close](HD_COMPONENT_DENSE.md)
for this ensemble. Rare routes create large actual low-weight subspaces.
The [routing gate](HD_ROUTING_CONDITION.md) controls one family of such
events, not every low-weight vector. No full K20 certificate was obtained.
The parked proof checkpoint passed 158 high-distance helper tests.

Do not resume accumulator micro-optimization as the default next step.
The active search is the blockwise SPIN construction linked at the top.
