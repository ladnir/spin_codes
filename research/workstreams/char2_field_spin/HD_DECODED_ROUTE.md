# Predecode the certified pair route

The 3.4-ms goal remains open. Predecoding the existing route improves the
certified full encoder from about 3.81 to 3.68 ms in a matched comparison.
It changes setup storage, not the sampled map or its distance proof.

## Exact representation change

The reference route stores each pair base in 15 bits and its independent
swap bit in bit 15. Every four-record load decodes two entries and constructs
the eight-bit SIMD shuffle mask during encoding.

The new setup representation stores 32 clean 16-bit bases followed by
16 shuffle-mask bytes for each 64-record inner packet. It occupies 80 bytes
instead of 64. At K20, hot route storage increases from 2 to 2.5 MiB.
The original sampled table remains available to the independent oracle.

For encoded pair entries a and b, setup computes bases a&32767 and b&32767.
Their mask is 15 times the first swap bit plus 240 times the second.
Encoding still loads the same two contiguous record pairs and applies the
same masked shuffle. The odd-region boundary pair still uses the copied
first record after the region. Setup makes no additional random draws.

`HDGROUP_DECODED_ROUTE=1` enables the isolated implementation in
`predecoded_route/`. The default remains zero. Production code is unchanged.

The outer, route distribution, persistent state recurrence, and refresh
keys are unchanged for every setup realization. Therefore the
[complete certificate](HD_SHIFTED_CLOSED_CERTIFICATE.md) still applies:
at K=2^20 and rate 1/2, the ideal sampled scalar GF256 code has distance
greater than 35%, except with probability at most 2^-55.092738508223.
Using its same scalar map on all sixteen byte components preserves this
record-Hamming distance. This does not certify the deterministic MT seed
used for repeatable performance measurements.

## Native evidence

Peach CPU2, GCC15, precomputed setup, 128-bit records, streaming output,
501 timed calls per run. Input generation, allocation, and setup are
excluded; the entire outer, routed inner, and output stores are included.
Each experiment holds both shared benchmark locks. All runs are serial.

| Implementation | Full run medians, ms | Median of run medians |
|---|---|---:|
| Matched reference | 3.783304, 3.837895, 3.812769 | 3.812769 |
| Predecoded masks | 3.682034, 3.755071, 3.635908 | 3.682034 |
| Direct-source first butterfly | 4.420804, 4.456831, 4.440610 | 4.440610 |
| Direct-source plus predecode | 4.357716, 4.350552, 4.353508 | 4.353508 |

Predecoding reduces matched median time by 3.43%. Separate inner-only
runs give 1.759545 ms for the reference and 1.601360 ms for predecoding.
Separate phase measurements are not an exact decomposition of full latency.
The older frozen 3.798160-ms checkpoint remains preserved.

All four builds pass the three-route 30,720-direction input/state basis
check, all pair-swap masks including wrap, small-size dense checks, and the
K20 dense oracle. Their certified-profile checksum is 7701295210980077712.
The predecoded variant also passes ASAN/UBSAN small-size and basis checks.
The independent source review found no sampling or layout change.

The certified cover was rechecked with its pinned manifest at 384-bit
aggregation precision. This authenticated the previously validated replay
cache and freshly checked complete coverage and directed aggregation;
it did not rerun every cached coefficient calculation.
All 11,465,830,401 feasible profiles are covered; the same 55.092738508223-bit
ceiling results.

## Rejected and follow-up implementation probes

The direct-source outer variant removes an intermediate scratch pass but
regresses materially. Assembly shows unchanged RS16 leaves and no new vector
spills. It changes sequential source staging into interleaved low/high loads.
The corrected traffic estimate is 56 MiB of logical L1 traffic removed at
K20, not 64 MiB. The dynamic CALL count does not decrease: the original
parent used a tail jump for its second leaf.

A second variant batches sixteen low-half loads before processing the high
half. Assembly confirms the intended batching without new spills. It
removes the large regression but does not produce a repeatable win:
3.717791, 3.733821, 3.786600 ms against predecoded-only controls of
3.763236, 3.649163, 3.708894 ms. Retain the original outer schedule.
The median of run medians is 3.733821 versus 3.708894 ms; the individual
rounds overlap and change order. Do not attribute the initial regression
to a uniquely identified hardware bottleneck from assembly alone.

Pre-scaling the clean route indices by eight permits byte addressing with
hardware scale two. It removes explicit address shifts without enlarging
the table, but limits this isolated 16-bit-index experiment to K<=2^20.
It does not improve full encoding: medians 3.730154, 3.779717, 3.696041 ms,
against 3.687695, 3.731527, 3.688105 ms for predecoding alone. Retain the
unscaled route; do not introduce the extra size restriction.

## Reproduction binding

The initial matched sources are frozen under `predecoded_route/initial/`.
`build.sh`, `bench.sh`, and saved logs describe the isolated variants.
The remote build directory is `/tmp/hd-exact-map-20261009-TTd1Wm`.
The initial binaries have these SHA256 values:

| Binary | SHA256 |
|---|---|
| Reference | `802015442c2cf8a3d664cdef601947f00ba9981632e8de3b08d6764d56983344` |
| Predecoded | `b1d27813dca7052929c3f118e25d42e2e3cdad83ced9f03d44ec1c0cd7a44db9` |
| Direct-source | `88ad3d0c2fc1ddaf727cb4451b83d5ae5179346a197a662957659889e62cbbe1` |
| Combined | `6cc82066c1fd97a85a1da44e5f1f44651084dcae1d3d192d94f38bbb7d771823` |

`hashes.log` pins the initial headers and driver. Later flags are recorded
separately; do not identify a variant only by its command-line arguments.
No library profile, paper, commit, or remote repository changed.
