# One-round mixer performance

The isolated one-round transvection encoder is about 4% faster at K=2^20.
It has **no full distance certificate**, so it does not replace the supported
encoder. These timings use the original selected A/B maps, not the later
no-constant map in `NO_CONSTANT_MAP.json`.

The final balanced-image map has since been certified and measured separately.
See [BALANCED_RESULT.md](BALANCED_RESULT.md) for that operating point and its
replay status. This document retains the original-map comparison as a control.

## Workload and measurements

Peach, Ryzen 9 7950X, CPU 15, GCC 15.2, release `-O3 -DNDEBUG -march=znver4`.
Quarter-rate BCH [128,32,32], t=128, s=19; complete 128-way bitsliced transposed
encoding in-place. Layout is packed24, tile size 4096. Setup is excluded.

| Implementation | Three 101-trial run medians (ms) | Median of medians |
|---|---|---:|
| Original certified mixer | 17.594, 17.707, 17.598 | 17.598 ms |
| New mixer, dense matrix rows | 17.643, 17.336, 17.724 | 17.643 ms |
| New mixer, sparse rank-one update | 17.069, 16.881, 16.802 | 16.881 ms |
| New mixer, fixed masked update | 16.838, 17.020, 16.895 | 16.895 ms |

Sparse and masked updates differ by only 0.014 ms; these runs do not establish
a meaningful winner between them. Both improve on the baseline in every run.
They retain 25,427,976 bytes of setup, versus 27,656,200 for the baseline:
a saving of 2.125 MiB. Workspace remains 72 MiB.

All three new implementations sample the same matrices and produce identical
output hashes. Their output intentionally differs from the original field
mixer. The dense-row implementation is a control: the gain comes from applying
the rank-one structure, not merely changing the sampled matrices.

The rebuilt baseline's binary hash matches the earlier baseline exactly.
[Raw timings and test logs](measurements/mixer/) and
[MIXER_PERFORMANCE.json](MIXER_PERFORMANCE.json) retain the evidence.
Every benchmark ran serially under the executable's shared lock and process
guard. Attempts rejected by the lock produced no timings. A host outage cleared
the earlier temporary builds; these measurements followed a clean rebuild and
new correctness checks.

## Implementation and verification

`generate_mixer.py` generates isolated copies of the implementation. It leaves
all certified sources unchanged. The dense oracle constructs transpose rows
from raw u/v; it does not use grouped hot-path masks or rank-one kernels.
All three variants pass the original quarter-rate test suite at K=2^16, 2^18
and 2^20, including alternate layouts/tiles, independent setup, boundary inputs,
linearity, in-place suffix preservation and setup compaction.

`MixerInner.h` reuses the emission table to compute the rank-one dot product.
The sparse variant updates selected state coordinates. The masked variant uses
a fixed unrolled SIMD sequence. Neither allocates or dispatches dynamically in
the hot path.

## Reproduction

```text
python workstreams/inner_design/generate_mixer.py
cmake -S workstreams/inner_design -B /tmp/inner-mixer-build -DCMAKE_BUILD_TYPE=Release
cmake --build /tmp/inner-mixer-build -j 3
ctest --test-dir /tmp/inner-mixer-build -R mixer --output-on-failure
python workstreams/inner_design/summarize_mixer.py
```

`run_mixer_bench.sh` records the serial comparison when supplied a root containing
the `baseline/` and `build/` build directories. The executables acquire the lock
internally. Do not hold a second lock on the same file while launching them.

Before adoption, close the dense proof and measure the final map/kernel together.
The present 4% result is not a performance claim for an unimplemented redesign.
