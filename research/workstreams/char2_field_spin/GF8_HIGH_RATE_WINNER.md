# Frozen high-rate winner

The user accepted the 3.2 ms K20 implementation on 2026-10-05. This ends the
current high-rate optimization round. Preserve this checkpoint while the
active search moves to high distance at rate 1/2. It remains research code;
the production library and paper are unchanged.

| Property | Selected checkpoint |
|---|---|
| Outer | GRS[64,40] over GF256, polynomial `0x11b` |
| Inner | Full-Cauchy four-input/four-state map, retaining three outputs |
| Routing | One independently sampled row permutation per outer coordinate |
| Requested K | 1,048,576 |
| Padded K / output length | 1,049,600 / 1,259,520 |
| Nominal / requested rate | 83.3333% / 83.2520% |
| Certified distance | At least 10% nonzero symbols |
| Whole-code setup-failure margin | 57 bits at the stated K20 profile |
| Payload | 128-bit records containing sixteen parallel GF256 lanes |
| Selected transpose | `RateSpeed4RoutedPlan<40>` and `rs4TransposeRouted16` |
| Precomputed encoding time | About 3.2 ms on Peach's Ryzen 9 7950X |

The [construction and certificate](GF8_SPEED_DESIGN.md) specify the setup
distribution and cover every nonzero message. The [matched native study](GF8_SPEED_ROUTED.md)
records 3.184--3.201 ms; the [follow-up audit](GF8_SPEED_PROFILE.md) measured
3.205 ms and found no material improvement among sixteen further schedules.
Restoring input before every timed call gives about 3.33--3.38 ms. Setup,
allocation, and restoration are excluded; the full encoding operation and
output stores are included. This is not an application-throughput claim.

## Immutable local source snapshot

[The source archive](../../.proof-archives/gf256-high-rate-k20-winner-20261005-source.zip)
contains the complete experimental CMake project, its source generators,
the proof scripts and tests, and the existing research notes. Each of its
459 source files has a SHA-256 entry in `MANIFEST.sha256`. Every archived
entry was read back and checked against that manifest. The archive excludes
raw measurements, numerical dumps, caches, binaries, and build directories.

Archive size: 1,562,217 bytes. Archive SHA-256:

```text
8ee1164660a486048316cfee33ca8e950a878e10c5e045feedfdf37c4614706
```

The archive is local and ignored by Git; it has not been published or
committed. Its directory tree preserves repository-relative source paths.
The selected implementation's three scheduling headers have these hashes:

```text
f9e0caa790955f090567a6e477022d62287c62c9c5014d20a1b0a079813f8b78  RateSpeed4RouteNext.h
8661a046241568f2fb0f2b92575c428fcb7d1b35e3339f0f9ca65470b1b85d12  RateSpeed4InnerNext.h
ca91210d38b86fac0c896075e2919a411e33f118a9ca139fe9c79a95aa8fe9aa  RateSpeed4OuterNext.h
```

## Reproduction

The original native validation used GCC 15.2 with AVX512/GFNI, pinned to
CPU 2 with boost disabled. Release and ASAN/UBSAN checks passed for the
selected entry point, including natural sizes, aliases, alignment offsets,
and wide-route fallback. Never run these benchmarks concurrently.

```text
cmake -S spin/experiments/char2_field_spin -B build/field -DCMAKE_BUILD_TYPE=Release
cmake --build build/field --target spin_char2_rate_speed4_next -j2
build/field/spin_char2_rate_speed4_next --check
build/field/spin_char2_rate_speed4_next --large-check
taskset -c 2 build/field/spin_char2_rate_speed4_next 1048576 151 selected
taskset -c 2 build/field/spin_char2_rate_speed4_next 1048576 151 selected-fresh
python research/workstreams/char2_field_spin/theory/gf8_speed_certificate.py
```

The proof replay depends on NumPy, SciPy, and python-flint. The environment
recorded at freezing has Python 3.14.4, NumPy 2.4.4, SciPy 1.18.0, and
python-flint 0.9.0. The source includes all numerical witnesses needed by
the selected certificate; it does not read an external proof dump.

New high-distance experiments must use separate files. Do not reinterpret
this checkpoint's distance or timing as evidence for those constructions.
