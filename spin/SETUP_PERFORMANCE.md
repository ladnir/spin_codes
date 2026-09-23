# Setup baseline: exact-map optimization

The first setup optimization pass preserves descriptor version 1 and its exact
seed-to-code mapping. It does not use a permutation bank, cyclic shifts, or
reused transvections. The encoder circuits and parameters are unchanged.

## What changed

- Generate SplitMix words eight at a time. Each consumer receives exactly the
  previous stream, including words consumed by rejection sampling. This also
  batches the generation of IMT transvection masks.
- Precompute a reciprocal and rejection threshold for each shuffle bound.
  Multiplication and one correction replace division inside the sampling loop,
  while producing the same remainder as the previous sampler.
- Construct only the routing-index width retained by the public plan. Previously,
  setup built both packed and 32-bit versions, then discarded one.
- Add `Code::prepare_workspace` to rebind compatible 128-bit scratch to a new
  plan. It does not allocate or clear the buffers. Other shapes and wide
  workspaces are recreated as needed.

The public plan still supports both forward and transpose encoding. Removing
forward-only tables from a dedicated transpose plan is not part of this pass.
Neither the underlying random generator nor the protocol seed policy changes.

## Measurements

Ryzen 9 7950X, CPU 15, Release, GCC 15.2, `SPIN_TUNE=znver4`, AVX-512 dispatch.
Baseline: public commit `b5a81b1d7fc98626419af184ce354cc2d4cf1b06`.
The candidate is the setup change described above. Each encoder handles one
128-bit record stream in place. K=2^16 uses T64S12R2; larger points use T128S19.

All times below include constructing a new seeded plan, preparing scratch, and
the first transpose call. They exclude input preparation, PPRF expansion,
communication, base OTs, hashing, and the separate byte-choice encoder.

| K | Old fresh setup + encode | New fresh setup + encode | New changed-seed batch, reused scratch |
|---:|---:|---:|---:|
| 2^16 | 4.031 ms | 2.941 ms | 2.189 ms |
| 2^18 | 19.360 ms | 14.174 ms | 8.444 ms |
| 2^20 | 81.383 ms | 63.015 ms | 39.550 ms |

Fresh-call latency falls by approximately 27%, 27%, and 23%. With scratch
reuse, the measured totals are approximately 1.84x, 2.29x, and 2.06x faster
than the old rebuild-every-batch path. Reuse is an explicit library API, not
yet enabled in the pinned libOTe integration. It does not eliminate first-use
workspace allocation.

Fresh plan construction alone:

| K | Old | New |
|---:|---:|---:|
| 2^16 | 3.030 ms | 1.663 ms |
| 2^18 | 12.950 ms | 7.879 ms |
| 2^20 | 51.550 ms | 33.104 ms |

Warm encode-only medians were 0.369/1.526/10.884 ms for the baseline and
0.368/1.496/10.859 ms for the fresh candidate. These measurements do not
establish an encoder-kernel speedup; setup is the intended change.
In the reuse experiment, warm encode-only medians were 0.368/1.506/11.092 ms.

The reuse experiment retains the preceding plan through its scratch handle
until rebinding. Rebinding can therefore include destruction of the preceding
plan, although it does not allocate scratch. The different allocation lifecycle
also changes allocator/cache behavior; do not derive its total by subtracting
the fresh workspace time from a fresh-call total.

Each reported value is the median of three process medians. Each process uses
nine changed-seed setups, excludes its first two samples, and measures one first
encode plus three warm encodes per setup. Process order reverses in the middle
repetition. All benchmarks ran serially under the shared encoder benchmark locks.
Component medians need not sum to the median total.

Every baseline/candidate/reuse process produced these matching output checksums:

| K | FNV regression checksum (not a cryptographic hash) |
|---:|---|
| 2^16 | `f19c7bd9550b4f9b` |
| 2^18 | `0f45c155c345f35c` |
| 2^20 | `0762faeec8ee1a0f` |

## Reproduction and checks

```sh
cmake -S spin -B out/setup -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_TUNE=znver4 -DSPIN_BUILD_BENCHMARKS=ON
cmake --build out/setup -j2
ctest --test-dir out/setup --output-on-failure
# Run only when no other benchmark is active.
out/setup/spin_setup_bench
out/setup/spin_setup_bench --reuse
```

For a matched comparison, build the old library in `baseline/build`, then compile
the current `bench/setup.cpp` against it without defining `SPIN_BENCH_REUSE`.
Name that executable `baseline/bench`. Build the candidate benchmark in
`candidate/build`, then run `bench/compare_setup.sh COMPARISON_ROOT CPU`.
That script acquires the three shared benchmark locks and alternates process
order. Raw CSV files remain outside version control.

`tests/setup.cpp` compares the batched stream with the previous scalar generator,
tests reciprocal remainders at arithmetic boundaries, forces a rejected draw,
compares complete shuffles, and compares 100,000 transvection pairs per state
size. Existing descriptor-v1 known answers remain unchanged. API tests compare
dirty rebound scratch with fresh scratch, reject its former owner, and cover
resizing and moved-from workspaces. The importer reapplies the reviewed setup
overlay and updates the kernel generation record.

Release API, sampler, and known-answer tests passed with GCC 13.3 on the local
AVX2 host, GCC 15.2 on the Ryzen AVX-512 host, and MSVC 19.50 on Windows.
The complete AVX2 API, sampler, and known-answer suite also passed ASan/UBSan.

## Next step

Setup still dominates the K=2^18 end-to-end encoder cost. Next, remove
forward-only preparation for transpose consumers and profile routing generation
separately from routing-table construction. Wire scratch reuse into libOTe when
updating its dependency. Only then compare a smaller permutation family against
this baseline; that alternative changes the sampling model and needs separate
analysis.
