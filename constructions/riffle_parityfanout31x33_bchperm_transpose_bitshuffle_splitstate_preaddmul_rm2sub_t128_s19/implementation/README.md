# Structured SPIN (B=256, t=128, s=19) for libOTe

This directory is the buildable, non-frozen cleanup of the exact frozen
Structured SPIN implementation. Its canonical layout is a libOTe source
overlay; the root CMake project is a standalone verification fallback. It does
not replace or modify `frozen_source/`.

## libOTe layout

- `libOTe/Tools/RiffleCode/StructuredSpinB256T128S19.h` is the public header.
- `StructuredSpinB256T128S19.cpp` compiles the AVX2/VPCLMUL kernel into
  `libOTe`; ISA flags do not leak to downstream targets.
- The public object and `Workspace` use compile-time-checked fixed inline
  storage. There is no PIMPL allocation, virtual dispatch, callback wrapper,
  or per-call allocation.
- `generated/Rm2SubS19Tables.h` contains immutable generated RM2Sub data.
- `RiffleExactPermFieldCheckpoint.h` is the byte-identical frozen production
  kernel. `RiffleRm2SubS19.h` retains the frozen hot body.
- `libOTe_Tests/` contains a native `CLP` test adapter and the independently
  staged correctness oracle.
- `integration/libote-root.patch` adds the implementation and test source to
  upstream libOTe with source-local ISA flags.
- `vendor/cryptoTools/Common/` exists only for the standalone fallback.

See `integration/README.md` for overlay and `oc::libOTe` linking instructions.

The setup object owns the permutation, packed routing schedules, fanout
schedules, and RM2Sub coefficient schedules. `Workspace` owns the per-call
40 MiB scratch area. Neither checked nor unchecked encode allocates.

## Standalone build on Peach

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_ARCH=znver4
cmake --build build
ctest --test-dir build --output-on-failure
./build/structured_spin_correctness
```

The public-header smoke target deliberately compiles without AVX2 or VPCLMUL
flags. The implementation, correctness oracle, and benchmark use C++20,
`-O3`, `-march=znver4`, `-mavx2`, `-mvpclmulqdq`, and `-mpclmul`.

The correctness output must end with `correctness=PASS`; the canonical output
checksum must be `0x0af28e9c011dfa1a`, and the historical aggregate must be
`0x95c9d722a9539fef`.

## Benchmark

Never run this concurrently with another benchmark. Immediately before each
run, inspect executable paths and abort if any benchmark is active:

```bash
busy=0
for p in /proc/[0-9]*; do
  e=$(readlink "$p/exe" 2>/dev/null || true)
  b=${e##*/}
  case "$b" in
    *Bench*|*bench*|*Benchmark*|*benchmark*)
      echo "ACTIVE_BENCHMARK:$p:$e"
      busy=1
      ;;
  esac
done
test "$busy" -eq 0
SPIN_BENCH_CPU=15 ./build/structured_spin_benchmark 21
```

`SPIN_BENCH_CPU` defaults to 15. The trial count must be an odd integer at
least three.

## Interface

Include:

```cpp
#include <libOTe/Tools/RiffleCode/StructuredSpinB256T128S19.h>
```

Call `init` with the permutation seed, setup-coefficient seed, and RM2Sub
coefficient seed. The checked `dualEncodeTo` validates initialization and
fixed input/output sizes once. `dualEncodeUnchecked` is for callers that
already enforce those invariants.

`retainOracleSchedules` defaults to false. Tests set it to true only so the
independent staged oracle can inspect the inverse route and fanout schedules
through the separate test-access header.

See `SOURCE_MANIFEST.json` for exact source hashes and
`workstreams/implementation_cleanup/` for correctness, integration,
performance, and operation receipts.
