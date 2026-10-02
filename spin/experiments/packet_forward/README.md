# RS packet forward encoding

This experiment developed fast forward encoding for the exact RS packet map
in SPIN revision `36394342820da99291cc736165ff13129a2ca359`. It preserves the
construction, sampled matrices, seeds and coordinate order. The kernels are now
promoted into `src/packet/PacketForward*.cpp` and exposed through both
`Code::forward_bytes` and `PacketCode::forward_bytes`, including their typed
`forward` wrappers. This directory retains the exploratory benchmark and its
original measurements; it links the production kernels rather than copies.

## Public API promotion

The public API was compared directly against the saved experimental executable,
with the same CPU, frequency, inputs, two seeds, 501 calls, and both execution
orders. All processes ran serially under the four shared benchmark locks.
Median of four process medians:

| K | Saved experiment | Public forward | Public transpose |
|---|---:|---:|---:|
| 2^16 | 0.167122 ms | 0.1705235 ms | 0.210733 ms |
| 2^18 | 0.6914005 ms | 0.669454 ms | 0.840213 ms |
| 2^20 | 3.3363085 ms | 3.2759655 ms | 3.237153 ms |

Public and experimental forward checksums match for every size and seed.
The public path retains the experimental performance; these small differences
do not establish an intrinsic speedup from the API promotion. Setup/allocation
are excluded. AVX-512 setup now also retains the forward coefficients, state
updates and inverse route, accounted for in `setup_bytes()`. Workspace size
is unchanged. The portable backend needs no additional prepared forward tables.

The production `spin_packet_forward_api_bench` target uses public headers and
`Code::forward_bytes`. It is available with `SPIN_BUILD_BENCHMARKS=ON`, without
enabling experiments:

```sh
spin_packet_forward_api_bench 1048576 1 501 huge
```

Public API tests cover all four permitted cache-line offsets, dense/zero/sparse
inputs, sizes on both sides of the layout threshold, the full K20 case, input
preservation and output guards, memory policies, dirty workspace reuse, seed
rebinding, overlap/size/alignment errors, moved/copied handles, concurrent calls
with separate workspaces, and allocation-free encoding. The six Linux release
tests, the 36 literal comparisons and adjoint checks, the Windows/MSVC packet
test, and a public-header-only installed-package consumer pass. The expanded
packet tests also pass with the entire library instrumented by ASan/UBSan
(including leak detection), and in a separate AVX-512-disabled release build.

Promotion receipts are in
`output/spin-forward-promotion-20261002/paired-20261002-010643/` at the repository
root. The original experiment receipts below remain separate.

## Measurements

Serial paired runs on Peach, Ryzen 9 7950X CPU15 at 4.5 GHz, boost disabled,
GCC 15.2, Release, `SPIN_TUNE=znver4`:

| K (128-bit records) | Public transpose | Selected forward | Forward / transpose |
|---|---:|---:|---:|
| 2^16 | 0.2139445 ms | 0.167969 ms | 0.785x |
| 2^18 | 0.834624 ms | 0.688726 ms | 0.825x |
| 2^20 | 3.2855805 ms | 3.329849 ms | 1.013x |

Each entry is the median of four process medians: seeds 1 and 17, each in
both execution orders. Each process excludes five warmups and measures 501
calls. All 24 benchmark processes ran sequentially under the shared locks.
These are precomputed encoding times; setup and allocation are excluded.
Normal allocation is used at K16/K18, and PreferHugePages at K20, matching
the earlier transpose comparison. Huge-page advice does not guarantee huge pages.

The transpose executable is the unchanged public `spin_packet_bench`. It
updates its buffer in place, as in the earlier published measurement. Forward
reuses a separate fixed message and output; it consumes K records and produces
2K, while transpose consumes 2K and produces K. Both use data-independent maps.
Forward is 17–21% faster at the smaller sizes and 1.35% slower at K20.
The K20 result is therefore approximately transpose-speed, not a measured win.
These are encoder timings, not FLOCK or PCS measurements.

Raw CSV summaries, all individual forward timing samples, resource receipts,
commands, binary hashes and source hashes are retained in
`output/spin-packet-forward-20261001/paired-20261002-005148/` at the repository
root. The timestamp is the host's UTC directory name; the experiment date is
October 1 in the user's timezone. `source.tgz` in the parent directory is the
exact measured source snapshot, and `measurements.tgz` preserves the receipts.

## Implementation

`tools/generate_packet_forward.py` transposes the retained factorized RS circuit and reverses its
byte/bit packing permutations. Each GF(2^8) multiplier is replaced by its
binary adjoint, prepared as a GFNI matrix during setup. The GF(2^32) tower
factorization remains intact. This avoids evaluating a dense 32-by-32 map
for every outer symbol.

The inner uses the same unrolled expansion and feedback schedule. Setup
transposes the sampled GL20 matrices and converts them into the existing
packed state basis. Forward processes epochs in increasing order. The
expansion columns are degree-at-most-two evaluations on six variables, so
the expansion matrix is self-orthogonal. Consequently, feeding back the
expanded output gives the same syndrome as feeding back its raw input.
This permits reuse of the fused output/feedback schedule without changing
the recurrence. The literal oracle independently computes feedback from raw input.

`forwardSelected` chooses between two memory layouts before encoding:

- Through K=2^18, the outer writes contiguous grouped scratch; the inner
  gathers four-record packets and streams contiguous output.
- Above that threshold, the outer scatters complete cache lines through a
  prepared inverse route using non-temporal stores. The inner then processes
  the output in place, with sequential reads, cached writes and write
  prefetches eight epochs ahead. An explicit fence orders the outer stores
  before those reads.

Both paths use fixed-width SIMD schedules and preallocated storage. They
allocate no memory, create no threads and perform no runtime callback dispatch
during encoding. Public calls require disjoint message/output buffers and
16-byte alignment. `spin::Buffer` supplies 64-byte alignment, enabling the
measured direct-output path at large sizes. Other output alignments use owned
scratch. Portable encoding also uses that workspace without allocating.

## Validation and reproduction

Validation covers 36 combinations of six lengths (including non-power-of-two
K=768 and the full K20), two distinct seed pairs, and dense/zero/sparse input.
Each compares the gathered, scattered, in-place and selected paths against
the independent literal `forwardScalar`, including recycled dirty scratch.
Additional dot-product checks verify the forward/transpose adjoint identity.
The separate-buffer paths exercise 16-byte output alignment and suffix guards.
All six existing library tests pass serially.

The original experiment translation units also passed AddressSanitizer and
UndefinedBehaviorSanitizer with leak detection. In that check, the pre-existing
`libspin.a` is not instrumented; the new forward kernels, setup and harness are.
Validation logs are retained beside the measurement archive.

```sh
python3 spin/experiments/packet_forward/generate.py
cmake -S spin -B out/forward -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BUILD_EXPERIMENTS=ON -DSPIN_BUILD_BENCHMARKS=ON -DSPIN_TUNE=znver4
cmake --build out/forward --target spin_packet_forward_bench spin_packet_bench -j4
out/forward/spin_packet_forward_bench verify
```

The experimental target requires Linux/GNU-compatible x86 compilation and an
AVX-512F/VL/BW/DQ, VBMI and GFNI CPU. Run measurements only one at a time:

```sh
flock /tmp/prindal-addition-encoder-benchmark.lock \
flock /tmp/bare-spin-benchmark.lock \
flock /tmp/hypercat-benchmark.lock \
flock /tmp/hypercat-global-benchmark.lock \
taskset -c 15 out/forward/spin_packet_forward_bench selected 1048576 1 501 huge
```

`compare.py`, `run_remote.py` and `check_remote.py` retain the exact Peach
commands with `/tmp/spin-forward-20261001` as their working directory.
`SPIN_FORWARD_RAW=/path/to/trials.csv` saves every measured forward call.
The generator leaves the production transpose sources untouched.

Next: integrate the public forward API into FLOCK and measure end-to-end.
