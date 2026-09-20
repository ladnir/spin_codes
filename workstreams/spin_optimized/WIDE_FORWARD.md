# Wide forward SPIN: equal-work comparison

This records the initial S19 comparison. The subsequent
[K16 extension and K20 tile tuning](WIDE_TUNING.md) improve the K20 result and
add support for the selected two-round K16 inner.

Hypercat's 256-bit forward kernel is worth retaining as an optional batched path.
On our Ryzen host, it improves encoding-plus-column-assembly throughput by
1.30–1.36x at K=2^14 through 2^18 when the input already has the wide layout.
At K=2^20, the improvement drops to 1.10x. Packing separate input planes nearly
eliminates that last gain. The 512-bit kernel does not beat 256-bit consistently.

These are half-rate BCH [256,128] measurements with the fixed (t,s,r)=(128,19,1)
IMT inner. They do not compare against the newer K16-specific T64S12R2 inner,
which is covered by the follow-up report. Nothing here changes
the selected parameters, the transposed implementation, or the submitted paper.

## Workload and results

One trial encodes **16 independent 128-bit planes**, or 2048 binary messages,
under one shared setup. Every binary message has length K and encodes to 2K bits.
The 128-bit path makes 16 calls, the 256-bit path makes eight, and the 512-bit
path makes four. Every path produces the same column-major output, checked
byte-for-byte against independent scalar column assembly after each trial.

The columns contain whole 128-bit words. This is a row-to-column layout change,
not an input bit transpose. We include this assembly in every total below.
The 128-bit control uses Hypercat's existing AVX-512 streaming 4x4 assembly.
The wide cases use the assembly from its `spin_wide.rs` comparison.

Total milliseconds for the same 16-plane workload, with input already in each
kernel's native layout:

| K | 128-bit | 256-bit | 512-bit | 256-bit throughput gain |
|---:|---:|---:|---:|---:|
| 2^14 | 2.278 | 1.680 | 1.974 | 1.36x |
| 2^16 | 10.296 | 7.880 | 8.956 | 1.31x |
| 2^18 | 46.569 | 35.938 | 39.902 | 1.30x |
| 2^20 | 202.461 | 184.890 | 200.226 | 1.10x |

If input instead arrives as 16 separate planes, the wide path must interleave
two or four planes before encoding. The following totals include that work:

| K | 128-bit control | Pack + 256-bit + columns | Gain | Pack + 512-bit + columns | Gain |
|---:|---:|---:|---:|---:|---:|
| 2^14 | 2.278 | 1.753 | 1.30x | 1.756 | 1.30x |
| 2^16 | 10.296 | 8.627 | 1.19x | 9.418 | 1.09x |
| 2^18 | 46.569 | 40.570 | 1.15x | 45.692 | 1.02x |
| 2^20 | 202.461 | 199.747 | 1.01x | 221.248 | 0.92x |

Both input cases are measured directly. Packing also warms the input cache;
subtracting its duration from the second table would not reproduce the first.

The gain mostly comes from encoding, not column assembly. At K=2^18, native
128-bit and 256-bit encoding take 39.16 and 28.85 ms; assembly takes 7.38 and
7.06 ms. At K=2^20, those encoding times are 172.18 and 156.06 ms, and assembly
takes 30.27 and 28.94 ms. Phase medians need not sum to the median total.

The wider records also enlarge the encoder's scratch space. At K=2^20,
the 128/256/512-bit workspaces occupy 40/80/160 MiB. At K=2^18 they occupy
9/18/36 MiB. These figures exclude setup, inputs, outputs, and assembly buffers.
The smaller gain at K=2^20 is consistent with greater memory pressure;
these timings alone do not identify which cache or bandwidth limit dominates.

## Measurement scope

Measurements use Peach's Ryzen 9 7950X, GCC 15.2, `-O3 -mtune=znver4`, and CPU 15.
There are three independent processes per K, each with seven trials after one
warmup of each case. Case order alternates. Tables report pooled medians over
21 samples. No benchmarks run concurrently; both shared benchmark locks are held.

At K=2^20, the three run medians range from 179.30 to 189.54 ms for native 256-bit,
versus 201.45 to 202.56 ms for 128-bit. With input packing, the 256-bit medians
range from 196.81 to 206.74 ms. Treat that packed-input result as a tie, not a win.

Setup, allocations, page initialization, and output checking are outside the
timers. Each width reuses one workspace across its calls. Native inputs are
prepacked outside the timed region; the separate-plane cases pay for packing
inside it. Output assembly uses streaming stores and fences before stopping
the timer. This measures encoding plus assembly, not a full PCS commitment:
hashing, folding, and any protocol-level layout changes are absent.

The routing and coefficient seeds are 1 and 2. The test uses the current tile
defaults unchanged. These measurements cover four representative lengths,
not the proposed full K=2^12 through 2^26 range. K=2^14 is a performance sample,
not an additional distance certificate.

## Integration and replay

`wide.py` copies the upstream wide kernels into the generated build tree.
It regenerates the S19 outer circuit, then lifts the integrated forward
recurrence and K16 map to the requested vector width. The original XOR circuits
and unrolled recurrence remain intact; map dispatch occurs once per call.
`WIDE_IMPORT.json` records source hashes; `WIDE_MANIFEST.json` records the
generator inputs and generated header hash. The generated tree also retains
the upstream `libote-LICENSE` for the imported column-assembly code.

```sh
cmake -S workstreams/spin_optimized -B out/spin-wide \
  -DCMAKE_BUILD_TYPE=Release -DSPIN_TUNE=znver4 \
  -DSPIN_BCH_AVX512=ON -DSPIN_BUILD_WIDE=ON -DSPIN_BUILD_BENCHMARK=ON \
  -DSPIN_BIDIRECTIONAL_SOURCE=/path/to/hypercat/hypercat/native/spin
cmake --build out/spin-wide -j2
ctest --test-dir out/spin-wide --output-on-failure -j1
out/spin-wide/spin_wide_benchmark 18 7
```

Link `spin_wide`, which brings in `spin_half_bidirectional`, and include `Wide.h`.
Do not also link the standalone transpose library; it defines the same `Spin`.

```cpp
bare_spin::Spin code(bare_spin::Configuration::T128S19, 18);
code.compact(bare_spin::Layout::Packed24);
bare_spin::WideWorkspace workspace(code, 2); // 256-bit logical elements
// input has 2*K blocks; output has 2*(2*K) blocks.
// Coordinate i occupies input[2*i], input[2*i+1].
workspace.forward(input, output); // spans of 128-bit blocks, disjoint buffers
```

The setup must outlive the workspace and remain unmoved and unmodified.
Workspaces are not shared between concurrent callers. Constructors check CPU
features before entering ISA-specific code, including workspace allocation.
The public wrapper and 256-bit kernel stay at the integration's existing AVX2
baseline. The separate 512-bit object requires AVX-512F/VL/BW/DQ; no LTO crosses
that boundary. The comparison's column-assembly kernels require AVX-512 too.

The width is a batch size, not a larger extension field. Each lane applies the
same binary generator to an independent message. Reusing the same setup does
not resample or modify the code. No new distance argument is needed for this
width lifting; the existing parameter restrictions still apply.

`wide_check.sh` reproduces the serial sweep in the Peach snapshot layout.
`wide_validate.sh` checks the masked-AVX-512 fallback and sanitizer builds.
`wide_summary.py` summarizes `measurements/wide/m*-r*.csv` and run variability.
Raw trials, logs, hashes, and archives remain ignored rather than committed.

Validation passed: all eight release tests; lane-by-lane equality with the
128-bit encoder at K=2^14 and 2^16; compacted and uncompressed setup;
default and two-row tiles; 16-byte-aligned inputs; malformed lengths and overlap
rejection; AVX-512-masked fallback; AddressSanitizer and UndefinedBehaviorSanitizer
checks of both wide encoders and all three column-assembly paths. The assembly
test also checks every timed case against the scalar reference through K=2^20.
Object-code inspection found AVX-512 instructions only in the designated objects.

## Next step

Keep 256-bit as an explicit option when callers already batch two planes.
Do not make 512-bit the default or automatically repack large planar inputs.
The K16 extension and K20 tile sweep are now complete; see [WIDE_TUNING.md](WIDE_TUNING.md).
Next measure the full Brakedown commitment pipeline before choosing an
application-level size/width dispatch policy.
