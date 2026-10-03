# K16 paired-s15 profile

`Parameters::PacketRsT64S15K16` selects the exact code from the frozen
`k16-rs16-paired15-v1` research checkpoint. It accepts only K=65,536 and maps
K message records to 2K encoded records. The 128-, 256-, and 512-bit forward
interfaces apply the same binary map independently to each bit of a record.
Transpose and in-place transpose accept 128-bit records.

```cpp
spin::Code code({65536, spin::Parameters::PacketRsT64S15K16, 1, 2});
auto workspace = code.make_workspace(spin::Width::Bits256);
// message has 65536*32 bytes; encoded has 131072*32 bytes.
code.forward_bytes(message, encoded, workspace);
// Use ordinary stores when the result will be read immediately.
code.forward_bytes(message, encoded, workspace, spin::OutputStores::Cached);
```

The public descriptor uses parameter identifier 6. The prior s20 packet
profile retains identifier 5 and its existing seeded map. Callers and proof
systems must bind the selected descriptor; the two codes are not interchangeable.
Full-mode `PreparedEncoder` supports this profile's transpose and setup refresh.
Banked setup, generic transpose, and bit-packed forward are unsupported.

Each outer group encodes 128 binary coordinates into 256 coordinates using
four GF(16) RS[16,8] rows and independent nonzero GF(2^16) symbol multipliers.
Four-coordinate packets are routed through 64 regions. The inner uses
64-coordinate steps, the frozen paired 15-coordinate expansion and its
transpose feedback, and invertible 15-by-15 state updates. State starts at
zero, persists across regions, and is discarded after the last step.

`route_seed` drives the original packet shuffles and the outer-symbol stream;
`inner_seed` drives the GL15 stream. Setup preserves the frozen research seed
derivation and coordinate order. The original ideal-randomness distance result
and its scope are recorded in the [research checkpoint](experiments/k16_codesign_100us/checkpoint/README.md).
Promotion adds an implementation of the same map, not a new distance claim.

`Automatic` selects AVX-512/VBMI/GFNI when available and the baseline SSE2
implementation otherwise. Setup and portable translation units do not require
AVX-512. Forward loads read interleaved records directly. Each lane has an outer
scratch plane, and inner output is interleaved one 64-coordinate step at a time.
There are no whole-message layout conversions or allocations inside encoding.
Foreign input/output buffers require 16-byte alignment and must not overlap.
The default `OutputStores::Streaming` preference streams final output when the
fast backend and 64-byte alignment permit it, with a cached fallback otherwise.
Pass `OutputStores::Cached` to force ordinary output stores at every alignment.
The same optional argument is available on typed `forward` calls. It is selected
per call, so one code and workspace can serve both large batches and immediately
consumed results. Policies produce identical bytes and do not enter descriptors.
Output is ready to read on return; this profile's outer scratch stays cached.
Streaming avoids the output traffic regression observed when encoding complete
FLOCK matrices, while cached stores can suit a reused standalone output buffer.
Workspace scratch is owned, aligned, and reusable across compatible setups.

The standalone `spin_paired15_test` checks scalar/SIMD agreement, all forward
widths, adjoint identities, alignment, guards, workspace reuse, and allocation
behavior. It also checks both output policies against the default, typed calls,
policy changes on reused scratch, and rejection of invalid policies before writes.
A repository checkout can additionally enable
`SPIN_BUILD_RESEARCH_ORACLE_TESTS=ON` to compare both directions against the
frozen research setup and scalar implementation. Ordinary builds and installed
consumers do not depend on research files or Python.

Build with `SPIN_BUILD_BENCHMARKS=ON` for the public-API comparison driver:

```sh
spin_paired15_bench paired15 forward 128 1 301 normal
spin_paired15_bench old forward 128 1 301 normal
```

The arguments select profile, direction, record width, setup seed, measured
calls, and memory policy. Setup and allocation are excluded. Run benchmarks
serially and compare the same direction, width, and memory policy in both orders.
The [promotion workstream](../research/workstreams/k16_forward_promotion/README.md)
records forward and FLOCK validation separately from the prior transpose timings.
