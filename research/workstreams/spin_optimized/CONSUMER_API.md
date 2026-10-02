# Consumer API: libOTe and Hypercat

Investigation date: 2026-09-20. The standalone core is now implemented in
[`spin/`](../../spin/README.md). The libOTe integration is implemented and tested;
Hypercat's migration remains planned.
The proposed contract below records the design investigation, not exact API spelling.
The optimized kernels and natural-length changes are committed at `3ced6dc0`.
The original investigation was read-only. Subsequent integration changes are in
`C:/Users/peter/repo/libOTe/libOTe/Tools/Spin/`, with a usage guide in `README.md`.

## Recommendation

Ship one self-contained C++ library, with an immutable code plan and reusable,
caller-owned workspaces. Expose forward and transpose explicitly. Keep the
existing specialized kernels behind that interface; do not rewrite their loops.
Rewrite both consumers against this library. Use its C++ interface directly in
libOTe and an exception-safe C interface with a Rust binding in Hypercat.

Existing consumer APIs are evidence of workloads, not compatibility constraints.
Do not preserve their class names, constructor signatures, workspace ownership,
or wrapper layers solely for backward compatibility. Remove duplicated encoder
implementations as each consumer migrates. Preserve the selected mathematical
maps and fast paths unless a change is separately justified and tested.

The main integration work is packaging, buffer contracts, and ownership.
Neither consumer needs a new encoding algorithm before using these kernels.

## Observed consumers

The local source snapshots inspected were:

| Checkout | Revision | Relevant files |
|---|---|---|
| `C:/Users/peter/repo/libOTe` | `2b686adbf8db1fba1d9d4a19a4c9e93411414e2b` | `libOTe/TwoChooseOne/Silent/SilentOtExt{Sender,Receiver}.cpp`, `libOTe/Vole/Silent/SilentVole{Sender,Receiver}.h`, `libOTe/Tools/CoeffCtx.h` |
| `C:/Users/peter/.codex/worktrees/spin-ot-perf/libOTe` | `910f88337fda607616bf2cd99ce2a0a974226472` | `libOTe/Tools/SpinPerf/SpinPerf.h`, `libOTe/TwoChooseOne/ConfigureCode.h`, `perf/spin_ot.cpp` |
| `C:/Users/peter/repo/hypercat` | `2582a9fb366d28750ef92afdd9b2dbe9538552c6` | `hypercat/src/spin.rs`, `hypercat/src/pcs/brakedown.rs`, `hypercat/src/pcs/brakedown/packed/{spin_layout,natural}.rs`, `hypercat/native/spin/ffi.cpp` |

These are local snapshots, not a claim about the latest remote branches.
The existing libOTe SPIN integration is on the separate performance branch,
not the inspected primary libOTe checkout.

### libOTe

Silent OT expands sparse correlations and then applies the transposed encoder.
The sender compresses a contiguous block buffer in place. The receiver has two
paths: choice bits packed into block low bits, or a separate byte vector.
The separate-vector path applies the same map to blocks and bytes.

The existing SPIN performance adapter already uses `encodeInplace` with prepared
setup and scratch. It fixes K=2^20, T128S19, and packed choices. It is useful
integration evidence, but is not a general-purpose adapter.

Required interface:

- An allocation-free transpose on `osuCrypto::block` buffers. The first K of
  2K elements receive the result; preserve the suffix as the existing API does.
- A characteristic-two generic transpose, particularly `uint8_t`, using the
  same realized code. Two calls with separate workspaces are sufficient initially;
  a fused block-plus-byte traversal is a later performance option.
- Explicit K and N queries, structural alignment checks, reusable scratch,
  and setup memory reporting. Keep original request size separate from code size.
- Reproducible setup across endpoints. Ordinary libOTe compression advances
  `mCodeSeed` after each batch; the SPIN performance branch instead reuses its
  prepared context. The adapter must choose the lifecycle explicitly.
- Synchronous encoding helpers. Protocol coroutines should own the plan and
  workspace across suspension, without putting SIMD scratch into coroutine locals.

Length adaptation belongs to libOTe. Its PPRF partition geometry is not generally
identical to SPIN's geometry. The performance branch fills uncovered coordinates
with zero and charges those coordinates against the effective distance.
Preserve that accounting when generalizing the adapter. Do not silently discard
input coordinates, change K, or reuse a certificate for a different code length.

The inspected VOLE interfaces also support general coefficient contexts.
Our generic XOR implementation is a characteristic-two map, not a replacement
for every prime-field or weighted-code backend. Start with binary OT; decide
additional VOLE profiles separately at the protocol layer.

### Hypercat

Brakedown's `LinearCode` trait requires `dimension`, `length`, `identifier`,
workspace allocation, and an allocation-free forward `encode` on B128 slices.
The code is `Send + Sync`; each active worker owns a separate `Send` workspace.
The verifier also encodes folded messages, so single-stream forward remains
important even when commitments use wide batches.

The SPIN wrapper additionally exposes transpose and bit-packed forward encoding.
Bit-packed forward is exercised by tests and the square-Brakedown example;
it is not a requirement of the main `LinearCode` trait.

Required interface:

- Disjoint K-to-2K forward buffers, fully overwritten output, and reusable scratch.
- Native 128-, 256-, and 512-bit records. For W 128-bit lanes, coordinate i and
  lane l occupy `input[i*W+l]`; output follows the same convention.
- No mandatory packing or row-to-column conversion inside SPIN. Hypercat already
  handles these layouts and can hash retained encoded batches directly.
- An immutable plan shareable across workers, with no global mutable scratch or
  internal thread pool. Rayon already controls parallelism in the ordinary PCS path.
- A stable description of the realized binary map for transcript binding.
- C status codes, opaque handles, checked sizes, and no exceptions across FFI.
- Preservation and testing of the existing bit-packed forward and transpose helpers.

Width has two meanings here. Within the encoder, wider records apply the same
binary map independently to more lanes. In `NaturalSpinPcs`, selecting the batch
width also chooses a polynomial variable partition and changes the PCS identifier.
SPIN must not silently change the caller's record layout when dispatching kernels.

Hypercat currently rejects non-power-of-two K and N in `Brakedown::new`, even if
the encoder accepts them. Generalizing the SPIN Rust wrapper does not generalize
the PCS protocol automatically. Keep those checks in their respective layers.

## Gaps identified in the previous package

The standalone core addresses these packaging, type, ownership, and dispatch
issues. Direction-specific setup pruning remains deferred; the initial plan
prepares both directions. libOTe now uses the package; Hypercat is not yet migrated.

1. **External source dependency.** Forward and wide builds copy Hypercat sources
   and patch them during configuration. A consumer cannot depend on a library
   that needs that consumer's own source tree. Import the required sources once,
   with their licenses and provenance, into this repository.
2. **Two incompatible classes.** `spin_half_transpose` and
   `spin_half_bidirectional` define alternative versions of `bare_spin::Spin`.
   They cannot be linked together. Consolidate the interface and implementation
   behind one target, while retaining direction-specific build options.
3. **Conflicting block definitions.** The transpose build exposes vendored
   cryptoTools headers; the forward import defines a minimal `osuCrypto::block`.
   Neither should be a public dependency of the new core. In particular, the
   minimal definition must not collide with libOTe's real cryptoTools type.
4. **Incomplete cross-direction convenience.** Generic transpose is currently
   generated only for the standalone target. A unified target must retain it.
   Generic forward can follow later; Hypercat's immediate requirement is B128
   and its wide representation, not arbitrary element types.
5. **Lifetime contract.** Our current `WideWorkspace` retains a reference to the
   code. Hypercat's old native wide workspace stores scratch without retaining
   the code, and Rust treats it as independent owned storage. Do not substitute
   the new workspace behind that FFI without fixing the lifetime contract.
6. **Portability and ISA.** Hypercat has MSVC SPIN build rules; our optimized
   package rejects MSVC. Our current baseline also requests x86-64-v3, PCLMUL,
   and VPCLMUL, beyond Hypercat's AVX2 availability check. Preserve a Windows
   path and centralize actual capability checks before replacing its backend.
   Keep ISA flags private to implementation objects, not consumer compile options.
7. **Identity and limits.** Hypercat's wrapper still fixes S19 and exponents
   14..20. Its identifier includes that version, exponent, and both seeds.
   Generalized parameters and actual K need a versioned descriptor. Route hashes
   alone are not full code identifiers.

## Proposed public contract

Use three small concepts:

- `CodeSpec`: actual message size K, a supported parameter-set identifier, setup
  seeds, and the version defining setup expansion and coordinate order.
- `ExecutionOptions`: backend preference, tile choice, and direction requirements.
  These affect performance and memory, not the mathematical map.
- `Code` plus a typed `Workspace`: setup is immutable after creation;
  scratch is owned by the caller and reused for every call.

The default library should finish setup and compaction before returning `Code`.
Users should not need to call `compact` correctly or select Packed24 versus
Indices32. Keep those choices internal or in an expert interface.

Illustrative call shape, not a promised spelling:

```cpp
auto code = spin::Code::create(spec, execution);
auto work = code.make_workspace128();
code.transpose_inplace(buffer_2k, work);   // libOTe
code.forward(input_k, output_2k, work);    // Hypercat's single stream

auto wide = code.make_forward_workspace256();
code.forward_wide(input_records, output_records, wide);
```

Use explicit `forward` and `transpose` names. Today's C++ `encode` means
transpose, whereas Hypercat's Rust `encode` means forward. Rewritten libOTe call
sites should use the explicit names too; no legacy `dualEncode` compatibility
shim is required for SPIN.

At the public boundary, check counts, overlap, workspace compatibility, and
capabilities once per call. Do not resize or allocate scratch during encoding.
Separate-buffer calls reject overlap; the dedicated transpose-in-place call
allows exact aliasing. Forward in-place expansion is not a v1 requirement.

Use scratch-only workspaces that do not retain raw plan references where possible.
If execution views must refer to plan storage, give their owner a stable lifetime
and make Rust retain that owner. Validate geometry, width, and execution layout;
matching K alone is not sufficient once different parameter sets are supported.

Give the core its own storage namespace. Adapters should accept consumer buffers
without whole-buffer conversion. Matching size and alignment alone do not justify
aliasing unrelated C++ struct types: use a defined byte-storage/load-store boundary
or a compile-time adapter with an explicit representation contract.

Expose canonical descriptor bytes so consumers can hash them with their existing
hash implementation. Include every map-defining choice; exclude tile sizes,
route packing, and ISA selection. Freeze cross-backend known-answer vectors.
This requirement concerns compatibility, not a new seed-policy discussion.

Natural lengths remain explicit: K is a positive multiple of 128t for the
currently supplied half-rate maps. Keep representation and allocation limits,
but add no certificate-size or benchmark-size caps. Offer a checked alignment
query; any rounding helper must return a proposed size rather than mutate the spec.
The caller selects parameters and handles certification.

## Implementation order and acceptance checks

1. **Standalone package and unified core.** Export one `spin::spin` CMake target
   and namespaced headers. Vendor the needed forward sources, preserve licenses,
   remove external-checkout requirements, and retain the selected fast paths.
   Generation may remain a maintainer tool; a release should build without Python
   or research workstream dependencies. Build flags must not leak ISA requirements.
2. **libOTe migration.** Replace the benchmark-only wrapper with direct library
   calls. Add zero-copy block transpose and byte fallback, explicit
   lengths, endpoint setup agreement, and reusable workspaces. Test packed and
   separate choices, non-power-of-two aligned lengths, and PPRF tail handling.
3. **Hypercat migration.** Add the C interface and rewrite the Rust binding's
   ownership, actual lengths, map selection, and identifiers. Adapt the PCS call
   sites as needed. Retain B128 and native wide workloads, the bit-packed helper,
   and Windows support, without retaining obsolete wrapper contracts.
4. **Compatibility and performance.** Check forward/transpose adjoint identities,
   byte/block equivalence, wide lane equivalence, invalid buffers, mismatched
   workspaces, and shared-plan calls with independent scratch. Check full Hypercat
   encoded outputs for unchanged specs and OT correlations in both choice
   layouts. Check complete Hypercat proofs under the new binding. If descriptor
   versioning changes the transcript, version it explicitly rather than requiring
   old proof bytes to match. Run allocation checks and serial performance regressions.

Do not add internal threading, a generic runtime callback layer, automatic
certificate selection, or protocol hashing to the core. A fused dual-stream
transpose and arbitrary-width forward fallback can be later additions.

Items 1 and 2 are implemented. libOTe uses optimized block transpose and a
compile-time coefficient-context fallback. Tests passed for all three maps,
both choice representations, minimum and non-power-of-two natural lengths,
PPRF tails, repeated batches, stationary noise, the malicious check, and hashed OT.
The libOTe caller selects only `MultType::Spin`. Its internal configuration picks
the inner and natural aligned size and uses BAA's 0.25 linear-attack tuning value.
This is a heuristic configuration parameter, not a claim of 25% minimum distance.
The 10% distance certificates at selected proof points remain separate evidence.
There is no certificate-range cap. The standalone library still exposes explicit
parameters. SPIN now derives both setup seeds from `mCodeSeed`, follows the usual
per-compression hash advancement, and rebuilds setup when the seed changes.
No Silent VOLE profile was added.

Next, measure one-party libOTe performance, separating expansion, compression,
and optional hashing, excluding setup and base OT. Then implement item 3.
