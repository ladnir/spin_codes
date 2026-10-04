# Small-outer s22 profile

`Parameters::PacketRsT64S22` provides the optimized small-outer packet code
through the ordinary `Code` interface. It maps K message records to 2K encoded
records. Its primary operation is precomputed transposed encoding of 128-bit
records, including in-place encoding.

```cpp
spin::Code code({262144, spin::Parameters::PacketRsT64S22, 1, 2});
auto work = code.make_workspace();
auto encoded = code.make_buffer();
// Populate all 2K records before encoding.
code.transpose_inplace_bytes(encoded.bytes(), work);
// The first K records hold the result; the remaining records are unchanged.
```

Setup is fixed by the immutable code. Encoding neither samples a new code nor
allocates scratch. Full-mode `PreparedEncoder` also supports the profile and
can rebuild its setup for new seeds while reusing compatible workspaces.
Banked heuristic setup, generic XOR transpose, and bit-packed forward are unsupported.

## Construction and certificate scope

Each outer group maps 128 binary coordinates to 256 coordinates using four
systematic GF(16) RS[16,8] rows. Independent nonzero GF(2^16) symbol multipliers
randomize the resulting symbols. Four-coordinate packets pass through 64
shuffled regions. The inner processes 64 coordinates per step, with the exact
22-coordinate expansion and feedback maps retained by the K18 certificate.
Invertible 22-by-22 updates mix its persistent state between steps.
The state starts at zero, continues across region boundaries, and is discarded
after the final step.

The profile retains the research setup's seed derivation and coordinate order.
The route seed determines the packet shuffles and outer-symbol stream; the
inner seed determines the GL22 stream. The descriptor uses parameter identifier
7. Identifiers 5 and 6, their seeded maps, and all default selections remain unchanged.

The complete certificates bound failure over the specified independent setup
ensemble. They certify the following whole-code distances, not the observed
distance of a particular seeded implementation:

| K | Minimum distance | Relative distance | Setup-failure margin |
|---|---:|---:|---:|
| 2^18 | 52,429 | At least 10% | 60.7509985457 bits |
| 2^20 | 209,716 | At least 10% | 58.7503767402 bits |

The [certified baseline index](../research/CERTIFIED_BASELINES.md) records the
proof sources and receipt identities. Promotion preserves these maps; it does
not extend the certificates to other lengths. The larger-outer s20 profile
remains available with its 68.103757-bit K20 certificate.

## Natural sizes and execution

The profile accepts positive multiples of 128. The largest representable K is
132,152,832, as limited by 32-bit byte offsets into padded scratch. Use
`valid_message_size` rather than duplicating this bound. Available memory may
impose a smaller practical limit. No power-of-two restriction or rounding applies.
The inner runs over the flattened route, so region boundaries need not coincide
with step boundaries. `tile_rows` must be zero or four.

`Automatic` selects AVX-512/VBMI/GFNI when available and portable SSE2 otherwise.
The fast transpose carries packed state between steps. It uses cached route and
outer stores through K=2^18, and streaming route stores plus unrolled outer mixing
above that size. Streaming final-output stores require 64-byte alignment;
16-byte-aligned outputs use the cached fallback. Output is readable on return.
These execution choices do not change the code or its descriptor.

Forward encoding accepts 128-, 256-, and 512-bit records through an allocation-free
scalar implementation. It is a correctness fallback, not a promoted fast forward
kernel. Each record lane uses the same binary map. The forward output-store
preference is accepted, but this fallback uses cached stores.

## Validation and reproduction

The 2026-10-03 promotion campaign measured the public API against the selected
research kernels in the same executable:

| K | Public API | Research kernel | Public/research difference |
|---|---:|---:|---:|
| 2^18 | 0.465073 ms | 0.463941 ms | +0.24% |
| 2^20 | 2.955466 ms | 2.943754 ms | +0.40% |

These are medians of four process medians: seeds 113 and 211, each in both
comparison orders. Every process used five warmups and 1001 timed in-place calls.
The Ryzen 7950X run used core 15, GCC 15.2, `SPIN_TUNE=znver4`, and explicit
large-page preference. All benchmarks ran serially. Setup and allocation were
excluded. Public and research checksums agree for every matched run.
The small differences are within the observed process-to-process variation;
promotion shows no material performance regression.

The larger-outer s20 public API took 3.360091 ms in this campaign, giving the
new K20 profile 12.04% less latency. That comparison retains different
certificates: 58.75 bits for s22 versus 68.10 bits for the larger outer.
K18's packed-state optimization is modest: a separate reversed-order holdout
measured 465.674 us versus 468.490 us for the original cached kernel.
Streaming stores and unrolled outer mixing were slower at K18 and were not selected.

`spin_certified22_test` exercises the public interface. A repository checkout can
add `SPIN_BUILD_RESEARCH_ORACLE_TESTS=ON` to compare setup and encoding with the
unchanged research implementation. Installed consumers need only the library and
its public headers, not research sources.

Validation on 2026-10-03 passed all ten Linux Release CTests, including both
research-oracle targets. The final s22 tests additionally cover all four legal
output alignments above the cached/streaming boundary and allocation-free
workspace rebinding between s15 and s22. Installed Linux and Windows consumers passed.
A separate build with AVX-512 disabled passed both s22 and s15 API suites.
The Windows/MSVC build passed its eight runnable CTests; both research-oracle
tests skipped because the host lacks their required ISA. Linux AddressSanitizer
and UndefinedBehaviorSanitizer also passed the s22 API suite, with the new
profile and public dispatch instrumented. Unchanged linked kernels were not
rebuilt with sanitizers. An independent source review found no blocking issue.

The serial comparison driver in [k20_small_outer](experiments/k20_small_outer)
accepts `public`, `small-packed`, and `small-nt-packed-unroll` modes. They share
the same buffers, warmups, timing loop, and input initialization. Setup and
allocation are excluded. The [optimization report](experiments/k20_small_outer/OPTIMIZATION.md)
retains the original controls and measured K20 improvements.

The measured source archive is local at `tmp/k20-small-promotion-v1.tar.gz`,
SHA256 `0291394cffe99bff2189ad408d171bb99d643b4d16961dca00da81c0fb4a139b`.
The later test-only overlay is `tmp/k20-small-promotion-tests-v2.tar.gz`,
SHA256 `359d214f2b76e84124a42dff63454b4f2f31dd382f4bd0afb11954b02a8de293`.
Raw timing summaries remain ignored:

| File under `tmp/` | SHA256 |
|---|---|
| `s22-public-screen-v1.txt` | `e969d8ad0545bee90ede75292d5939f84e83af3b88b1b9d168a384d2c3e5ca3e` |
| `s22-public-holdout-v1.txt` | `332c75b54e1e6460c2cf9e5038c92e13bee218dfa9c4837aa2db2b65cf07be4c` |
| `k18-cached-holdout-v1.txt` | `dcdf68320186067c319a067b16fe37603ac302c08e11b4503107a8f86a91bbb8` |
