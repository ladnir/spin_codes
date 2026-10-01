# Precomputed Packet Encoder at Smaller Message Sizes

This report follows the exact-map optimizations committed in `710f7be8`.
It studies K=2^16 and K=2^18 with the same t64/s16 packet construction.
The performance goal is a practical tuning plateau, not a claim of global
optimality. No production default or paper result changes in this work.

## Result and Selected Recipe

The selected precomputed encoder measures **0.222 ms at K=2^16** and
**0.899 ms at K=2^18**. These sizes use one shared recipe, not separate
handwritten encoders. Multiplying K by four increases latency by 4.04 times.

| Matched implementation | K=2^16 | K=2^18 |
|---|---:|---:|
| Retained packet recipe, streaming stores, original outer | 0.302314 ms | 1.101825 ms |
| Prior VBMI conversion win, still streaming stores and original outer | 0.293743 ms | 1.065927 ms |
| VBMI packet stream, cached stores, original outer | 0.231286 ms | 0.934573 ms |
| Same inner, byte/GFNI outer | 0.227354 ms | 0.916599 ms |
| Same inner, byte/GFNI outer plus factored GL32 | **0.222225 ms** | **0.898586 ms** |

The final row reduces latency by 26.5% and 18.4% relative to the retained
packet recipe. Relative to the already improved VBMI recipe, the further
reductions are 24.3% and 15.7%. These controls implement the same construction;
they are not comparisons with the paper's older parameter sets.

Use `InnerRecipe::StreamVbmi`, cached routing stores, compact coefficients,
and the outer generated with `--tile-mode 5 --mode 2`. In the probe this is
`tile5layout2:6`. No prefetch, pruned evaluation, or early update is required.
The eight process medians span 0.221553--0.223016 ms at K=2^16 and
0.892795--0.901330 ms at K=2^18.

Early updates yield 0.222084/0.899312 ms. Pruned evaluation yields
0.222680/0.895439 ms. Combining both yields 0.222891/0.895961 ms.
These differences are below 0.4%, do not win consistently at both sizes,
and overlap observed variation. We retain the simpler shared recipe.

A separate matched confirmation repeats the selected recipe at
0.222170/0.898762 ms. It also carries both prior VBMI recipes without the
store/layout changes: 0.293052/1.065017 ms for materialized packets and
0.293037/1.063389 ms for the packet stream. The new result is therefore
not an artifact of choosing the slower of the two prior winners.
The checked non-VBMI fallback measures approximately 0.245/0.987 ms
with cached stores and the original outer; it is not the optimization target.

Separate phase measurements explain where the larger gains occur:

| K and recipe | Inner plus routing | Outer GL32/BCH |
|---|---:|---:|
| 2^16, retained packet | 0.08955 ms | 0.21335 ms |
| 2^16, VBMI stream and cached stores | 0.07279 ms | 0.16039 ms |
| 2^16, selected combined outer | 0.07279 ms | 0.15173 ms |
| 2^18, retained packet | 0.35678 ms | 0.74822 ms |
| 2^18, VBMI stream and cached stores | 0.29531 ms | 0.64291 ms |
| 2^18, selected combined outer | 0.29438 ms | 0.60750 ms |

These are averages of four process phase means, using two seeds, both
orders, five warmups, and 101 calls. They include warmups and timer overhead;
they need not sum to the separately measured whole-call medians. In particular,
cached routing also benefits the following outer phase, not just the store loop.

## Workload and Measurement

The rate-one-half forward code has K input bits and N=2K output bits.
Its transpose consumes N 128-bit XOR elements and writes K elements
in place. We measure complete precomputed encodes, including routing and
the outer GL32/BCH stage. Setup, allocation, validation, and input filling
are outside the timed region. The input allocation is reused between calls.
These are steady-state timings, not cold-start latency or end-to-end OT.

All native measurements run serially on Peach's Ryzen 7950X, pinned to
CPU15 under the three shared benchmark locks. The build uses GCC 15.2.0,
`-O3 -DNDEBUG -std=gnu++20 -mtune=znver4 -fno-strict-aliasing`,
AVX2, AVX512F/VL/BW, GFNI, and AVX512VBMI. The fallback build excludes VBMI
in every translation unit. Hardware selection remains the caller's job.

Screening uses two setup seeds. Final confirmation uses seeds 1, 17, 43,
and 91. For each seed, the runner visits candidates in both forward and
reverse order. Each process performs five warmups and 301 timed calls.
The reported aggregate is the median of the eight process medians.
The main comparison binaries share one driver object; only the linked outer
implementation changes between outer variants. The no-VBMI fallback is
compiled separately and is not part of this main table.

## Changes and Their Scope

The dominant small-size change is cached routing stores. Each emitted
packet fills one aligned 64-byte destination. The outer stage then reads
the routed scratch. This immediate reuse plausibly favors cached stores
at these sizes; the timings do not isolate cache misses or write allocation.
Input plus scratch occupy approximately 4, 16, and 64 MiB at K=2^16,
2^18, and 2^20. The larger case retains its separately measured store policy.

The inner uses VBMI state conversion and streams four packets at a time.
The outer combines byte/GFNI conversion with factored GL32 mixing.
The two operations act on different coordinates: packing permutes payload
bits, while GL32 mixes row lanes independently for every payload bit.
Their composition therefore preserves the map. The generator verifies the
packing permutation, its inverse, and the GL32 identities on complete bases.

The reusable [inner kernel](PacketInnerKernel.h) takes original-coordinate
reverse matrices M^T. Preparation converts them once to P M^T P^-1.
The expansion is A P^-1 and feedback is P A^T; feedback is not the
transpose of the transformed expansion. The kernel has no allocation,
virtual dispatch, or runtime recipe selection inside its hot loop.
The emitter controls cached versus streaming stores and supplies any
required streaming-store fence before the outer reads scratch.
Feedback can consume emitted packets because this fixed map satisfies A^T A=0.

## Explored Tuning Space

The search includes:

- All 27 prior inner variants, carried to both smaller sizes.
- Four inner recipes crossed with cached/streaming stores and compact/expanded coefficients.
- Cached-store prefetch leads of 8, 16, and 32 packets for those recipes.
- Original, single-output, and byte/GFNI BCH schedules, plus factored GL32 mixing.
- Early state updates, pruned packet evaluation, and both together.
- The combined byte/GFNI outer layout and factored GL32 mixer.

Prefetching did not improve the first cached-store screen. Expanded
coefficients did not offer a consistent advantage across both sizes.
The last schedules are compared again with the combined outer before
selecting the shared implementation. Those comparisons support stopping
this exact-map tuning campaign: no remaining tested scheduling variant
offers a material, repeatable gain over the selected shared recipe.
The experimental alternatives remain
available for reproducibility; they are not separate public encoder APIs.

## Correctness and Certificates

[Boundary tests](PacketInnerBoundaryTests.h) cover lengths 0, 64, 128,
and 192. The independent oracle implements the original-coordinate
recurrence with scalar 64-bit XORs. Fixtures use distinct, nonsymmetric,
full-rank matrices. Tests include every input coordinate and every payload
bit at physical-step boundaries. They check output order, input immutability,
guards, empty calls, malformed sizes, matrix widths, and rank validation.

Complete-encoder checks compare with both the retained implementation
and a scalar inner followed by the explicit route, scalar GL32, and BCH.
They also check the untouched input suffix and scratch canaries. The
VBMI conversions have separate complete 2,048-bit pack and unpack tests.
All checks remain active in release builds.

The final build passes 216 complete configurations: 46 base modes, five
combined-outer modes, and three no-VBMI modes, each at both sizes with two
setup seeds. Every check-mode process also runs the boundary/API tests.
Both clean-source build scripts completed successfully. Generation and
compilation leave the source archive unchanged.

The selected map and setup distribution are unchanged. The existing
greater-than-10% distance certificate with more than 55.08 bits of
setup-failure margin is for K=2^20. These timings do not transfer that
numerical certificate to K=2^16 or K=2^18. All 521 files in the original
certificate manifest still match their recorded hashes. The whole receipt
remains `698e849a142641f6af3885f72a4ae5b2f285f0100a17051640c6a97bc0409c3f`.

## Reproduction

[inner_size_prepare.sh](inner_size_prepare.sh) builds from a clean checkout
or source archive, without retained remote objects. It extracts authenticated
reference headers, regenerates the outer implementations, and records
source hashes, commands, dependency files, symbols, and disassembly.
[inner_size_extra.sh](inner_size_extra.sh) adds the factored layouts and
the no-VBMI fallback. Both scripts require the source repository and an
existing build directory. The preparation script requires that directory
to be empty and never deletes its contents.

```sh
repo=/path/to/spin_codes
build=$(mktemp -d /tmp/spin-size-XXXXXX)
workstream="$repo/research/workstreams/permutation_locality"
bash "$workstream/inner_size_prepare.sh" "$repo" "$build"
bash "$workstream/inner_size_extra.sh" "$repo" "$build"
for k in 16 18; do
  SPIN_SIZE_VARIANT=tile5layout2 bash "$workstream/inner_size_run.sh" "$build" check "$k" 6
  bash "$workstream/inner_size_compare.sh" "$build" "$k" base:0 base:1 base:6 tile5:6 tile5layout2:6
done
```

The `inner_size_run.sh` build/link phases retain the historical object-path
interface for old campaigns. Use `inner_size_prepare.sh` for clean builds;
the check and timing phases work with the binaries that it produces.
The compare runner selects `base:MODE` or `BINARY_SUFFIX:MODE` and holds
all three locks throughout a matched run. Never launch benchmarks concurrently.
Raw measurements, binaries, and disassembly remain in ignored scratch storage.

The final clean build and runs are under `/tmp/spin-inner-packet-YMwLe7/final`,
archived locally in `tmp/inner-packet/size-final-records/`. Main confirmations
are `compare-pelb3a` and `compare-sW05SV`; phase runs are `compare-woXuTf`
and `compare-LfAxcn`. The earlier broad sweep is retained separately in
`tmp/inner-packet/size-second/`.
The independent control confirmation is `compare-Q8nsNQ` / `compare-Os4NJb`.

| Measured artifact | SHA256 |
|---|---|
| Reusable inner header | `5c3d32f42507fb76fdb537834954f67aea0dc3c2176a10df01cdbb7f72dd6b4f` |
| Probe source | `2f0500dbcbed6099e0aa29c7c4f4b20b3b54f1bde11e2fbf606a5798dd05fb3d` |
| Outer layout generator | `729d6133911e3e94bc7774c4c56b69f6f3a84ba58087d0fc6151ac8207ac4fbf` |
| Selected executable | `1458007d97c5d212a2c0d3f1066bbd298215af9114162f379b79bff3f8d245d0` |

Next, promote the measured shared recipe into the consumer-facing encoder
with size/ISA selection outside the hot loop. Keep the K=2^20 route policy
until it is remeasured. Numerical distance claims at the smaller lengths
still require their own certificate replay.
The current header requires consistent VBMI macros across translation units;
mixed-ISA library integration needs isolated target wrappers before deployment.
