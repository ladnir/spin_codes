# Precomputed Packet Encoder at Smaller Message Sizes

This report follows the exact-map optimizations committed in `710f7be8`.
It studies K=2^16 and K=2^18 with the same t64/s16 packet construction.
The performance goal is a practical tuning plateau, not a claim of global
optimality. No production default or paper result changes in this work.

The structural follow-ups below supersede the original selected recipe.
The latest exact-map kernel measures **0.205008/0.826656 ms** with 64-byte-aligned
input/output at K=2^16/2^18. Prepared huge-page buffers reach
**0.203269/0.817720 ms**. The final section separates kernel and allocation gains.
Earlier results remain measured controls, not an optimization ceiling.

## Initial Size-Tuning Campaign

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
The last schedules were compared again with the combined outer before
selecting the shared implementation. None of those tested scheduling variants
offered a material, repeatable gain. This conclusion was limited to that
search space; the structural follow-up found further exact-map improvements.
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

## Structural Follow-up: Compose Maps Before Tuning Loops

The previous search retained the boundaries between state conversion, matrix
multiplication, and output layout. The follow-up instead composes those maps.
That combined candidate measures **0.21196/0.85148 ms** at K=2^16/2^18.
It preserves the original message coordinates and every sampled setup map.

| Matched candidate | K=2^16 | K=2^18 |
|---|---:|---:|
| Previous shared recipe | 0.222621 ms | 0.894390 ms |
| Composed inner and fused BCH output layout | 0.214170 ms | 0.859234 ms |
| Also absorb the BCH bit transpose into its products | **0.211956 ms** | **0.851484 ms** |

Each entry is the median of eight process medians: four setup seeds, both
candidate orders, five warmups, and 301 timed calls. The combined change
reduces latency by 4.79% and 4.80% against its matched control. The recipe
is `transpose:4` in the follow-up harness. This remains a research implementation;
public library defaults and paper numbers are unchanged.

An earlier matched round isolates the inner change: 0.223192/0.896189 ms
for the control and 0.219916/0.881326 ms with only the composed inner.
The selected final candidate spans 0.211615--0.212597 ms and
0.849987--0.855167 ms across its eight process medians.

### Three Exact Compositions

[InnerComposed.h](InnerComposed.h) combines the state update with the inverse
bit transpose. GFNI accepts both an input vector and a binary matrix. Swapping
their roles, with the required coefficient-byte reversal, computes the
transposed matrix product directly. The selected variant holds the state in
four wide vectors and streams emitted packets. It replaces four GFNI
instructions per steady-state epoch with four wide XORs. Setup still uses
32 bytes per epoch; it only changes coefficient orientation.

[outer_unpack_codegen.py](outer_unpack_codegen.py) combines three output
layout stages into one byte-permutation network. After eight unchanged bit
transposes, 24 two-source byte permutations replace 48 byte/lane shuffles per
eight output coordinates. Both transformations preserve the exact binary map.
Their individual timings and combined timings are measured separately.

[outer_transpose_codegen.py](outer_transpose_codegen.py) then moves the BCH
bit transpose through its matrix products. GL32 preparation writes the
systematic half in transposed form; the dense half keeps its existing form.
Swapped GFNI operands produce BCH contributions directly in output form.
The parity coordinate becomes one selected byte per matrix lane, handled
with a byte shuffle and XOR-AND. This removes 128 final GFNI instructions
and replaces 128 parity GFNI instructions per four-row tile. Setup coefficient
bytes and all output coordinates remain unchanged.

### Alternatives That Did Not Win

The following results are screening evidence, not selected implementations.
Screens use two setup seeds, both orders, and 61 calls per process.

| Structural alternative | K=2^16 | K=2^18 | Outcome |
|---|---:|---:|---|
| Smaller BCH plane batches | 0.224--0.263 ms | 0.903--1.054 ms | Extra passes outweigh lower register pressure |
| Pack data during scattered routing | 0.359 ms | 1.623 ms | Partial stores lose badly |
| Sequential inner stores, random outer reads | 0.235 ms | 0.939 ms | Slower without prefetch |
| Same gather route, best screened prefetch | 0.221 ms | 0.896 ms | Recovers the control, not the fused winner |
| Polynomial BCH, best screened leaf size | about 0.240 ms | about 0.963 ms | Slower despite fewer GFNI operations |

The routing experiment changes where the existing permutation is materialized,
not its distribution. Full cached packet stores remain the selected route.
Streaming stores followed by random reads were substantially slower still.

The polynomial experiment changes BCH message coordinates. An independent
row-space calculation finds an invertible binary matrix U with
`B_raw = U B_old`. Hence every fixed downstream mixer, route, and inner
produces the same forward code image. The transpose output is transformed
by the corresponding block-diagonal U, so this is not a byte-compatible
replacement. This candidate is retained only as an algorithmic experiment.
Its reduced GFNI count does not compensate for the additional XOR network,
intermediates, and register pressure in this implementation.

### Validation and Replay

[creative_probe.cpp](creative_probe.cpp) checks all nine initial candidates
at K=2^14, 2^16, and 2^18 with two setup seeds. It compares full outputs with
the retained encoder and an independent scalar recurrence, route, GL32 map,
and literal BCH matrix. The tests also cover the untouched input suffix,
scratch guards, padding, and every coordinate for one through four epochs.
Native helper tests cover complete state, feedback, and routed-packing bit
bases. The fused output network has a symbolic 4,096-bit basis check and
native full-encoder checks. The operand-swap generator additionally checks
all 4,096 bilinear input/matrix basis pairs and complete BCH output forms.
The final mixed-layout outer also passes [creative_basis.cpp](creative_basis.cpp):
all 131,072 physical input bits for each of its Original, Aligned, and Compact
exports, under two setup seeds. Those 786,432 native calls use scalar GL32
and literal BCH rows as the oracle. They check unaligned data/output,
guards, and input/coefficient immutability.

The raw polynomial candidate separately passes all 131,072 physical input-bit
bases for each of two setup seeds. Its full-encoder checks use the different
raw message basis explicitly. These checks are not distance certificates.
All 521 original certificate source hashes remain unchanged.

[creative_run.sh](creative_run.sh) builds and runs the main candidates.
[creative_more_run.sh](creative_more_run.sh) adds gather-prefetch variants,
the fused output network, and the mixed-layout BCH products. Both reuse the
clean control build described above.
[raw_outer_run.sh](raw_outer_run.sh) builds and checks the separately labeled
polynomial candidate. Every runner holds the shared benchmark locks.
After building the final outer, `transpose-basis 1 17` replays its exhaustive
outer checks. These are correctness checks, not performance measurements.

The campaign is retained at `/tmp/spin-inner-packet-8ze5YN` and locally under
`tmp/inner-packet/creative-final/`. Final matched confirmations are
`compare-xIgBuU` and `compare-FcJhe8`; the earlier inner/output-layout round is
`compare-BbCBFf` and `compare-MTHB7s`. Raw measurements remain ignored.

| Final artifact | SHA256 |
|---|---|
| Composed inner | `3fd6f90b28302aebb2f438aac68288c84eba96fbc624489ae765d280295b6cce` |
| Fused output generator | `64520ba7b61c9ae31dc0228fc4a1db968e5d0a443bb69dd2b900a32a6cc0d326` |
| Mixed-layout BCH generator | `3ecd30a8997969b3e466c3dae6250907416bfd4495532a81048f03d71c1a5776` |
| Final measured executable | `6765e18f34325885926d265fd71792bcad2da9adf12e62f765a0ebe8c6988c1d` |

This campaign left a concrete question: can a polynomial-native representation
make BCH cheaper after accounting for conversion costs? The next campaign
tests that question against this exact-map winner.

## Final Structural Search: Wide Emission and Buffer Geometry

This campaign pursued concrete alternatives until the remaining proposals
lacked a plausible implementation-level saving. Three independent reviews
covered polynomial BCH arithmetic, inner-map composition, and memory layout.
All benchmarks ran serially on the same Ryzen 7950X core, with the preceding
exact-map winner as the control. This is a practical local plateau, not a
claim of global optimality.

### Confirmed Results

| Configuration | K=2^16 | K=2^18 |
|---|---:|---:|
| Previous exact-map winner, original buffers | 0.212627 ms | 0.856432 ms |
| Previous kernel, 64-byte-aligned input/output | 0.208820 ms | 0.845366 ms |
| Wide emission, original buffers | 0.208289 ms | 0.840758 ms |
| Wide emission, 64-byte-aligned input/output | **0.205008 ms** | **0.826656 ms** |
| Wide emission, prepared 2 MiB-aligned huge-page buffers | **0.203269 ms** | **0.817720 ms** |

Each result is the median of eight process medians: seeds 1, 17, 43, and 91,
both candidate orders, five warmups, and 301 timed calls. The aligned kernel
reduces latency by 3.58%/3.48% relative to its matched original-buffer control.
Huge-page preparation raises those reductions to 4.40%/4.52%.
Changing K by a factor of four changes the aligned latency by 4.03 times.
Both sizes use the same algorithm and parameters.

The aligned result spans 0.204351--0.205684 ms and 0.824257--0.843833 ms across
the eight process medians. The huge-page result spans 0.203069--0.203639 ms
and 0.815961--0.820390 ms. Setup, allocation, page preparation, and correctness
checks are excluded from these precomputed-encoding timings.

[further_mapWide.h](further_mapWide.h) is the selected new kernel component.
Use `transposeInnerFurtherWide<false,false>` with the existing composed setup,
cached full-packet routing stores, compact GL32 coefficients, and the outer
from [outer_transpose_codegen.py](outer_transpose_codegen.py).
The inner keeps four state vectors wide through coefficient preparation.
It reuses lane broadcasts across masked packet assembly and the quadratic
coefficient circuit. This removes the extract-to-128-bit/reinsert-to-512-bit
boundary: 39 source-level operations replace 61 in that prefix.
The feedback and state-update recipes remain unchanged.

The kernel still accepts ordinary 16-byte-aligned elements. Alignment is a
performance recommendation, not a new correctness requirement. A separate
offset screen found that offsets 0 and 32 modulo 64 performed similarly;
offsets 16 and 48 were slower on this processor.

[further_memory_codegen.py](further_memory_codegen.py) isolates buffer geometry.
Its two allocation modes reserve equal owner sizes but align their interior
views differently. Both advise only pages wholly owned by their allocation.
The huge-page result requires successful Linux page collapse for input and
scratch; all confirmation runs reported success. At K=2^16 the comparable
64-byte-aligned allocation could not collapse a full huge page. At K=2^18,
both allocation modes could collapse pages, and the residual difference was
small: 0.819458 ms versus 0.817720 ms. This is not a portable kernel guarantee.

Separate phase measurements support the decomposition. At K=2^18, the
original-buffer control spends about 0.283 ms in inner/routing and 0.578 ms
in GL32/BCH. Wide emission changes these to 0.263/0.577 ms. Alignment and
huge-page preparation then reduce both phases, reaching about 0.254/0.568 ms.
These are phase means from separate runs, not sums of the table's medians.

### Rejected Alternatives and Stopping Decision

Screening uses two seeds, both candidate orders, and 61 calls per process.
Numbers from different screens are not treated as paired comparisons.

| Alternative | K=2^16 | K=2^18 | Decision |
|---|---:|---:|---|
| Fully unroll the current BCH input-pair loop | 0.216876 ms | 0.872103 ms | Slower than its 0.211290/0.857336 ms control |
| Reverse disjoint BCH output-tile order | 0.213309 ms | 0.858007 ms | No shared win |
| Exact Winograd BCH, four-plane products | 0.221709 ms | 0.892201 ms | Slower than 0.211375/0.854610 ms control |
| Exact Winograd BCH, two-plane products | 0.225296 ms | 0.908977 ms | Lower register pressure does not recover the cost |
| Raw-basis carryless-multiply BCH | 0.249165 ms | 1.001535 ms | Slower than its 0.212593/0.861333 ms exact control |
| Wide feedback in addition to wide emission | 0.211875 ms | 0.854453 ms | Slower than 0.208284/0.839971 ms wide-emission control |

[further_outer_winograd_codegen.py](further_outer_winograd_codegen.py) reduces
dense BCH GFNI operations from 2,048 to 1,792 per tile. Extra additions,
coefficient broadcasts, and temporary traffic outweigh that reduction in
both tested schedules. Both compiled helpers use stack temporaries.
This result does not justify deeper factorization.

[outer_clmul_codegen.py](outer_clmul_codegen.py) computes the Q/P polynomial
middle products with native carryless multiplication. It preserves the code
image but changes message coordinates, as established above. Both modes in
its comparison use the same current composed inner. The raw outer needs
768 carryless multiplications, 400 arithmetic GFNI operations, and 128 output
GFNI operations per tile, besides unchanged GL32 work. It also needs 1,152
conversion shuffles and 384 final output byte permutations. The complete
encoder loses despite the smaller dense-arithmetic count. Keeping GL32 in
the polynomial layout is not a free fusion: independently sampled matrices
act on different coordinate groups within each proposed polynomial word.

[further_feedbackWide.h](further_feedbackWide.h) batches low-lane reductions
and keeps seven syndrome coordinates wide. Its boundary drops from 98 to 82
source-level operations, but all three tested allocation configurations lose
to their wide-emission counterparts. A separate full-process counter check
records fewer retired instructions but slightly more cycles. Thus fewer
source operations alone do not predict a faster instruction schedule.
Pruned and incremental emission schedules, separately and together, also
lose to the simpler grouped evaluator after the new wide prefix.

No concrete high-confidence exact-map candidate remains from these reviews.
Retain wide emission and the existing feedback, not the unsuccessful rewrites.
Further basis resynthesis or instruction rescheduling may still help, but
there is no demonstrated large saving to pursue there. A larger advance
likely requires a construction or routing-layout change that reduces
conversion or scattered-memory work; such a change needs its own proof.
The next engineering task is to promote the shared winner into the library,
with optional owned-buffer preparation outside encoding, then test consumer
integration. Public defaults and paper numbers remain unchanged here.

### Validation and Reproduction

The wide-emission generator verifies all 44 coefficient lane forms and the
complete 64-input/16-state packet and feedback-moment maps. Its native helper
test exhausts all 2,048 state-bit bases. Selected variants pass full scalar
encoder comparisons at K=2^16/2^18 under two seeds and four input patterns,
plus one-through-four-epoch coordinate bases, suffix checks, and guards.
The ordinary-buffer wide kernel also passes the K=2^14 full check.

The rejected feedback candidate passes all 8,192 physical input-bit bases
against literal feedback masks and the full K=2^14 scalar encoder check.
Its symbolic verification executes the actual shuffle immediates and masks.
Winograd passes all 8,192 symbolic dense-input bit forms and native full
encoder checks. The carryless-multiply candidate passes 131,072 physical
outer input-bit bases under each of two seeds, including unaligned buffers
and guards. Its full-encoder checks use its raw message basis explicitly.
These implementation tests do not add numerical distance certificates.
The original certificate receipt and all 521 source hashes remain unchanged.

The follow-up runners consume the same retained build directories described
above: `ROOT`, `PREVIOUS_CREATIVE_BUILD`, `CONTROL_BUILD`, and `SOURCE_REPO`.
Place each runner's new sources in ROOT; imports resolve through ROOT,
the previous creative overlay, and the source workstream. Build stages are:

- [further_control_run.sh](further_control_run.sh): alignment, BCH unrolling,
  and output-order controls.
- [further_candidates_run.sh](further_candidates_run.sh): Winograd, wide
  emission, and separately labeled raw CLMUL with a common composed inner.
- [further_fine_run.sh](further_fine_run.sh): emission schedules, alignment,
  and equally sized owned-buffer geometry variants.
- [further_feedback_run.sh](further_feedback_run.sh): rejected wide-feedback
  variants, retained for replay.

Use `inner_size_compare.sh ROOT 16 base:4 wide:4 widealigned:4 memoryhuge:4`
and repeat with exponent 18. Here ROOT/inner-size is the preceding measured
exact-map winner. `SPIN_SIZE_VARIANT=widealigned inner_size_run.sh ROOT check
16 4` selects its full untimed checks; invoke the shell scripts with bash.
`inner-size-clmul basis 1` replays the raw outer's exhaustive native basis.
Build scripts and comparison runners hold all three shared locks.

The campaign is `/tmp/spin-inner-packet-nU76i3`, archived locally in ignored
`tmp/inner-packet/further-final/`. Confirmations are `compare-GNF5MC` and
`compare-gYybg8`; phase runs are `compare-h0cfzD` and `compare-uaf0BR`.
The report retains results and source-level ideas; raw measurements and
compiled experiments remain outside version control.

| Selected artifact | SHA256 |
|---|---|
| Wide-emission header | `f4daab3c89ff7a5fb24f8542f9670d54881a8dc900b2d9b7be18d3d85a1b13f7` |
| Ordinary-buffer executable | `b169a0b12471421fb1406a3d2d763af4e721fd87697fd305e9d34048f7da8096` |
| 64-byte-aligned executable | `357b86fb0637382fba0fe6b9610b7c6f9e368921aef0eac4bf5d6c48067564bc` |
| Huge-page executable | `3bc2180c0901bc07cd4f6f1568756364de5ce9019198801ea6f9e81b5640ae39` |
