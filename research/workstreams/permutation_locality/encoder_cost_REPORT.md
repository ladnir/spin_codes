# Where the Packed Encoder Spends Time

2026-09-30. **The precomputed t64/s16 encoder now measures about 4.77 ms.**
Packet expansion and feedback fusion first reduced latency from 5.320 ms
to 4.885 ms. A follow-up campaign adds another approximately 2% through
VBMI state conversion. Four-packet streaming gives a similar result when
combined with that conversion. Packed GL32/BCH processing remains about
3.1 ms. Both campaigns and unsuccessful alternatives are preserved below.
No code distribution or distance claim changes.

The follow-up [small-size campaign](inner_size_REPORT.md) reaches 0.222 ms
at K=2^16 and 0.899 ms at K=2^18. It adds cached routing and a combined
outer layout. The measurements below remain the original K=2^20 campaign.

All measurements use K=2^20, rate one-half, 128-bit XOR elements, and
Peach's Ryzen 7950X on CPU15. The forward code maps K bits to 2K bits;
the timed transpose consumes 2K elements and writes K elements in place.
Setup and allocation are excluded from latency measurements. Benchmarks
ran serially under the three shared locks. No certified source or
production default was changed.

## Full-Encoder Phase Timers

The original executables already have coarse phase timers. Two setup
seeds, each with three warmups and 31 measured calls, give:

| Encoder | Inner + routing mean | Packed GL32/BCH mean |
|---|---:|---:|
| t=64, s=16, uniform GL16 | 2.2064--2.2074 ms | 3.1407--3.1410 ms |
| t=128, s=19, four transvections | 2.1663--2.1881 ms | 3.1382--3.1654 ms |

These phase means include warmups; they do not exactly decompose the
separately reported whole-call medians. Both configurations have the same
outer phase. Shrinking the state has not reduced the inner/routing total.

## Separating the Outer Phase

Cycle sampling of the unchanged t64 executable used 5,001 measured calls,
seed 17, and 1,999 samples/second. No samples were lost. The self-sample
shares were 39.54% in the fused inner/routing function, 25.40% in the
GL32 preparation function, and 31.10% across the eight BCH output kernels.
Setup and checks account for most remaining samples; they are excluded
when normalizing the encoder shares below.

| Work | Approximate encoder share | Approximate time |
|---|---:|---:|
| Inner arithmetic, state conversion, and routing | 41% | 2.21 ms |
| Input packing and GL32 mixing | 26% | 1.41 ms |
| BCH output computation, unpacking, and output stores | 32% | 1.73 ms |

The latter two times split the measured outer mean using the sampled
ratio 25.40:31.10. They are estimates, not separately timed stages.
Instruction-level samples can skid and do not assign every stall to its
cause. In particular, this experiment does not isolate state-update
latency from the rest of the fused inner/routing loop.

The source explains why the small-state GFNI update is not an obvious
dominant cost. At this size, nominal GFNI512 counts are:

| Component | GFNI512 operations per complete call |
|---|---:|
| GL32 local mixing | 2,097,152 |
| BCH dense/parity arithmetic | 4,456,448 |
| 16-bit state matrix updates | 262,128 |

The outer has about 25 times as many GFNI operations as the state matrix
updates. Counts are not cycle estimates. The outer also has packing,
cross-lane shuffles, logical operations, loads, and stores. The t64 inner
doubles the frequency of state unpacking, feedback packing, and output
table construction compared with t128.

## Cache-Footprint Experiment

[encoder_cost_probe.cpp](encoder_cost_probe.cpp) links the exact
`PackedCoeff1.o` used by both measured encoders. Every diagnostic call
invokes this kernel 2,048 times. Only the number of distinct input,
coefficient, and output tiles changes; the arithmetic kernel is unchanged.
Four tiles are checked against scalar GL32 mixing and the retained BCH
reference before timing.

The experiment uses seeds 1 and 17, both footprint orders, three warmups,
and 101 measured calls per case. The following values are medians of the
four case medians. Logical footprint includes padded input, coefficients,
and output, but excludes the 16-KiB local packed buffer and code.

| Distinct tiles | Logical footprint | Time for 2,048 kernel calls |
|---|---:|---:|
| 1 | 28.1 KiB | 2.5036 ms |
| 16 | 449 KiB | 2.5067 ms |
| 64 | 1.75 MiB | 2.5347 ms |
| 256 | 7.02 MiB | 2.5394 ms |
| 1,024 | 28.1 MiB | 2.6478 ms |
| 2,048 | 56.1 MiB | 3.3356 ms |

About 75% of isolated bulk time remains with a small reused working set.
That remainder includes arithmetic, layout operations, and cache traffic;
it is not a pure arithmetic lower bound. The roughly 0.83-ms footprint
penalty is real in this diagnostic, but is not directly subtractable from
the full encoder. Its independently allocated buffers and preceding
operations differ: the actual outer phase measures about 3.14 ms.

A longer paired counter run used 1,001 measured calls per order:

| Whole diagnostic process | One tile | 2,048 tiles |
|---|---:|---:|
| Retired instructions | 47.2932 billion | 47.2933 billion |
| Cycles | 24.2694 billion | 31.7654 billion |
| Generic `cache-misses:u` events | 160,191 | 214,018,642 |

These counters include identical setup and reference checks, unlike the
latency table. Almost unchanged instruction counts, more cycles, and far
more cache-miss events support a footprint-related memory penalty. The
generic event is not converted into DRAM bytes or attributed to one
specific cache level.

## Next Optimization Targets

1. Reduce work and local-buffer traffic in the shared packed outer kernel.
   Both GL32 preparation and BCH output computation are substantial.
   Consider scheduling or fusion while preserving the exact maps.
2. Check whether the current outer-buffer layout can retain more locality.
   The full-footprint penalty is worth pursuing, but eliminating it alone
   cannot remove most outer cost.
3. Profile the fused inner/routing loop separately before another inner
   redesign. Its 2.21 ms includes scatter stores, expansion, feedback,
   packing, and state updates; it is not a measurement of GFNI alone.

For scale, a 20% reduction in the 3.14-ms outer phase would save about
0.63 ms, or 12% end to end. This is an opportunity estimate, not an
achieved improvement. No optimization was implemented in this diagnostic.

## Exact-Map Optimization Follow-Up

A subsequent campaign tested nine implementation variants of the same outer
map. None demonstrated a useful full-encoder improvement. Keep the retained
kernel: the current t64/s16 encoder remains approximately **5.3 ms** at K=2^20.
This result does not establish an implementation lower bound.

The candidates change instructions or scheduling, not the GL32 distribution,
BCH matrix, routing, or inner. The three new generators retain the reference
entry points and replace only the compact-coefficient kernel:

- [outer_prefetch_codegen.py](outer_prefetch_codegen.py): output write/read
  prefetching, optionally with input and coefficient lookahead.
- [outer_schedule_codegen.py](outer_schedule_codegen.py): force-inline the
  existing BCH tiles, or use four-output/four-plane tiles with packed writeback.
- [outer_layout_codegen.py](outer_layout_codegen.py): ternary reduction or
  factored GL32 mixing. Factoring reduces each plane from three lane rotations
  to two; it still uses four GFNI operations. The inverse rotation of odd
  coefficients is absorbed into their existing compact expansion.

Initial screening used two setup seeds and 31 measured calls per process.
Values below are medians of the two process medians. Controls belong to their
own screening campaign; small differences between campaigns are not speedups.

| Change | Full encoder | Matched retained control |
|---|---:|---:|
| Output write prefetch | 5.407 ms | 5.327 ms |
| Output write plus input/coefficient prefetch | 5.928 ms | 5.327 ms |
| Output read prefetch | 5.447 ms | 5.327 ms |
| Force-inline BCH helpers | 5.328 ms | 5.319 ms |
| Four-output/four-plane BCH, paired inputs | 5.582 ms | 5.319 ms |
| Four-output/four-plane BCH, single inputs | 5.807 ms | 5.319 ms |
| Ternary-only GL32 reduction | 5.319 ms | 5.302 ms |
| Factored GL32, one plane at a time | 5.291 ms | 5.302 ms |
| Factored GL32, four-plane scheduling | 5.318 ms | 5.302 ms |

Factoring received longer confirmation: two seeds, two execution orders,
and 101 measured calls per process. In the first campaign, retained and
factored medians were 5.3187 and 5.3311 ms. In a second campaign they were
5.3381 and 5.3244 ms. A freshly compiled, unchanged-kernel control in that
second campaign measured 5.3268 ms. Thus the candidate does not show a
consistent end-to-end advantage over unchanged code.

The isolated full-footprint outer test suggested a small reduction, from
3.3422 to 3.3055 ms. This is a weaker screen: one seed, two passes, and
31 calls per pass. It does not override the complete-encoder results.
The wider BCH schedules reduce source loads but add coefficient broadcasts
and a writeback pass. Their measured regressions illustrate why counting
loads alone is insufficient.

All executables link the same `Driver.o`; there is no runtime kernel dispatch.
The retained full executable is byte-identical to the preceding campaign.
Disassembly confirms identical inner/routing instructions and addresses in
the retained, freshly rebuilt, and factored binaries. The freshly rebuilt
baseline also has identical compact-outer instructions and placement.
Output write prefetch emitted actual `prefetchw` instructions.

The factored kernel passed all 131,072 physical input-basis vectors for
each of two coefficient seeds, against scalar GL32 mixing and reference BCH.
It also passed misaligned-output and canary checks, and complete-encoder
reference/adjoint/suffix checks at K=2^14 and K=2^20 with both seeds.
These C++ comparisons use explicit exceptions and remain enabled under
`-DNDEBUG`. Generator checks also verify the GL32 factoring symbolically
and on coefficient/input bases.

No candidate is promoted to the production or certified implementation.
All 521 source/input hashes in the t64 whole-certificate receipt remain
unchanged, as does the receipt itself. The 10% distance and greater-than-55.08-bit
certificate therefore retains its original implementation and evidence.

These results motivated a fusion experiment: produce the one-use systematic
BCH contributions directly in the output accumulators, while retaining
packed storage for the reused dense inputs. The next section records its
implementation and measurements.

The serial runner is [outer_opt_run.sh](outer_opt_run.sh). Like the original
diagnostic, it uses retained Peach sources and objects, not a standalone
library build. Its campaign is `/tmp/spin-outer-opt-hq0sio`; generated sources,
objects, executables, and every timing/check log are preserved locally in
ignored `tmp/outer-opt/records/`. No raw data is proposed for commit.

| Follow-up artifact | SHA256 |
|---|---|
| Shared driver object | `f8eb873339b0e69683d95a2753362ae84244de35868a3cabcf89e147bda47a6d` |
| Factored outer object | `db785a6c7989188a2d17c4bd706bf169ec8c56850a1141ce49bea6deeec52caf` |
| Factored full executable | `f2bf27fa1d7696d52c86b10afafc205cacb9d73fb91b6e009bd50c71ba2d100d` |

## Systematic-Fusion Follow-Up

The fusion is correct but slower. Its best tested version measures
**5.7414 ms versus 5.3248 ms** for the retained complete encoder, a 7.8%
regression. Keep the retained implementation. This conclusion applies to
the tested schedules and compiler, not every possible fused implementation.

[outer_fusion_codegen.py](outer_fusion_codegen.py) retains packed storage
for the 16 dense input groups and group 15, which supplies parity coordinate
127. It prepares systematic groups 0--14 directly in BCH accumulators.
Group 15 is mixed once and retained unmasked; only its systematic contribution
is masked by `0x7f`. Every GL32 group is still evaluated exactly once.
The matrices, randomness, input/output layouts, and encoder map are unchanged.

At source level, the two-output schedule removes 120 vector stores and
120 vector loads per four-row tile: 15 KiB of combined local-buffer traffic.
Its explicit packed buffer falls from 16 KiB to 8.5 KiB. The GFNI operation
count does not fall. The single-output variant also doubles dense-source
loads, so it tests register pressure rather than offering the same traffic saving.

All six default-compiler variants passed exhaustive correctness checks.
Serial screening used two seeds and 31 measured calls per process:

| Mode | Schedule | Full encoder | Matched control |
|---|---|---:|---:|
| 20 | Two output groups, direct systematic preparation | 5.978 ms | 5.308 ms |
| 21 | One output group | 5.954 ms | 5.308 ms |
| 22 | Shared parity pass, immediate cached stores | 6.044 ms | 5.329 ms |
| 23 | Mode 22 with one runtime-indexed direct-call helper | 5.944 ms | 5.329 ms |
| 24 | Mode 23 with explicit three-term ternary reduction | 5.907 ms | 5.330 ms |
| 25 | Mode 23 with two-rotation GL32 factoring | 5.852 ms | 5.330 ms |

The runtime-indexed helper has no indirect dispatch or allocation. It
reduces duplicated instruction bodies, while keeping the arithmetic loops
fixed-width. This variant was tested because code footprint grew alongside
register pressure in the fully specialized helpers.

Assembly inspection explains a substantial part of the initial regression.
Mode 20 materializes many GFNI terms before reducing them. Its helpers add
178 spill/reload pairs per tile. An initial preparation schedule also delayed
all eight stores until all eight mixed planes were ready, unlike the retained
immediate-store ordering. That added 102 further spill/reload pairs. Together,
these add 35 KiB of local stack traffic, exceeding the nominal 15-KiB saving.
The dense BCH accumulation loop itself remains spill-free in this version.

Modes 22--25 restore immediate stores for cached groups. Register pressure
still remains in their fused preparation. GCC reports 2,184 bytes of stack
for ordinary mode-22 helpers and 456 bytes for the shared mode-23 helper;
the retained BCH helpers use only eight bytes. These are compiler stack-size
reports, not a measurement of memory latency or a complete causal attribution.

Three outer-only compiler-control builds did not produce a win. Disabling
tree reassociation in mode 20 gave 5.918 ms. Disabling RTL instruction
scheduling gave 5.805 ms in mode 22 and 5.747 ms in mode 23. The latter
received confirmation with two seeds, both execution orders, and 101 measured
calls per process. Its median of four process medians is the 5.7414-ms result
above; the matched retained median is 5.3248 ms.

Phase timers localize that regression. The retained outer phase measured
3.118--3.129 ms, versus 3.540--3.544 ms for mode 23 with scheduling disabled.
Inner/routing means stayed approximately 2.2 ms. As before, phase means
include warmups and do not exactly decompose whole-call medians.

Every default variant passed all 131,072 physical basis vectors for each
of two coefficient seeds, plus offset-output and canary checks. The selected
compiler-control binary passed the same checks. Complete-encoder checks
cover K=2^20 for all timed variants; additional checks at K=2^14 and K=2^20
cover mode 25 and the selected compiler-control binary with both seeds.
No optimized candidate was installed. All 521 original proof/input hashes
and the whole-certificate receipt remain unchanged.

The campaign is `/tmp/spin-outer-fusion-BRzsCE`; all generated sources, logs,
objects, binaries, and inspected disassembly are archived in ignored
`tmp/outer-fusion/records/`. The runner now accepts optional compiler flags
after the build tag and records the exact compilation command. These flags
affect only the candidate outer object; every executable links the same driver.

| Fusion artifact | SHA256 |
|---|---|
| Final generator | `3f9af7a5f5e73abfe7dae2d9d5477dbf9994c51a9618448d256f7cb699e5cf35` |
| Selected compiler-control outer object | `5b1ff563f8d77d0649facf3c127223f854db39dcb319cef485c69f5f56e408b8` |
| Selected compiler-control full executable | `38126e8c850f1d12c230b7d832091a4a57dfb0245716e069f2606f01d87df7cd` |

The next recommended target is the unchanged 2.2-ms inner/routing phase:
separate emission-table construction, feedback/state conversion, and routing
costs before choosing another optimization. The present experiments do not
justify replacing the retained outer or changing the code distribution.

## Packet Expansion and Feedback Fusion

The winning candidate is mode 5 of
[inner_packet_probe.cpp](inner_packet_probe.cpp). It generates complete
four-element packets from the fixed quadratic expansion, changes the
internal state basis, and computes feedback from those emitted packets.
GFNI still applies the random state maps. The candidate remains isolated
research code; it does not replace a production default.

Final confirmation used two setup seeds, both execution orders, three
warmups, and 101 timed calls per process. The table reports the median
of the four process medians, not the fastest observed sample.

| Mode | Configuration | Full encoder |
|---|---|---:|
| 0 | Retained implementation | 5.319982 ms |
| 5 | Packet expansion, new basis, emitted-packet feedback | **4.885253 ms** |
| 9 | Mode 5 with the pruned quadratic evaluator | 4.896204 ms |
| 16 | Mode 5 with feedback finished/packed eight coordinates at a time | 4.927011 ms |
| 17 | Pruned evaluator and grouped feedback packing | 4.950400 ms |

Mode 5 gives 1.089x throughput at unchanged parameters. The difference
between ordinary and pruned evaluation is small; fewer source XORs did
not establish an additional speedup. Grouped packing reduced the compiler's
frame-size report from 256 to 192 bytes, but did not improve latency.

Separate phase measurements, with 31 calls and the same two seeds, give:

| Phase | Retained | Mode 5 |
|---|---:|---:|
| Inner and routing | 2.190--2.201 ms | 1.750--1.750 ms |
| Packed GL32/BCH | 3.126--3.126 ms | 3.119--3.136 ms |

Phase means include warmups and do not exactly sum to whole-call medians.
The approximately 0.45-ms improvement is localized to inner/routing.

### Exact Transformations

Let A be the retained 64-by-16 binary expansion. Let x be one 64-element
input step and r its entering reverse state. Each element is a 128-bit
XOR value. For the original sampled forward matrix M, the reverse step is

    y = x + A r,
    r_previous = M^T r + A^T x.

The fixed invertible 16-by-16 matrix P changes only the state coordinates:
z=P r. The candidate therefore uses

    y = x + (A P^-1) z,
    z_previous = (P M^T P^-1) z + (P A^T) x.

Setup computes the conjugated matrices from the original sampled matrices.
It does not resample them or introduce extra random choices. In general,
the transformed feedback is not the transpose of the transformed expansion;
both maps must be changed consistently.

The authenticated selected map satisfies A^T A=0. Consequently,

    (P A^T) y = (P A^T) x.

Feedback can consume emitted packets already in registers. It need not
reload the raw input or retain a separate raw-input array. This equality
holds for every input and state, not only on average over setup.

Each expansion output is quadratic in its six coordinate-index bits.
The two low bits select a lane within a packet. Evaluating the remaining
four bits in ZMM registers produces all sixteen packets directly.
[inner_packet_codegen.py](inner_packet_codegen.py) authenticates the
selected map and verifies both expansion circuits, both feedback maps,
the basis inverse, and A^T A=0 with exact binary linear forms.

The state basis reduces coefficient preparation from 29 to 12 source
XORs. Its feedback finisher uses 32 XORs versus the original 36. The
split finisher uses 31 XORs, but its additional scheduling change did not
help measured performance. These are circuit counts, not cycle predictions.

### Other Candidates and Assembly Findings

Seventeen candidates were implemented in addition to the retained control.
They vary the state basis, ordinary/pruned expansion, feedback source,
feedback layout, state representation, and grouped feedback packing.

Corrected expansion-only screens measured about 5.18--5.23 ms. Packet
feedback from raw input improved further, but emitted-packet feedback gave
the largest gain. An early expansion-only implementation called a 1-KiB
`memcpy` on every step. That call forced SIMD state spills. The final
implementation retains raw packets while emitting outputs instead; the
early slower timings are not evidence against packet expansion itself.

Modes 12--15 retain sixteen ordinary state words and replace GFNI updates
with four nibble tables. They eliminate state packing/unpacking, but build
1 KiB of tables and perform 64 indexed table reads per update. Screens
measured approximately 5.11--5.26 ms, slower than the selected candidate.
Those variants are retained as tested alternatives, not promoted.

The first inspected mode-5 assembly keeps emitted packets in ZMM registers
through their routing stores and immediately reuses them for feedback.
It has no steady-state ZMM spills. Some XMM intermediates still spill.
The final compiler frame reports are 1,664 bytes for the retained kernel,
256 for mode 5, 192 for grouped packing, and 1,472 for ordinary-word updates.
Frame size alone is not a performance metric.

Generic dense GFNI expansion/feedback and two-step matrix composition were
also considered analytically, not benchmarked. Dense maps add packing
costs. Two-step composition still needs the intermediate state to emit its
output, so it cannot simply remove one state transition. Neither offers
the same direct reuse as the quadratic packet circuit.

### Validation and Records

The final executable passes all 72 combinations of modes 0--17, K=2^14
and K=2^20, and setup seeds 1 and 17. Each complete-encoder check uses
four input patterns and compares against the independent scalar inner,
explicit routing, scalar GL32, and retained BCH implementation. Separate
comparisons cover the entire in-place buffer, unchanged suffix, and scratch
canaries. Candidate checks also cover:

- all 2,048 expansion input-bit bases;
- all 8,192 feedback input-bit bases, including grouped pack/unpack;
- every coordinate of a four-step inner sequence, exercising two actual
  state updates plus both boundary cases;
- all 65,536 row masks and 4,096 state/feedback input-bit bases for the
  ordinary-word update helper.

These checks remain active with `NDEBUG`. Exact map equivalence preserves
the [existing certificate](packed_mixing/s16_closure/CLOSURE.md): greater
than 10% distance with more than 55.08 bits of setup-failure margin at
K=2^20. All 521 original manifest files remain byte-identical. The whole
certificate receipt retains SHA256
`698e849a142641f6af3885f72a4ae5b2f285f0100a17051640c6a97bc0409c3f`.

[inner_packet_run.sh](inner_packet_run.sh) serializes build, check, screen,
confirmation, and phase-timing commands under the three shared locks.
The compiler is GCC 15.2.0 with the retained AVX2/AVX-512/GFNI flags.
The same `PackedCoeff1.o` outer object is linked for every variant.
Generate `InnerPacketMaps.h` and the authenticated `T64Reference.h` using
`inner_packet_codegen.py --output ... --reference-output ...` before building.
The runner relies on the retained Peach objects; it is not a standalone
library build.

The campaign is `/tmp/spin-inner-packet-DJlYTv`. Raw records are archived
locally in ignored `tmp/inner-packet/records/`. Final confirmation is
`confirm-d0u76v`, phase timing is `profile-I4W8Kr`, and the final check
summary contains 72 passing invocations. `first-disasm.txt` predates the
raw-copy correction; `final-disasm.txt` describes the final binary.

| Final artifact | SHA256 |
|---|---|
| Candidate probe | `1f796772fc380519fecec1dd45418fcbaeab701f56e9fcd71d3bdd9faaaeb733` |
| Generated packet maps | `842f3e726d65a9719bd7677dfdf81d2152b39778475873ac157559d8613bb84f` |
| Complete executable | `36997879c69fcece8e13029da1d8256aa5c748c149a55b34d14e642f358b1d05` |

At this stage, the recommended next experiment targeted conversion and
spill costs within the remaining 1.75-ms inner/routing phase. The following
campaign implements that experiment without editing the prior winner.

## State Conversion and Four-Packet Streaming

The follow-up tests 27 compile-time variants of the same encoder. All use
the retained outer object and the same sampled route and state matrices
for a given seed. None changes the selected fixed maps. The experiments
remain isolated in [inner_fusion_probe.cpp](inner_fusion_probe.cpp).

### Confirmed Timings

The final comparison uses setup seeds 1, 17, 43, and 91. Each seed has
two process runs in opposite mode orders, with three warmups and 101
measured calls per process. The table gives the median of eight process
medians; the range shows those medians, not individual-call percentiles.

| Mode | Exact implementation | Median | Process-median range | Latency reduction |
|---|---|---:|---:|---:|
| 0 | Prior packet-fusion winner | 4.8754 ms | 4.8496--4.9115 ms | reference |
| 6 | Four-packet streaming, existing conversion | 4.8419 ms | 4.8251--4.8637 ms | 0.69% |
| 17 | VBMI conversion, existing packet schedule | 4.7779 ms | 4.7387--4.8050 ms | 2.00% |
| 19 | VBMI conversion and grouped four-packet streaming | 4.7706 ms | 4.7330--4.7902 ms | 2.15% |

Both VBMI candidates beat the matched baseline in every confirmation run.
Their 0.0073-ms aggregate difference does not establish a reliable advantage
for streaming: their order reverses across seeds and repetitions. Prefer
mode 17 for initial integration because it changes only the conversion
helpers. Retain mode 19 as a scaling candidate. Relative to the earlier
5.320-ms campaign, 4.77 ms is approximately 10% lower latency; that comparison
spans campaigns, unlike the matched 2% result above.

Separate phase runs at seeds 1 and 17 give inner/routing means of
1.743--1.779 ms for mode 0, 1.636--1.659 ms for mode 17, and 1.655--1.665 ms
for mode 19. The unchanged outer phase ranges from 3.092 to 3.181 ms.
These measurements include warmups and instrumentation; they are not an
exact decomposition of the confirmation medians.

### Changes and Rejected Alternatives

[InnerWideState.h](InnerWideState.h) replaces narrow byte-transpose
networks with AVX-512 word permutations and byte shuffles. It preserves
the existing packed format. The combined conversion barely improves the
screen; replacing only one direction loses performance.

[InnerVbmiState.h](InnerVbmiState.h) folds each word-permutation/byte-shuffle
pair into one VBMI byte permutation. Each direction uses four byte
permutations and four GFNI operations. Boundary lane inserts and extracts
remain. This is the conversion used by modes 17 and 19. It requires a
VBMI-capable build and processor; it is not a replacement for the fallback
on unsupported machines.

[InnerPacketStream.h](InnerPacketStream.h) expands, routes, and consumes
four output packets at a time. Both schedules preserve descending packet
store order. The incremental schedule immediately accumulates eleven
high-coordinate moments. The grouped schedule retains partial moments and
reduces them through a smaller XOR tree. The high moments then feed the
unchanged transformed feedback map. Neither schedule needs sixteen live
completed output packets.

[InnerDirectFeedback.h](InnerDirectFeedback.h) instead composes the
feedback map into byte-matrix accumulation, followed by four GFNI
operations. It never constructs sixteen syndrome words. However, its
23 broadcasts and 46 byte shuffles outweigh that benefit in these tests.
The best combined screens remain around 5.01--5.05 ms. Do not promote it.

Early state-update schedules compute the matrix product immediately after
unpacking the old state. Expansion still reads the old coordinate words;
the new packed state receives feedback afterward. This is exactly the
same recurrence, but loses four fused ternary operations. The tested
early schedules do not improve the finalists. Pruning three expansion
XORs likewise gives no winning combination.

| Modes | Tested change |
|---|---|
| 1--3 | Wide word-permutation conversion: both directions, pack only, unpack only |
| 4--5 | Direct packed feedback, with either unpack implementation |
| 6--10 | Incremental/grouped streaming; pruned expansion; wide conversion combinations |
| 11--14 | Streaming plus direct feedback and wide-unpack combinations |
| 15--16 | Early update with or without streaming |
| 17--20 | VBMI conversion alone, incremental/grouped streaming, early update |
| 21--23 | Pruned early streaming; VBMI pack only; VBMI unpack only |
| 24--26 | Grouped VBMI streaming with pruning, early update, or both |

Assembly inspection explains why source operation counts were insufficient.
In the second build, the baseline has 336 bytes of XMM stack traffic per
steady epoch. Word-permutation conversion increases that traffic to 432
bytes; incremental streaming reduces it to 144 bytes. Mode 19 reduces it
to 32 bytes and replaces 112 XMM unpack instructions with eight VBMI byte
permutations. Its GFNI count remains sixteen. It also loads eight 64-byte
permutation controls per epoch. None of these loops calls a helper or
spills YMM/ZMM state. These are instruction counts, not a single-cause
explanation of the timings.

### Screening Different Fixed Maps

[inner_map_cost_screen.py](inner_map_cost_screen.py) separately screens
new constructions. It exactly enumerates the 65,536 expansion states,
checks packet ranks, and counts feedback-kernel supports through weight
four. MacWilliams coefficients independently verify those support counts.
The feedback-kernel minimum below is the minimum weight of a nonzero
64- or 128-bit input annihilated by the feedback map, not the SPIN distance.

| Fixed map | Minimum expansion weight | Feedback-kernel minimum | Weight-four kernel words |
|---|---:|---:|---:|
| Selected t64/s16 | 16 | 6 | 0 |
| Nested t64/s16 | 16 | 4 | 16 |
| Nested t128/s16 | 48 | 4 | 320 |

All three maps have rank sixteen, zero feedback-times-expansion product,
distinct nonzero columns, rank-four aligned packets, and rank-seven
two-packet unions. Those local rank properties alone therefore miss a
material difference in cancellation behavior.

The tested t64 monomial maps simplify the circuits but introduce 192--560
weight-four kernel words; the tested t128 monomial map has 3,168. Some
also reduce two-packet ranks. Omitting
the monomial involving the two low coordinate bits makes every aligned
packet rank three. These candidates need new analysis; the screen does
not prove that their complete SPIN codes have poor distance.

A bounded search tests 135 single-monomial mutations of the selected map;
69 pass the initial cancellation and rank filters. Six low-cost finalists
receive exact state-spectrum enumeration. One candidate toggles
`x0*x5` in generator row 12. It retains both minima and at least six active
output packets, but increases weight-sixteen expansion states from 20 to
36. Its small native-basis XOR saving is not comparable to the already
basis-optimized kernel. It is a proof-screen candidate, not a measured
performance improvement or a replacement certificate. Changing physical
step size to 128 also changes refresh frequency and requires new analysis.

### Validation and Reproduction

The exact-code candidates preserve the conjugated state representation
and the identity allowing feedback from emitted packets. Initialization
uses raw packets with zero state. The final physical epoch emits output
without an unused feedback calculation or state update.

Both wide conversion helpers have separate exhaustive pack and unpack
checks on 2,048-bit bases. The VBMI checks additionally compare with an
explicit scalar byte-layout oracle. Direct feedback is checked on all
8,192 input bits. The streaming generator checks every output and moment
as a symbolic linear form in all input and state coordinates. Complete
tests compare against the scalar inner, explicit route, scalar GL32 and
retained BCH implementation, including output suffix and scratch canaries.
Every candidate also exercises all 256 input coordinates across four
physical epochs. Checks remain active under `NDEBUG`. The final executable
passes all 108 combinations of modes 0--26, K=2^14 and K=2^20, and setup
seeds 1 and 17.

The selected construction retains its greater-than-10% distance certificate
and more than 55.08 bits of setup-failure margin at K=2^20. All 521 original
manifest files were rehashed and remain unchanged, as does the original
whole-certificate receipt. The new-map screen has no whole-code certificate.

[inner_fusion_prepare.py](inner_fusion_prepare.py) authenticates and extracts
the prior winner into `PacketReference.h`. Generate `T64Reference.h` with
the earlier packet generator. The new streaming and direct-feedback
headers have separate generators; the conversion headers are handwritten
and exhaustively checked. [inner_fusion_run.sh](inner_fusion_run.sh)
provides serial build, check, screen, confirm and profile phases under
all three benchmark locks. Final builds add `-mavx512vbmi` to the retained
GCC 15.2.0 flags. This remains a research build linked to retained Peach
objects, not a standalone consumer API.

Campaigns are `/tmp/spin-inner-packet-ihP4SR` (first fifteen variants),
`/tmp/spin-inner-packet-q8T5o5` (VBMI and early update), and
`/tmp/spin-inner-packet-qw5uUw` (final 27-mode executable). Final confirmation
is `confirm-DuOriZ`; phase timing is `profile-Uou8H6`; complete validation
is `check-yqoWCJ`. The final campaign is archived locally under ignored
`tmp/inner-packet/fusion-final/`. Raw measurements and disassembly are not
versioned experiment data.

| Final artifact | SHA256 |
|---|---|
| Candidate probe | `2b6e197f555ba3fec597ca0b77905751cde0754e72d753ee520a50daf86ded92` |
| VBMI conversions | `c13dba64599bac5c67159ea03c03211ec9176f8548acb0908ddc834088c84e94` |
| Complete executable | `4ce4ce451a523fb90c0c8030cbd9412b66cf28c6340d0b54741c7536244cb767` |

The [small-size follow-up](inner_size_REPORT.md) extracts the reusable
kernel and completes the K=2^16/K=2^18 tuning campaign. It retains the existing
conversion as a checked non-VBMI fallback. At K=2^20, larger gains require attacking
the approximately 3.1-ms outer phase or changing the construction; another
small reduction in inner source XOR counts is unlikely to dominate.

## Reproduction and Provenance

[encoder_cost_run.sh](encoder_cost_run.sh) supplies serialized `build`,
`coarse`, `footprint`, `sample`, and `counters` phases. It is a research
runner for the retained Peach objects, not a standalone library build.
The original campaign is `/tmp/spin-cost-profile-LaoICZ`; all collected
records are copied to ignored `tmp/encoder-cost/records/` locally.

| Artifact | SHA256 |
|---|---|
| Diagnostic source | `884beb7841f6c6c515e75355127f2a0bc2cb4b0ebec6664cfdcf8feeeee5dbff` |
| Diagnostic executable | `6b8fedb5182da60299be4eed4a68465ae0314619b3e23519a63dac6b7f11965e` |
| Exact linked outer object | `c61633de0d685f7e229c5882c0db2f55c244fd68eda69a5cf10ad6e8a25b1abd` |
| Unchanged t64 executable | `efd97950c52ba6df3d7f45ed757043f6428d45c3f33f9865b4e580e7d3ea3af7` |
| Unchanged S19 executable | `3b61d526a74314e23ea14880f1dede7f4d9a212cc0ebd43b54fe7770d52230a3` |

The [t64 closure](packed_mixing/s16_closure/CLOSURE.md) and
[retained S19 closure](packed_mixing/R4_CLOSURE.md) are unchanged.
