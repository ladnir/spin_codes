# Packed GL32/BCH Implementation Tuning

Follow-up: [the GFNI inner experiment](gfni_inner_REPORT.md) lowers the same
certified encoder to 5.334629 ms. The measurements below remain the matched
record of the preceding BCH optimization.

2026-09-30. Wide output stores and compact coefficient storage reduce the
certified R4 encoder to **5.854159 ms**, versus **6.156269 ms** for the
retained BCH kernel in the same comparison executable: **4.91% less time**.
The GL32 coefficient table shrinks from **16 MiB to 8 MiB**. A later
tile-comparison batch measures the retained candidate at 5.889110 ms;
**about 5.85--5.90 ms** is the practical summary of this campaign.

The workload is the complete precomputed transposed encoder at K=2^20,
rate one-half, with 128-bit XOR elements on Peach's Ryzen 7950X. Setup,
allocation, input filling, reference checks, and checksums are outside
the timed interval. The four-update inner retains its exact-map scalar
fusion from [the previous optimization](packed_fused_r4_REPORT.md).

No distribution or code parameter changes. The
[>10% distance / >52.05-bit certificate](packed_mixing/R4_CLOSURE.md)
therefore gives the same whole-code setup-failure bound under the unchanged
ideal setup distribution. It does not individually certify the benchmark
seeds. Production defaults and paper claims remain unchanged.

## Retained Changes

The existing BCH kernel extracts and stores each 128-bit row element
separately. The new output stage rearranges four adjacent column vectors
into contiguous row vectors, then writes them with 512-bit stores. It
retains the two-output-group, paired-input arithmetic schedule, systematic
coordinates, and exceptional parity coordinate 127.

Each GL32 coefficient previously appeared twice, once for each 64-bit
half of a row lane. Setup now retains one copy in a 64-byte-aligned table.
The kernel loads four distinct coefficients and expands them in registers
using the fixed indices [0,0,1,1,2,2,3,3]. This is a representation change
to the original sampled matrices, not resampling or a restricted family.

The hot path has no new allocation or indirect dispatch. Kernel selection
is static; loops and fixed-width operations remain explicit. The complete
comparison executable allocates all coefficient forms outside timing for
every mode. Thus 8 MiB describes the chosen representation, not all memory
allocated by the test harness. Merely aligning the expanded table did not
give a consistent gain.

## Matched Measurements

GCC 15.2, pinned CPU 15, the retained optimization/ISA flags. Each process
uses three warmups and 101 measured calls. Both process orders are tested
for seeds 1 and 17. All compilation and benchmarks run serially under the
three shared benchmark locks. No samples or processes are discarded.
There is no input reset between timed calls, matching the prior harness.

The primary matched batch is `coeff-1-confirm-bRhThg`:

| Seed | Order | Retained BCH (ms) | Wide stores (ms) | Wide stores + compact coefficients (ms) |
|---|---|---:|---:|---:|
| 1 | retained, wide, compact | 6.160612 | 6.082456 | 5.847717 |
| 1 | compact, wide, retained | 6.168366 | 6.086272 | 5.884456 |
| 17 | retained, wide, compact | 6.151926 | 6.088036 | 5.848338 |
| 17 | compact, wide, retained | 6.134473 | 6.078678 | 5.859980 |
| | Median of process medians | **6.156269** | **6.084364** | **5.854159** |

The primary saving is 0.302110 ms. The historical 6.212096 ms measurement
used the previous executable; it is retained, but is not the denominator
for the matched 4.91% claim here. The rerun's lower baseline is not credited
to this optimization.

A separate 31-call instrumented batch, `coeff-1-profile-F8xc2I`, gives:

| Kernel | Inner + route phase means (ms) | Packed GL32/BCH phase means (ms) |
|---|---:|---:|
| Retained BCH | 2.738071--2.744678 | 3.383521--3.392282 |
| Wide stores | 2.739700--2.743303 | 3.339793--3.366469 |
| Wide stores + compact coefficients | 2.711260--2.714679 | 3.165324--3.198130 |

These include warmups and are not a decomposition of the headline medians.
The inner algorithm is unchanged; its smaller timing variation is not
claimed as a separate algorithmic improvement.

## Alternatives Tested

- Fully unrolling BCH input pairs did not beat the selected bounded loop.
- Three-output-group tiles were slower, approximately 6.65--6.68 ms with
  expanded coefficients in the initial screen.
- One-output-group tiles with compact coefficients were competitive. In
  `coeff-1-tile-compare-A4xQiX`, their median is 5.858843 ms, versus
  5.889110 ms for two-output-group tiles. The 0.030267 ms difference is
  small relative to the variation across this campaign. We retain the
  original two-group arithmetic schedule rather than claim a robust win.
- A new byte/GFNI packing layout passes exact symbolic and compiled checks
  but loses in the screen: 6.364362 and 6.363670 ms with compact
  coefficients. Fewer source-level packing operations did not improve
  full-encoder throughput. It remains an isolated negative result.

The initial BCH-only variants are selected by `SPIN_PACKED_TUNE`: 0 is
retained, 1 wide stores, 2 unrolled pairs, 3 one-output-group tiles,
4 three-output-group tiles, and 5 byte/GFNI packing. Their CSV `mode` field
identifies the inner mode, not this compile-time selection; log filenames
record the BCH variant. The coefficient driver reports both its storage
mode and the selected BCH tile mode explicitly.

## Correctness and Scope

For seeds 1 and 17, compiled checks pass:

1. All 131,072 physical input basis vectors of a four-row tile, against
   the independent scalar GL32 action and production BCH reference.
2. Dense, sparse, and alternate-dense full inputs at K=2^14 and K=2^20,
   with original sequential/dense R4 references, inner adjoint checks,
   full output equality, and the untouched suffix.
3. Deliberately 16-byte-offset output pointers and prefix/suffix canaries.
4. Agreement of all 88 timed-process checksums, grouped by seed and call
   count, including the slower candidates.

Generator checks cover BCH submatrices, transpose basis vectors, output
lane order, tile coverage, and the alternate layout's invertible payload
permutation. Independent source review found no blocker in coefficient
expansion, alignment, parity handling, reference lifetimes, or tail stores.
No new sanitizer campaign was run.

The original proof receipt remains unchanged, with SHA256
`637a1ac049a7aa81441705e5a4fbeebf79e78fa80e63fe54f93cfe6dad05a8f0`.
No numerical proof replay is needed for this exact-map implementation change.

## Sources and Reproduction

- [BCH schedule generator](packed_bch_tune_codegen.py) and
  [serial schedule runner](packed_bch_tune_run.sh).
- [Coefficient-layout generator](packed_coeff_codegen.py),
  [comparison driver](packed_coeff_r4.cpp), and
  [serial coefficient runner](packed_coeff_run.sh).

From the prepared remote build directory, the selected run is:

```sh
bash packed_coeff_run.sh /tmp/spin-bch-tune-fN0NiP build 1
bash packed_coeff_run.sh /tmp/spin-bch-tune-fN0NiP check 1
bash packed_coeff_run.sh /tmp/spin-bch-tune-fN0NiP confirm 1
# Explicit selected implementation: tile mode 1, coefficient mode 2.
taskset -c 15 /tmp/spin-bch-tune-fN0NiP/coeff-1 20 1 101 2 0
```

The runner names the retained reference headers and object files. It is a
research runner for this machine, not an independently packaged library.
The generated tile-mode-1 source is byte-identical before and after adding
the optional tile-mode-5 generator support.

| Source or executable | SHA256 |
|---|---|
| `packed_bch_tune_codegen.py` | `d58daea63388b6bb9cf0a6b3d7f3c7014d1f7b04cff3bf21c34e22369dab26eb` |
| `packed_coeff_codegen.py` | `2526f7caf159fc4b209a89c9063633799c9b4eea74a9d66b84a8cb235182ffbe` |
| `packed_coeff_r4.cpp` | `543f1b65c92e168248e0dc36cb896babe9d252c3529f7bb14d6d1e8913c316f6` |
| Generated `PackedCoeff1.cpp` | `b3020fdfc4c7ba8fc6a8d87306cfd04d7c9084afe6258c7d0a96f549ca0178b3` |
| Selected `coeff-1` executable | `335638fb505d87b4191bebdb4e0181ef245df87991ca5f5d73e4b26aa9c4fdd0` |

All measurement logs, generated sources, the selected executable, and compiler stack/symbol reports
are retained locally under ignored `tmp/packed-hill/bch-tune-performance/`.
No raw experiment data is committed. The retained R2 and sequential/fused
R4 sources are unchanged.

Next: integrate the compact-table and wide-store kernel into the reusable
precomputed path, retaining this exact-map baseline and scalar oracle.
Further large gains likely need a new approach to the remaining roughly
3.2 ms packed GL32/BCH work or 2.7 ms inner/routing work; this campaign does
not establish which one has the larger attainable improvement.
