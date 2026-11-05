# Native implementation of the near-half-distance window

Research checkpoint, 2026-10-21. This implements the
[>49.1% / 40-bit checkpoint](NEAR_HALF_WINDOW.md) at K=2^20, rate 1/2,
over GF(2^128). The fixed RS[256,128] outer and single stratified route
are unchanged. Memory is 1024 symbols; each recursive step emits 256.

The exact Toeplitz control and the explicit coordinate variant below are
both retained. The latter changes the feedback distribution, but preserves
the inverse-input distribution used by the whole-code certificate.
Neither implementation is a library profile or a replacement for the
earlier fast small-alphabet encoder.

## Exact Toeplitz implementation

`native/ToeplitzWindow1024.h` computes the original matrix
T[i,j]=d[1023+i-j] from its 1279 sampled diagonals. Split its input into
four 256-symbol chunks. Each square product has block form

    [ A B ] [lo]     [u+w]
    [ C A ] [hi]  =  [u+v],
    u=A(lo+hi),  v=(A+C)lo,  w=(A+B)hi.

Each block is a 128-square Toeplitz map. A leaf becomes a polynomial
middle product: multiply degrees at most 254 and 127, then extract
coefficients 127 through 254. Evaluation at all 256 GF256 points computes
the product modulo Z^256-Z. Its wrapped coefficients have degrees at most
126, so this extraction is exact.

The rolling state caches eight 128-to-256 transforms. Appending one batch
requires only two new transforms. The four contributions and the shared
Karatsuba terms combine before interpolation, leaving two inverse transforms.
This control uses twelve full-field products per output.

`native/ToeplitzWindow1024Nine.h` reduces that count to nine. Regard the
eight state halves as coefficients of a degree-seven polynomial. The
two-output block map is the adjoint of multiplication by a linear polynomial.
Evaluate that multiplication at eight GF256 points and retain its leading
coefficient separately. The resulting nine leaf keys are linear combinations
of the same raw diagonals, not independently sampled spectral coefficients.
An eight-point additive transform performs the small block mixing.

Both implementations use the sparse-byte representation
GF256[Y]/(Y^16+Y^5+Y^3+0x32). Deferred arithmetic XORs unreduced products
before field reduction. The optional raw-scaling path applies GF256
scalars before reduction. For each 16-bit slot, reduction is
lo + 0x1b*hi; bytewise GF256 scaling commutes with this map and with the
extension-field reduction. This optimization changes no matrix or output.

## Explicit coordinate variant and certificate transfer

Polynomial-basis conversion accounts for much of the initial runtime.
The optional `SPIN_OR_WINDOW_COORDINATES=1` construction removes those
runtime conversions. It is not the identical Toeplitz encoder.

Let P_j be the existing novel basis, with degree j and nonzero leading
coefficient. Define V to convert 128 novel coefficients into monomial
coefficients. Define U to convert the top 128 monomial coefficients of a
degree-below-256 polynomial into its top 128 novel coefficients. Both maps
are invertible. The latter conversion is independent of the lower 128
coefficients, because lower-degree basis elements cannot affect higher terms.

During setup, multiply every prepared leaf evaluation by its evaluation
point. This replaces its diagonal polynomial D by ZD. Setup still evaluates
monomial polynomials. Runtime evaluates the state in the novel basis and
interpolates the output in that basis, extracting coefficients 128 through
255. Product degree is at most 382, so reduction modulo Z^256-Z aliases
only into degrees at most 127. Each leaf therefore computes U T_leaf V.
Both block decompositions consequently implement

    T'_t = diag(U,U) T_t diag(V,V,V,V,V,V,V,V).

The proof transfers through the inverse-input distribution, not through
Hamming-weight invariance. Fix a candidate output word y. Its preceding
raw-output window S_t(y) is fixed. The right-hand invertible map preserves
whether that whole window is zero. For a nonzero window, T_t V S_t(y) is
jointly uniform over F^256; the left-hand invertible map preserves that
distribution. Independent raw diagonal tapes preserve independence across
batches. A zero state still gives zero feedback.

Thus the same activation constraints, forced-zero tails, affine inverse
space, output count, and finite cover apply. The coordinate variant retains
the >49.1% symbol-distance / 40-bit setup-failure guarantee at the certified
size. Its output words differ from those of the exact control.

## Correctness and measurement

The research harness is `native/near_half_native.cpp`. Tests cover all
1024 state-basis inputs, dense independent schoolbook field multiplication,
diagonal endpoints, and the full 16-step recursive encoder. The coordinate
oracle constructs polynomial bases from root products and performs direct
triangular division; it does not use FFT conversion routines. The tests also
check regional bijections, the large-region routing branch, ring wraparound,
the RS outer, contiguous in-place encoding, and equality of arithmetic variants.
Raw scaling is checked on all 65536 coefficient-slot values and all 512
binary directions of an unreduced field product. UBSan checks are separate.

Benchmarks are single-threaded on Peach's Ryzen 9 7950X, pinned to CPU 2,
using GCC 15.2.0, `-O3 -std=c++20 -march=znver4`. Setup and allocation are
outside encoding time; outer encoding, routing, and the recursive inner are
inside. Each timing uses three warmups followed by the stated median sample
count. Runs are strictly serial. Fixed MT19937 fixtures reproduce arithmetic
and timings; they are not a production cryptographic setup sampler.

Selected full-encoder medians, milliseconds. The shared host shows run-to-run
variation; repeated medians below are not confidence intervals.

| Implementation | Full encoder | Samples |
|---|---:|---:|
| Exact twelve-leaf control | 154.243 | 31 |
| Exact twelve-leaf, deferred reduction | 135.417 | 31 |
| Exact nine-leaf, deferred reduction | 130.209 | 31 |
| Nine-leaf coordinate variant, deferred reduction | 75.103 | 31 |
| Same, raw-scaling optimization | 70.904; 74.303 on repeat | 31 each |

Earlier seven-sample medians were 151.884, 135.793, and 123.981 ms for the
first three rows. The fastest variant's final phase run measured outer
6.890 ms, standalone route 5.816 ms, routed inner 66.578 ms, and full encoder
74.303 ms. Routing is fused with the inner in the full encoder; do not add
the standalone routing measurement again. UBSan passes for the exact
twelve-leaf control and the nine-leaf coordinate/raw-scaling variant.

The exact implementations agree on checksum 11208093802884930491. The
coordinate implementations, including the twelve-leaf control, agree on
checksum 8004028696846130412. These are regression receipts, not proofs.

The prepared key tape drops from 384 MiB to 288 MiB. The state-transform
cache uses 32 KiB; nine-leaf scratch uses 8.5 KiB. The sampled diagonal count
remains 1279 per batch. A selected nine-leaf setup took about 278 ms,
including allocation and initialization; route sampling took 4.4 ms and
coefficient generation plus preparation took 179 ms.

Isolated hot-stream timings normalized to the encoder's operation counts
give 57.63 ms for two forward/two inverse monomial transforms per batch,
versus 7.13 ms without those conversions. They are not additive phase timings.
Reading the complete 288-MiB prepared tape alone took 6.88 ms. The nine-leaf
encoder still performs 18,874,368 full-field products per codeword. The
remaining gap to the earlier 16.6-ms RS reference is substantial; neither
the proof closure nor this approximately twofold optimization establishes
a competitive high-distance encoder.

## Reproduction and next gate

Use this native directory plus the include directory
`spin/experiments/char2_field_spin/field128_lift`. No generated tape is saved.
Build outside the source tree. For the fastest measured variant:

```text
g++ -O3 -std=c++20 -march=znver4 -I spin/experiments/char2_field_spin/field128_lift -DSPIN_OR_WINDOW_NINE=1 -DSPIN_OR_WINDOW_COORDINATES=1 -DSPIN_OR_WINDOW_RAW_SCALE=1 research/workstreams/one_route_coupled_scan/native/near_half_native.cpp -o /tmp/near_half_native
/tmp/near_half_native --check
taskset -c 2 /tmp/near_half_native --bench 20 31 all
```

Omit `SPIN_OR_WINDOW_COORDINATES` for the exact Toeplitz distribution.
Omit `SPIN_OR_WINDOW_NINE` as well for the twelve-leaf control. `full` selects
ordinary field reduction; `deferred` selects deferred reduction. `micro`
and `keys` are separate diagnostics, not end-to-end timings. Compile with
`-O1 -g1 -fsanitize=undefined -fno-sanitize-recover=all` for UBSan checks.
The harness restricts benchmark sizes to K12..K20; only K20 has the claimed
49.1% certificate. This is a research-driver bound, not a proposed library cap.

Root owns all builds, checks, and benchmarks. Runs use external limits of
120 seconds for compilation, 180 for correctness, and 90 for benchmarking,
with at most 2 GiB virtual memory. More than 20 GiB remains free on both
machines. Only source and selected receipts are retained locally.

Final source SHA-256 values match the remote build inputs:

| File in native/ | SHA-256 |
|---|---|
| near_half_native.cpp | 0c628ae34c1b2b3d8155b89b0bb83b847f7e9188ee3be6300b975e8d39d973cd |
| ToeplitzWindow1024.h | d0fa39bc0fe7ce8eb9c595e9ddc8f848274cef689df527bdc8be84bfab061ed8 |
| ToeplitzWindow1024Nine.h | ea48f4d55599dde809234fde2aeaaef1bb30169c7266085c333ad9b666157800 |

All six transform/field dependency headers also match local sources. The
four certificate sources retain the hashes in the near-half proof note.
The independent final static review found no correction to the coordinate
oracle, block factorization, raw scaling, or certificate-transfer argument.
No numerical process remains running; both machines have over 60 GiB free.

Next, reduce full-field arithmetic while preserving the actual inverse-law
proof. The current measurements do not justify further distance hill climbing
or assuming another small kernel change will close the performance gap.

## Phase breakdown and the 3--5 ms target (2026-10-21)

The measurement-only driver `native/near_half_phases.cpp` reuses the frozen
harness and setup. `native/NearHalfPhaseKernel.h` splits the fastest nine-leaf
coordinate kernel into its pointwise stage and inverse transforms. Neither
file changes the retained encoder or the certified distribution.

At K=2^20, the retained encoder takes 74.319 ms. The same computation split
into untimed phase functions takes 74.332 ms. Instrumenting every batch takes
75.255 ms. These are 15-sample medians after three warmups, with alternating
order, on the same pinned Ryzen core and compiler configuration. The sampled
core frequency was 4.5 GHz. Setup and allocation remain excluded.

| Instrumented phase | Median ms |
|---|---:|
| RS outer, including stores in regional layout | 7.788 |
| State combinations, full-field products, reduction, and key reads | 55.280 |
| Inverse transforms and feedback extraction | 3.650 |
| Routed input reads, feedback XOR, and output stores | 4.835 |
| Append output to state and compute forward transforms | 3.470 |
| Remaining loop/timer work | 0.203 |
| Instrumented full encoder | 75.255 |

Phase medians need not sum exactly to the total median. Instrumentation adds
about 0.94 ms; some timer overhead is charged to adjacent phases, not solely
to the residual row. Full outputs match the frozen encoder byte-for-byte at
K12 and K20, including both routing branches. The K20 checksum remains
8004028696846130412.

Two separate diagnostics locate the dominant costs:

- With a fixed state and streamed prepared keys, the pointwise stage takes
  55.765 ms. Reusing one prepared key takes 49.474 ms. This is a cache
  diagnostic, not a valid replacement code or certificate. The key, state,
  and scratch together exceed L1. The difference is not an additive memory
  charge, but the remaining cost rules out key traffic as the sole bottleneck.
- The identical four-block RS arithmetic takes 6.984 ms with regional stores
  and 3.547 ms with contiguous four-block tiles. After reindexing, all output
  values agree. This is an isolated outer comparison: tiled storage scatters
  the later reads for each region across the full workspace. The end-to-end
  routing tradeoff remains unmeasured. Controls without the FFT take 5.862
  and 2.338 ms respectively; subtracting them is not an exact FFT cost.

Reading the 288-MiB key tape alone takes 6.870 ms. It is already larger than
the whole target budget on this core. Reducing tape traffic alone is also
insufficient: the cache-hot pointwise diagnostic still takes about 49 ms.

The current implementation performs nine full-field products per output.
For the unchanged 256-by-1024 Toeplitz family, a bilinear-rank lower bound is
1279 products per batch, or about five per output. The independent diagonal
matrices are linearly independent, giving this tensor-flattening bound; fixed
invertible coordinate changes preserve it. This is not a lower bound for
other code families, binary circuits, or hardware instruction counts.
It nevertheless limits the benefit of another factorization of this map.

The target therefore requires more than optimizing the existing factorization.
An illustrative engineering budget is 1--1.5 ms for the outer, 0.5--1 ms for
routing/output, and 1--2 ms for the inner. These are targets, not predictions.
Even the isolated tiled outer has not met its share of this budget.

There are two distinct next steps. First, improve the exact RS implementation
and its layout, measuring the whole routed encoder after each change. Its
four-block transform lacks the register-level small leaves used by the
single-block path. A blocked output layout may also retain some regional
locality without the current widely separated stores.

Second, reduce the work demanded by the inner's proof interface. The present
proof gets a jointly uniform 256-symbol feedback vector for every fixed
nonzero 1024-symbol state. This strong local property makes the inverse
enumeration tractable, but it also motivates substantial full-field arithmetic
and independent prepared material. A faster family may need a weaker property
that retains the outer's constraints on reachable values and cancellation.
Simply shrinking the state or reusing keys does not establish that property.
The revised proof must still cover all messages for one sampled setup.

The near-half certificate remains the reference, not a 3--5 ms performance
candidate. No smaller-state or sparse-feedback replacement is certified here.
The primary interface remains GF(2^128)-linear; a merely binary-linear mixer
would be a separate application choice.

Build this diagnostic with the same include path and optimization flags as
above, replacing the source with `native/near_half_phases.cpp`. Its three
variant macros are set in the source. Run `taskset -c 2 near_half_phases 20 15`.
Root ran the K12 and K20 diagnostics serially with explicit time and memory
bounds. No key tapes or raw timing dumps were retained.

## Exact small-RS optimization (2026-10-21)

The earlier 7--8 ms outer timing was not an inherent cost of RS[256,128].
It included poorly scheduled regional stores and an incompletely optimized
small transform. `native/RsOuterOptimized.h` now provides an exact replacement;
the frozen outer and transform headers remain unchanged as controls.

The selected `encodeRegional` path uses the full-RS implementation techniques:

1. Fuse source loads with the first two FFT levels. The known-zero upper
   coefficient half needs neither initialization nor a general butterfly.
2. Compute the upper levels for 128 outer blocks in a cache-sized tile.
3. Hold the final three levels in eight explicit SIMD registers. Fuse those
   layers with stores to the existing coordinate-major output. This avoids
   writing and rereading completed leaves in scratch, and limits the active
   output streams to eight at a time.

The field representation, novel input basis, evaluation points, generator
matrix, output order, and sampled permutation are unchanged. Scratch is
512 KiB on the ordinary stack, with no hot-path allocation. This is a K20
research selection, not a small-size dispatch policy or a library promotion.

The separate `routeRegional` path batches sixteen symbols and first reads
each contiguous region into cache with four independent SIMD chains. Actual
loads outperform a burst of nonbinding prefetch hints here. This reads the
same region an additional time; it does not change the permutation or perform
another global permutation. Both methods require separate, nonoverlapping
input/output buffers and 64-byte-aligned output. Block counts are positive
multiples of four; partial final tiles and route batches are supported.

Two independent 31-sample confirmation runs, after three warmups, gave:

| Measurement, verified huge-page buffers | Frozen control ms | Selected ms |
|---|---:|---:|
| Outer alone, first run | 8.297 | 3.303 |
| Outer plus routing, first run | 12.372 | 5.923 |
| Outer alone, second run | 8.258 | 3.245 |
| Outer plus routing, second run | 12.234 | 5.930 |

A matched ordinary-page run gave 7.951 versus 3.554 ms for the outer,
and 12.778 versus 6.519 ms for outer plus routing. Huge pages are therefore
not the explanation for the improvement. Each comparison uses identical
buffer backing for control and candidate. The driver checks `/proc/self/smaps`:
the input/work/output mappings have 16/32/32 MiB of huge-page backing when
requested, and none in the ordinary-page control. The route labels retain
the original allocation and values.

All timings use one pinned core at the sampled 4.5-GHz setting, GCC 15.2.0,
and `-O3 -std=c++20 -march=znver4`. Setup, allocation, page initialization,
checksums, and validation are excluded. Case order alternates between trials.
These ranges are repeated process medians, not confidence intervals.

Rejected alternatives remain reproducible in the diagnostic driver. Contiguous
four-block output is fast in isolation but loses after routing. Non-temporal
outer stores also lose in the combined measurement. Register leaves alone
do not fix the original store pattern. A single-block FFT with a blocked
transpose improves the combined path, but loses to fused final-layer stores.

The driver checks every outer and routed output against the frozen map at
K12 and K20. Kernel checks cover each input position in four blocks, all 128
binary field directions, dense inputs, and an independent root-product RS
generator matrix. Additional block counts 4, 12, 36, 132, and 260 test partial
tiles and route tails. Input preservation and trailing guards pass. UBSan
passes at K12 and K20; its timings are not performance measurements.

Outer-plus-routing materializes a complete permuted 32-MiB buffer. The full
encoder instead fuses routing with feedback, so its total speedup cannot be
obtained by subtracting these combined timings. `fullOptimized` in the driver
integrates only the exact outer replacement with the unchanged fused inner.
The original implementation remains directly callable for comparison.

The final-source replay measures outer 8.312 versus 3.329 ms and outer plus
routing 12.366 versus 5.952 ms, with verified huge pages and 31 samples.
Ordinary pages give 8.184 versus 3.585 ms and 13.602 versus 6.517 ms.
The matched full-encoder comparison, with three warmups and nine samples,
measures **74.406 versus 69.588 ms**. Every full output agrees byte-for-byte,
and the checksum remains 8004028696846130412. This full path uses the
original fused routing, not the new standalone cache-fill route.

Build `native/rs_outer_bench.cpp` using the same include path and release
flags as above. The CLI is `rs_outer_bench logK repeats huge [final|full]`.
Omitting the last argument runs all screening variants. `final` selects the
control and finalists; `full` also checks and times the full encoder with
the outer replaced. Run processes serially. For example:

```text
taskset -c 2 /tmp/rs_outer_bench 20 31 1 full
taskset -c 2 /tmp/rs_outer_bench 20 31 0 final
```

Selected source SHA-256 values:

| File in native/ | SHA-256 |
|---|---|
| RsOuterOptimized.h | a1c15859700fbec68681baedbe0ea1caf076a99b5dcf927edb6708d68bde0829 |
| rs_outer_bench.cpp | 085318ff167f26dcc394f4179f1aa1b92543f82c4b98825e97584809f85c9595 |

The construction and certificate do not change. No library or paper was
modified. Next, carry the same implementation audit into the feedback kernel
before treating its current timing as a structural lower bound.
