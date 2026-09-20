# Transposed binary encoder comparison

Checkpoint scope: code and documentation only. The new JSON results and raw
measurement records referenced below remain local, not committed. The table
and report retain their measured values; rerunning the measurement scripts
is necessary to obtain a fresh data bundle from this checkout alone.

This is the reproducible comparison for the paper's
`sec:transposed-comparison` and `tab:transposed-comparison`.
It measures SPIN, chosen-block BAA, Expand--Convolute (not
Expand--Accumulate), and **binary RAA only**. Forward-encoder measurements
are separate future work.

## Results

Median latency in milliseconds; each cell is the median of three process
medians. Each process has three warmups and 101 timed calls. Inputs are
initialized once and never copied or reset between calls.

| Encoder | Rate | Distance / margin (largest K except *) | K near 2^16 | K near 2^18 | K near 2^20 |
|---|---:|---:|---:|---:|---:|
| SPIN, BCH-256/RM2Sub (128,19) | 1/2 | 0.10 / 50.448 bits | 0.561 | 2.323 | 11.104 |
| Chosen Golay BAA | 1/2 | 0.10 / 46.386 bits | 0.595 | 2.781 | 23.997 |
| Chosen RM BAA | 1/2 | 0.10 / 47.093 bits | 0.866 | 3.807 | 27.660 |
| Expand--Convolute (33,25) | 1/2 | 0.05 / 20 bits | 3.189 | 13.128 | 181.825 |
| Binary RAA, repetition 3 | 1/3 | unassigned | 2.834 | 14.171 | 87.643 |
| Binary RAA, repetition 4 | 1/4 | 0.19 / ~13 bits (reference size)* | 1.009 | 6.643 | 45.176 |

The [retained JSON](results_no_reset_20260911.json) contains all 5454 samples,
process medians, output hashes, compiler commands, hardware information,
and measured source/binary hashes. There are 54 serial processes and
18 configuration/length cells. No two benchmark processes were launched
concurrently. The machine is Peach, Ryzen 9 7950X, Linux, GCC 15.2.0, CPU 15.

BAA uses K=65,472; 262,080; 1,048,512 and N=2K. Each K is 64 below the
corresponding power of two. This matches the companion paper's common
Golay/RM tiling at the largest size. The other encoders use K=2^m.
RAA receives R*K blocks; the other encoders receive 2K blocks.
One block is 128 parallel binary instances, not one packed binary codeword.

At the largest sizes, the BAA/SPIN latency ratios are 2.16 for Golay and
2.49 for RM. Expand--Convolute/SPIN is 16.37, but EC has the lower guarantee
shown in the table. These ratios describe these
implementations on this host, not the fastest possible implementation of
each code family. In particular, RAA uses the Feistel backend shared with
the chosen-BAA comparison, not an independently tuned RAA library.
Repetition 3 needs cycle walking on a non-power-of-two domain; repetition 4
does not. A precomputed-permutation RAA variant is a worthwhile follow-up.

## Parameters and distance context

- SPIN: unchanged sources from [bare_bch_rm2sub](../bare_bch_rm2sub/README.md),
  selected T128S19, grouped nibble kernel, packed 24-bit routing, default
  tiles (256 rows at m=16/18, 2048 at m=20), route/coefficient seeds 1/2.
- BAA: upstream `ChosenBlockBaa` with `Golay24` or `Rm32`, XOR circuits on,
  one thread, four-round domain Feistel permutations. No upstream kernel edits.
- Expand--Convolute: upstream `ExConvCodeOld::dualEncode<block>(input, output)`,
  the legacy nonsystematic implementation, with weight 33 and memory 25.
  Its API argument is **24 random taps**, followed by one fixed tap at
  offset 25. It applies one nonwrapping convolution pass and the original
  expander. No upstream kernel edits. Both buffers are preallocated; the
  input is modified in place and the result is written to the output.
  We do not use the allocating one-buffer convenience overload.
- RAA: [BinaryRaa.h](BinaryRaa.h), with two independent permutation seeds,
  binary prefix scans, and contiguous repetition blocks of size 3 or 4.
  It uses eight-way scan/scatter batching and preallocated workspaces.
  No diagonal weights or GF(2^128) multiplications are used.

At K=2^20, the SPIN finite certificate gives 50.448203 bits at relative
distance 0.10. The companion chosen-BAA manuscript reports full binary
first-moment bounds of 46.386 bits (Golay) and 47.093 bits (RM) at
N=2,097,024 and the same distance target. Its stated bounds include a
one-bit numerical rounding allowance. Source: `gen-BAA` revision
[`bb18ff8b7d20fd66efec79235accbc9e68c7ddb3`](https://github.com/ladnir/gen-BAA/tree/bb18ff8b7d20fd66efec79235accbc9e68c7ddb3),
binary instantiations and finite-certificate appendix. We read those
results; this comparison does not independently replay the BAA certificate.

The EC row uses the older conservative profile, explicitly accepted for
comparison at its **lower relative distance 0.05 and margin 20 bits**.
The source is the table in [the earlier BA evaluation](../../BA_paper/eval.tex)
at K=2^20. This is not a 40-bit, distance-0.10 EC claim. The paper table's
distance/margin column applies only to the largest measured length, except
for the explicitly marked RAA reference below.

These are ensemble bounds, not an individual setup's measured distance.
The BAA timing implementation uses the companion paper's pseudorandom
permutations. We do not transfer the largest-size BAA/EC margins to smaller sizes.

*For repetition-4 RAA, [Blaze](https://eprint.iacr.org/2024/1609), Section 3.1,
reports relative distance 0.19 with approximately 13 bits of margin **without
setup testing at K=2^22** (N=2^24). This is a published reference bound,
not a bound evaluated at our timed K=2^20. One-time testing of low-weight
messages can improve the margin to about 40 bits: Figure 2 reports 41.5 bits
at K=2^22 with near-quadratic setup work. These tests are outside encoding
time and are not implemented or run by our benchmark. Repetition-3 remains
unassigned; Blaze's stated theorem covers even repetition factors at least 4.
The authors provide a [calculation notebook](https://github.com/raa-code-analysis/raa-code-analysis).
Evaluating the no-test bound at our exact sizes remains a follow-up; do not
extrapolate the 13-bit or 41.5-bit reference values to those sizes.

## Build and reproduce

Requires Linux, Python 3, CMake 3.20+, C++20, and an x86 CPU with the
instructions explicitly enabled by CMake. The published CPU target is
`znver4`; change it only for a new, separately labeled host measurement.
Run from the repository root:

```sh
git clone https://github.com/ladnir/libOTe.git out/comparison-libote
git -C out/comparison-libote checkout fa7bfc3c8a0178fc99b0b666d5f31fc2a9277579
git -C out/comparison-libote submodule update --init cryptoTools
cmake -S workstreams/transposed_comparison -B out/transposed-comparison \
  -DLIBOTE_SOURCE="$PWD/out/comparison-libote" -DCMAKE_BUILD_TYPE=Release \
  -DCOMPARISON_ARCH=znver4
cmake --build out/transposed-comparison --target transposed_benchmark -j 4
ctest --test-dir out/transposed-comparison --output-on-failure
python3 workstreams/transposed_comparison/run.py \
  out/transposed-comparison/transposed_benchmark out/comparison-fresh.json
```

cryptoTools is pinned by the submodule to
`643b3ca0b57ee4ec5e3388992c5906e26f3ca72d`. Its build fetches libdivide
v5.2.0. No full libOTe protocol build or networking dependencies are needed.
Dependency clones and build products belong in ignored `out/`, not Git.

**Run benchmarks one at a time.** The harness takes
`/tmp/bare-spin-benchmark.lock`, checks `/proc` for competing benchmark
executables, and pins CPU 15. These guards supplement, rather than replace,
coordination with other users of the host. `run.py` launches each process
synchronously and refuses to overwrite a result file. An interrupted run
has `status=incomplete` and must not be reported as the complete comparison.

To validate retained records and regenerate/check the paper table without
running a benchmark:

```sh
python -B workstreams/transposed_comparison/report.py
python -B workstreams/transposed_comparison/report.py --check
```

## Timing and correctness boundary

The harness allocates and initializes native buffers once. After warmup,
it calls each encoder repeatedly without a copy, reset, or input generation.
EC operates on its current, modified input contents. The other interfaces
preserve their input. There is no extra input-snapshot buffer.
The timed interval includes the encoder's full online work, including
on-demand permutation and coefficient generation. It excludes setup,
caller buffer allocation, process checks, and output hashing. Allocations
inside an upstream encoder call, including its PRNG construction, remain
inside the timed interval. The process guard runs
before warmup and after the batch, not between encoding calls.
BAA/RAA do not include a forward pass. RAA's forward routine is a correctness
oracle only.

The checks compare SPIN with its dense oracle and BAA with its upstream
reference, compare legacy EC with a separate scalar traversal of the
random/fixed taps and expander indices on 128-bit blocks, and test RAA's bilinear transpose identity on
boundary and non-power-of-two lengths. These adapter checks complement
the source implementations' tests; they are not a new distance proof.

The historical [SPIN performance record](../bare_bch_rm2sub/PERFORMANCE.json)
is unchanged. Comparisons use the fresh SPIN row from this campaign, not
the historical 11.259 ms. Neither record is a cold-cache or end-to-end protocol
measurement. Source hashes in the result are byte-sensitive; Git line-ending
conversion can cause a hash check to fail without changing the code.

## Superseded input-reset experiment

[results_20260911.json](results_20260911.json) retains the first campaign
for provenance only. It copied input before every call, scanned `/proc`
between calls, and measured a different EC implementation at `(7,16)`.
Its 38.199/41.026 ms BAA values and 24.546 ms EC value are superseded.
Do not use its EC row as a proved baseline or mix its timings into the
current table. Its hashes describe the earlier harness, not current sources.

The two campaigns differ in buffer/cache policy and scheduling overhead;
they are not a controlled experiment isolating memcpy alone. The no-reset
BAA timing loop now follows the companion benchmark's repeated-input policy.
The companion's 22.810 ms Golay and 24.278 ms RM results used Intel
i9-14900HX/MSVC rather than Ryzen 7950X/GCC.

Next useful step: compare precomputed and on-demand RAA routing, retaining
its explicit certification status. Forward-encoder timings remain a
separate extension. The older EC point is intentionally included at its
weaker distance and margin, not treated as security-matched to SPIN.
