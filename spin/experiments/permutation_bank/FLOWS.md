# Two full transposed-encoder flows

The initial measurements below are retained for comparison. The subsequent
[optimization report](OPTIMIZATION.md) reduces the fresh 16-entry flow to
2.211 ms and tracks the still-open 10% overhead target.

`Flows.h` implements both flows for K=2^18, BCH [256,128], T128S19, and one
stream of 128-bit elements. These are experimental entry points, not installed
public APIs. Production descriptors, default sampling, and libOTe remain
unchanged. Other sizes and forward encoding are outside this experiment.

## Fully precomputed code

`PrecomputedCode(routeSeed, maskSeed)` samples the original shuffle family
and prepares all routing schedules, IMT masks, and scratch once. Each later
call performs only encoding. The code, route, and masks remain fixed.
This uses the existing tuned transposed kernel, not the slower compact
materialized-route experiment.

```cpp
using namespace spin::experimental::bank;
PrecomputedCode code(17, 29);           // One-time setup, original family.
std::vector<kernel::block> values(CodeSize);
// Fill values with the next input, then:
code.encode(values.data());            // Reuse code for subsequent inputs.
```

CodeSize=2*K. Encoding is in-place; the first K entries contain the result.
The wrapper owns its mutable workspace, so do not call it concurrently.
The constructor also accepts a fully materialized bank-derived route for
same-map comparisons. The bank is no longer needed after materialization.

## Precomputed bank, fresh parameters

`BankFlow<16>` or `BankFlow<64>` owns persistent packed inverse tables.
Each `fresh(routeSeed, maskSeed)` call samples new bank indices, XOR offsets,
additions, rotations, and IMT masks. It does not build a full route. During
encoding, the inner emits values and routing computes batches of 16
destination addresses for immediate use.

```cpp
BankFlow<16> bank(913);                 // One-time reusable bank.
std::vector<kernel::block> scratch(CodeSize); // Reuse across fresh codes.
std::vector<kernel::block> values(CodeSize);
// For each fresh code, choose the agreed seeds and fill values:
auto code = bank.fresh(routeSeed, maskSeed);
code.encode(values.data(), scratch);   // No full route array.
```

The returned code borrows the bank; the bank must outlive it. Moving/copying
the owner and creating a code from a temporary owner are disabled. Concurrent
encodes need separate input and scratch buffers. Encoder arithmetic is shared
with the compact experimental kernel through a compile-time route type.
There is no virtual dispatch or encode-time allocation.

This is the heuristic family evaluated in `README.md`, not the uniform
permutation family of the existing distance certificates. Correctness tests
and family screens do not transfer those certificates.

## Full-encoder measurements

Ryzen 7950X, CPU 15, GCC 15.2, `-O3 -mavx2 -mtune=znver4`, runtime AVX-512
BCH. Each call encodes N=524,288 128-bit input elements into K=262,144 output
elements.

| Flow | Per-code setup | Encode | Total per call |
|---|---:|---:|---:|
| Fully precomputed, original family | None | 1.486 ms | 1.486 ms |
| Fully precomputed, 16-entry bank-derived code | None | 1.475 ms | 1.475 ms |
| Fully precomputed, 64-entry bank-derived code | None | 1.495 ms | 1.495 ms |
| Persistent 16-entry bank, fresh offsets and masks | 0.0378 ms | 3.554 ms | 3.592 ms |
| Persistent 64-entry bank, fresh offsets and masks | 0.0381 ms | 3.678 ms | 3.716 ms |

Fixing the bank-derived code restores approximately the same throughput as
fixing an original-family code. On-demand address computation and its
interaction with data movement account for most of the fresh-flow overhead.

The original-family fixed flow took about 13.56 ms once to prepare its code
and workspace. Bank-derived fixed flows took about 10.52/10.53 ms for route
materialization, schedules, masks, and workspace, plus bank construction.
These costs are not paid on each encode.

For fresh flows, bank construction took about 0.10/0.38 ms for 16/64 entries.
Allocating and initializing reusable scratch took 1.80/1.83 ms once. Bank
modes validate the encoder before timing, warming allocator and generation
code; these are not cold process-startup measurements. Earlier cold
address-only experiments measured higher bank construction times. Fresh setup
includes allocation of the small parameter and mask vectors; scratch is
retained across calls.

Storage, excluding the input/output buffer:

| Flow | Persistent bank | Per-code plan | Scratch |
|---|---:|---:|---:|
| Fully precomputed | Not needed after setup | 11,075,596 bytes | 9 MiB |
| Fresh, 16 entries | 68 KiB | 59 KiB | 8 MiB |
| Fresh, 64 entries | 272 KiB | 59 KiB | 8 MiB |

The existing full plan retains schedules beyond those needed by the compact
on-demand plan; the storage difference is not just a route array. The fresh
flow never materializes the 2 MiB route during timed setup or encoding.

## Measurement and validation boundaries

Timings are medians of five process medians. Flow order reverses on alternating
passes. Each process performs 13 encodes, discards two, and initializes a
separate input for every encode. Fixed flows retain one code; fresh flows
advance both instance and mask seeds each call. Input initialization and
output checksumming are outside timing. Fresh totals include per-code plan
construction and destruction. Component medians need not sum exactly.

Process-median totals ranged from 1.473--1.491 ms for the original fixed flow,
3.584--3.683 ms for fresh 16-entry encoding, and 3.596--3.824 ms for fresh
64-entry encoding. These measure the linear map only: no base OT, PPRF,
communication, hashing, or separate receiver choice-byte encoding. They are
not full PCG or Silent OT timings.

Validation compares complete in-place buffers between both flows for each
bank size, three code seeds, and zero/random inputs. Checks passed on AVX2
locally and AVX-512 on Ryzen. Focused ASan/UBSan checks instrument new wrappers,
routing, and inline inner code while linking existing release kernels; this
is not a full sanitizer rebuild of the reference library. Existing API and
descriptor known-answer tests also pass. Production sampling is unchanged.

## Reproduce

```sh
cmake -S spin -B out/flows -DCMAKE_BUILD_TYPE=Release \
  -DSPIN_BUILD_EXPERIMENTS=ON -DSPIN_TUNE=znver4
cmake --build out/flows --target spin_bank_flows -j2
ctest --test-dir out/flows -R '^spin_bank_flows$' --output-on-failure
out/flows/spin_bank_flows fixed-exact
out/flows/spin_bank_flows fresh16
out/flows/spin_bank_flows fresh64
```

The serial comparison expects executables under `COMPARISON_ROOT/build`:

```sh
bash spin/experiments/permutation_bank/compare_flows.sh COMPARISON_ROOT 15
python3 spin/experiments/permutation_bank/summarize_flows.py COMPARISON_ROOT/results-bank-flows
```

Raw data stays outside source control. Next optimize fused on-demand routing
and stores, not parameter or mask generation. Retain the original-family
fully precomputed flow as the fixed-code baseline.
