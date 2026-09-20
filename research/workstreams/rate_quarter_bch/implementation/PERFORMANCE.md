# Quarter-rate transposed encoder performance

The fixed BCH-derived [128,32,32] outer with RM2Sub `(t,s)=(128,19)` takes **17.996 ms** at `K=2^20`.
This is the complete in-place transposed encoder, not an isolated inner kernel.
Each 128-bit block represents 128 parallel binary instances. The operation reads `N=4K` blocks and replaces the first `K` blocks.

| Message size K | Encoded size N | In-place time | Three run medians (ms) | Tile rows |
|---:|---:|---:|---|---:|
| 2^16 | 2^18 | 0.764 ms | 0.757, 0.764, 0.775 | 256 |
| 2^18 | 2^20 | 3.849 ms | 3.829, 3.900, 3.849 | 2048 |
| 2^20 | 2^22 | 17.996 ms | 17.939, 17.996, 18.054 | 4096 |

The same build measures the existing rate-half [256,128] / t128_s19 encoder at **11.214 ms** for `K=2^20`.
The quarter-rate encoder takes **1.60x** as long at the same message size, while processing twice as many encoded coordinates.
The rate-half measurement is a throughput reference, not an equal-distance comparison.

## Distance certificates

The implemented outer spans exactly the certified fixed [128,32,32] code. Its systematic message basis does not change its spectrum.
The inner reuses the exact selected columns, field multiplication, and optimized kernels from the existing t128_s19 implementation.
For the setup distribution in [the outward certificate](../SMALLER_OUTWARD_CERTIFICATE.md), the two operating points at `K=2^20` are:

| Relative distance target | Certified setup-failure probability | Margin diagnostic |
|---:|---:|---:|
| 16.5% | < 2^-40 | 41.083485 bits |
| 19% | < 2^-30 | 30.052513 bits |

Both rows use the same encoder and have the same runtime. The performance rows at smaller K do not add new distance certificates.

## Optimizations and memory

The generated outer uses 664 XORs per row, versus 1,568 for dense evaluation. AVX2 evaluates two independent rows together.
The implementation retains the shared unrolled inner, optimized nibble grouping, pruned zeta transform, and two-stage routing.
There are no online allocations or runtime callback wrappers. Outer selection dispatches once before the specialized hot loop.

A bounded sweep compared packed24 and indices32 routing, with 128 through 32,768 rows per tile at K=2^20.
A follow-up in-place sweep selected packed24 with 4,096 rows. Representative screening medians were:

| K=2^20, in-place | Time |
|---|---:|
| Packed24, 2,048 rows | 19.343 ms |
| Packed24, 4,096 rows | 18.221 ms |
| Packed24, 8,192 rows | 19.100 ms |
| Indices32, 4,096 rows | 19.587 ms |

These screening runs used 31 trials. The headline result uses the separate three-run measurement above.
Smaller tile sweeps selected 256 rows at K=2^16 and 2,048 rows at K=2^18. These choices are now the defaults.

At K=2^20, the in-place buffer occupies 64 MiB; reusable workspace occupies 72 MiB.
Retained setup occupies 26.375 MiB after discarding diagnostics and unused routing.
Setup itself took a median 81.6 ms, outside the online interval.

## Measurement and verification

Measurements used Peach, Ryzen 9 7950X, one thread pinned to CPU 15, GCC 15.2.0, and the release flags in PERFORMANCE.json.
The recorded machine state had frequency boost disabled. Runs share the existing benchmark lock and reject other visible benchmark executables.
Each timing run contains three warmups and 101 measured calls. The summary is the median of the three run medians.
Setup, allocation, and input initialization are outside the timer. In-place calls reuse the buffer without resetting or copying it between calls.
The timed entry point is `encodeUnchecked` with identical input and output pointers; the checked `encodeInplace` API adds argument validation.
Thus later calls see the preceding output prefix; the encoder uses the same fixed operation schedule for every input value.
One separate out-of-place check measured 17.701 ms at K=2^20. It is not the in-place headline measurement.

Release and AddressSanitizer/UndefinedBehaviorSanitizer tests pass for both outers. The quarter-rate tests cover:

- Every outer coordinate impulse in both AVX2 lanes.
- Full independent dense oracles at K=2^16, 2^18, and 2^20.
- Both routing layouts, alternate tile sizes, a second setup, boundary impulses, zero input, and linearity.
- In-place equality, unchanged buffer suffix, compaction, and rejected invalid API arguments.

Three Python identity tests verify source hashes, equality of the outer row spaces, and the exact selected inner.
The existing 28 quarter-rate proof tests also pass. The code generator verifies its XOR circuit symbolically.

See [README.md](README.md) for build commands and the API. [receipts](receipts/) retains compact timing rows, test logs, environment data, and binary hashes.
Rebuild this summary with `python workstreams/rate_quarter_bch/implementation/report.py`.
