# Native GF16-coset outer and exact store fusion

One mode-5-only store expansion changes the initial negative into a native
win without changing the encoder map. The final matched median is
3.397943 ms, 6.04% below mode 4 and 16.45% below the original mode-5 build.
All checks pass; old modes retain byte-identical binaries. The distance
proof remains incomplete, so this is not a promoted construction.

The first mode-5 control is retained below. Its median was 4.044269 ms:
12.64% slower than mode 4,
15.88% slower than XOR2, and 7.22% slower than certified Toeplitz.
All correctness checks pass. This first control does not isolate the
cause of the regression, and no whole-code certificate is claimed.

## Structural change

`HDGROUP_XOR_MIXER=5` selects sixteen RS[16,8] constituents in each
128-input, 256-output outer group. The fixed column matrix is
`H=I+11^T+uu^T` over AES GF256, with the ordered labels

    06 07 0a 0b 56 57 5a 5b b6 b7 ba bb e6 e7 ea eb.

These are `06+GF16` in the basis `(01,0c,50,b0)`. The implementation uses
the explicit eight-fixed-product, 62-binary-XOR circuit from
[the algebraic proposal](HD_TWO_SUM_GF16_PROPOSAL.md). Independent
nonzero D1/D2 sampling and its scalar gauge are unchanged. Group width 16
is enforced at compile time. The paired route, independent swaps, h4
inner, interpolation refresh, output layout, and payload sizes remain
unchanged. The default remains mode 0.

Relative to mode 4, the algebraic count per equal 128-input layout drops
from 1064 to 896 field products and from 1864 to 1648 binary XORs. The
normalized diagonal count rises from 480 to 496; the saving is in the
shorter RS transforms and fixed mixer products. These are circuit counts,
not measured instruction counts. Packed diagonal keys increase by 16 bytes
per source group, or 128 KiB at K20; no new mixer keys are introduced.

This control deliberately retains the existing Group16 FFT schedule.
`HDGROUP_FUSE_LAST=1` applies only when `Group<=8`, so the mode-5 profile
reports `fused_last,0`. Its producer materializes all eight final RS values
per constituent before gathering a sixteen-vector column. Register
pressure and this different schedule are possible costs, not established
causes of the measured regression. No tile, fusion, or marker sweep was
performed for this control.

## Correctness and compatibility

The new fixture checks all 8192 input bit directions across sixteen ZMM
registers, 16 random columns, and the involution identity. Its scalar
oracle uses a literal label array and dense coefficients
`delta(i,j)+1+u_i*u_j`; it does not reuse the SIMD bit-sum circuit.

Mode 5 passes `--check` and `--large-check`, including the 30720 inner
packet input/state directions, pair swaps and wrap, K=8192 and 24576,
K20, rank-2/rank-4 controls, guards, padding, output/scratch offsets,
cached output, full aliasing, and unaligned input. No new sanitizer run
is included. The matched Toeplitz, mode-2, and mode-4 controls also pass
both commands. Rebuilt modes 0, 2, 3, and 4 are byte-identical to their
frozen binaries from the preceding experiment.

The shorter constituent has minimum support nine rather than seventeen.
Its outer comparison measure and support-floor proof must be recomputed.
The coordinating r17 sparse FLOAT gate is promising but is not a directed
whole-code proof; dense cohorts and other source weights remain obligations.

## Matched complete timing

All runs used Peach CPU 2, GCC 15.2, 501 calls, and the driver's three
warmups. Setup, allocation, and input generation were outside the timer.
Process inspection found no active benchmark before the series. Both
`hypercat-global-benchmark.lock` and `spin-char2-benchmark.lock` were held.
All twelve runs were serial, in the table's order. They began at
2026-10-09 13:15:05 UTC. Times are microseconds.

| Round | Build | Median | p10 | p90 |
|---:|---|---:|---:|---:|
| 1 | Mode 4 | 3586.705 | 3569.524 | 3608.876 |
| 1 | Mode 2 | 3473.343 | 3453.306 | 3494.353 |
| 1 | Toeplitz | 3771.921 | 3655.504 | 3800.625 |
| 1 | Mode 5 | 4046.073 | 4030.323 | 4070.558 |
| 2 | Mode 4 | 3590.282 | 3571.768 | 3610.519 |
| 2 | Mode 2 | 3495.916 | 3478.343 | 3517.726 |
| 2 | Toeplitz | 3765.379 | 3741.925 | 3793.240 |
| 2 | Mode 5 | 4040.994 | 4023.751 | 4065.378 |
| 3 | Mode 4 | 3592.877 | 3574.422 | 3611.301 |
| 3 | Mode 2 | 3489.974 | 3469.306 | 3516.204 |
| 3 | Toeplitz | 3797.338 | 3739.711 | 3824.108 |
| 3 | Mode 5 | 4044.269 | 4025.895 | 4069.597 |

The medians across rounds are 3590.282, 3489.974, 3771.921, and 4044.269
microseconds for modes 4, 2, Toeplitz, and mode 5, respectively.

## Reproduction and pins

Remote scratch: `/tmp/hd-two-sum-g16-20261009-3VElrn`. Use the exact common
GCC flags and library from [HD_TWO_SUM_NATIVE.md](HD_TWO_SUM_NATIVE.md#reproduction-and-pins).
Add `-DHDGROUP_XOR_MIXER=5` for the candidate. Its command is

    taskset -c 2 ./gf16-two-sum --bench 501 full 20 nt 16 4

The other commands use group 8 for modes 2/4 and group 4 for Toeplitz.
All use rank 4, independent pair swaps, and streaming final output.
Run checks first, acquire both locks, and never benchmark concurrently.

| File/build | SHA256 | K20 checksum |
|---|---|---:|
| Initial `HdGroupedOuter.h`, saved as `HdGroupedOuter.initial.h` | `cb29372e59b8c3ccf1f6e6c2292c2dcc550a548dff96d3aaf7e71d6ff3b92b83` | — |
| `hd_grouped_forward.cpp` | `2223df3595805e54fe548ce9beac96915f758b639794387736467facb1fb837e` | — |
| Unchanged `HdGroupedInner.h` | `4550e45454cdfb03b570384ba9a540e4cb2d39d27aa6340b10b9ac79ea5911d1` | — |
| Library | `3d3c7654037e31e8cea8e5a4e59b32a6b53a663ee89ac650f591030db4d02a38` | — |
| Toeplitz | `0b31e23d66948665135485da064868ad769d007bb5cc8a2c4c4570698f1638bc` | 7701295210980077712 |
| Mode 2 | `2c0d38b79383f54b343ec9f40997bdee0e590314d95fee10774c4648a819a595` | 2343136256414281289 |
| Mode 3, rebuilt but not timed | `700a22c2f560667e21dbbe21cc22c7af9883bd44efd306300002e289cf506222` | — |
| Mode 4 | `c1e5bb9804e6b3b7f0462c9633af62d9ef8c17ae5cacbbf56b842dccf2c5a01c` | 1506686965355192349 |
| Mode 5 | `566e48cb274b8dd9372ee8601d3cee2cc3d17099680c392b71c1caf836a0bdca` | 5924045154917714530 |

The remote logs and scripts are pinned as follows:

- `matched-three.log`: `ff075b06fda00faeeab699f7fb273e71122798829abacd8bf606d61c0ca2ac53`.
- `build-check.log`: `5a9da5f6409227678e00c0b5f439ccfec525259d013147839b3df2e48e1e3757`.
- `hd-two-sum-g16-build.sh`: `626f0a547c6bdb7ef7a7f9cb3102056d159ee4b951ba14ef684abff5cc27cec5`.
- `hd-two-sum-g16-bench.sh`: `c4f88016f3d69d0c0a39d529ba4727bbf5b1bae34a0a8f0e83e47cbaaaada582`.

## Localized diagnosis and one exact change

After the initial result was frozen, one authorized phase comparison used
501 calls per entry, serially under the same locks and CPU affinity.

| Build | Phase | Median, us | p10 | p90 |
|---|---|---:|---:|---:|
| Mode 4 | Outer | 1781.055 | 1764.453 | 1791.073 |
| Initial mode 5 | Outer | 2232.678 | 2224.903 | 2249.650 |
| Mode 4 | Inner | 1728.185 | 1719.469 | 1738.275 |
| Initial mode 5 | Inner | 1749.456 | 1733.896 | 1758.873 |

This comparison places the large regression in the outer phase. These
isolated times have different cache histories from the complete call;
they must not be added to explain its runtime.

Disassembly identifies a concrete additional boundary. GCC outlines the
sixteen-row output-store `c8Unroll` closure from `mixStore`. The initial
coset writes all sixteen mixed ZMM vectors to stack offsets `0xc0..0x480`,
then calls the helper. The helper reloads every vector, applies D2, and
stores it. It also repeatedly dereferences the closure and recomputes row
addresses. Mode 4 has no corresponding outlined output-store helper.

There are 2048 four-group batches and sixteen output columns per batch.
Consequently, this interface makes 32768 calls and materializes 32768
one-KiB mixed columns. Its vector stores and reloads total 64 MiB of
architectural stack traffic. This is not measured DRAM traffic, and it
does not include the unchanged RS workspace.

The single follow-up expands the sixteen stores explicitly inside the
existing always-inlined `mixStore`, only for mode 5. Row r still applies
`diagonalProduct(column[r], after+(16*j+r)*4)` and writes
`planes+(16*j+r)*stride`, in the same order. Cached versus streaming
selection remains compile-time. No matrix, key, route, FFT, or field
operation changes. The output-store source now exposes a register-to-store
path instead of requiring the mixed column to cross a function boundary.

The final disassembly has no output-helper call and no direct rsp/rbp-based
vector references in either hot coset. The 8 KiB RS workspace is still
present and accessed through other base registers. The frame decreases
from `0x2500` to `0x2080` bytes. The old coset plus helper contain 23+16
GFNI instructions per column body; the inline body contains the same 39.
Thus the observed removal is materialization and address/call overhead,
not an arithmetic reduction. Both offsets retain `fused_last,0`.

All mixer, packet, natural-size, and K20 checks pass again, with the same
mode-5 checksum `5924045154917714530`. Rebuilding modes 0, 2, 3, and 4
again gives their exact frozen binaries. The source change also received
an independent read-only review.

## Matched store-fusion result

The follow-up began at 2026-10-09 13:22:20 UTC. Each entry again uses 501
calls, CPU 2, both locks, and serial execution. No second tuning change
was made. Times are microseconds.

| Round | Build | Median | p10 | p90 |
|---:|---|---:|---:|---:|
| 1 | Mode 4 | 3629.475 | 3610.640 | 3652.839 |
| 1 | Initial mode 5 | 4079.625 | 4069.736 | 4104.411 |
| 1 | Inline-store mode 5 | 3436.394 | 3424.523 | 3456.833 |
| 2 | Mode 4 | 3616.210 | 3603.787 | 3638.372 |
| 2 | Initial mode 5 | 4044.460 | 4029.843 | 4067.252 |
| 2 | Inline-store mode 5 | 3397.943 | 3387.042 | 3419.082 |
| 3 | Mode 4 | 3607.435 | 3588.569 | 3630.857 |
| 3 | Initial mode 5 | 4067.052 | 4053.938 | 4093.651 |
| 3 | Inline-store mode 5 | 3395.970 | 3382.965 | 3420.415 |

The corresponding medians across rounds are 3616.210, 4067.052, and
3397.943 microseconds. Comparisons against the earlier mode-2 or Toeplitz
times would cross series; they are not reported as new matched percentages.

The final header is
`f9d4aa8e934533ef2dd0c0c0ab39d5d1f3b64a3d14697e0dabcd1ab53856f962`.
The driver remains `2223df3595805e54fe548ce9beac96915f758b639794387736467facb1fb837e`.
The final `gf16-inline` binary is
`239278430a28e0a1f6b5e563c64cb0eb7e8f9a03e3de34dfb06f5912f18d28d5`.
It uses the same mode-5 build flags and benchmark arguments. Both initial
source snapshots and all initial binaries remain in the remote directory.

Additional remote log SHA256 values:

- `phase-control.log`: `8f964ecf67fd55b7acb935cfc6d45427237b2e49802a99ac6f3ad401c41e4040`.
- `asm-summary.log`: `0a1ac3f83988b1e62a8f7b1cc282da5408d52162d61a239e9904ac20f2eaea2f`.
- `inline-asm-summary.log`: `5454416a5582e39ad4977bd9b1d92f4d048158e651135b204fe9f76ced76956f`.
- `inline-build-check.log`: `dedb16b29e97cf168eba82403d5389d4dc5adcd6908debdc2b5da170cc4e6dd6`.
- `inline-matched-three.log`: `2140f98a908b298388db373e9a8e51f588ea65af9d335ed5969859652d0b3816`.

The next priority is the changed RS[16,8] distance proof, especially dense
cohorts and source weights beyond the sparse r17 screen. Keep the initial
negative and the exact-map implementation win distinct. No further native
sweep or production promotion is justified by the sparse screen alone.
