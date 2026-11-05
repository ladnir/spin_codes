# Native two-sum outer controls

The default-off two-sum outer retains a modest speed advantage over the
certified Toeplitz encoder, but is slower than the unproved XOR2 candidate.
The affine variant has the simpler sparse-support law. Neither new variant
has a whole-code 35% distance certificate with 40 failure bits.

On 2026-10-09, three matched serial runs gave these median complete times:

| Comparison | Certified Toeplitz | XOR2, mode 2 | New candidate | Candidate versus Toeplitz |
|---|---:|---:|---:|---:|
| Mode 3 | 3.662937 ms | 3.465099 ms | 3.536632 ms | -3.45% |
| Mode 4 | 3.667536 ms | 3.471000 ms | 3.573931 ms | -2.55% |

Mode 3 is 2.06% slower than its matched XOR2 control. Mode 4 is 2.97%
slower than its matched XOR2 control. These are two separate matched
comparisons, not a phase-time decomposition or a parameter sweep.

## Construction and implementation

Both candidates retain eight RS[32,16] constituents per 128-input group.
For each output column, let x contain the eight RS values after the
existing nonzero input diagonal. The new fixed matrix is

    H = I + 11^T + uu^T,

over AES GF256. `HDGROUP_XOR_MIXER=3` selects u=(0,...,7), and mode 4
selects u=(8,...,15). The existing nonzero output diagonal follows H.
Both use the unchanged gauge-normalized D1/D2 sampler and packed keys.
The paired route, independent fair swaps, h4 inner, interpolation refresh,
output layout, and all payload dimensions remain unchanged.

The explicit mode-3 circuit computes S=sum(x_i), T=sum(i*x_i), and
y_i=x_i+S+i*T. Binary subset sums compute T with two fixed field products;
two further products generate 2T and 4T. The complete circuit has four
fixed GFNI products and 28 binary XORs per eight vectors. Mode 4 computes
T=T_low+8S and common=S+8T, then emits x_i+common+i*T. It adds two fixed
products and two XORs. These are source-circuit counts; GCC may combine
binary XORs into ternary-logic instructions.

The hot path uses explicit SIMD expressions and compile-time selection.
There is no runtime matrix dispatch or additional key table. Modes 3/4
require group width 8. The default remains mode 0, and mode 2 is unchanged.
The final source rebuilds the pre-edit Toeplitz and mode-2 binaries exactly.
It also rebuilds the initial mode-3 binary exactly after adding mode 4.

The independent algebra and support-law analysis is in
[HD_TWO_SUM_OUTER.md](HD_TWO_SUM_OUTER.md). Invertibility and a better
column branch number do not transfer the old whole-code receipts. The
stateful inner bound must use the new joint support law.

## Checks and timing protocol

All four builds pass `--check` and `--large-check`. The new mixer fixture
checks all 4096 input bit directions across eight ZMM registers, 16 random
columns, and H squared against the identity. A separate scalar oracle
evaluates every matrix coefficient as delta(i,j)+1+u_i*u_j.

The existing checks exercise 30720 packet input/state directions, all pair
swap masks, wrapping pairs, K=8192 and 24576, K20, rank-2/rank-4 controls,
output and scratch offsets, padding, guards, cached output, full aliasing,
and unaligned input. Partial-overlap checks belong to the disabled coset
pipeline experiment and are not claimed here. No new sanitizer run was
performed in this experiment.

Peach process checks found no active benchmark before either comparison.
Every timed process ran serially on CPU 2, holding both
`/tmp/hypercat-global-benchmark.lock` and `/tmp/spin-char2-benchmark.lock`.
Each run used 501 calls after the driver's three warmups. Setup, allocation,
and input generation were outside the timed calls. No tuning followed.

Raw run results below are in microseconds. The first series ran from
12:54:18 UTC; the affine series ran from 12:56:57 UTC.

| Series | Round | Build | Median | p10 | p90 |
|---|---:|---|---:|---:|---:|
| Mode 3 | 1 | Toeplitz | 3663.138 | 3636.037 | 3690.630 |
| Mode 3 | 1 | XOR2 | 3462.113 | 3448.988 | 3481.419 |
| Mode 3 | 1 | Mode 3 | 3536.632 | 3518.217 | 3557.761 |
| Mode 3 | 2 | Toeplitz | 3648.140 | 3628.243 | 3678.917 |
| Mode 3 | 2 | XOR2 | 3465.099 | 3445.031 | 3486.278 |
| Mode 3 | 2 | Mode 3 | 3526.012 | 3510.894 | 3548.123 |
| Mode 3 | 3 | Toeplitz | 3662.937 | 3637.240 | 3690.038 |
| Mode 3 | 3 | XOR2 | 3470.509 | 3457.224 | 3493.571 |
| Mode 3 | 3 | Mode 3 | 3554.676 | 3534.588 | 3579.793 |
| Mode 4 | 1 | Toeplitz | 3670.511 | 3647.839 | 3692.723 |
| Mode 4 | 1 | XOR2 | 3447.556 | 3431.235 | 3471.380 |
| Mode 4 | 1 | Mode 4 | 3573.931 | 3556.959 | 3596.644 |
| Mode 4 | 2 | Toeplitz | 3664.079 | 3636.448 | 3688.625 |
| Mode 4 | 2 | XOR2 | 3481.699 | 3463.896 | 3499.903 |
| Mode 4 | 2 | Mode 4 | 3571.737 | 3555.727 | 3593.568 |
| Mode 4 | 3 | Toeplitz | 3667.536 | 3643.611 | 3695.618 |
| Mode 4 | 3 | XOR2 | 3471.000 | 3452.455 | 3495.074 |
| Mode 4 | 3 | Mode 4 | 3574.783 | 3556.539 | 3594.630 |

## Reproduction and pins

The isolated remote directory is `/tmp/hd-two-sum-20261009-2DAzj6`.
Its initial sources were copied from the previous coset experiment;
only `HdGroupedOuter.h` and `hd_grouped_forward.cpp` were replaced.
Use GCC 15.2 and the existing `libspin_char2_field.a` with these flags:

    -std=gnu++20 -O3 -DNDEBUG -mavx2 -mavx512f -mavx512bw
    -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto
    -DSPIN_CHAR2_CLMUL_REDUCTION=0 -DSPIN_CHAR2_SCALAR_LOOKUP=0
    -DHDGROUP_SUBFIELD=1 -DHDGROUP_FUSE_LAST=1
    -DHDGROUP_BYTE_KEYS=1 -DHD_SUMRANK_BYTE_KEYS=1
    -DHDGROUP_INNER=5 -DHDGROUP_FUSE_REFRESH=1
    -DHDGROUP_NATIVE_OUTPUT=1 -DHDGROUP_PACKET_ROUTE=1
    -DHDGROUP_ROUTE_WIDTH=2 -DHDGROUP_SHIFTED_PAIRS=2
    -DHDGROUP_PAIR_SWAPS=1

Add `-DHDGROUP_TOEPLITZ=1` for the certified control. Otherwise add exactly
one of `-DHDGROUP_XOR_MIXER=2`, `=3`, or `=4`. Leave the refresh and coset
pipeline experiments disabled. Compile `hd_grouped_forward.cpp` against
the library, then run the checks before timing:

    taskset -c 2 ./toeplitz --bench 501 full 20 nt 4 4
    taskset -c 2 ./xor2 --bench 501 full 20 nt 8 4
    taskset -c 2 ./two-sum --bench 501 full 20 nt 8 4
    taskset -c 2 ./affine-two-sum --bench 501 full 20 nt 8 4

Run only one benchmark at a time. The two series above each used the
control, XOR2, and one candidate, in that order, for three rounds.

| File/build | SHA256 | K20 output checksum |
|---|---|---:|
| Final `HdGroupedOuter.h` | `89311e8a92efc047c02eb13a5bb2097db24bf0e4be117b9d64eda7ed371ea9fe` | — |
| Final `hd_grouped_forward.cpp` | `1cefa022b2414e788e798b70ce21d45be4d1b5a600848e800348e0565bb4c5b7` | — |
| Unchanged `HdGroupedInner.h` | `4550e45454cdfb03b570384ba9a540e4cb2d39d27aa6340b10b9ac79ea5911d1` | — |
| Linked library | `3d3c7654037e31e8cea8e5a4e59b32a6b53a663ee89ac650f591030db4d02a38` | — |
| Toeplitz | `0b31e23d66948665135485da064868ad769d007bb5cc8a2c4c4570698f1638bc` | 7701295210980077712 |
| Mode 2 | `2c0d38b79383f54b343ec9f40997bdee0e590314d95fee10774c4648a819a595` | 2343136256414281289 |
| Mode 3 | `700a22c2f560667e21dbbe21cc22c7af9883bd44efd306300002e289cf506222` | 7504277393477830387 |
| Mode 4 | `c1e5bb9804e6b3b7f0462c9633af62d9ef8c17ae5cacbbf56b842dccf2c5a01c` | 1506686965355192349 |

The immediate pre-edit mode-2 source was rebuilt independently and matched
`2c0d38...`. The older `/tmp/hd-xor2-swaps-20261009-inAmNF` binary has a
different hash; it predates intervening source changes and was not the
matched control here. The initial mode-3 source snapshot remains in the
new directory as `HdGroupedOuter.mode3.h` and `hd_grouped_forward.mode3.cpp`.

Remote log SHA256 values:

- `matched-three.log`: `5259bbc6d1c93f7360a3e243ae077f130def91c7369702b938c1972407a70c1c`.
- `affine-matched-three.log`: `6e0be3b184d37a7651f21c7ed7e5197f356257449b56bb93e2b3b29aeb980cc1`.
- `affine-build-check.log`: `714a82de36181639c5bd91c5996ce7124f55c1fa49404f18bafefa55d981836c`.

Next, test whether mode 4's label-independent one-/two-row support law
permits broad cohort bounds without the XOR2 pair-multiplicity strata.
Retain the dense-cohort obligations and the whole-code failure budget.
Do not promote either candidate or begin another native sweep before that
proof question is resolved.
