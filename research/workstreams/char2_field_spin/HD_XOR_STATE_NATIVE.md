# Native diagonal-XOR-diagonal refresh experiment

The cheaper state refresh saves about 20% in an isolated inner measurement,
but only about 1--5% in complete encoding. It is a research candidate, not
a certified replacement for the current inner. The complete result does
not justify a large proof search solely on the operation count.

## Change and validation

Keep the 64-symbol packet, 16-symbol state, four RS[16,4] emission blocks,
feedback transpose, route, and independently sampled nonzero diagonal keys.
Replace the fixed interpolation map in the refresh by

    Q_next = D2 (I + 11^T) D1 (Q + CX).

On four striped ZMM registers, the new fixed map uses nine XORs and two
lane shuffles instead of 20 fixed GFNI products. Both random diagonals remain.
The default-off `HDGROUP_XOR_REFRESH=1` flag is confined to the research
header and driver. The existing h2 control and production code are unchanged.
The [support analysis](HD_XOR_STATE_REFRESH.md) describes the new setup law;
the old distance certificate does not transfer automatically.

Validation covers 30,720 packet basis directions, all pair-swap masks and
wrapping loads, small and natural lengths, alias/alignment cases, and K20
comparisons against an independent dense oracle. Both the Toeplitz-outer
and XOR2-outer builds pass. An ASAN/UBSAN build of the combined candidate
also passes `--check`; its linked static library is not instrumented.

## Matched complete encoding

All runs use K=2^20 input records of 128 bits, rate 1/2, h4, independent
pair swaps, native output layout, and streaming output stores. Setup and
allocation are outside the timed interval. Each entry below is the median
of 501 complete calls. Three matched rounds reverse the order in round two.
Runs are serial on Peach CPU 2 under the shared benchmark locks.

| Outer and refresh | Three run medians, ms | Median, ms |
|---|---|---:|
| Toeplitz, interpolation control | 3.773865, 3.655404, 3.764006 | 3.764006 |
| Toeplitz, XOR refresh | 3.699576, 3.570645, 3.566958 | 3.570645 |
| XOR2, interpolation control | 3.538776, 3.485596, 3.529789 | 3.529789 |
| XOR2, XOR refresh | 3.486348, 3.515452, 3.393865 | 3.486348 |

The median reductions are 5.14% with the certified Toeplitz outer and 1.23%
with the candidate XOR2 outer. One matched XOR2 round regresses. Shared-host
variation prevents treating the combined result as a robust large gain.
The earlier 3.456162-ms XOR2 result used the unchanged interpolation refresh;
do not compare its best older run to a new control.

Separate warmed phase measurements confirm that the intended work disappeared:

| Phase | Interpolation, ms | XOR refresh, ms |
|---|---:|---:|
| Toeplitz routed inner, first pair | 1.735379 | 1.402418 |
| Toeplitz routed inner, second pair | 1.722255 | 1.403450 |
| XOR2 routed inner | 1.748775 | 1.398952 |
| XOR2 outer | 1.680077 | 1.675709 |

These phases have different cache histories. They do not add to complete
encoding time. The gap motivates investigation of the outer-to-inner data
flow; it does not by itself identify a bottleneck.

## Reproduction and pins

Use GCC 15.2, `-std=gnu++20 -O3 -DNDEBUG -mtune=znver4 -fno-lto`, AVX2,
AVX512F/BW/VL/DQ, and GFNI. Common definitions are:

    -DSPIN_CHAR2_CLMUL_REDUCTION=0 -DSPIN_CHAR2_SCALAR_LOOKUP=0
    -DHDGROUP_SUBFIELD=1 -DHDGROUP_FUSE_LAST=1
    -DHDGROUP_BYTE_KEYS=1 -DHD_SUMRANK_BYTE_KEYS=1
    -DHDGROUP_INNER=5 -DHDGROUP_FUSE_REFRESH=1
    -DHDGROUP_NATIVE_OUTPUT=1 -DHDGROUP_PACKET_ROUTE=1
    -DHDGROUP_ROUTE_WIDTH=2 -DHDGROUP_SHIFTED_PAIRS=2
    -DHDGROUP_PAIR_SWAPS=1

Select `-DHDGROUP_TOEPLITZ=1` or `-DHDGROUP_XOR_MIXER=2`. Add
`-DHDGROUP_XOR_REFRESH=1` only to the corresponding candidate. After
`--check` and `--large-check`, run `--bench 501 full 20 nt 4 4` for Toeplitz
or `--bench 501 full 20 nt 8 4` for XOR2. Replace `full` by `inner` or
`outer` only for phase diagnostics. Never run two benchmarks concurrently.

Source SHA256 values:

- `HdGroupedOuter.h`: `175dd5c03a68ee442edbb4653eebf3f95ffb09ff39bcaa08addca4a402480035`.
- `HdGroupedInner.h`: `4550e45454cdfb03b570384ba9a540e4cb2d39d27aa6340b10b9ac79ea5911d1`.
- `hd_grouped_forward.cpp`: `f32660c25dde195aa69e5fd2216a6e01dd3cfc6c319d485e521cee1eaf2a9eec`.
- `libspin_char2_field.a`: `3d3c7654037e31e8cea8e5a4e59b32a6b53a663ee89ac650f591030db4d02a38`.

| Build | Binary SHA256 | K20 output checksum |
|---|---|---:|
| Toeplitz control | `0b31e23d66948665135485da064868ad769d007bb5cc8a2c4c4570698f1638bc` | 7701295210980077712 |
| Toeplitz XOR refresh | `54cf184e6760f52e703852c36d5e73e63ed835cb771d36b1c73a9ceb11c27d7e` | 1892448702060305324 |
| XOR2 control | `2c0d38b79383f54b343ec9f40997bdee0e590314d95fee10774c4648a819a595` | 2343136256414281289 |
| XOR2 XOR refresh | `e1be241611fbf1e4e691834864300da302da92eb69ec4bb7a76fe608ce940e64` | 1565471324984973961 |

Keep this prototype and its exact local proof tools available, but prioritize
the cheaper outer's unfinished certificate and proof-preserving data-flow
changes before expanding the XOR refresh to a full regional proof.
