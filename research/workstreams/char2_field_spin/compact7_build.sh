#!/usr/bin/env bash
set -euo pipefail
cd /tmp/hd-compact7-20261009-uot619
common=(-std=gnu++20 -O3 -DNDEBUG -mavx2 -mavx512f -mavx512bw -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto
 -DSPIN_CHAR2_CLMUL_REDUCTION=0 -DSPIN_CHAR2_SCALAR_LOOKUP=0
 -DHDGROUP_SUBFIELD=1 -DHDGROUP_FUSE_LAST=1 -DHDGROUP_BYTE_KEYS=1 -DHD_SUMRANK_BYTE_KEYS=1
 -DHDGROUP_INNER=5 -DHDGROUP_FUSE_REFRESH=1 -DHDGROUP_NATIVE_OUTPUT=1
 -DHDGROUP_PACKET_ROUTE=1 -DHDGROUP_ROUTE_WIDTH=2 -DHDGROUP_SHIFTED_PAIRS=2 -DHDGROUP_PAIR_SWAPS=1
 -DHDGROUP_TOEPLITZ=1)
g++-15 "${common[@]}" hd_grouped_forward.cpp libspin_char2_field.a -o control
g++-15 "${common[@]}" -DHDGROUP_TOEPLITZ_COMPACT7=1 hd_grouped_forward.cpp libspin_char2_field.a -o compact7
cmp control /tmp/hd-two-sum-g16-20261009-3VElrn/inline-toeplitz
taskset -c 2 ./compact7 --check
taskset -c 2 ./compact7 --large-check
sha256sum HdGroupedOuter.h hd_grouped_forward.cpp control compact7
