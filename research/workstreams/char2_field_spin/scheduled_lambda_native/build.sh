#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"

# Local scratch build only. The caller supplies the same support headers and
# static library and unchanged HdGroupedInner.h from bit_syndrome_native;
# this script runs no checks or benchmarks. Both variants use the same route,
# keys and RS transform. The candidate changes only the outer column map.
common=(-std=gnu++20 -O3 -DNDEBUG -mavx2 -mavx512f -mavx512bw -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto
  -DSPIN_CHAR2_CLMUL_REDUCTION=0 -DSPIN_CHAR2_SCALAR_LOOKUP=0
  -DHDGROUP_SUBFIELD=1 -DHDGROUP_FUSE_LAST=1 -DHDGROUP_BYTE_KEYS=1 -DHD_SUMRANK_BYTE_KEYS=1
  -DHDGROUP_TOEPLITZ=0 -DHDGROUP_XOR_MIXER=5 -DHDGROUP_BIT_SYNDROMES=1
  -DHDGROUP_INNER=5 -DHDGROUP_FUSE_REFRESH=1 -DHDGROUP_NATIVE_OUTPUT=1
  -DHDGROUP_PACKET_ROUTE=1 -DHDGROUP_ROUTE_WIDTH=2 -DHDGROUP_SHIFTED_PAIRS=2 -DHDGROUP_PAIR_SWAPS=1
  -DHDGROUP_DECODED_ROUTE=1 -DHDGROUP_ROUTE_ADDRESS_SHIFT=0 -DHDGROUP_SOURCE_PREFIX=0
  -DHDGROUP_CROSS_LAYERS=0)

if [[ "${1:-release}" == sanitize ]]; then
  g++-15 "${common[@]}" -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -DHDGROUP_SCHEDULED_LAMBDA=1 hd_grouped_forward.cpp libspin_char2_field.a -o sanitize
elif [[ "${1:-release}" == natural ]]; then
  g++-15 "${common[@]}" -UHDGROUP_NATIVE_OUTPUT -DHDGROUP_NATIVE_OUTPUT=0 -DHDGROUP_SCHEDULED_LAMBDA=1 hd_grouped_forward.cpp libspin_char2_field.a -o lambda-natural
else
  g++-15 "${common[@]}" -DHDGROUP_SCHEDULED_LAMBDA=0 hd_grouped_forward.cpp libspin_char2_field.a -o bit
  g++-15 "${common[@]}" -DHDGROUP_SCHEDULED_LAMBDA=1 hd_grouped_forward.cpp libspin_char2_field.a -o lambda
fi
