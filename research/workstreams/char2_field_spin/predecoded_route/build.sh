#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
common=(-std=gnu++20 -O3 -DNDEBUG -mavx2 -mavx512f -mavx512bw -mavx512vl -mavx512dq -mgfni -mtune=znver4 -fno-lto
  -DSPIN_CHAR2_CLMUL_REDUCTION=0 -DSPIN_CHAR2_SCALAR_LOOKUP=0
  -DHDGROUP_SUBFIELD=1 -DHDGROUP_FUSE_LAST=1 -DHDGROUP_BYTE_KEYS=1 -DHD_SUMRANK_BYTE_KEYS=1
  -DHDGROUP_INNER=5 -DHDGROUP_FUSE_REFRESH=1 -DHDGROUP_NATIVE_OUTPUT=1
  -DHDGROUP_PACKET_ROUTE=1 -DHDGROUP_ROUTE_WIDTH=2 -DHDGROUP_SHIFTED_PAIRS=2 -DHDGROUP_PAIR_SWAPS=1
  -DHDGROUP_TOEPLITZ=1)
case "${1:-route}" in
  route)
    g++-15 "${common[@]}" hd_grouped_forward.cpp libspin_char2_field.a -o baseline
    g++-15 "${common[@]}" -DHDGROUP_DECODED_ROUTE=1 hd_grouped_forward.cpp libspin_char2_field.a -o decoded
    ;;
  source)
    g++-15 "${common[@]}" -DHDGROUP_SOURCE_PREFIX=1 hd_grouped_forward.cpp libspin_char2_field.a -o source
    g++-15 "${common[@]}" -DHDGROUP_SOURCE_PREFIX=1 -DHDGROUP_DECODED_ROUTE=1 hd_grouped_forward.cpp libspin_char2_field.a -o both
    ;;
  sanitize)
    g++-15 "${common[@]}" -O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer -DHDGROUP_DECODED_ROUTE=1 hd_grouped_forward.cpp libspin_char2_field.a -o sanitize
    ;;
  address)
    g++-15 "${common[@]}" -DHDGROUP_DECODED_ROUTE=1 -DHDGROUP_ROUTE_ADDRESS_SHIFT=3 hd_grouped_forward.cpp libspin_char2_field.a -o address
    ;;
  batched)
    g++-15 "${common[@]}" -DHDGROUP_SOURCE_PREFIX=2 hd_grouped_forward.cpp libspin_char2_field.a -o batched
    g++-15 "${common[@]}" -DHDGROUP_SOURCE_PREFIX=2 -DHDGROUP_DECODED_ROUTE=1 hd_grouped_forward.cpp libspin_char2_field.a -o batched-decoded
    ;;
esac
