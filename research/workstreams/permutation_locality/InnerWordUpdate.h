#pragma once

#include "Inner.h"
#include <cstdint>
#include <immintrin.h>

namespace spin::research::packet_inner {
namespace word_update_detail {

template<unsigned Row>
static SPIN_FORCEINLINE void updateRow(__m128i* state,
    const std::uint16_t* rows, const __m128i* feedback,
    const __m128i table[][16]) {
    const unsigned mask = rows[Row];
    const auto first = _mm_ternarylogic_epi64(
        table[0][mask & 15], table[1][(mask >> 4) & 15],
        table[2][(mask >> 8) & 15], 0x96);
    state[Row] = _mm_ternarylogic_epi64(
        first, table[3][mask >> 12], feedback[Row], 0x96);
}

} // namespace word_update_detail

// state'[r] = XOR { state[c] : bit c of rows[r] is set } XOR feedback[r].
// Rows are the actual desired map, already transposed/conjugated by setup
// when required. This function does not change matrix or state coordinates.
// State, rows, and feedback must be disjoint, with ordinary element alignment.
// Requires AVX-512F+VL. The 1-KiB table uses only 16-byte alignment; invoke
// from an ordinary function, not a coroutine with scratch in its frame.
static SPIN_FORCEINLINE void wordUpdate(__m128i* __restrict state,
    const std::uint16_t* __restrict rows,
    const __m128i* __restrict feedback) {
    __m128i table[4][16];
    // Capture every old state word before the first state overwrite. Each
    // nibble table needs 11 XORs; the four tables need 44 XORs in total.
    spin::detail::kernel::tables<16>(state, table);
    using word_update_detail::updateRow;
    updateRow<0>(state, rows, feedback, table);
    updateRow<1>(state, rows, feedback, table);
    updateRow<2>(state, rows, feedback, table);
    updateRow<3>(state, rows, feedback, table);
    updateRow<4>(state, rows, feedback, table);
    updateRow<5>(state, rows, feedback, table);
    updateRow<6>(state, rows, feedback, table);
    updateRow<7>(state, rows, feedback, table);
    updateRow<8>(state, rows, feedback, table);
    updateRow<9>(state, rows, feedback, table);
    updateRow<10>(state, rows, feedback, table);
    updateRow<11>(state, rows, feedback, table);
    updateRow<12>(state, rows, feedback, table);
    updateRow<13>(state, rows, feedback, table);
    updateRow<14>(state, rows, feedback, table);
    updateRow<15>(state, rows, feedback, table);
}

// Optional setup-only validation. No assertions are disabled by NDEBUG.
// First check every 16-bit row mask, then all 4096 state/feedback input-bit
// bases for a fixed dense matrix. The latter checks both 64-bit halves and
// every bit position of each opaque 128-bit element.
[[nodiscard]] inline bool wordUpdateSelfCheck() {
    __m128i state[16], feedback[16];
    std::uint16_t rows[16];
    std::uint64_t words[2], actual[2];
    for (unsigned base = 0; base < 65536; base += 16) {
        for (unsigned row = 0; row < 16; ++row) {
            rows[row] = std::uint16_t(base + row);
            words[0] = std::uint64_t(1) << row;
            words[1] = std::uint64_t(1) << (32 + row);
            state[row] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(words));
            words[0] = 0x0123456789abcdefULL ^ row;
            words[1] = 0xfedcba9876543210ULL ^ row;
            feedback[row] = _mm_loadu_si128(reinterpret_cast<const __m128i*>(words));
        }
        wordUpdate(state, rows, feedback);
        for (unsigned row = 0; row < 16; ++row) {
            _mm_storeu_si128(reinterpret_cast<__m128i*>(actual), state[row]);
            if (actual[0] != (std::uint64_t(rows[row]) ^ 0x0123456789abcdefULL ^ row) ||
                actual[1] != ((std::uint64_t(rows[row]) << 32) ^ 0xfedcba9876543210ULL ^ row))
                return false;
        }
    }

    constexpr std::uint16_t dense[16] = {
        0x0000, 0xffff, 0xa5a5, 0x5a5a, 0x1234, 0x8001, 0x8421, 0x1248,
        0x1111, 0x2222, 0x4444, 0x8888, 0x0f0f, 0xf0f0, 0x00ff, 0xff00
    };
    for (unsigned input = 0; input < 32; ++input) {
        for (unsigned bit = 0; bit < 128; ++bit) {
            for (unsigned row = 0; row < 16; ++row) {
                state[row] = _mm_setzero_si128();
                feedback[row] = _mm_setzero_si128();
            }
            const auto value = std::uint64_t(1) << (bit & 63);
            words[0] = bit < 64 ? value : 0;
            words[1] = bit >= 64 ? value : 0;
            const auto basis = _mm_loadu_si128(reinterpret_cast<const __m128i*>(words));
            if (input < 16) state[input] = basis;
            else feedback[input - 16] = basis;
            wordUpdate(state, dense, feedback);
            for (unsigned row = 0; row < 16; ++row) {
                const bool selected = input < 16
                    ? ((dense[row] >> input) & 1) != 0 : row == input - 16;
                const auto expected = selected ? value : 0;
                _mm_storeu_si128(reinterpret_cast<__m128i*>(actual), state[row]);
                if (actual[0] != (bit < 64 ? expected : 0) ||
                    actual[1] != (bit >= 64 ? expected : 0)) return false;
            }
        }
    }
    return true;
}

} // namespace spin::research::packet_inner
