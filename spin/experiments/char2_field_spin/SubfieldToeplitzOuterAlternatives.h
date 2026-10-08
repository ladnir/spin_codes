#pragma once
#include "Subfield16.h"
#include <type_traits>
#include <utility>

#if defined(_MSC_VER)
#define SPIN_TOEPLITZ_ALT_INLINE __forceinline
#else
#define SPIN_TOEPLITZ_ALT_INLINE inline __attribute__((always_inline))
#endif

namespace spin::research::char2 {

// Exact-map alternatives using the original 243-leaf coefficient layout.
// Packed Karatsuba retains all products and XORs, but broadcasts each group
// of three leaf coefficients once. Dense2 instead recovers the 2x2 matrix
// from those leaves and trades one data XOR for one additional GFNI product.
template<unsigned Size>
inline constexpr unsigned sfToeplitzAlternativeLeaves = 3 * sfToeplitzAlternativeLeaves<Size / 2>;
template<>
inline constexpr unsigned sfToeplitzAlternativeLeaves<1> = 1;

template<class Function, std::size_t... Index>
SPIN_TOEPLITZ_ALT_INLINE void sfToeplitzAlternativeUnroll(std::index_sequence<Index...>, Function&& function) {
    (function(std::integral_constant<std::size_t, Index>{}), ...);
}

SPIN_TOEPLITZ_ALT_INLINE __m512i sfToeplitzBroadcastPackedFour(std::uint32_t coefficients) {
    const auto packet = _mm512_broadcast_i32x4(_mm_cvtsi32_si128(static_cast<int>(coefficients)));
    const auto selectors = _mm512_setr_epi64(0, 0, 0x0101010101010101LL, 0x0101010101010101LL,
        0x0202020202020202LL, 0x0202020202020202LL, 0x0303030303030303LL, 0x0303030303030303LL);
    return _mm512_shuffle_epi8(packet, selectors);
}

template<unsigned Size, bool Dense2 = false, unsigned Offset = 0,
         unsigned TotalLeaves = sfToeplitzAlternativeLeaves<Size>>
SPIN_TOEPLITZ_ALT_INLINE void subfieldToeplitzApplyAlternative4(const std::uint8_t* leaves,
                                                               const __m512i* input, __m512i* output) {
    if constexpr (Size == 1) {
        output[0] = _mm512_gf2p8mul_epi8(input[0], sfBroadcastFour(leaves + 4 * Offset));
    } else if constexpr (Size == 2) {
        if constexpr (Dense2) {
            std::uint32_t a, b, c;
            std::memcpy(&a, leaves + 4 * Offset, 4);
            std::memcpy(&b, leaves + 4 * (Offset + 1), 4);
            std::memcpy(&c, leaves + 4 * (Offset + 2), 4);
            const auto diagonal = sfToeplitzBroadcastPackedFour(a);
            const auto lower = sfToeplitzBroadcastPackedFour(b ^ a);
            const auto upper = sfToeplitzBroadcastPackedFour(c ^ a);
            output[0] = _mm512_xor_si512(_mm512_gf2p8mul_epi8(input[0], diagonal),
                                         _mm512_gf2p8mul_epi8(input[1], upper));
            output[1] = _mm512_xor_si512(_mm512_gf2p8mul_epi8(input[0], lower),
                                         _mm512_gf2p8mul_epi8(input[1], diagonal));
        } else {
            __m128i compact;
            if constexpr (Offset + 3 == TotalLeaves)
                compact = _mm_maskz_loadu_epi32(0x7, leaves + 4 * Offset);
            else {
                static_assert(Offset + 4 <= TotalLeaves);
                compact = _mm_loadu_si128(reinterpret_cast<const __m128i*>(leaves + 4 * Offset));
            }
            // The final node reads precisely twelve bytes. Other nodes may
            // read the next node's first coefficient, but never select it.
            const auto packet = _mm512_broadcast_i32x4(compact);
            const auto a = _mm512_gf2p8mul_epi8(_mm512_xor_si512(input[0], input[1]), sf16ExpandRow<0>(packet));
            const auto b = _mm512_gf2p8mul_epi8(input[0], sf16ExpandRow<1>(packet));
            const auto c = _mm512_gf2p8mul_epi8(input[1], sf16ExpandRow<2>(packet));
            output[0] = _mm512_xor_si512(a, c);
            output[1] = _mm512_xor_si512(a, b);
        }
    } else {
        constexpr unsigned half = Size / 2, count = sfToeplitzAlternativeLeaves<half>;
        __m512i temporary[half];
        sfToeplitzAlternativeUnroll(std::make_index_sequence<half>{}, [&](auto i) {
            temporary[i] = _mm512_xor_si512(input[i], input[i + half]);
        });
        subfieldToeplitzApplyAlternative4<half, Dense2, Offset, TotalLeaves>(leaves, temporary, output);
        subfieldToeplitzApplyAlternative4<half, Dense2, Offset + count, TotalLeaves>(leaves, input, output + half);
        subfieldToeplitzApplyAlternative4<half, Dense2, Offset + 2 * count, TotalLeaves>(leaves, input + half, temporary);
        sfToeplitzAlternativeUnroll(std::make_index_sequence<half>{}, [&](auto i) {
            output[i + half] = _mm512_xor_si512(output[i], output[i + half]);
            output[i] = _mm512_xor_si512(output[i], temporary[i]);
        });
    }
}

} // namespace spin::research::char2

#undef SPIN_TOEPLITZ_ALT_INLINE
