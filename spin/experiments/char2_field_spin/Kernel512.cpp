#include "Gfni.h"
#include <stdexcept>

namespace spin::research::char2 {

void multiplyGfniPrepared(const std::uint8_t* input, const Coeff4* coefficients,
                          std::uint8_t* output, std::size_t count) {
    const auto complete = count / 4;
    for (std::size_t i = 0; i != complete; ++i)
        _mm512_storeu_si512(output + 64 * i,
            smallmul4(_mm512_loadu_si512(input + 64 * i), coefficients[i]));
    for (std::size_t i = 4 * complete; i != count; ++i) {
        const auto& c = coefficients[i / 4];
        _mm_storeu_si128(reinterpret_cast<__m128i*>(output + 16 * i),
            smallmul128(_mm_loadu_si128(reinterpret_cast<const __m128i*>(input + 16 * i)),
                        Coeff128{c.low[2 * (i % 4)], c.high[2 * (i % 4)]}));
    }
}

void multiplyGfni512(const std::uint8_t* input, const std::uint8_t* coefficients,
                     std::uint8_t* output, std::size_t count, unsigned bits) {
    if (!validBits(bits)) throw std::invalid_argument("unsupported field size");
    if (count == 0) return;
    std::size_t i = 0;
    if (bits == 128) {
        for (; i + 4 <= count; i += 4)
            _mm512_storeu_si512(output + 16 * i,
                smallmul4(_mm512_loadu_si512(input + 16 * i), coefficients + i));
    } else if (bits == 64 || bits == 32) {
        for (; i + 8 <= count; i += 8) {
            const auto packed = _mm_loadl_epi64(reinterpret_cast<const __m128i*>(coefficients + i));
            const auto indices = _mm256_cvtepu8_epi32(packed);
            const auto low = _mm512_i32gather_epi64(indices, affineTable.low.data(), 8);
            const auto high = _mm512_i32gather_epi64(indices, affineTable.high.data(), 8);
            if (bits == 64)
                _mm512_storeu_si512(output + 8 * i,
                    smallmul8x64(_mm512_loadu_si512(input + 8 * i), low, high));
            else {
                const auto values = _mm512_maskz_expandloadu_epi32(0x5555, input + 4 * i);
                _mm512_mask_compressstoreu_epi32(output + 4 * i, 0x5555,
                    smallmul8x32Expanded(values, low, high));
            }
        }
    } else if (bits == 8) {
        for (; i + 64 <= count; i += 64)
            _mm512_storeu_si512(output + i, _mm512_gf2p8mul_epi8(
                _mm512_loadu_si512(input + i), _mm512_loadu_si512(coefficients + i)));
    } else throw std::invalid_argument("unsupported field size");
    multiplyPortable(input + i * (bits / 8), coefficients + i,
                     output + i * (bits / 8), count - i, bits);
}

} // namespace spin::research::char2
