#include "Gfni.h"

namespace spin::research::char2 {
void multiplyGfni128(const std::uint8_t* input, const std::uint8_t* coefficients,
                     std::uint8_t* output, std::size_t count) {
    for (std::size_t i = 0; i != count; ++i)
        _mm_storeu_si128(reinterpret_cast<__m128i*>(output + 16 * i),
            smallmul128(_mm_loadu_si128(reinterpret_cast<const __m128i*>(input + 16 * i)),
                        coefficients[i]));
}

void chainGfni128(std::uint8_t* state, const std::uint8_t* coefficients,
                  std::size_t count) {
    auto value = _mm_loadu_si128(reinterpret_cast<const __m128i*>(state));
    for (std::size_t i = 0; i != count; ++i) value = smallmul128(value, coefficients[i]);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(state), value);
}

void chainGfni128Prepared(std::uint8_t* state, const Coeff128* coefficients,
                          std::size_t count) {
    auto value = _mm_loadu_si128(reinterpret_cast<const __m128i*>(state));
    for (std::size_t i = 0; i != count; ++i) value = smallmul128(value, coefficients[i]);
    _mm_storeu_si128(reinterpret_cast<__m128i*>(state), value);
}
} // namespace spin::research::char2
