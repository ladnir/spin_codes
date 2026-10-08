#include "SmallField.h"
#include <cstring>
#include <stdexcept>

#if defined(_MSC_VER)
#include <intrin.h>
#endif

namespace spin::research::char2 {

void multiplyOracle(const std::uint8_t* input, std::uint8_t coefficient,
                    std::uint8_t* output, unsigned bits) {
    if (!validBits(bits)) throw std::invalid_argument("unsupported field size");
    // Unreduced product, followed by bitwise polynomial long division. This is
    // independent of the GFNI affine decomposition and optimized shift loop.
    std::array<std::uint8_t, 17> product{};
    for (unsigned inputBit = 0; inputBit < bits; ++inputBit) {
        if (!((input[inputBit / 8] >> (inputBit % 8)) & 1U)) continue;
        for (unsigned coefficientBit = 0; coefficientBit < 8; ++coefficientBit)
            if ((coefficient >> coefficientBit) & 1U) {
                const auto bit = inputBit + coefficientBit;
                product[bit / 8] ^= std::uint8_t(1U << (bit % 8));
            }
    }
    for (unsigned position = bits + 7; position-- > bits;) {
        if (!((product[position / 8] >> (position % 8)) & 1U)) continue;
        product[position / 8] ^= std::uint8_t(1U << (position % 8));
        for (unsigned bit = 0; bit != 8; ++bit)
            if ((reduction(bits) >> bit) & 1U) {
                const auto target = position - bits + bit;
                product[target / 8] ^= std::uint8_t(1U << (target % 8));
            }
    }
    std::memcpy(output, product.data(), bits / 8);
}

template<unsigned Bytes, std::uint8_t Modulus>
void portable(const std::uint8_t* input, const std::uint8_t* coefficients,
              std::uint8_t* output, std::size_t count) {
    for (std::size_t record = 0; record != count; ++record) {
        std::array<std::uint8_t, Bytes> current{}, result{};
        std::memcpy(current.data(), input + record * Bytes, Bytes);
        const auto c = coefficients[record];
        for (unsigned bit = 0; bit != 8; ++bit) {
            const auto mask = std::uint8_t(0U - ((c >> bit) & 1U));
            for (unsigned j = 0; j != Bytes; ++j) result[j] ^= current[j] & mask;
            const auto overflow = current[Bytes - 1] >> 7;
            for (unsigned j = Bytes - 1; j != 0; --j)
                current[j] = std::uint8_t((current[j] << 1) | (current[j - 1] >> 7));
            current[0] = std::uint8_t((current[0] << 1) ^ ((0U - overflow) & Modulus));
        }
        std::memcpy(output + record * Bytes, result.data(), Bytes);
    }
}

void multiplyPortable(const std::uint8_t* input, const std::uint8_t* coefficients,
                      std::uint8_t* output, std::size_t count, unsigned bits) {
    switch (bits) {
    case 8: portable<1, 0x1b>(input, coefficients, output, count); break;
    case 32: portable<4, 0x8d>(input, coefficients, output, count); break;
    case 64: portable<8, 0x1b>(input, coefficients, output, count); break;
    case 128: portable<16, 0x87>(input, coefficients, output, count); break;
    default: throw std::invalid_argument("unsupported field size");
    }
}

#if defined(_MSC_VER) && (defined(_M_X64) || defined(_M_IX86))
static bool features(bool wide) noexcept {
    int cpu[4]{};
    __cpuid(cpu, 0);
    if (cpu[0] < 7) return false;
    __cpuidex(cpu, 1, 0);
    if ((cpu[2] & (1 << 27)) == 0) return false;
#if SPIN_CHAR2_CLMUL_REDUCTION
    if ((cpu[2] & (1 << 1)) == 0) return false;
#endif
    const auto xcr0 = _xgetbv(0);
    if ((xcr0 & 0x6) != 0x6) return false;
    __cpuidex(cpu, 7, 0);
    if ((cpu[2] & (1 << 8)) == 0 || (cpu[1] & (1 << 5)) == 0) return false;
#if SPIN_CHAR2_CLMUL_REDUCTION
    if (wide && (cpu[2] & (1 << 10)) == 0) return false;
#endif
    return !wide || ((xcr0 & 0xe6) == 0xe6 &&
        (cpu[1] & (1 << 16)) && (cpu[1] & (1 << 30)) && (cpu[1] & (1 << 31)));
}
bool gfni128Available() noexcept { return features(false); }
bool gfni512Available() noexcept { return features(true); }
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
bool gfni128Available() noexcept {
    __builtin_cpu_init();
#if SPIN_CHAR2_CLMUL_REDUCTION
    if (!__builtin_cpu_supports("pclmul")) return false;
#endif
    return __builtin_cpu_supports("gfni") && __builtin_cpu_supports("avx2");
}
bool gfni512Available() noexcept {
    __builtin_cpu_init();
#if SPIN_CHAR2_CLMUL_REDUCTION
    if (!__builtin_cpu_supports("vpclmulqdq")) return false;
#endif
    return gfni128Available() && __builtin_cpu_supports("avx512f") &&
           __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("avx512vl");
}
#else
bool gfni128Available() noexcept { return false; }
bool gfni512Available() noexcept { return false; }
#endif

} // namespace spin::research::char2
