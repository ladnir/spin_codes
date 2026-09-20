#pragma once

#include <cstdint>
#include <immintrin.h>

#if defined(_MSC_VER)
#define OC_FORCEINLINE __forceinline
#else
#define OC_FORCEINLINE inline __attribute__((always_inline))
#endif

namespace osuCrypto
{
    using u8 = std::uint8_t;
    using u16 = std::uint16_t;
    using u32 = std::uint32_t;
    using u64 = std::uint64_t;

    struct alignas(16) block
    {
        __m128i mData;

        block() noexcept : mData(_mm_setzero_si128()) {}
        block(u64 high, u64 low) noexcept
            : mData(_mm_set_epi64x(
                static_cast<long long>(high), static_cast<long long>(low))) {}
        explicit block(__m128i value) noexcept : mData(value) {}

        block& operator^=(block other) noexcept
        {
            mData = _mm_xor_si128(mData, other.mData);
            return *this;
        }
    };

    OC_FORCEINLINE block operator^(block left, block right) noexcept
    {
        return block(_mm_xor_si128(left.mData, right.mData));
    }

    inline const block ZeroBlock{};
}
