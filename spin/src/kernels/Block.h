#pragma once
#include <immintrin.h>
#include <cstdint>
#include "../Cpu.h"
#ifdef _MSC_VER
#define SPIN_FORCEINLINE __forceinline
#define SPIN_NOINLINE __declspec(noinline)
#define SPIN_MAY_ALIAS
#else
#define SPIN_FORCEINLINE inline __attribute__((always_inline))
#define SPIN_NOINLINE __attribute__((noinline))
#define SPIN_MAY_ALIAS __attribute__((may_alias))
#endif
namespace spin::detail::storage {
// Private SIMD storage. GCC/Clang may_alias plus private -fno-strict-aliasing
// permits the documented byte-view boundary without type-punning in callers.
struct alignas(16) SPIN_MAY_ALIAS block {
    __m128i mData;
    block()=default;
    block(__m128i v):mData(v) {}
    block(std::uint64_t hi,std::uint64_t lo):mData(_mm_set_epi64x(hi,lo)) {}
    block& operator^=(block b) {mData=_mm_xor_si128(mData,b.mData);return *this;}
};
static_assert(sizeof(block)==16 && alignof(block)==16);
}
