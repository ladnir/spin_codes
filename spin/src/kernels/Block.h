#pragma once
#include <cstdint>
#include <cstring>
#ifdef _MSC_VER
#include <cstdlib>
#endif
#include "../Cpu.h"
#if SPIN_X86
#include <immintrin.h>
#endif
#ifdef _MSC_VER
#define SPIN_FORCEINLINE __forceinline
#define SPIN_NOINLINE __declspec(noinline)
#define SPIN_MAY_ALIAS
#elif defined(__GNUC__) || defined(__clang__)
#define SPIN_FORCEINLINE inline __attribute__((always_inline))
#define SPIN_NOINLINE __attribute__((noinline))
#define SPIN_MAY_ALIAS __attribute__((may_alias))
#else
#define SPIN_FORCEINLINE inline
#define SPIN_NOINLINE
#define SPIN_MAY_ALIAS
#endif
#if defined(_MSC_VER) || defined(__GNUC__) || defined(__clang__)
#define SPIN_RESTRICT __restrict
#else
#define SPIN_RESTRICT
#endif
#if defined(__GNUC__) || defined(__clang__)
#define SPIN_UNROLL8 _Pragma("GCC unroll 8")
#else
#define SPIN_UNROLL8
#endif
// Prefetch hint (no-op where unavailable and in SPIN_PORTABLE builds).
#if SPIN_PORTABLE_BUILD
#define SPIN_PREFETCH(p,write) ((void)0)
#elif defined(__GNUC__) || defined(__clang__)
#define SPIN_PREFETCH(p,write) __builtin_prefetch((p),(write),3)
#elif defined(_MSC_VER) && SPIN_X86
#define SPIN_PREFETCH(p,write) _mm_prefetch(reinterpret_cast<const char*>(p),_MM_HINT_T0)
#else
#define SPIN_PREFETCH(p,write) ((void)0)
#endif
namespace spin::detail {
inline std::uint64_t byteswap64(std::uint64_t x) noexcept {
#if defined(_MSC_VER) && !defined(__clang__)
    return _byteswap_uint64(x);
#elif defined(__GNUC__) || defined(__clang__)
    return __builtin_bswap64(x);
#else
    x=((x&0x00ff00ff00ff00ffULL)<<8)|((x>>8)&0x00ff00ff00ff00ffULL);
    x=((x&0x0000ffff0000ffffULL)<<16)|((x>>16)&0x0000ffff0000ffffULL);
    return (x<<32)|(x>>32);
#endif
}
}
namespace spin::detail::storage {
// The common 128-bit record: one bit per OT instance, little-endian bytes.
// Portable code uses only the members below (construction from two 64-bit
// words, word access, XOR, zero), which every architecture provides. On
// x86-64 the representation is __m128i (SSE2 is baseline) and ISA kernels
// may use mData directly; elsewhere it is two 64-bit words.
// GCC/Clang may_alias plus private -fno-strict-aliasing permits the
// documented byte-view boundary without type-punning in callers.
struct alignas(16) SPIN_MAY_ALIAS block {
#if SPIN_X86
    __m128i mData;
    block()=default;
    block(__m128i v):mData(v) {}
    block(std::uint64_t hi,std::uint64_t lo):mData(_mm_set_epi64x(std::int64_t(hi),std::int64_t(lo))) {}
    block& operator^=(block b) {mData=_mm_xor_si128(mData,b.mData);return *this;}
    friend block operator^(block a,block b) {return block(_mm_xor_si128(a.mData,b.mData));}
    static block zero() {return block(_mm_setzero_si128());}
#else
    std::uint64_t mWords[2];
    block()=default;
    block(std::uint64_t hi,std::uint64_t lo):mWords{lo,hi} {}
    block& operator^=(block b) {mWords[0]^=b.mWords[0];mWords[1]^=b.mWords[1];return *this;}
    friend block operator^(block a,block b) {a^=b;return a;}
    static block zero() {return block(0,0);}
#endif
    std::uint64_t lo() const {std::uint64_t w;std::memcpy(&w,this,8);return w;}
    std::uint64_t hi() const {std::uint64_t w;std::memcpy(&w,reinterpret_cast<const char*>(this)+8,8);return w;}
};
static_assert(sizeof(block)==16 && alignof(block)==16);
}
