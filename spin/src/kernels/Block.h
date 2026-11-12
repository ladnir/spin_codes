#pragma once
#include <cstdint>
#include <cstring>
#if defined(__x86_64__) || defined(_M_X64)
#define SPIN_X86 1
#include <immintrin.h>
#else
#define SPIN_X86 0
#endif
#if defined(__aarch64__) || defined(_M_ARM64)
#define SPIN_AARCH64 1
#else
#define SPIN_AARCH64 0
#endif
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
