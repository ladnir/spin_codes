#pragma once
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#if defined(__linux__) && !defined(SPIN_WORKSPACE_ROUTING_PORTABLE_TEST)
#include <cerrno>
#include <sys/mman.h>
#include <unistd.h>
#endif

#ifndef SPIN_WORKSPACE_ROUTING_OPT
#define SPIN_WORKSPACE_ROUTING_OPT 0
#endif

namespace spin::detail::kernel::workspace_routing {
// A compile-time tag selects a separate kernel instantiation. The original map
// instantiation stays unchanged for small messages and disabled builds.
template<class Base> struct TunedMap:Base {static constexpr bool tunedWorkspaceRouting=true;};
#if SPIN_WORKSPACE_ROUTING_OPT && defined(__linux__) && !defined(SPIN_WORKSPACE_ROUTING_PORTABLE_TEST)
inline constexpr bool compiled=true;
#else
inline constexpr bool compiled=false;
#endif

inline constexpr bool eligible(bool quarter,std::size_t messageBlocks) noexcept {
    return compiled && quarter && messageBlocks>=(std::size_t{1}<<18);
}

struct PageRange { std::uintptr_t begin=0;std::size_t bytes=0; };
// Never advise pages shared with another allocation, including allocator metadata.
inline constexpr PageRange wholePages(std::uintptr_t address,std::size_t bytes,std::size_t page) noexcept {
    constexpr auto max=std::numeric_limits<std::uintptr_t>::max();
    if(!page || (page&(page-1)) || bytes>max-address || address>max-(page-1)) return {};
    const auto begin=(address+page-1)&~std::uintptr_t(page-1);
    const auto end=(address+bytes)&~std::uintptr_t(page-1);
    return end>begin?PageRange{begin,end-begin}:PageRange{};
}

struct PageAdvice {std::size_t bytes=0;bool hinted=false,collapsed=false;};

// Setup only. Failure leaves the allocation usable; no logging or global changes.
inline PageAdvice adviseOwned(void* memory,std::size_t bytes) noexcept {
    PageAdvice result;
#if defined(__linux__) && !defined(SPIN_WORKSPACE_ROUTING_PORTABLE_TEST)
    if constexpr(compiled) {
        const int savedErrno=errno;
        const long page=sysconf(_SC_PAGESIZE);
        if(page>0) {
            const auto range=wholePages(reinterpret_cast<std::uintptr_t>(memory),bytes,static_cast<std::size_t>(page));
            result.bytes=range.bytes;
            if(range.bytes) {
#ifndef SPIN_WORKSPACE_ROUTING_ADVICE_FAILURE_TEST
                auto* begin=reinterpret_cast<void*>(range.begin);
                result.hinted=madvise(begin,range.bytes,MADV_HUGEPAGE)==0;
#ifdef MADV_COLLAPSE
                result.collapsed=madvise(begin,range.bytes,MADV_COLLAPSE)==0;
#endif
#endif
            }
        }
        errno=savedErrno;
    }
#else
    (void)memory;(void)bytes;
#endif
    return result;
}

// Same scalar scatter and lookahead as the measured winner. The write hint is
// selected at compile time; there is no per-coordinate policy branch.
template<bool Packed,class Block>
inline void scatterWrite(const Block* values,Block* tile,const std::uint8_t* offsets24,
                         const std::uint32_t* offsets32,std::size_t base,std::size_t count) noexcept {
    auto offset=[&](std::size_t position) {
        if constexpr(Packed) {
            std::uint32_t value;
            std::memcpy(&value,offsets24+3*position,4);
            return value&0xffffffU;
        } else return offsets32[position];
    };
    for(std::size_t j=0;j<count;++j) {
#if defined(__GNUC__) || defined(__clang__)
        if(j+32<count) __builtin_prefetch(tile+offset(base+j+32),1,3);
#endif
        tile[offset(base+j)]=values[base+j];
    }
}
}
