#pragma once
#include "Full.h"
#include <cstdlib>
#if defined(__linux__)
#include <sys/mman.h>
#endif

namespace spin::research::char2 {
// Benchmark allocation only. Setup and page faults stay outside timed calls.
struct BenchPayload {
    Block* p{};
    std::size_t count, allocated;
    explicit BenchPayload(std::size_t n, bool preferHuge = true) : count(n) {
        const bool huge = preferHuge && n * sizeof(Block) >= 16 * 1024 * 1024;
        const std::size_t alignment = huge ? 2 * 1024 * 1024 : 64;
        allocated = (n * sizeof(Block) + alignment - 1) & ~(alignment - 1);
#if defined(_MSC_VER)
        p = static_cast<Block*>(_aligned_malloc(allocated, alignment));
#else
        p = static_cast<Block*>(std::aligned_alloc(alignment, allocated));
#endif
        if (!p) throw std::bad_alloc();
#if defined(__linux__)
        if (huge) madvise(p, allocated, MADV_HUGEPAGE);
#endif
        std::memset(p, 0, allocated);
    }
    ~BenchPayload() {
#if defined(_MSC_VER)
        _aligned_free(p);
#else
        std::free(p);
#endif
    }
    BenchPayload(const BenchPayload&) = delete;
    BenchPayload& operator=(const BenchPayload&) = delete;
};
}
