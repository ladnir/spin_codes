#pragma once
#include "Spin.h"
#include <span>

namespace bare_spin {
bool wideAvailable(unsigned lanes) noexcept;
// Forward E only, S19 or R2, with packed or 32-bit routes. lanes=2/4 means 256/512 bits.
// Element i occupies [i*lanes, (i+1)*lanes) in both spans: streams are interleaved.
// Each 128-bit lane contains 128 independent binary messages, not field arithmetic.
// code must outlive this workspace and may not be moved or mutated while in use.
// A workspace is not thread-safe; use one per concurrent caller.
class WideWorkspace {
public:
    WideWorkspace(const Spin& code, unsigned lanes);
    ~WideWorkspace();
    WideWorkspace(const WideWorkspace&)=delete;
    WideWorkspace& operator=(const WideWorkspace&)=delete;
    // Exact sizes: K*lanes input blocks and 2K*lanes output blocks. No overlap.
    // Only block (16-byte) alignment is needed. Reuses scratch; no per-call allocation.
    void forward(std::span<const block> input, std::span<block> output);
    std::size_t bytes() const noexcept;
private:
    const Spin& mCode;
    unsigned mLanes;
    void* mWork;
};
}
