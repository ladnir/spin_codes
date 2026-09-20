#pragma once
#include "Spin.h"
#include <span>

namespace bare_spin {
bool wideAvailable(unsigned lanes) noexcept;
// Forward E only, S19 or K16 R2 with Packed24. Element i occupies lanes consecutive blocks.
// Each lane is an independent binary plane, NOT extension-field multiplication.
// code must outlive this workspace and may not be moved or mutated while in use.
// A workspace is not thread-safe; use one per concurrent caller.
class WideWorkspace {
public:
    WideWorkspace(const Spin& code, unsigned lanes);
    ~WideWorkspace();
    WideWorkspace(const WideWorkspace&)=delete;
    WideWorkspace& operator=(const WideWorkspace&)=delete;
    void forward(std::span<const block> input, std::span<block> output);
    std::size_t bytes() const noexcept;
private:
    const Spin& mCode;
    unsigned mLanes;
    void* mWork;
};
}
