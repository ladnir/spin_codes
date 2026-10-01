#include <spin/Code.h>
#include "kernels/WorkspaceRouting.h"
#include <cstring>
#include <limits>
#include <new>
#include <stdexcept>

namespace spin::detail {
struct BufferAllocation {
    std::byte* data = nullptr;
    std::size_t logical = 0, allocated = 0, alignment = 64;

    BufferAllocation(std::size_t bytes, MemoryPolicy policy) : logical(bytes) {
        if(policy != MemoryPolicy::Normal && policy != MemoryPolicy::PreferHugePages)
            throw std::invalid_argument("SPIN unknown memory policy");
        if(policy == MemoryPolicy::PreferHugePages) alignment = 2 * 1024 * 1024;
        if(bytes > std::numeric_limits<std::size_t>::max() - (alignment - 1))
            throw std::bad_alloc();
        allocated = (bytes + alignment - 1) & ~(alignment - 1);
        data = static_cast<std::byte*>(::operator new(allocated, std::align_val_t(alignment)));
        std::memset(data, 0, allocated);
        // Only pages wholly owned by this allocation are advised. This is a
        // setup hint, never an encoding operation or a system-wide setting.
        if(policy == MemoryPolicy::PreferHugePages)
            kernel::workspace_routing::adviseOwned(data, allocated);
    }
    ~BufferAllocation() { ::operator delete(data, std::align_val_t(alignment)); }
    BufferAllocation(const BufferAllocation&) = delete;
    BufferAllocation& operator=(const BufferAllocation&) = delete;
};
}

namespace spin {
Buffer::Buffer(std::size_t bytes, MemoryPolicy policy) {
    if(policy != MemoryPolicy::Normal && policy != MemoryPolicy::PreferHugePages)
        throw std::invalid_argument("SPIN unknown memory policy");
    const std::size_t alignment = policy == MemoryPolicy::PreferHugePages ? 2 * 1024 * 1024 : 64;
    if(bytes > std::numeric_limits<std::size_t>::max() - (alignment - 1))
        throw std::bad_alloc();
    if(bytes) allocation_ = std::make_unique<detail::BufferAllocation>(bytes, policy);
}
Buffer::~Buffer() = default;
Buffer::Buffer(Buffer&&) noexcept = default;
Buffer& Buffer::operator=(Buffer&&) noexcept = default;
std::span<std::byte> Buffer::bytes() noexcept {
    return allocation_ ? std::span<std::byte>(allocation_->data, allocation_->logical)
                       : std::span<std::byte>{};
}
std::span<const std::byte> Buffer::bytes() const noexcept {
    return allocation_ ? std::span<const std::byte>(allocation_->data, allocation_->logical)
                       : std::span<const std::byte>{};
}
std::size_t Buffer::allocation_bytes() const noexcept {
    return allocation_ ? allocation_->allocated : 0;
}
}
