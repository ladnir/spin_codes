#pragma once
#include <spin/Code.h>
#include "kernels/WorkspaceRouting.h"

namespace spin::detail {
inline constexpr std::size_t automaticHugePageMinBytes = 16 * 1024 * 1024;

// Resolve per allocation at setup, never in an encoding kernel. Keep the
// requested policy on workspaces so resizing reevaluates Automatic.
inline MemoryPolicy resolveMemoryPolicy(std::size_t bytes, MemoryPolicy policy) {
    switch(policy) {
    case MemoryPolicy::Automatic:
        return kernel::workspace_routing::compiled && bytes >= automaticHugePageMinBytes
            ? MemoryPolicy::PreferHugePages : MemoryPolicy::Normal;
    case MemoryPolicy::Normal:
    case MemoryPolicy::PreferHugePages:
        return policy;
    default:
        throw std::invalid_argument("SPIN unknown memory policy");
    }
}
}
