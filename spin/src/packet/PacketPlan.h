#pragma once
#include "../kernels/Block.h"
#include "PacketLargeTypes.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::detail::packet {
using Block=storage::block;
inline constexpr std::size_t groupStride=516;
bool validMessageSize(std::size_t k) noexcept;

// The 16+4 state coordinates are packed as two full byteplanes and one
// four-bit plane. Layout exactly matches the measured fused refresh kernel.
struct BorderRow {
    large::Dense16Row base;
    std::uint64_t upper[2]{},lower[2]{},corner=0;
};
static_assert(sizeof(BorderRow)==72);

// Fixed rate-1/2 RS packet construction: eight GF16 RS[16,8] rows, independent
// adjoint GF32 multipliers, 128 packet regions, and the t64/s20 inner.
// All counts and route offsets are in 128-bit records, not bytes.
struct Plan {
    explicit Plan(std::size_t k,std::uint64_t seed);
    Plan(std::size_t k,std::uint64_t routeSeed,std::uint64_t innerSeed);
    Plan(const Plan&)=delete;
    Plan& operator=(const Plan&)=delete;
    Plan(Plan&&)=delete;
    Plan& operator=(Plan&&)=delete;
    std::size_t k=0,n=0,groups=0;
    std::vector<std::uint32_t> route;
    std::vector<std::array<std::uint32_t,20>> reverseMatrices;
    std::vector<std::array<std::uint32_t,32>> outerMatrices;
    std::vector<std::uint64_t> compactStorage;
    std::vector<BorderRow> updates;

    const std::uint64_t* compactCoefficients() const noexcept;
    std::size_t scratchBlocks() const noexcept { return groups*groupStride; }
    std::size_t setupBytes() const noexcept;
};

// Allocation-free. The whole input is consumed before any output is written,
// so output may equal input. Scratch is disjoint, 64-byte aligned, and holds
// scratchBlocks() records. Input/output need only 16-byte alignment.
void transposeScalar(const Block* input,Block* output,Block* scratch,const Plan&);
void transposeFast(const Block* input,Block* output,Block* scratch,const Plan&);

// Private stage boundaries used by the retained fast kernel.
void reverseRoute(const Block* input,Block* scratch,const Plan&);
void outerFast(const Block* scratch,Block* output,const Plan&);
// Exact same map and instruction schedule, cached stores for small working sets.
void reverseRouteCached(const Block* input,Block* scratch,const Plan&);
void outerFastCached(const Block* scratch,Block* output,const Plan&);

// Allocation-free portable forward; same buffer contract as the fast path.
void forwardScalar(const Block* message,Block* encoded,Block* scratch,const Plan&);

// Literal independent forward oracle, used by package correctness tests.
// Unlike the encoder entry points above, it allocates scratch storage.
void forwardScalar(const Block* message,Block* encoded,const Plan&);
}
