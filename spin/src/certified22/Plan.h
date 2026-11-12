#pragma once
#include "../kernels/Block.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::detail::certified22 {
using Block = storage::block;
inline constexpr std::size_t groupStride = 260;

bool validMessageSize(std::size_t k) noexcept;

// Automatic uses the AVX-512 setup when the CPU has the packet ISA, else the
// AVX2 setup when available, else the portable fast setup; all paths
// produce bit-identical tables. Scalar (the literal oracle), Avx2 and Portable
// force a path for differential tests (Avx2 requires cpu_avx2()).
enum class SetupPath { Automatic, Scalar, Avx2, Portable };

// Four systematic GF(16) RS[16,8] rows per 128-to-256-coordinate group,
// native GF(2^16) symbol maps, four-bit routing through 64 regions, and
// the exact t64/s22 inner. State runs continuously across region boundaries.
struct Plan {
    Plan(std::size_t k, std::uint64_t routeSeed, std::uint64_t innerSeed,
         bool withOuter = true, SetupPath path = SetupPath::Automatic);
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;
    std::size_t k = 0, n = 0, groups = 0;
    std::vector<std::uint32_t> route, routeBytes;
    std::vector<std::array<std::uint32_t,22>> reverseMatrices;
    std::vector<std::uint64_t> updates;
    // Forward-direction updates in the same fast-circuit basis and packed
    // 3x3 GFNI block layout: basis * R^T * inverseBasis per epoch.
    std::vector<std::uint64_t> forwardUpdates;
    std::vector<std::array<std::uint16_t,16>> outerMatrices16;
    std::vector<std::array<std::uint8_t,3>> outerField;
    std::size_t scratchBlocks() const noexcept { return groups * groupStride; }
    std::size_t setupBytes() const noexcept;
};

// Scratch is 64-byte aligned and disjoint from input/output. No entry point
// allocates. Exact input/output aliasing is supported at 16-byte alignment.
// Forward lanes=1,2,4 means interleaved 128/256/512-bit records; scratch
// consists of lanes contiguous planes of scratchBlocks() 128-bit blocks.
void transposeScalar(const Block*, Block*, Block* scratch, const Plan&);
void forwardScalar(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes);
void transposeFast(const Block*, Block*, Block* scratch, const Plan&);

// Internal stage entry points; also used for differential kernel checks.
// withOuter=false retains identical route/inner streams without RS16 tables.
void reverseRouteScalar(const Block*, Block* scratch, const Plan&);
// plane is the lane-plane stride in blocks; 0 selects scratchBlocks().
void forwardInnerScalar(const Block* scratch, Block*, const Plan&, unsigned lanes, std::size_t plane = 0);
void reverseRouteFast(const Block*, Block* scratch, const Plan&);
// Forward inner: gathers routed packets from lanes contiguous scratch planes
// and writes 2K interleaved records. Streaming requires 64-byte output.
void forwardInnerFast(const Block* scratch, Block* output, const Plan&, unsigned lanes, bool stream,
                      std::size_t plane = 0);
// AVX-512 exact inner setup (SetupFast.cpp); requires cpu_packet512().
void prepareInnerFast(Plan&, std::uint64_t innerSeed);
// AVX2 exact inner setup (SetupAvx2.cpp), bit-identical; requires cpu_avx2().
// Portable fast exact inner setup (SetupPortable.cpp), bit-identical.
void prepareInnerPortable(Plan&, std::uint64_t innerSeed);
void prepareInnerAvx2(Plan&, std::uint64_t innerSeed);
// Words first..first+count-1 (count a multiple of 8) of kernel::setup::Words(seed).
void fillWords(std::uint64_t seed, std::uint64_t first, std::uint64_t* out, std::size_t count) noexcept;
void outerFastCached(const Block* scratch, Block*, const Plan&);
void outerFastCachedUnrolled(const Block* scratch, Block*, const Plan&);
void outerFastStreaming(const Block* scratch, Block*, const Plan&);
}
