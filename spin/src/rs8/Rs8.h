#pragma once
#include "../paired15/Paired15.h"
#include "../certified22/Plan.h"
#include <optional>

namespace spin::detail::rs8 {
using Block = storage::block;
inline constexpr std::size_t groupStride = 260;
static_assert(groupStride == paired15::groupStride);
static_assert(groupStride == certified22::groupStride);

bool validMessageSize(std::size_t k) noexcept;

// Eight systematic GF(16) RS[8,4] rows per 128-to-256-coordinate group,
// independent nonzero GF(2^32) reverse symbol maps, and four-bit routing
// through 64 regions, with the exact t64/s22 inner at every size (the map
// certified in Lean at K = 2^14..2^26). State crosses region boundaries.
struct Plan {
    Plan(std::size_t k, std::uint64_t routeSeed, std::uint64_t innerSeed,
         certified22::SetupPath path = certified22::SetupPath::Automatic);
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;
    std::size_t k = 0, n = 0, groups = 0;
    std::optional<paired15::Plan> paired15Plan;
    std::optional<certified22::Plan> certified22Plan;
    std::vector<std::array<std::uint8_t,9>> field;
    // The sampled nonzero GF(2^32) scalar per symbol; the scalar oracles derive
    // the literal reverse map's columns from it on demand.
    std::vector<std::uint32_t> scalars;
    // GFNI matrices of the nine adjoint tower byte multiplications per column.
    std::vector<std::array<std::uint64_t,9>> adjoint;
    std::size_t scratchBlocks() const noexcept { return groups * groupStride; }
    std::size_t setupBytes() const noexcept;
};

// Scratch is 64-byte aligned and disjoint from input/output. No entry point
// allocates. Exact input/output aliasing is supported at 16-byte alignment.
// Forward lanes=1,2,4 selects interleaved 128/256/512-bit records; scratch
// consists of lanes contiguous planes of scratchBlocks() 128-bit blocks.
void transposeScalar(const Block*, Block*, Block* scratch, const Plan&);
void forwardScalar(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes);
void transposeFast(const Block*, Block*, Block* scratch, const Plan&);
// AVX-512 forward: same lanes/scratch contract as forwardScalar. Streaming
// output stores are used only when requested and output is 64-byte aligned.
void forwardFast(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes, bool stream);
void outerForwardFast(const Block*, Block* scratch, const Plan&, unsigned lanes);
// Portable kernels (Portable.cpp; plain C++, any 64-bit CPU): same contracts
// as the scalar entry points.
void transposePortable(const Block*, Block*, Block* scratch, const Plan&);
void forwardPortable(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes);
// AVX2 kernels (Avx2.cpp; no AVX-512, VBMI or GFNI): same contracts as the
// scalar entry points, require cpu_avx2().
void transposeAvx2(const Block*, Block*, Block* scratch, const Plan&);
void forwardAvx2(const Block*, Block*, Block* scratch, const Plan&, unsigned lanes, bool stream);
void reverseRouteAvx2(const Block*, Block* scratch, const Plan&);
void outerAvx2(const Block* scratch, Block*, const Plan&);
void outerForwardAvx2(const Block*, Block* scratch, const Plan&, unsigned lanes);
void forwardInnerAvx2(const Block* scratch, Block*, const Plan&, unsigned lanes);

// Internal stage entry points and independent literal outer oracles. Fast
// outer stores stream above K18 when output is 64-byte aligned.
void outerFast(const Block* scratch, Block*, const Plan&);
void outerCached(const Block* scratch, Block*, const Plan&);
void outerScalar(const Block* scratch, Block*, const Plan&);
void outerForwardScalar(const Block*, Block* scratch, const Plan&, unsigned lanes);
bool selfTest() noexcept;
// Rewrites field (and scalars when requested) in place for a new outer seed,
// exactly as the constructor would; adjoint is left stale. No allocation.
void refreshOuterTranspose(Plan&, std::uint64_t seed, bool scalars) noexcept;
// Literal per-symbol definitions (bit-serial field arithmetic), for oracles.
std::array<std::uint32_t,32> reverseColumns(std::uint32_t scalar) noexcept;
std::array<std::uint64_t,9> literalAdjoint(std::uint32_t scalar) noexcept;
}
