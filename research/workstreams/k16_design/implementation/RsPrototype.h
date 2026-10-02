#pragma once
#include "../../../../spin/src/kernels/Block.h"
#include "../../../../spin/src/packet/PacketLargeTypes.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::research::rs {
using Block = detail::storage::block;
inline constexpr std::size_t groupStride = 260;
enum class Variant { Rs8Gf256, Rs16Gf16 };
enum class InnerKernel { Cached, RetainedStreaming };

// A new research ensemble, not PacketT64S16's seed-compatible setup. Four
// parallel systematic RS[8,4] rows use GF256 modulus 0x11b and points 0..7.
// Each group independently shuffles all64 packet columns; each of64 regions
// independently shuffles all groups. Each of8 symbols/group receives an
// independent GL32, followed by the selected t64/s16 inner and GL16 updates.
// The RS16 variant uses four parallel RS[16,8] rows over GF16 modulus0x13,
// points0..15, and independent GL16 at each of16 four-row symbols instead.
// The existing deterministic Words stream implements reproducible sampling;
// distributional statements use the ideal independent-word setup model.
struct Plan {
    Plan(std::size_t messageBits, std::uint64_t routeSeed, std::uint64_t innerSeed,
         Variant selected = Variant::Rs8Gf256, InnerKernel kernel = InnerKernel::Cached);
    Plan(const Plan&) = delete;
    Plan& operator=(const Plan&) = delete;
    Plan(Plan&&) = default;
    Plan& operator=(Plan&&) = default;
    std::size_t k = 0, n = 0, groups = 0;
    Variant variant = Variant::Rs8Gf256;
    std::vector<std::uint32_t> route;
    // These are binary transpose matrices; forwardScalar explicitly applies
    // their adjoints. Keeping both scalar and compact forms aids validation.
    std::vector<std::array<std::uint16_t,16>> reverseMatrices;
    std::vector<std::array<std::uint32_t,32>> outerMatrices;
    std::vector<std::array<std::uint16_t,16>> outerMatrices16;
    std::vector<std::uint64_t> compactStorage, composedUpdates;
    std::vector<detail::packet::large::Dense16Row> denseUpdates;
    InnerKernel innerKernel = InnerKernel::Cached;
    const std::uint64_t* compactCoefficients() const noexcept;
    std::size_t scratchBlocks() const noexcept { return groups * groupStride; }
    std::size_t setupBytes() const noexcept;
};

// Fast entry points allocate nothing. Scratch must be64-byte aligned and
// disjoint from input/output, holding plan.scratchBlocks() records. Input and
// output need only16-byte alignment; output may equal input. The input has
// plan.n records, output plan.k. No flush symbols are emitted.
// The Fast entry points require a Cached plan; the Large entry points require
// RetainedStreaming. outerFast and the scalar oracles accept either plan.
void reverseRouteFast(const Block* input, Block* scratch, const Plan&);
void outerFast(const Block* scratch, Block* output, const Plan&);
void transposeFast(const Block* input, Block* output, Block* scratch, const Plan&);
// Explicit research candidate; same binary map, complete retained mode-19
// inner plus full-cache-line streaming stores. Does not change the default.
void reverseRouteLarge(const Block* input, Block* scratch, const Plan&);
void transposeLarge(const Block* input, Block* output, Block* scratch, const Plan&);
void transposeScalar(const Block* input, Block* output, Block* scratch, const Plan&);

// Independent literal forward map for binary-adjoint checks. Allocates its
// temporary route buffer; this validation-only function is not benchmarked.
void forwardScalar(const Block* message, Block* encoded, const Plan&);
}
