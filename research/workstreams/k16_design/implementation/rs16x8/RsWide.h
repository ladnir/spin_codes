#pragma once
#include "../RsPrototype.h"

namespace spin::research::rswide {
using Block = rs::Block;
inline constexpr std::size_t groupStride = 516;
enum class Randomizer { DenseGl32, TowerField32, TowerByte32 };
enum class InnerRandomizer { DenseGl16, TowerByte16 };
struct Plan {
    Plan(std::size_t messageBits, std::uint64_t routeSeed, std::uint64_t innerSeed,
         Randomizer selected = Randomizer::DenseGl32,
         InnerRandomizer innerSelected = InnerRandomizer::DenseGl16);
    std::size_t k = 0, n = 0, groups = 0;
    Randomizer randomizer = Randomizer::DenseGl32;
    InnerRandomizer innerRandomizer = InnerRandomizer::DenseGl16;
    std::vector<std::uint32_t> route;
    std::vector<std::array<std::uint16_t,16>> reverseMatrices;
    std::vector<detail::packet::large::Dense16Row> denseUpdates;
    std::vector<std::array<std::uint8_t,3>> fieldUpdates;
    std::vector<std::array<std::uint32_t,32>> outerMatrices;
    std::vector<std::uint64_t> compactStorage;
    const std::uint64_t* compactCoefficients() const noexcept;
    std::size_t coefficientsPerSymbol() const noexcept {
        return randomizer == Randomizer::DenseGl32 ? 16 : 9;
    }
    std::size_t coefficientBytesPerSymbol() const noexcept {
        return randomizer == Randomizer::TowerByte32 ? 9 : 8*coefficientsPerSymbol();
    }
    std::size_t scratchBlocks() const noexcept { return groups * groupStride; }
};
// Eight GF16 RS[16,8] rows, independent GL32 on each aligned eight-nibble
// symbol; independent group-column and regional permutations, 128 regions.
// Retains the original t64/s16 streaming inner, without a new proof claim.
// The two four-row halves use precisely the retained fast16 packed layout.
void reverseRoute(const Block*, Block*, const Plan&);
void outerFast(const Block*, Block*, const Plan&);
void transposeFast(const Block*, Block*, Block*, const Plan&);
void transposeScalar(const Block*, Block*, Block*, const Plan&);
void forwardScalar(const Block*, Block*, const Plan&);
void outerScalar(const Block*, Block*, const Plan&);
}
