#include "PacketPlan.h"
#include "../kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spin::detail::packet {
namespace {
using kernel::setup::Divisor;
using kernel::setup::Words;

// Basis used only by the fused kernel. These rows and their inverse are over
// GF(2); the fallback retains the original physical state coordinates.
constexpr std::array<std::uint16_t, 16> basis{
    0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0x3500,
    0x1180, 0x5980, 0x2c80, 0xcf80, 0x9400, 0xe400, 0xe00, 0xd980};
constexpr std::array<std::uint16_t, 16> inverseBasis{
    0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40, 0xa500,
    0x6980, 0x1a00, 0x5f80, 0x580, 0xcd80, 0xfb00, 0x680, 0x8200};

template<class Word, std::size_t Width>
bool fullRank(std::array<Word, Width> rows) noexcept {
    unsigned rank = 0;
    for(unsigned column = 0; column < Width; ++column) {
        unsigned pivot = rank;
        while(pivot < Width && !(rows[pivot] & (Word{1} << column))) ++pivot;
        if(pivot == Width) return false;
        std::swap(rows[pivot], rows[rank]);
        for(unsigned j = rank + 1; j < Width; ++j)
            if(rows[j] & (Word{1} << column)) rows[j] ^= rows[rank];
        ++rank;
    }
    return true;
}

template<class Container>
void shuffle(Container& values, Words& words, const std::vector<Divisor>& divisors) {
    std::iota(values.begin(), values.end(), 0U);
    for(std::size_t i = values.size(); i > 1; --i)
        std::swap(values[i - 1], values[divisors[i].sample(words, i)]);
}

void prepareRoute(Plan& plan, std::uint64_t seed) {
    const auto rows = plan.n / 256;
    const auto groups = rows / 4;
    Words words(seed);
    std::vector<Divisor> divisors(std::max(rows, std::size_t{256}) + 1);
    for(std::size_t i = 2; i < divisors.size(); ++i) divisors[i] = Divisor(i);
    // All four rows in a tile share one independent column permutation.
    std::vector<std::array<std::uint32_t, 256>> columns(groups);
    for(auto& permutation : columns) shuffle(permutation, words, divisors);
    std::vector<std::uint32_t> destinations(groups);
    plan.route.resize(plan.n / 4);
    for(unsigned region = 0; region < 256; ++region) {
        shuffle(destinations, words, divisors);
        for(std::size_t group = 0; group < groups; ++group) {
            // The measured construction's descriptor stream also drew a GF16
            // multiplier at this point. GL32 replaced that operation, but keep
            // its (possibly rejected) draw to preserve exact seed compatibility.
            (void)divisors[15].sample(words, 15);
            plan.route[region * groups + destinations[group]] =
                static_cast<std::uint32_t>(group * tileStride + 4 * columns[group][region]);
        }
    }
}

void prepareUpdates(Plan& plan, std::uint64_t seed) {
    Words words(seed ^ 0x3f625a92ULL);
    plan.reverseMatrices.resize(plan.n / 64);
    plan.composedUpdates.resize(4 * (plan.n / 64));
    for(std::size_t epoch = 0; epoch < plan.n / 64; ++epoch) {
        std::array<std::uint16_t, 16> forward{};
        do {
            for(auto& row : forward) row = static_cast<std::uint16_t>(words());
        } while(!fullRank(forward));
        auto& reverse = plan.reverseMatrices[epoch];
        for(unsigned i = 0; i < 16; ++i)
            for(unsigned j = 0; j < 16; ++j)
                reverse[j] |= static_cast<std::uint16_t>(((forward[i] >> j) & 1U) << i);

        std::array<std::uint16_t, 16> transformed{};
        for(unsigned output = 0; output < 16; ++output) {
            unsigned originalMask = 0;
            for(unsigned mask = basis[output]; mask; mask &= mask - 1)
                originalMask ^= reverse[std::countr_zero(mask)];
            for(unsigned mask = originalMask; mask; mask &= mask - 1)
                transformed[output] ^= inverseBasis[std::countr_zero(mask)];
        }
        // GFNI data operands: ascending output rows in each 8x8 submatrix.
        for(unsigned out = 0; out < 2; ++out)
            for(unsigned in = 0; in < 2; ++in)
                for(unsigned j = 0; j < 8; ++j)
                    plan.composedUpdates[4 * epoch + 2 * out + in] |=
                        std::uint64_t((transformed[8 * out + j] >> (8 * in)) & 255U) << (8 * j);
    }
}

void prepareOuter(Plan& plan, std::uint64_t seed) {
    Words words(seed ^ 0x75a1dc09ULL);
    const auto tiles = plan.n / 1024;
    plan.outerMatrices.resize(32 * tiles);
    plan.compactStorage.resize(512 * tiles + 7);
    auto* coefficients = const_cast<std::uint64_t*>(plan.compactCoefficients());
    for(std::size_t tile = 0; tile < tiles; ++tile)
        for(unsigned group = 0; group < 32; ++group) {
            auto& matrix = plan.outerMatrices[32 * tile + group];
            do {
                for(auto& row : matrix) row = static_cast<std::uint32_t>(words());
            } while(!fullRank(matrix));
            for(unsigned diagonal = 0; diagonal < 4; ++diagonal)
                for(unsigned lane = 0; lane < 4; ++lane) {
                    std::uint64_t value = 0;
                    for(unsigned j = 0; j < 8; ++j)
                        value |= std::uint64_t((matrix[8 * lane + j] >> (8 * ((lane + diagonal) % 4))) & 255U)
                                 << (8 * (7 - j));
                    coefficients[512 * tile + 16 * group + 4 * diagonal + lane] = value;
                }
        }
}
}

bool validMessageSize(std::size_t k) noexcept {
    // The padded route representation, not a performance/certificate policy.
    constexpr auto maxTiles = std::numeric_limits<std::uint32_t>::max() / tileStride;
    return k && k % 512 == 0 && k / 512 <= maxTiles;
}

Plan::Plan(std::size_t k, std::uint64_t seed) : Plan(k, seed, seed) {}

Plan::Plan(std::size_t k, std::uint64_t routeSeed, std::uint64_t innerSeed) {
    if(!validMessageSize(k))
        throw std::invalid_argument("packet SPIN K must be a positive multiple of 512 within the 32-bit padded routing range");
    n = 2 * k;
    prepareRoute(*this, routeSeed);
    prepareUpdates(*this, innerSeed);
    prepareOuter(*this, routeSeed);
}

const std::uint64_t* Plan::compactCoefficients() const noexcept {
    const auto address = reinterpret_cast<std::uintptr_t>(compactStorage.data());
    const auto offset = ((64 - (address & 63U)) & 63U) / sizeof(std::uint64_t);
    return compactStorage.data() + offset;
}

std::size_t Plan::setupBytes() const noexcept {
    return route.capacity() * sizeof(route[0])
         + reverseMatrices.capacity() * sizeof(reverseMatrices[0])
         + outerMatrices.capacity() * sizeof(outerMatrices[0])
         + compactStorage.capacity() * sizeof(compactStorage[0])
         + composedUpdates.capacity() * sizeof(composedUpdates[0]);
}
}
