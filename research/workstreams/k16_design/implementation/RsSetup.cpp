#include "RsPrototype.h"
#include "../../../../spin/src/kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spin::research::rs {
namespace {
using detail::kernel::setup::Divisor;
using detail::kernel::setup::Words;
constexpr std::array<std::uint16_t,16> basis{
    0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x3500,
    0x1180,0x5980,0x2c80,0xcf80,0x9400,0xe400,0xe00,0xd980};
constexpr std::array<std::uint16_t,16> inverseBasis{
    0x1,0x2,0x4,0x8,0x10,0x20,0x40,0xa500,
    0x6980,0x1a00,0x5f80,0x580,0xcd80,0xfb00,0x680,0x8200};

template<class Word, std::size_t Width>
bool fullRank(std::array<Word,Width> rows) {
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

template<class Word, std::size_t Width>
std::array<Word,Width> randomTranspose(Words& words) {
    std::array<Word,Width> forward{}, reverse{};
    do {
        for(auto& row : forward) row = static_cast<Word>(words());
    } while(!fullRank(forward));
    for(unsigned i = 0; i < Width; ++i)
        for(unsigned j = 0; j < Width; ++j)
            reverse[j] |= Word(((forward[i] >> j) & 1U) << i);
    return reverse;
}

template<class Container>
void shuffle(Container& values, Words& words, const std::vector<Divisor>& divisors) {
    std::iota(values.begin(), values.end(), 0U);
    for(std::size_t i = values.size(); i > 1; --i)
        std::swap(values[i - 1], values[divisors[i].sample(words, i)]);
}

void prepareRoute(Plan& plan, std::uint64_t seed) {
    Words words(seed);
    std::vector<Divisor> divisors(std::max(plan.groups, std::size_t{64}) + 1);
    for(std::size_t i = 2; i < divisors.size(); ++i) divisors[i] = Divisor(i);
    std::vector<std::array<std::uint32_t,64>> columns(plan.groups);
    for(auto& permutation : columns) shuffle(permutation, words, divisors);
    std::vector<std::uint32_t> destinations(plan.groups);
    plan.route.resize(plan.n / 4);
    for(unsigned region = 0; region < 64; ++region) {
        shuffle(destinations, words, divisors);
        for(std::size_t group = 0; group < plan.groups; ++group)
            plan.route[region * plan.groups + destinations[group]] =
                static_cast<std::uint32_t>(group * groupStride + 4 * columns[group][region]);
    }
}

void prepareUpdates(Plan& plan, std::uint64_t seed) {
    Words words(seed ^ 0x3f625a92ULL);
    plan.reverseMatrices.resize(plan.n / 64);
    const bool large = plan.innerKernel == InnerKernel::RetainedStreaming;
    if(large) plan.denseUpdates.resize(plan.n / 64);
    else plan.composedUpdates.resize(4 * (plan.n / 64));
    for(std::size_t epoch = 0; epoch < plan.n / 64; ++epoch) {
        const auto reverse = randomTranspose<std::uint16_t,16>(words);
        plan.reverseMatrices[epoch] = reverse;
        std::array<std::uint16_t,16> transformed{};
        for(unsigned output = 0; output < 16; ++output) {
            unsigned originalMask = 0;
            for(unsigned mask = basis[output]; mask; mask &= mask - 1)
                originalMask ^= reverse[std::countr_zero(mask)];
            for(unsigned mask = originalMask; mask; mask &= mask - 1)
                transformed[output] ^= inverseBasis[std::countr_zero(mask)];
        }
        for(unsigned out = 0; out < 2; ++out)
            for(unsigned in = 0; in < 2; ++in)
                for(unsigned j = 0; j < 8; ++j) {
                    const auto byte = std::uint64_t((transformed[8 * out + j] >> (8 * in)) & 255U);
                    if(large) plan.denseUpdates[epoch].matrix[2 * out + in] |= byte << (8 * (7 - j));
                    else plan.composedUpdates[4 * epoch + 2 * out + in] |= byte << (8 * j);
                }
    }
}

void prepareOuter(Plan& plan, std::uint64_t seed) {
    Words words(seed ^ 0x75a1dc09ULL);
    plan.outerMatrices.resize(8 * plan.groups);
    plan.compactStorage.resize(128 * plan.groups + 7);
    auto* coefficients = const_cast<std::uint64_t*>(plan.compactCoefficients());
    for(std::size_t group = 0; group < plan.groups; ++group)
        for(unsigned symbol = 0; symbol < 8; ++symbol) {
            auto& matrix = plan.outerMatrices[8 * group + symbol];
            matrix = randomTranspose<std::uint32_t,32>(words);
            for(unsigned diagonal = 0; diagonal < 4; ++diagonal)
                for(unsigned lane = 0; lane < 4; ++lane) {
                    std::uint64_t value = 0;
                    for(unsigned j = 0; j < 8; ++j)
                        value |= std::uint64_t((matrix[8 * lane + j] >> (8 * ((lane + diagonal) % 4))) & 255U)
                                 << (8 * (7 - j));
                    coefficients[128 * group + 16 * symbol + 4 * diagonal + lane] = value;
                }
        }
}

void prepareOuter16(Plan& plan, std::uint64_t seed) {
    Words words(seed ^ 0x75a1dc09ULL);
    plan.outerMatrices16.resize(16 * plan.groups);
    plan.compactStorage.resize(64 * plan.groups + 7);
    auto* coefficients = const_cast<std::uint64_t*>(plan.compactCoefficients());
    for(std::size_t group = 0; group < plan.groups; ++group)
        for(unsigned symbol = 0; symbol < 16; ++symbol) {
            auto& matrix = plan.outerMatrices16[16 * group + symbol];
            matrix = randomTranspose<std::uint16_t,16>(words);
            // Coordinates0..7 concatenate rows0/1;8..15 concatenate rows2/3.
            // GFNI's high matrix byte supplies the low output bit.
            for(unsigned out = 0; out < 2; ++out)
                for(unsigned in = 0; in < 2; ++in) {
                    std::uint64_t value = 0;
                    for(unsigned j = 0; j < 8; ++j)
                        value |= std::uint64_t((matrix[8 * out + j] >> (8 * in)) & 255U) << (8 * (7 - j));
                    coefficients[64 * group + 4 * symbol + 2 * out + in] = value;
                }
        }
}
}

Plan::Plan(std::size_t messageBits, std::uint64_t routeSeed, std::uint64_t innerSeed,
           Variant selected, InnerKernel kernel) {
    constexpr auto maxGroups = std::numeric_limits<std::uint32_t>::max() / groupStride;
    // Whole physical inner epochs in each region; K4096 and K65536 are valid.
    if(!messageBits || messageBits % 2048 || messageBits / 128 > maxGroups)
        throw std::invalid_argument("RS prototype K must be a positive multiple of2048 within the32-bit padded route range");
    k = messageBits;
    n = 2 * k;
    groups = k / 128;
    if(selected != Variant::Rs8Gf256 && selected != Variant::Rs16Gf16)
        throw std::invalid_argument("unknown RS research variant");
    variant = selected;
    if(kernel != InnerKernel::Cached && kernel != InnerKernel::RetainedStreaming)
        throw std::invalid_argument("unknown RS research inner kernel");
    innerKernel = kernel;
    prepareRoute(*this, routeSeed);
    prepareUpdates(*this, innerSeed);
    if(variant == Variant::Rs8Gf256) prepareOuter(*this, routeSeed);
    else prepareOuter16(*this, routeSeed);
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
         + outerMatrices16.capacity() * sizeof(outerMatrices16[0])
         + compactStorage.capacity() * sizeof(compactStorage[0])
         + composedUpdates.capacity() * sizeof(composedUpdates[0])
         + denseUpdates.capacity() * sizeof(denseUpdates[0]);
}
}
