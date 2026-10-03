#include "Paired15.h"
#include "../kernels/SetupRandom.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>

namespace spin::detail::paired15 {
namespace {
using kernel::setup::Words;
using kernel::setup::Divisor;
constexpr unsigned permutation[16] = {0,1,2,7,3,4,5,6,8,15,10,12,11,13,14,9};
using Matrix = std::array<std::uint16_t,16>;

bool fullRank15(Matrix rows) {
    for(unsigned pivot = 0; pivot < 15; ++pivot) {
        unsigned selected = pivot;
        while(selected < 15 && !((rows[selected] >> pivot) & 1U)) ++selected;
        if(selected == 15) return false;
        std::swap(rows[pivot], rows[selected]);
        for(unsigned row = pivot + 1; row < 15; ++row)
            if((rows[row] >> pivot) & 1U) rows[row] ^= rows[pivot];
    }
    return true;
}
template<class Container>
void shuffle(Container& values, Words& words, const std::vector<Divisor>& divisors) {
    std::iota(values.begin(), values.end(), 0U);
    for(std::size_t i = values.size(); i > 1; --i)
        std::swap(values[i - 1], values[divisors[i].sample(words, i)]);
}
void prepareRoute(Plan& p, std::uint64_t seed) {
    Words words(seed);
    std::vector<Divisor> divisors(std::max(p.groups, std::size_t{64}) + 1);
    for(std::size_t i = 2; i < divisors.size(); ++i) divisors[i] = Divisor(i);
    std::vector<std::array<std::uint32_t,64>> columns(p.groups);
    for(auto& column : columns) shuffle(column, words, divisors);
    std::vector<std::uint32_t> destinations(p.groups);
    p.route.resize(p.n / 4);
    p.routeBytes.resize(p.n / 4);
    for(unsigned region = 0; region < 64; ++region) {
        shuffle(destinations, words, divisors);
        for(std::size_t group = 0; group < p.groups; ++group) {
            const auto i = region * p.groups + destinations[group];
            const auto offset = std::uint32_t(group * groupStride + 4 * columns[group][region]);
            p.route[i] = offset;
            p.routeBytes[i] = 16 * offset;
        }
    }
}
void compact(const Matrix& rows, std::uint64_t* reverse, std::uint64_t* forward) {
    Matrix r{}, f{};
    for(unsigned i = 0; i < 16; ++i)
        for(unsigned j = 0; j < 16; ++j) {
            r[i] |= std::uint16_t(((rows[permutation[i]] >> permutation[j]) & 1U) << j);
            f[i] |= std::uint16_t(((rows[permutation[j]] >> permutation[i]) & 1U) << j);
        }
    for(unsigned o = 0; o < 2; ++o)
        for(unsigned i = 0; i < 2; ++i) {
            std::uint64_t a = 0, b = 0;
            for(unsigned j = 0; j < 8; ++j) {
                a |= std::uint64_t((r[8*o+j] >> (8*i)) & 255U) << (8*j);
                b |= std::uint64_t((f[8*o+j] >> (8*i)) & 255U) << (8*j);
            }
            reverse[2*o+i] = a;
            forward[2*o+i] = b;
        }
}
void prepareInner(Plan& p, std::uint64_t seed) {
    Words words(seed ^ 0x3f625a92ULL);
    p.reverseMatrices.resize(p.n / 64);
    p.updates.resize(4 * p.reverseMatrices.size());
    p.forwardUpdates.resize(p.updates.size());
    for(std::size_t epoch = 0; epoch < p.reverseMatrices.size(); ++epoch) {
        auto& rows = p.reverseMatrices[epoch];
        do {
            for(unsigned row = 0; row < 15; ++row) rows[row] = std::uint16_t(words() & 0x7fffU);
            rows[15] = 0x8000U;
        } while(!fullRank15(rows));
        compact(rows, p.updates.data() + 4*epoch, p.forwardUpdates.data() + 4*epoch);
    }
}
std::uint8_t multiply8(std::uint8_t a, std::uint8_t b) {
    std::uint8_t result = 0;
    for(unsigned bit = 0; bit < 8; ++bit) {
        if(b & 1U) result ^= a;
        a = std::uint8_t((unsigned(a) << 1) ^ ((a & 128U) ? 0x11bU : 0U));
        b >>= 1;
    }
    return result;
}
std::uint16_t multiply16(std::uint16_t a, std::uint16_t b) {
    const auto a0 = std::uint8_t(a), a1 = std::uint8_t(a >> 8);
    const auto b0 = std::uint8_t(b), b1 = std::uint8_t(b >> 8);
    const auto lo = multiply8(a0,b0) ^ multiply8(0x20,multiply8(a1,b1));
    const auto hi = multiply8(a0,b1) ^ multiply8(a1,b0) ^ multiply8(a1,b1);
    return std::uint16_t(lo | (unsigned(hi) << 8));
}
std::uint64_t adjoint8(std::uint8_t scalar) {
    std::uint64_t matrix = 0;
    for(unsigned bit = 0; bit < 8; ++bit)
        matrix |= std::uint64_t(multiply8(scalar, std::uint8_t(1U << bit))) << (8*(7-bit));
    return matrix;
}
void prepareOuter(Plan& p, std::uint64_t seed) {
    Words words(seed ^ 0x75a1dc09ULL);
    p.outerMatrices16.resize(16*p.groups);
    p.outerField.resize(16*p.groups);
    p.outerAdjoint.resize(16*p.groups);
    for(std::size_t symbol = 0; symbol < 16*p.groups; ++symbol) {
        std::uint16_t scalar;
        do { scalar = std::uint16_t(words()); } while(!scalar);
        const auto c0 = std::uint8_t(scalar), c1 = std::uint8_t(scalar >> 8);
        p.outerField[symbol] = {c0, multiply8(0x20,c1), std::uint8_t(c0 ^ c1)};
        for(unsigned j = 0; j < 3; ++j)
            p.outerAdjoint[symbol][j] = adjoint8(p.outerField[symbol][j]);
        auto& rows = p.outerMatrices16[symbol];
        for(unsigned column = 0; column < 16; ++column) {
            const auto image = multiply16(scalar, std::uint16_t(1U << column));
            const auto physical = 4*(column%4) + column/4;
            for(unsigned row = 0; row < 16; ++row)
                rows[row] |= std::uint16_t(((image >> row) & 1U) << physical);
        }
    }
}
}
Plan::Plan(std::size_t size, std::uint64_t routeSeed, std::uint64_t innerSeed)
    : k(size), n(2*size), groups(size/128) {
    if(k != 65536) throw std::invalid_argument("paired-s15 profile requires K=65536");
    prepareRoute(*this, routeSeed);
    prepareInner(*this, innerSeed);
    prepareOuter(*this, routeSeed);
}
std::size_t Plan::setupBytes() const noexcept {
    return route.capacity()*sizeof(route[0]) + routeBytes.capacity()*sizeof(routeBytes[0])
        + reverseMatrices.capacity()*sizeof(reverseMatrices[0])
        + outerMatrices16.capacity()*sizeof(outerMatrices16[0])
        + updates.capacity()*sizeof(updates[0]) + forwardUpdates.capacity()*sizeof(forwardUpdates[0])
        + outerField.capacity()*sizeof(outerField[0]) + outerAdjoint.capacity()*sizeof(outerAdjoint[0]);
}
}
