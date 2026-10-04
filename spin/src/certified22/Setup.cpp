#include "Plan.h"
#include "productionMap.h"
#include "../kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spin::detail::certified22 {
namespace {
using kernel::setup::Words;
using kernel::setup::Divisor;
using Matrix = std::array<std::uint32_t,22>;

bool fullRank(Matrix rows) {
    for(unsigned column=0; column<22; ++column) {
        unsigned pivot=column;
        while(pivot<22 && !((rows[pivot]>>column)&1U)) ++pivot;
        if(pivot==22) return false;
        std::swap(rows[pivot],rows[column]);
        for(unsigned row=column+1; row<22; ++row)
            if((rows[row]>>column)&1U) rows[row]^=rows[column];
    }
    return true;
}
Matrix multiply(const Matrix& left,const Matrix& right) {
    Matrix result{};
    for(unsigned row=0; row<22; ++row)
        for(auto mask=left[row]; mask; mask&=mask-1)
            result[row]^=right[std::countr_zero(mask)];
    return result;
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

void prepareInner(Plan& p,std::uint64_t seed) {
    const auto epochs=p.n/64;
    p.reverseMatrices.resize(epochs);
    p.updates.resize(9*epochs);
    Words words(seed^0x3f625a92ULL);
    for(std::size_t epoch=0; epoch<epochs; ++epoch) {
        auto& literal=p.reverseMatrices[epoch];
        // Exact uniform-GL22 rejection sampler from the certified setup.
        do {
            for(auto& row:literal) row=static_cast<std::uint32_t>(words())&0x3fffffU;
        } while(!fullRank(literal));
        const auto transformed=multiply(multiply(map::basis,literal),map::inverseBasis);
        for(unsigned out=0; out<3; ++out)
            for(unsigned in=0; in<3; ++in) {
                std::uint64_t packed=0;
                for(unsigned j=0; j<8; ++j) {
                    const auto row=8*out+j;
                    if(row<22) packed|=std::uint64_t((transformed[row]>>(8*in))&255U)<<(8*j);
                }
                p.updates[9*epoch+3*out+in]=packed;
            }
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

void prepareOuter(Plan& p, std::uint64_t seed) {
    Words words(seed ^ 0x75a1dc09ULL);
    p.outerMatrices16.resize(16*p.groups);
    p.outerField.resize(16*p.groups);
    for(std::size_t symbol = 0; symbol < 16*p.groups; ++symbol) {
        std::uint16_t scalar;
        do { scalar = std::uint16_t(words()); } while(!scalar);
        const auto c0 = std::uint8_t(scalar), c1 = std::uint8_t(scalar >> 8);
        p.outerField[symbol] = {c0, multiply8(0x20,c1), std::uint8_t(c0 ^ c1)};
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
bool validMessageSize(std::size_t k) noexcept {
    if(!k || k%128) return false;
    const auto groups=k/128;
    // Each route points to a complete 64-byte packet. Also reserve workspace
    // for the widest supported records without overflowing size_t.
    constexpr auto maxByteRoute=std::numeric_limits<std::uint32_t>::max();
    constexpr auto maxSize=std::numeric_limits<std::size_t>::max();
    return groups<=maxByteRoute/(groupStride*sizeof(Block))
        && groups<=maxSize/(groupStride*sizeof(Block)*4)
        && k<=maxSize/(2*sizeof(Block)*4);
}
Plan::Plan(std::size_t size,std::uint64_t routeSeed,std::uint64_t innerSeed,bool withOuter) {
    if(!validMessageSize(size))
        throw std::invalid_argument("certified-s22 K must be a positive multiple of 128 within the byte-route and workspace limits");
    k=size;
    n=2*size;
    groups=size/128;
    prepareRoute(*this,routeSeed);
    prepareInner(*this,innerSeed);
    if(withOuter) prepareOuter(*this,routeSeed);
}
std::size_t Plan::setupBytes() const noexcept {
    return route.capacity()*sizeof(route[0])+routeBytes.capacity()*sizeof(routeBytes[0])
        +reverseMatrices.capacity()*sizeof(reverseMatrices[0])
        +updates.capacity()*sizeof(updates[0])
        +outerMatrices16.capacity()*sizeof(outerMatrices16[0])
        +outerField.capacity()*sizeof(outerField[0]);
}
}
