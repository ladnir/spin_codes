#include "Packet8Wide24.h"
#include "../../src/kernels/SetupRandom.h"
#include "../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"
#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spin::research::packet8wide24 {
namespace {
using detail::kernel::setup::Words;
using detail::kernel::setup::Divisor;

template<class Container>
void shuffle(Container& values, Words& random, const std::vector<Divisor>& divisors) {
    std::iota(values.begin(), values.end(), 0U);
    for(std::size_t i=values.size(); i>1; --i)
        std::swap(values[i-1], values[divisors[i].sample(random,i)]);
}
}

Plan::Plan(std::size_t messageBits, std::uint64_t seed) {
    // Routes hold block offsets, not byte offsets. Every group starts on a
    // 64-byte boundary because groupStride is a multiple of four blocks.
    constexpr auto maxRouteGroups=std::numeric_limits<std::uint32_t>::max()/groupStride;
    constexpr auto maxStorageGroups=std::numeric_limits<std::size_t>::max()/(groupStride*sizeof(Block));
    constexpr auto maxGroups=std::min(std::size_t(maxRouteGroups),maxStorageGroups);
    if(!messageBits || messageBits%2048 || messageBits/256>maxGroups)
        throw std::invalid_argument("packet8wide24 K must be a positive multiple of2048 within the routing/storage range");
    k=messageBits; n=2*k; groups=k/256;

    Words outer(seed^0x7d656cf648f0ab33ULL);
    outerRows.resize(16*groups);
    outerCoefficients.resize(16*groups);
    for(std::size_t symbol=0; symbol<16*groups; ++symbol) {
        const auto scalar=rs::tower32::sampleNonzero(outer);
        outerRows[symbol]=rs::tower32byte::multiplyRows(scalar);
        outerCoefficients[symbol]=rs::tower32byte::coefficients(scalar);
    }

    Words routing(seed^0xa13df7412678bb08ULL);
    std::vector<Divisor> divisors(std::max(groups,std::size_t(64))+1);
    for(std::size_t i=2; i<divisors.size(); ++i) divisors[i]=Divisor(i);
    std::vector<std::array<std::uint32_t,64>> columns(groups);
    for(auto& permutation:columns) shuffle(permutation,routing,divisors);
    std::vector<std::uint32_t> positions(groups);
    route.resize(n/8);
    for(unsigned region=0; region<64; ++region) {
        shuffle(positions,routing,divisors);
        for(std::size_t group=0; group<groups; ++group)
            route[region*groups+positions[group]]=std::uint32_t(group*groupStride+8*columns[group][region]);
    }

    Words inner(seed^0xb862dec09a47c103ULL);
    updates.resize(n/64);
    for(auto& update:updates) {
        do { update.scalar=std::uint32_t(inner())&0xffffffU; } while(!update.scalar);
        const auto r0=std::uint8_t(update.scalar);
        const auto r1=std::uint8_t(update.scalar>>8);
        const auto r2=std::uint8_t(update.scalar>>16);
        update.coefficients={r0,r1,r2,std::uint8_t(r0^r1),std::uint8_t(r0^r2),std::uint8_t(r1^r2)};
    }
}
}
