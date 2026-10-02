#include "PacketPlan.h"
#include "PacketRsField.h"
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
// Change of basis for the retained 16-coordinate core, extended by identity
// on the four appended coordinates. Scalar fallback uses the literal basis.
constexpr std::uint16_t basis[16]={0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x3500,0x1180,0x5980,0x2c80,0xcf80,0x9400,0xe400,0xe00,0xd980};
constexpr std::uint16_t inverseBasis[16]={0x1,0x2,0x4,0x8,0x10,0x20,0x40,0xa500,0x6980,0x1a00,0x5f80,0x580,0xcd80,0xfb00,0x680,0x8200};

bool fullRank(std::array<std::uint32_t,20> rows) {
    for(unsigned c=0;c<20;++c) {
        unsigned p=c;while(p<20 && !(rows[p]&(1U<<c)))++p;
        if(p==20)return false;std::swap(rows[p],rows[c]);
        for(unsigned j=c+1;j<20;++j)if(rows[j]&(1U<<c))rows[j]^=rows[c];
    }
    return true;
}

template<class Container>
void shuffle(Container& values,Words& words,const std::vector<Divisor>& divisors) {
    std::iota(values.begin(),values.end(),0U);
    for(std::size_t i=values.size();i>1;--i)
        std::swap(values[i-1],values[divisors[i].sample(words,i)]);
}

void prepareRoute(Plan& plan,std::uint64_t seed) {
    Words words(seed);
    std::vector<Divisor> divisors(std::max(plan.groups,std::size_t{128})+1);
    for(std::size_t i=2;i<divisors.size();++i)divisors[i]=Divisor(i);
    std::vector<std::array<std::uint32_t,128>> columns(plan.groups);
    for(auto& permutation:columns)shuffle(permutation,words,divisors);
    std::vector<std::uint32_t> destinations(plan.groups);
    plan.route.resize(plan.n/4);
    for(unsigned region=0;region<128;++region) {
        shuffle(destinations,words,divisors);
        for(std::size_t group=0;group<plan.groups;++group)
            plan.route[region*plan.groups+destinations[group]]=
                static_cast<std::uint32_t>(group*groupStride+4*columns[group][region]);
    }
}

void prepareUpdates(Plan& plan,std::uint64_t seed) {
    Words words(seed^0x3f625a92ULL);
    plan.reverseMatrices.resize(plan.n/64);plan.updates.resize(plan.n/64);
    for(std::size_t epoch=0;epoch<plan.n/64;++epoch) {
        std::array<std::uint32_t,20> forward{},reverse{},transformed{};
        do {for(unsigned i=0;i<20;++i)forward[i]=std::uint32_t(words())&((1U<<20)-1);} while(!fullRank(forward));
        for(unsigned i=0;i<20;++i)for(unsigned j=0;j<20;++j)reverse[j]|=((forward[i]>>j)&1U)<<i;
        plan.reverseMatrices[epoch]=reverse;
        // diag(B,I) * M^T * diag(B^-1,I), preserving the measured update.
        for(unsigned out=0;out<20;++out) {
            unsigned mask=0;
            if(out<16)for(unsigned selected=basis[out];selected;selected&=selected-1)mask^=reverse[std::countr_zero(selected)];
            else mask=reverse[out];
            transformed[out]=mask&~65535U;
            for(unsigned selected=mask&65535U;selected;selected&=selected-1)transformed[out]^=inverseBasis[std::countr_zero(selected)];
        }
        auto& update=plan.updates[epoch];
        for(unsigned out=0;out<2;++out)for(unsigned in=0;in<2;++in)for(unsigned j=0;j<8;++j)
            update.base.matrix[2*out+in]|=std::uint64_t((transformed[8*out+j]>>(8*in))&255U)<<(8*(7-j));
        for(unsigned out=0;out<2;++out)for(unsigned j=0;j<8;++j)
            update.upper[out]|=std::uint64_t(transformed[8*out+j]>>16)<<(8*(7-j));
        for(unsigned j=0;j<4;++j) {
            update.lower[0]|=std::uint64_t(transformed[16+j]&255U)<<(8*(7-j));
            update.lower[1]|=std::uint64_t((transformed[16+j]>>8)&255U)<<(8*(7-j));
            update.corner|=std::uint64_t(transformed[16+j]>>16)<<(8*(7-j));
        }
    }
}

void prepareOuter(Plan& plan,std::uint64_t seed) {
    Words words(seed^0x75a1dc09ULL);
    plan.outerMatrices.resize(16*plan.groups);
    plan.compactStorage.resize((144*plan.groups+7)/8+7);
    auto* coefficients=reinterpret_cast<std::uint8_t*>(const_cast<std::uint64_t*>(plan.compactCoefficients()));
    for(std::size_t group=0;group<plan.groups;++group)
        for(unsigned symbol=0;symbol<16;++symbol) {
            const auto scalar=tower32::sampleNonzero(words);
            plan.outerMatrices[16*group+symbol]=tower32byte::multiplyRows(scalar);
            const auto packed=tower32byte::coefficients(scalar);
            std::copy(packed.begin(),packed.end(),coefficients+144*group+9*symbol);
        }
}
}

bool validMessageSize(std::size_t k) noexcept {
    // Only structural/representation limits: N is automatically divisible by
    // 64, and region boundaries need not coincide with inner epochs. Padded
    // packet offsets fit uint32_t; every vector and scratch byte count fits
    // size_t. There is no performance threshold or certificate-policy limit.
    constexpr auto maxRouteGroups=std::numeric_limits<std::uint32_t>::max()/groupStride;
    constexpr auto maxByteGroups=std::numeric_limits<std::size_t>::max()/(groupStride*sizeof(Block));
    return k && k%256==0 && k/256<=maxRouteGroups && k/256<=maxByteGroups;
}

Plan::Plan(std::size_t messageBits,std::uint64_t seed):Plan(messageBits,seed,seed) {}
Plan::Plan(std::size_t messageBits,std::uint64_t routeSeed,std::uint64_t innerSeed) {
    if(!validMessageSize(messageBits))
        throw std::invalid_argument("RS packet SPIN K must be a positive multiple of 256 within padded routing and byte-size limits");
    k=messageBits;n=2*k;groups=k/256;
    prepareRoute(*this,routeSeed);
    prepareUpdates(*this,innerSeed);
    prepareOuter(*this,routeSeed);
}

const std::uint64_t* Plan::compactCoefficients() const noexcept {
    const auto address=reinterpret_cast<std::uintptr_t>(compactStorage.data());
    return compactStorage.data()+(((64-(address&63U))&63U)/8);
}

std::size_t Plan::setupBytes() const noexcept {
    return route.capacity()*sizeof(route[0])
        +reverseMatrices.capacity()*sizeof(reverseMatrices[0])
        +outerMatrices.capacity()*sizeof(outerMatrices[0])
        +compactStorage.capacity()*sizeof(compactStorage[0])
        +updates.capacity()*sizeof(updates[0]);
}
}
