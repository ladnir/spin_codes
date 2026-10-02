#include "RsWide64.h"
#include "Tower64ByteRandomizer.h"
#include "../../../../../spin/src/kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spin::research::rswide64 {
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
    for(unsigned column=0; column<Width; ++column) {
        unsigned pivot=column;
        while(pivot<Width && !(rows[pivot]&(Word{1}<<column))) ++pivot;
        if(pivot==Width) return false;
        std::swap(rows[pivot],rows[column]);
        for(unsigned j=column+1;j<Width;++j)
            if(rows[j]&(Word{1}<<column))rows[j]^=rows[column];
    }
    return true;
}
template<class Word, std::size_t Width>
std::array<Word,Width> randomTranspose(Words& words) {
    std::array<Word,Width> forward{}, reverse{};
    do { for(auto& row:forward)row=static_cast<Word>(words()); } while(!fullRank(forward));
    for(unsigned i=0;i<Width;++i)
        for(unsigned j=0;j<Width;++j)
            reverse[j]|=Word(((forward[i]>>j)&1U)<<i);
    return reverse;
}
template<class Container>
void shuffle(Container& values, Words& words, const std::vector<Divisor>& divisors) {
    std::iota(values.begin(),values.end(),0U);
    for(std::size_t i=values.size();i>1;--i)
        std::swap(values[i-1],values[divisors[i].sample(words,i)]);
}
}

Plan::Plan(std::size_t messageBits,std::uint64_t routeSeed,std::uint64_t innerSeed) {
    constexpr auto maxGroups=std::numeric_limits<std::uint32_t>::max()/groupStride;
    if(!messageBits || messageBits%8192 || messageBits/512>maxGroups)
        throw std::invalid_argument("RS16x16 K must be a positive multiple of 8192 within 32-bit routing range");
    k=messageBits;n=2*k;groups=k/512;
    Words routeWords(routeSeed);
    std::vector<Divisor> divisors(std::max(groups,std::size_t{256})+1);
    for(std::size_t i=2;i<divisors.size();++i)divisors[i]=Divisor(i);
    std::vector<std::array<std::uint32_t,256>> columns(groups);
    for(auto& permutation:columns)shuffle(permutation,routeWords,divisors);
    std::vector<std::uint32_t> destinations(groups);
    route.resize(n/4);
    for(unsigned region=0;region<256;++region) {
        shuffle(destinations,routeWords,divisors);
        for(std::size_t group=0;group<groups;++group)
            route[region*groups+destinations[group]]=static_cast<std::uint32_t>(group*groupStride+4*columns[group][region]);
    }
    Words innerWords(innerSeed^0x3f625a92ULL);
    reverseMatrices.resize(n/64);denseUpdates.resize(n/64);
    for(std::size_t epoch=0;epoch<n/64;++epoch) {
        const auto reverse=randomTranspose<std::uint16_t,16>(innerWords);
        reverseMatrices[epoch]=reverse;
        std::array<std::uint16_t,16> transformed{};
        for(unsigned output=0;output<16;++output) {
            unsigned originalMask=0;
            for(unsigned mask=basis[output];mask;mask&=mask-1)originalMask^=reverse[std::countr_zero(mask)];
            for(unsigned mask=originalMask;mask;mask&=mask-1)transformed[output]^=inverseBasis[std::countr_zero(mask)];
        }
        for(unsigned out=0;out<2;++out)
            for(unsigned in=0;in<2;++in)
                for(unsigned j=0;j<8;++j)
                    denseUpdates[epoch].matrix[2*out+in]|=
                        std::uint64_t((transformed[8*out+j]>>(8*in))&255U)<<(8*(7-j));
    }
    Words outerWords(routeSeed^0x75a1dc09ULL);
    outerMatrices.resize(16*groups);coefficientBytes.resize(432*groups);
    for(std::size_t group=0;group<groups;++group)
        for(unsigned symbol=0;symbol<16;++symbol) {
            const auto scalar=rs::tower64::sampleNonzero(outerWords);
            outerMatrices[16*group+symbol]=rs::tower64byte::multiplyRows(scalar);
            const auto packed=rs::tower64byte::coefficients(scalar);
            std::copy(packed.begin(),packed.end(),coefficientBytes.data()+432*group+27*symbol);
        }
}
}
