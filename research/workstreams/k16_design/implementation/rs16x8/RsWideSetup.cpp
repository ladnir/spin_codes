#include "RsWide.h"
#include "Tower32Randomizer.h"
#include "Tower32ByteRandomizer.h"
#include "../../../../../spin/src/kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace spin::research::rswide {
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

Plan::Plan(std::size_t messageBits,std::uint64_t routeSeed,std::uint64_t innerSeed,Randomizer selected,InnerRandomizer innerSelected) {
    constexpr auto maxGroups=std::numeric_limits<std::uint32_t>::max()/groupStride;
    if(!messageBits || messageBits%4096 || messageBits/256>maxGroups)
        throw std::invalid_argument("RS16x8 K must be a positive multiple of4096 within32-bit routing range");
    if(selected!=Randomizer::DenseGl32 && selected!=Randomizer::TowerField32 && selected!=Randomizer::TowerByte32)
        throw std::invalid_argument("unknown outer randomizer");
    if(innerSelected!=InnerRandomizer::DenseGl16 && innerSelected!=InnerRandomizer::TowerByte16)
        throw std::invalid_argument("unknown inner randomizer");
    randomizer=selected;innerRandomizer=innerSelected;k=messageBits;n=2*k;groups=k/256;
    Words routeWords(routeSeed);
    std::vector<Divisor> divisors(std::max(groups,std::size_t{128})+1);
    for(std::size_t i=2;i<divisors.size();++i)divisors[i]=Divisor(i);
    std::vector<std::array<std::uint32_t,128>> columns(groups);
    for(auto& permutation:columns)shuffle(permutation,routeWords,divisors);
    std::vector<std::uint32_t> destinations(groups);
    route.resize(n/4);
    for(unsigned region=0;region<128;++region) {
        shuffle(destinations,routeWords,divisors);
        for(std::size_t group=0;group<groups;++group)
            route[region*groups+destinations[group]]=static_cast<std::uint32_t>(group*groupStride+4*columns[group][region]);
    }
    Words innerWords(innerSeed^0x3f625a92ULL);
    reverseMatrices.resize(n/64);denseUpdates.resize(n/64);
    if(innerRandomizer==InnerRandomizer::TowerByte16)fieldUpdates.resize(n/64);
    for(std::size_t epoch=0;epoch<n/64;++epoch) {
        std::array<std::uint16_t,16> reverse{};
        if(innerRandomizer==InnerRandomizer::DenseGl16) {
            reverse=randomTranspose<std::uint16_t,16>(innerWords);
        } else {
            std::uint16_t scalar;
            do {scalar=std::uint16_t(innerWords());} while(!scalar);
            rs::tower32byte::coefficients16(scalar,fieldUpdates[epoch].data());
            std::array<std::uint16_t,16> packedRows{};
            for(unsigned column=0;column<16;++column) {
                const auto image=rs::tower32::multiply16(scalar,std::uint16_t(1U<<column));
                for(unsigned row=0;row<16;++row)packedRows[row]|=std::uint16_t(((image>>row)&1U)<<column);
            }
            // The fast state is P times the literal state. Choose its reverse
            // update T=M_scalar; the literal reverse matrix is P^-1 T P.
            for(unsigned row=0;row<16;++row) {
                unsigned mask=0;
                for(unsigned bits=inverseBasis[row];bits;bits&=bits-1)mask^=packedRows[std::countr_zero(bits)];
                for(;mask;mask&=mask-1)reverse[row]^=basis[std::countr_zero(mask)];
            }
        }
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
    outerMatrices.resize(16*groups);compactStorage.resize((16*coefficientBytesPerSymbol()*groups+7)/8+7);
    auto* coefficients=const_cast<std::uint64_t*>(compactCoefficients());
    for(std::size_t group=0;group<groups;++group)
        for(unsigned symbol=0;symbol<16;++symbol) {
            auto& matrix=outerMatrices[16*group+symbol];
            if(randomizer==Randomizer::TowerByte32) {
                const auto scalar=rs::tower32::sampleNonzero(outerWords);
                matrix=rs::tower32byte::multiplyRows(scalar);
                const auto packed=rs::tower32byte::coefficients(scalar);
                std::copy(packed.begin(),packed.end(),reinterpret_cast<std::uint8_t*>(coefficients)+144*group+9*symbol);
                continue;
            }
            if(randomizer==Randomizer::TowerField32) {
                const auto scalar=rs::tower32::sampleNonzero(outerWords);
                matrix=rs::tower32::adjointRows(scalar);
                const auto packed=rs::tower32::adjointCoefficients(scalar);
                std::copy(packed.begin(),packed.end(),coefficients+144*group+9*symbol);
                continue;
            }
            matrix=randomTranspose<std::uint32_t,32>(outerWords);
            for(unsigned out=0;out<4;++out)
                for(unsigned in=0;in<4;++in) {
                    std::uint64_t value=0;
                    for(unsigned j=0;j<8;++j)
                        value|=std::uint64_t((matrix[8*out+j]>>(8*in))&255U)<<(8*(7-j));
                    coefficients[256*group+16*symbol+4*out+in]=value;
                }
        }
}
const std::uint64_t* Plan::compactCoefficients() const noexcept {
    const auto address=reinterpret_cast<std::uintptr_t>(compactStorage.data());
    return compactStorage.data()+(((64-(address&63U))&63U)/8);
}
}
