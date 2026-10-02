#include "RsBorder.h"
#include "../../../../../spin/src/kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <stdexcept>
namespace spin::research::rsborder {
namespace {
constexpr std::uint16_t basis[16]={0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x3500,0x1180,0x5980,0x2c80,0xcf80,0x9400,0xe400,0xe00,0xd980};
constexpr std::uint16_t inverseBasis[16]={0x1,0x2,0x4,0x8,0x10,0x20,0x40,0xa500,0x6980,0x1a00,0x5f80,0x580,0xcd80,0xfb00,0x680,0x8200};
bool fullRank(std::array<std::uint32_t,20> rows,unsigned bits) {
    for(unsigned c=0;c<bits;++c) {
        unsigned p=c;while(p<bits && !(rows[p]&(1U<<c)))++p;
        if(p==bits)return false;std::swap(rows[p],rows[c]);
        for(unsigned j=c+1;j<bits;++j)if(rows[j]&(1U<<c))rows[j]^=rows[c];
    }
    return true;
}
}
Plan::Plan(std::size_t k,std::uint64_t routeSeed,std::uint64_t innerSeed,unsigned bits,rswide::Randomizer randomizer)
    :outer(k,routeSeed,innerSeed,randomizer),stateBits(bits) {
    if(bits<16 || bits>20)throw std::invalid_argument("border prototype supports s16 through s20 only");
    if(bits==16)return;
    detail::kernel::setup::Words words(innerSeed^0x3f625a92ULL);
    reverseMatrices.resize(outer.n/64);updates.resize(outer.n/64);
    for(std::size_t epoch=0;epoch<outer.n/64;++epoch) {
        std::array<std::uint32_t,20> forward{},reverse{},transformed{};
        do {for(unsigned i=0;i<bits;++i)forward[i]=std::uint32_t(words())&((1U<<bits)-1);} while(!fullRank(forward,bits));
        for(unsigned i=0;i<bits;++i)for(unsigned j=0;j<bits;++j)reverse[j]|=((forward[i]>>j)&1U)<<i;
        reverseMatrices[epoch]=reverse;
        // diag(B,I) * M^T * diag(B^-1,I), retaining the old16-coordinate basis.
        for(unsigned out=0;out<bits;++out) {
            unsigned mask=0;
            if(out<16)for(unsigned selected=basis[out];selected;selected&=selected-1)mask^=reverse[std::countr_zero(selected)];
            else mask=reverse[out];
            transformed[out]=mask&~65535U;
            for(unsigned selected=mask&65535U;selected;selected&=selected-1)transformed[out]^=inverseBasis[std::countr_zero(selected)];
        }
        auto& update=updates[epoch];
        for(unsigned out=0;out<2;++out)for(unsigned in=0;in<2;++in)for(unsigned j=0;j<8;++j)
            update.base.matrix[2*out+in]|=std::uint64_t((transformed[8*out+j]>>(8*in))&255U)<<(8*(7-j));
        for(unsigned out=0;out<2;++out)for(unsigned j=0;j<8;++j)
            update.upper[out]|=std::uint64_t(transformed[8*out+j]>>16)<<(8*(7-j));
        for(unsigned j=0;j<bits-16;++j) {
            update.lower[0]|=std::uint64_t(transformed[16+j]&255U)<<(8*(7-j));
            update.lower[1]|=std::uint64_t((transformed[16+j]>>8)&255U)<<(8*(7-j));
            update.corner|=std::uint64_t(transformed[16+j]>>16)<<(8*(7-j));
        }
    }
}
}
