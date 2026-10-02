#include "PacketForward.h"
#include "PacketRsField.h"
#include <bit>
namespace spin::detail::packet {
ForwardPlan::ForwardPlan(const Plan& p):updates(p.n/64),outer(16*p.groups) {
    inverseRoute.resize(p.scratchBlocks()/4);
    for(std::size_t i=0;i<p.route.size();++i)inverseRoute[p.route[i]/4]=4*i;
    constexpr std::uint16_t basis[16]={0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x3500,0x1180,0x5980,0x2c80,0xcf80,0x9400,0xe400,0xe00,0xd980};
    constexpr std::uint16_t inverse[16]={0x1,0x2,0x4,0x8,0x10,0x20,0x40,0xa500,0x6980,0x1a00,0x5f80,0x580,0xcd80,0xfb00,0x680,0x8200};
    for(std::size_t e=0;e<updates.size();++e) {
        std::uint32_t rows[20]{},mapped[20]{};
        for(unsigned i=0;i<20;++i)for(unsigned j=0;j<20;++j)
            rows[i]|=((p.reverseMatrices[e][j]>>i)&1U)<<j;
        for(unsigned i=0;i<20;++i) {
            unsigned mask=0;
            for(auto m=i<16?unsigned(basis[i]):1U<<i;m;m&=m-1)mask^=rows[std::countr_zero(m)];
            mapped[i]=mask&~65535U;
            for(auto m=mask&65535U;m;m&=m-1)mapped[i]^=inverse[std::countr_zero(m)];
        }
        auto& u=updates[e];
        for(unsigned a=0;a<2;++a)for(unsigned b=0;b<2;++b)for(unsigned j=0;j<8;++j)
            u.base.matrix[2*a+b]|=std::uint64_t((mapped[8*a+j]>>(8*b))&255)<<(8*(7-j));
        for(unsigned a=0;a<2;++a)for(unsigned j=0;j<8;++j)
            u.upper[a]|=std::uint64_t(mapped[8*a+j]>>16)<<(8*(7-j));
        for(unsigned j=0;j<4;++j) {
            u.lower[0]|=std::uint64_t(mapped[16+j]&255)<<(8*(7-j));
            u.lower[1]|=std::uint64_t((mapped[16+j]>>8)&255)<<(8*(7-j));
            u.corner|=std::uint64_t(mapped[16+j]>>16)<<(8*(7-j));
        }
    }
    const auto* coeff=reinterpret_cast<const std::uint8_t*>(p.compactCoefficients());
    for(std::size_t s=0;s<outer.size();++s)for(unsigned j=0;j<9;++j) {
        // Row i of multiplication's transpose is column i of multiplication.
        for(unsigned i=0;i<8;++i)outer[s][j]|=std::uint64_t(tower32::multiply8(coeff[9*s+j],1U<<i))<<(8*(7-i));
    }
}
}
