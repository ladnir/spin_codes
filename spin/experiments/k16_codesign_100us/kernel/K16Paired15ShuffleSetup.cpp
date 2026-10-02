#include "K16Paired15Shuffle.h"
namespace spin::research::k16codesign {
void customizePaired15Shuffle(rs::Plan& plan,std::uint64_t seed,PairedTables& tables,PairedOptimizedTables& optimized) {
    customizePaired15(plan,seed,tables,optimized);
    constexpr unsigned permutation[16]={0,1,2,7,3,4,5,6,8,15,10,12,11,13,14,9};
    for(std::size_t epoch=0;epoch<plan.n/64;++epoch) {
        std::array<std::uint16_t,16> rows{};
        for(unsigned r=0;r<16;++r)
            for(unsigned c=0;c<16;++c)
                rows[r]|=std::uint16_t(((plan.reverseMatrices[epoch][permutation[r]]>>permutation[c])&1U)<<c);
        for(unsigned out=0;out<2;++out)
            for(unsigned in=0;in<2;++in) {
                std::uint64_t value=0;
                for(unsigned j=0;j<8;++j)value|=std::uint64_t((rows[8*out+j]>>(8*in))&255U)<<(8*j);
                tables.updates[4*epoch+2*out+in]=value;
            }
    }
}
}
