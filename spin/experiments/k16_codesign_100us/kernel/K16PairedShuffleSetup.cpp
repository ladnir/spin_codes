#include "K16PairedShuffle.h"
namespace spin::research::k16codesign {
namespace {
constexpr unsigned permutation[16]={0,1,2,7,3,4,5,6,8,10,11,13,12,14,15,9};
std::array<std::uint16_t,16> toFast(const std::array<std::uint16_t,16>& rows) {
    std::array<std::uint16_t,16> result{};
    for(unsigned r=0;r<16;++r)
        for(unsigned c=0;c<16;++c)result[r]|=std::uint16_t(((rows[permutation[r]]>>permutation[c])&1U)<<c);
    return result;
}
std::array<std::uint16_t,16> toLiteral(const std::array<std::uint16_t,16>& rows) {
    std::array<std::uint16_t,16> result{};
    for(unsigned r=0;r<16;++r)
        for(unsigned c=0;c<16;++c)result[permutation[r]]|=std::uint16_t(((rows[r]>>c)&1U)<<permutation[c]);
    return result;
}
}
void preparePairedShuffle(const rs::Plan& plan,PairedTables& tables,PairedOptimizedTables& opt) {
    preparePairedOptimized(plan,opt);tables.updates.resize(4*(plan.n/64));
    for(std::size_t epoch=0;epoch<plan.n/64;++epoch) {
        const auto rows=toFast(plan.reverseMatrices[epoch]);
        for(unsigned out=0;out<2;++out)
            for(unsigned in=0;in<2;++in) {
                std::uint64_t value=0;
                for(unsigned j=0;j<8;++j)value|=std::uint64_t((rows[8*out+j]>>(8*in))&255U)<<(8*j);
                tables.updates[4*epoch+2*out+in]=value;
            }
    }
}
void customizePairedField16Shuffle(rs::Plan& plan,std::uint64_t seed,PairedTables& tables,PairedOptimizedTables& opt) {
    customizePairedField16(plan,seed,tables,opt);
    for(auto& rows:plan.reverseMatrices)rows=toLiteral(rows);
}
}
