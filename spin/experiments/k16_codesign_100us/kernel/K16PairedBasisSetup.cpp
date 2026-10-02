#include "K16PairedBasis.h"
namespace spin::research::k16codesign {
namespace {
constexpr unsigned p(unsigned j) {return j==3?7:j==7?3:j;}
std::array<std::uint16_t,16> conjugate(const std::array<std::uint16_t,16>& rows) {
    std::array<std::uint16_t,16> result{};
    for(unsigned r=0;r<16;++r)
        for(unsigned c=0;c<16;++c)result[r]|=std::uint16_t(((rows[p(r)]>>p(c))&1U)<<c);
    return result;
}
}
void preparePairedBasis(const rs::Plan& plan,PairedTables& tables,PairedOptimizedTables& opt) {
    preparePairedOptimized(plan,opt);tables.updates.resize(4*(plan.n/64));
    for(std::size_t epoch=0;epoch<plan.n/64;++epoch) {
        const auto rows=conjugate(plan.reverseMatrices[epoch]);
        for(unsigned out=0;out<2;++out)
            for(unsigned in=0;in<2;++in) {
                std::uint64_t value=0;
                for(unsigned j=0;j<8;++j)value|=std::uint64_t((rows[8*out+j]>>(8*in))&255U)<<(8*j);
                tables.updates[4*epoch+2*out+in]=value;
            }
    }
}
void customizePairedField16Basis(rs::Plan& plan,std::uint64_t seed,PairedTables& tables,PairedOptimizedTables& opt) {
    // This prepares M in the fast basis. P^-1 M P below is the literal matrix.
    customizePairedField16(plan,seed,tables,opt);
    for(auto& rows:plan.reverseMatrices)rows=conjugate(rows);
}
}
