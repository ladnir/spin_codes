#include "Kernels.h"
#include "../../src/kernels/SetupRandom.h"
#include <algorithm>
#include <stdexcept>

namespace spin::research::rs8sizes {
namespace {
using Matrix = std::array<std::uint16_t,16>;
constexpr unsigned permutation[16] = {0,1,2,7,3,4,5,6,8,15,10,12,11,13,14,9};
bool fullRank15(Matrix rows) {
    for(unsigned pivot=0; pivot<15; ++pivot) {
        unsigned selected=pivot;
        while(selected<15 && !((rows[selected]>>pivot)&1U)) ++selected;
        if(selected==15) return false;
        std::swap(rows[pivot],rows[selected]);
        for(unsigned row=pivot+1; row<15; ++row)
            if((rows[row]>>pivot)&1U) rows[row]^=rows[pivot];
    }
    return true;
}
void compact(const Matrix& rows, std::uint64_t* reverse, std::uint64_t* forward) {
    Matrix r{},f{};
    for(unsigned i=0;i<16;++i) for(unsigned j=0;j<16;++j) {
        r[i]|=std::uint16_t(((rows[permutation[i]]>>permutation[j])&1U)<<j);
        f[i]|=std::uint16_t(((rows[permutation[j]]>>permutation[i])&1U)<<j);
    }
    for(unsigned o=0;o<2;++o) for(unsigned i=0;i<2;++i) {
        std::uint64_t a=0,b=0;
        for(unsigned j=0;j<8;++j) {
            a|=std::uint64_t((r[8*o+j]>>(8*i))&255U)<<(8*j);
            b|=std::uint64_t((f[8*o+j]>>(8*i))&255U)<<(8*j);
        }
        reverse[2*o+i]=a;forward[2*o+i]=b;
    }
}
std::uint8_t multiply8(std::uint8_t a,std::uint8_t b) {
    std::uint8_t result=0;
    for(unsigned i=0;i<8;++i) {
        if(b&1U) result^=a;
        a=std::uint8_t((unsigned(a)<<1)^((a&128U)?0x11bU:0U));b>>=1;
    }
    return result;
}
std::uint64_t adjoint8(std::uint8_t scalar) {
    std::uint64_t matrix=0;
    for(unsigned bit=0;bit<8;++bit)
        matrix|=std::uint64_t(multiply8(scalar,std::uint8_t(1U<<bit)))<<(8*(7-bit));
    return matrix;
}
}

void resizePaired15(detail::paired15::Plan& p,const detail::certified22::Plan& geometry,
                    std::uint64_t innerSeed) {
    static_assert(detail::paired15::groupStride==detail::certified22::groupStride);
    if(geometry.k<2048) throw std::invalid_argument("research paired15 requires K>=2048");
    // Both production setups have identical seeded routing and outer tables.
    // Retain the literal production Plan at K16, and check that correspondence.
    if(geometry.k==65536) {
        if(p.route!=geometry.route || p.routeBytes!=geometry.routeBytes ||
           p.outerMatrices16!=geometry.outerMatrices16 || p.outerField!=geometry.outerField)
            throw std::runtime_error("paired15/certified22 setup geometry differs");
        return;
    }
    p.k=geometry.k;p.n=geometry.n;p.groups=geometry.groups;
    p.route=geometry.route;p.routeBytes=geometry.routeBytes;
    p.outerMatrices16=geometry.outerMatrices16;p.outerField=geometry.outerField;
    p.outerAdjoint.resize(p.outerField.size());
    for(std::size_t i=0;i<p.outerField.size();++i)
        p.outerAdjoint[i]={adjoint8(p.outerField[i][0]),adjoint8(p.outerField[i][1]),
                          adjoint8(p.outerField[i][2])};
    p.reverseMatrices.resize(p.n/64);
    p.updates.resize(4*p.reverseMatrices.size());p.forwardUpdates.resize(p.updates.size());
    detail::kernel::setup::Words words(innerSeed^0x3f625a92ULL);
    for(std::size_t epoch=0;epoch<p.reverseMatrices.size();++epoch) {
        auto& rows=p.reverseMatrices[epoch];
        do {
            for(unsigned row=0;row<15;++row) rows[row]=std::uint16_t(words()&0x7fffU);
            rows[15]=0x8000U;
        } while(!fullRank15(rows));
        compact(rows,p.updates.data()+4*epoch,p.forwardUpdates.data()+4*epoch);
    }
}
}
