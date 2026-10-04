#include "K18Certified22.h"
#include "K18Map.h"
#include "../../src/kernels/SetupRandom.h"
#include <algorithm>
#include <bit>
#include <limits>
#include <stdexcept>

namespace spin::research::k18cert22 {
namespace {
using Matrix = std::array<std::uint32_t,22>;
bool fullRank(Matrix rows) {
    for(unsigned column=0; column<22; ++column) {
        unsigned pivot=column;
        while(pivot<22 && !((rows[pivot]>>column)&1U)) ++pivot;
        if(pivot==22) return false;
        std::swap(rows[pivot],rows[column]);
        for(unsigned row=column+1; row<22; ++row)
            if((rows[row]>>column)&1U) rows[row]^=rows[column];
    }
    return true;
}
Matrix multiply(const Matrix& left,const Matrix& right) {
    Matrix result{};
    for(unsigned row=0; row<22; ++row)
        for(auto mask=left[row]; mask; mask&=mask-1)
            result[row]^=right[std::countr_zero(mask)];
    return result;
}
}

void prepare(rs::Plan& plan,std::uint64_t seed,Tables& tables) {
    if(plan.variant!=rs::Variant::Rs16Gf16 || !plan.n || plan.n%64)
        throw std::invalid_argument("certified22 requires complete t64 RS16x8 steps");
    tables.routeBytes.resize(plan.route.size());
    for(std::size_t packet=0; packet<plan.route.size(); ++packet) {
        if(plan.route[packet]>std::numeric_limits<std::uint32_t>::max()/sizeof(Block))
            throw std::invalid_argument("certified22 route exceeds 32-bit byte offsets");
        tables.routeBytes[packet]=static_cast<std::uint32_t>(sizeof(Block)*plan.route[packet]);
    }
    const auto epochs=plan.n/64;
    tables.reverseMatrices.resize(epochs);
    tables.updates.assign(9*epochs,0);
    detail::kernel::setup::Words words(seed^0x3f625a92ULL);
    for(std::size_t epoch=0; epoch<epochs; ++epoch) {
        auto& literal=tables.reverseMatrices[epoch];
        // Rejection from all 22x22 binary matrices gives uniform GL22 in
        // the ideal independent-word model. No dummy state bits are sampled.
        do {
            for(auto& row:literal) row=static_cast<std::uint32_t>(words())&0x3fffffU;
        } while(!fullRank(literal));
        const auto transformed=multiply(multiply(map::basis,literal),map::inverseBasis);
        for(unsigned out=0; out<3; ++out)
            for(unsigned in=0; in<3; ++in) {
                std::uint64_t packed=0;
                for(unsigned j=0; j<8; ++j) {
                    const auto row=8*out+j;
                    if(row<22) packed|=std::uint64_t((transformed[row]>>(8*in))&255U)<<(8*j);
                }
                tables.updates[9*epoch+3*out+in]=packed;
            }
    }
}
}
