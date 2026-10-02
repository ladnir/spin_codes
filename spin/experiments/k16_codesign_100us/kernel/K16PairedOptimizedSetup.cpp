#include "K16PairedOptimized.h"
#include "../../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"
#include "../../../src/kernels/SetupRandom.h"
#include <limits>
#include <stdexcept>
namespace spin::research::k16codesign {
namespace {
std::uint64_t multiply8Rows(std::uint8_t scalar) {
    std::uint64_t rows=0;
    for(unsigned c=0;c<8;++c) {
        const auto image=rs::tower32::multiply8(scalar,std::uint8_t(1U<<c));
        for(unsigned r=0;r<8;++r)rows|=std::uint64_t((image>>r)&1U)<<(8*r+c);
    }
    return rows;
}
}
void preparePairedOptimized(const rs::Plan& plan,PairedOptimizedTables& tables) {
    tables.routeBytes.resize(plan.route.size());
    for(std::size_t packet=0;packet<plan.route.size();++packet) {
        if(plan.route[packet]>std::numeric_limits<std::uint32_t>::max()/16)
            throw std::invalid_argument("paired byte route exceeds 32-bit byte-offset range");
        tables.routeBytes[packet]=16*plan.route[packet];
    }
}
void customizePairedField16(rs::Plan& plan,std::uint64_t seed,PairedTables& dense,PairedOptimizedTables& tables) {
    preparePairedOptimized(plan,tables);
    detail::kernel::setup::Words words(seed^0x3f625a92ULL);
    tables.fieldUpdates.resize(3*(plan.n/64));
    for(std::size_t epoch=0;epoch<plan.n/64;++epoch) {
        std::uint16_t scalar;
        do {scalar=static_cast<std::uint16_t>(words());} while(!scalar);
        std::uint8_t coefficients[3];rs::tower32byte::coefficients16(scalar,coefficients);
        for(unsigned j=0;j<3;++j)tables.fieldUpdates[3*epoch+j]=multiply8Rows(coefficients[j]);
        auto& rows=plan.reverseMatrices[epoch];rows.fill(0);
        for(unsigned c=0;c<16;++c) {
            const auto image=rs::tower32::multiply16(scalar,std::uint16_t(1U<<c));
            for(unsigned r=0;r<16;++r)rows[r]|=std::uint16_t(((image>>r)&1U)<<c);
        }
    }
    preparePaired(plan,dense);
}
}
