#include "K16NativeOuter.h"
#include "../../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"
#include "../../../src/kernels/SetupRandom.h"
#include <stdexcept>

namespace spin::research::k16codesign {
namespace {
// Old coordinate c lives at physical input 4*(c%4)+c/4. This permutation is
// its own inverse. Apply to matrix columns, leaving output coordinates intact.
std::array<std::uint16_t,16> permuteColumns(const std::array<std::uint16_t,16>& rows) {
    std::array<std::uint16_t,16> result{};
    for(unsigned row=0;row<16;++row)
        for(unsigned c=0;c<16;++c)
            result[row]|=std::uint16_t(((rows[row]>>c)&1U)<<(4*(c%4)+c/4));
    return result;
}
void compact(const std::array<std::uint16_t,16>& rows,std::uint64_t* out) {
    for(unsigned o=0;o<2;++o)
        for(unsigned i=0;i<2;++i) {
            std::uint64_t value=0;
            for(unsigned j=0;j<8;++j)
                value|=std::uint64_t((rows[8*o+j]>>(8*i))&255U)<<(8*(7-j));
            out[2*o+i]=value;
        }
}
void requirePlan(const rs::Plan& plan) {
    if(plan.variant!=rs::Variant::Rs16Gf16)
        throw std::invalid_argument("native outer requires Rs16Gf16");
}
}
void prepareNativeGl16(const rs::Plan& plan,NativeOuterTables& tables) {
    requirePlan(plan);
    tables.gl.resize(64*plan.groups);
    for(std::size_t symbol=0;symbol<16*plan.groups;++symbol)
        compact(permuteColumns(plan.outerMatrices16[symbol]),tables.gl.data()+4*symbol);
}
void customizeNativeField16(rs::Plan& plan,std::uint64_t seed,NativeOuterTables& tables) {
    requirePlan(plan);
    detail::kernel::setup::Words words(seed^0x75a1dc09ULL);
    tables.field.resize(16*plan.groups);
    auto* oldCompact=const_cast<std::uint64_t*>(plan.compactCoefficients());
    for(std::size_t symbol=0;symbol<16*plan.groups;++symbol) {
        std::uint16_t scalar;
        do { scalar=static_cast<std::uint16_t>(words()); } while(!scalar);
        rs::tower32byte::coefficients16(scalar,tables.field[symbol].data());
        std::array<std::uint16_t,16> rows{};
        for(unsigned c=0;c<16;++c) {
            const auto image=rs::tower32::multiply16(scalar,std::uint16_t(1U<<c));
            for(unsigned r=0;r<16;++r)rows[r]|=std::uint16_t(((image>>r)&1U)<<c);
        }
        plan.outerMatrices16[symbol]=permuteColumns(rows);
        compact(plan.outerMatrices16[symbol],oldCompact+4*symbol);
    }
    prepareNativeGl16(plan,tables);
}
}
