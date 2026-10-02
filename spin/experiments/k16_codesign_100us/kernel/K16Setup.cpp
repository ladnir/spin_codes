#include "K16CodeSign.h"
#include "../../../../research/workstreams/k16_design/implementation/rs16x8/Tower32ByteRandomizer.h"
#include "../../../src/kernels/SetupRandom.h"
#include <bit>
#include <stdexcept>

namespace spin::research::k16codesign {
namespace {
using detail::kernel::setup::Words;
constexpr std::array<std::uint16_t,16> basis{
    0x1,0x2,0x4,0x8,0x10,0x20,0x40,0x3500,
    0x1180,0x5980,0x2c80,0xcf80,0x9400,0xe400,0xe00,0xd980};
constexpr std::array<std::uint16_t,16> inverseBasis{
    0x1,0x2,0x4,0x8,0x10,0x20,0x40,0xa500,
    0x6980,0x1a00,0x5f80,0x580,0xcd80,0xfb00,0x680,0x8200};

void requirePlan(const rs::Plan& plan) {
    if(plan.variant != rs::Variant::Rs16Gf16 ||
       plan.innerKernel != rs::InnerKernel::RetainedStreaming)
        throw std::invalid_argument("codesign requires Rs16Gf16/RetainedStreaming");
}
std::uint16_t nonzero(Words& words) {
    std::uint16_t value;
    do { value = static_cast<std::uint16_t>(words()); } while(!value);
    return value;
}
std::array<std::uint16_t,16> multiplyRows(std::uint16_t scalar) {
    std::array<std::uint16_t,16> rows{};
    for(unsigned column=0;column<16;++column) {
        const auto image=rs::tower32::multiply16(scalar,std::uint16_t(1U<<column));
        for(unsigned row=0;row<16;++row)
            rows[row]|=std::uint16_t(((image>>row)&1U)<<column);
    }
    return rows;
}
void compactRows(const std::array<std::uint16_t,16>& rows,std::uint64_t* output) {
    for(unsigned out=0;out<2;++out)
        for(unsigned in=0;in<2;++in) {
            std::uint64_t value=0;
            for(unsigned j=0;j<8;++j)
                value|=std::uint64_t((rows[8*out+j]>>(8*in))&255U)<<(8*(7-j));
            output[2*out+in]=value;
        }
}
}

void customizeOuterField16(rs::Plan& plan,std::uint64_t seed,Tables& tables) {
    if(plan.variant != rs::Variant::Rs16Gf16)
        throw std::invalid_argument("codesign outer requires Rs16Gf16");
    Words words(seed ^ 0x75a1dc09ULL);
    tables.outerField.resize(16*plan.groups);
    auto* coefficients=const_cast<std::uint64_t*>(plan.compactCoefficients());
    for(std::size_t symbol=0;symbol<16*plan.groups;++symbol) {
        const auto scalar=nonzero(words);
        rs::tower32byte::coefficients16(scalar,tables.outerField[symbol].data());
        // Forward symbol map is M_scalar^T; transpose uses ordinary multiply.
        plan.outerMatrices16[symbol]=multiplyRows(scalar);
        compactRows(plan.outerMatrices16[symbol],coefficients+4*symbol);
    }
}

void customizeInnerField16(rs::Plan& plan,std::uint64_t seed,Tables& tables) {
    requirePlan(plan);
    Words words(seed ^ 0x3f625a92ULL);
    tables.innerField.resize(plan.n/64);
    for(std::size_t epoch=0;epoch<plan.n/64;++epoch) {
        const auto scalar=nonzero(words);
        rs::tower32byte::coefficients16(scalar,tables.innerField[epoch].data());
        const auto packedRows=multiplyRows(scalar);
        auto& reverse=plan.reverseMatrices[epoch];
        reverse.fill(0);
        // Packed state is P times literal state; literal reverse is P^-1 M P.
        for(unsigned row=0;row<16;++row) {
            unsigned mask=0;
            for(unsigned bits=inverseBasis[row];bits;bits&=bits-1)
                mask^=packedRows[std::countr_zero(bits)];
            for(;mask;mask&=mask-1)reverse[row]^=basis[std::countr_zero(mask)];
        }
        compactRows(packedRows,plan.denseUpdates[epoch].matrix);
    }
}
}
