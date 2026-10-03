// Reuse the exact frozen mode52 shared-parity circuit in this translation
// unit, just as ScalarAdapters reuses the retained literal outer oracle.
// The frozen source is not otherwise compiled into this executable.
#include "../k16_codesign_100us/outer_variants/OuterVariants.cpp"
#include "Packet8.h"

namespace spin::research::packet8 {
namespace {
static SPIN_NOINLINE void group(const Block* __restrict input,Block* __restrict output,
    const std::array<std::uint8_t,3>* __restrict coefficients) {
    alignas(64) __m512i packed[64];
    for(unsigned symbol=0;symbol<16;++symbol) {
        const auto a=_mm512_load_si512(input+16*symbol);
        const auto b=_mm512_load_si512(input+16*symbol+4);
        const auto c=_mm512_load_si512(input+16*symbol+8);
        const auto d=_mm512_load_si512(input+16*symbol+12);
        const auto lo=rs::tower32byte::quadraticMultiply(a,c,coefficients[symbol].data());
        const auto hi=rs::tower32byte::quadraticMultiply(b,d,coefficients[symbol].data());
        packed[4*symbol]=lo.lo;packed[4*symbol+1]=hi.lo;
        packed[4*symbol+2]=lo.hi;packed[4*symbol+3]=hi.hi;
    }
    k16codesign::outervariants::finish<0,true>(packed,output);
    k16codesign::outervariants::finish<2,true>(packed,output);
}
}
void outerFastShared(const Block* scratch,Block* output,const Plan& plan) {
    for(std::size_t g=0;g<plan.outer.groups;++g)
        group(scratch+rs::groupStride*g,output+128*g,plan.native.field.data()+16*g);
}
}
