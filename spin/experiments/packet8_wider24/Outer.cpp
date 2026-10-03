// Instantiate the retained 15-GFNI shared-parity circuit in this translation
// unit, as packet8_codesign/OuterShared.cpp does. Do not compile the retained
// OuterVariants.cpp separately into the same executable.
#include "../k16_codesign_100us/outer_variants/OuterVariants.cpp"
#include "Packet8Wide24.h"

namespace spin::research::packet8wide24 {
namespace {
using Coefficients = std::array<std::uint8_t,9>;

static SPIN_FORCEINLINE void mixSymbol(const Block* input,__m512i* low,
    __m512i* high,const std::uint8_t* coefficients) {
    // Native packet bytes already contain contiguous physical coordinates.
    // Even vectors are the low payload half; odd vectors are the high half.
    // There is no nibble permutation or bit-plane packing at this boundary.
    auto v0=_mm512_load_si512(input);
    auto v1=_mm512_load_si512(input+4);
    auto v2=_mm512_load_si512(input+8);
    auto v3=_mm512_load_si512(input+12);
    auto v4=_mm512_load_si512(input+16);
    auto v5=_mm512_load_si512(input+20);
    auto v6=_mm512_load_si512(input+24);
    auto v7=_mm512_load_si512(input+28);
    rs::tower32byte::applyMultiply(v0,v2,v4,v6,coefficients);
    rs::tower32byte::applyMultiply(v1,v3,v5,v7,coefficients);
    // Each output byte holds two logical GF16 row symbols: its low and
    // high nibbles are coordinates c=4*lane+bit for consecutive RS lanes.
    low[0]=v0;low[1]=v1;low[2]=v2;low[3]=v3;
    high[0]=v4;high[1]=v5;high[2]=v6;high[3]=v7;
}

static SPIN_NOINLINE void outerGroup(const Block* __restrict input,
    Block* __restrict output,const Coefficients* __restrict coefficients) {
    // Fixed 8KiB aligned workspace; never placed in a coroutine frame.
    // Retain the symbol loop: fully unrolling sixteen field multiplications
    // would lengthen live ranges before the independent parity planes.
    alignas(64) __m512i low[64],high[64];
    for(unsigned symbol=0;symbol<16;++symbol)
        mixSymbol(input+32*symbol,low+4*symbol,high+4*symbol,
            coefficients[symbol].data());
    k16codesign::outervariants::finish<0,true>(low,output);
    k16codesign::outervariants::finish<2,true>(low,output);
    k16codesign::outervariants::finish<0,true>(high,output+128);
    k16codesign::outervariants::finish<2,true>(high,output+128);
}
}

void outerFast(const Block* scratch,Block* output,const Plan& plan) {
    static_assert(sizeof(Block)==16 && groupStride%4==0);
    for(std::size_t group=0;group<plan.groups;++group)
        outerGroup(scratch+groupStride*group,output+256*group,
            plan.outerCoefficients.data()+16*group);
}
}
