#pragma once
#include "Map64S12.h"
#include <type_traits>
namespace bare_spin {
struct Map64S12R2 : Map64S12 {};
struct Map64S12R2Transpose : Map64S12Transpose {};
template<class Map> inline constexpr bool isTwoRoundMap=
    std::is_same_v<Map,Map64S12R2> || std::is_same_v<Map,Map64S12R2Transpose>;
// Called explicitly twice, with the order reversed by the transpose.
OC_FORCEINLINE void imtStep(__m128i* state,u32 dotMask,u32 xorMask) {
    auto dot=_mm_setzero_si128();
    while(dotMask) {const auto j=std::countr_zero(dotMask);dotMask&=dotMask-1;dot=_mm_xor_si128(dot,state[j]);}
    while(xorMask) {const auto j=std::countr_zero(xorMask);xorMask&=xorMask-1;state[j]=_mm_xor_si128(state[j],dot);}
}
// Independent dense-reference row of (T2 T1)^T, not the optimized SIMD update.
inline u32 imtR2ReferenceRow(unsigned j,const u32* masks) {
    u32 row=1U<<j;
    for(unsigned r=0;r<2;++r) if(std::popcount(row&masks[2*r+1])&1) row^=masks[2*r];
    return row;
}
}
