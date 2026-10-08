#pragma once
#include "RateSpeed4RouteNext.h"
#include "RateSpeed4InnerMore.h"
#include <chrono>

namespace spin::research::char2 {
template<unsigned D=40>struct RateSpeed4ProfilePlan:RateSpeed4RoutedPlan<D>{
    std::vector<std::uint32_t> byteOffset;
    std::vector<std::uint16_t> qwordOffset;
    explicit RateSpeed4ProfilePlan(std::size_t k,std::uint64_t seed):RateSpeed4RoutedPlan<D>(k,seed){
        if(this->groups<=std::numeric_limits<std::uint32_t>::max()/sizeof(Block)){
            byteOffset.resize(this->inputN);
            for(std::size_t i=0;i!=this->inputN;++i)byteOffset[i]=static_cast<std::uint32_t>((this->foldedInverseRoute[i]-(i/this->groups)*this->groups)*sizeof(Block));
        }
        if(this->groups<=32768){
            qwordOffset.resize(this->inputN);
            for(std::size_t i=0;i!=this->inputN;++i)qwordOffset[i]=static_cast<std::uint16_t>(2*(this->foldedInverseRoute[i]-(i/this->groups)*this->groups));
        }
    }
};
struct Rs4PhaseTimes {double compute=0,route=0,outer=0;};
inline __m512i rs4GatherBalanced(const Block* input,const std::uint16_t* index){
    const auto low=_mm256_inserti128_si256(_mm256_castsi128_si256(reg(input[index[0]])),reg(input[index[1]]),1);
    const auto high=_mm256_inserti128_si256(_mm256_castsi128_si256(reg(input[index[2]])),reg(input[index[3]]),1);
    return _mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
}
inline __m512i rs4BroadcastFourDword(const std::uint8_t* coefficients){
    std::uint32_t packed;std::memcpy(&packed,coefficients,4);
    const auto repeated=_mm512_set1_epi32(static_cast<int>(packed));
    return sf16ExpandRow<0>(repeated);
}
template<bool Balanced>inline __m512i rs4GatherByteOffsets(const Block* input,const std::uint32_t* index){
    const auto* bytes=reinterpret_cast<const unsigned char*>(input);
    const auto a=_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+index[0]));
    const auto b=_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+index[1]));
    const auto c=_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+index[2]));
    const auto d=_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+index[3]));
    if constexpr(Balanced){
        const auto low=_mm256_inserti128_si256(_mm256_castsi128_si256(a),b,1);
        const auto high=_mm256_inserti128_si256(_mm256_castsi128_si256(c),d,1);
        return _mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
    }else{
        auto result=_mm512_castsi128_si512(a);result=_mm512_inserti32x4(result,b,1);result=_mm512_inserti32x4(result,c,2);return _mm512_inserti32x4(result,d,3);
    }
}
template<bool Balanced>inline __m512i rs4GatherQwordOffsets(const Block* input,const std::uint16_t* index){
    const auto* bytes=reinterpret_cast<const unsigned char*>(input);
    const auto a=_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+8*std::size_t(index[0])));
    const auto b=_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+8*std::size_t(index[1])));
    const auto c=_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+8*std::size_t(index[2])));
    const auto d=_mm_loadu_si128(reinterpret_cast<const __m128i*>(bytes+8*std::size_t(index[3])));
    if constexpr(Balanced){
        const auto low=_mm256_inserti128_si256(_mm256_castsi128_si256(a),b,1);
        const auto high=_mm256_inserti128_si256(_mm256_castsi128_si256(c),d,1);
        return _mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
    }else{
        auto result=_mm512_castsi128_si512(a);result=_mm512_inserti32x4(result,b,1);result=_mm512_inserti32x4(result,c,2);return _mm512_inserti32x4(result,d,3);
    }
}
template<bool Balanced>inline __m512i rs4GatherRecordLoads(const Block* input,const std::uint16_t* index){
    const auto a=_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+index[0]));
    const auto b=_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+index[1]));
    const auto c=_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+index[2]));
    const auto d=_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+index[3]));
    if constexpr(Balanced){
        const auto low=_mm256_inserti128_si256(_mm256_castsi128_si256(a),b,1);
        const auto high=_mm256_inserti128_si256(_mm256_castsi128_si256(c),d,1);
        return _mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
    }else{
        auto result=_mm512_castsi128_si512(a);result=_mm512_inserti32x4(result,b,1);result=_mm512_inserti32x4(result,c,2);return _mm512_inserti32x4(result,d,3);
    }
}

// Measurement clocks are compiled out of the ordinary encoder variants.
// Timed mode measures each compute/route region in its real execution order,
// retaining the recursive history and the producer-to-consumer cache state.
template<unsigned Variant,bool Timed=false,unsigned D>
inline void rs4TransposeProfile(const RateSpeed4ProfilePlan<D>& p,const Block* input,Block* output,
                                Block* planes,Block* temp,Rs4PhaseTimes* phases=nullptr){
    using Clock=std::chrono::steady_clock;
    auto history=_mm512_setzero_si512();
    for(unsigned r=64;r--;){
        const auto offset=r*p.groups;Clock::time_point start,mid;
        if constexpr(Timed)start=Clock::now();
        if constexpr(((Variant>>3)&3)==1)rs4InnerTransposeCompact(input+3*(offset/4),temp,history,p.stateDiagonal.data()+offset,p.groups/4);
        else if constexpr(((Variant>>3)&3)==2)rs4InnerTransposeRolling<false>(input+3*(offset/4),temp,history,p.stateDiagonal.data()+offset,p.groups/4);
        else if constexpr(((Variant>>3)&3)==3)rs4InnerTransposeRolling<true>(input+3*(offset/4),temp,history,p.stateDiagonal.data()+offset,p.groups/4);
        else for(auto step=(offset+p.groups)/4;step!=offset/4;){
            step-=16;const auto begin=4*step;
            rs4InnerTransposeSixteenPipeline<true,true>(input+12*(step/4),temp+begin-offset,history,p.stateDiagonal.data()+begin);
        }
        if constexpr(Timed){mid=Clock::now();phases->compute+=std::chrono::duration<double,std::micro>(mid-start).count();}
        auto* region=planes+r*p.planeStride;
        if((Variant&64)&&!p.qwordOffset.empty()){
            for(std::size_t j=0;j!=p.groups;j+=4){
                const auto value=rs4GatherQwordOffsets<bool(Variant&1)>(temp,p.qwordOffset.data()+offset+j);
                const auto diagonal=(Variant&2)?rs4BroadcastFourDword(p.routedDiagonal.data()+offset+j):sfBroadcastFour(p.routedDiagonal.data()+offset+j);
                _mm512_storeu_si512(region+j,_mm512_gf2p8mul_epi8(value,diagonal));
            }
        }else if((Variant&4)&&!p.byteOffset.empty()){
            for(std::size_t j=0;j!=p.groups;j+=4){
                const auto value=rs4GatherByteOffsets<bool(Variant&1)>(temp,p.byteOffset.data()+offset+j);
                const auto diagonal=(Variant&2)?rs4BroadcastFourDword(p.routedDiagonal.data()+offset+j):sfBroadcastFour(p.routedDiagonal.data()+offset+j);
                _mm512_storeu_si512(region+j,_mm512_gf2p8mul_epi8(value,diagonal));
            }
        }else if(!p.foldedInverseLocalRoute.empty()){
            const auto* indices=p.foldedInverseLocalRoute.data()+offset;
            for(std::size_t j=0;j!=p.groups;j+=4){
                const auto value=(Variant&128)?rs4GatherRecordLoads<bool(Variant&1)>(temp,indices+j):((Variant&1)?rs4GatherBalanced(temp,indices+j):c8GatherLocal4(temp,indices+j));
                const auto diagonal=(Variant&2)?rs4BroadcastFourDword(p.routedDiagonal.data()+offset+j):sfBroadcastFour(p.routedDiagonal.data()+offset+j);
                _mm512_storeu_si512(region+j,_mm512_gf2p8mul_epi8(value,diagonal));
            }
        }else{
            for(std::size_t j=0;j!=p.groups;j+=4){
                const auto value=sf16GatherRegion4(temp,p.foldedInverseRoute.data()+offset+j,static_cast<std::uint32_t>(offset));
                const auto diagonal=(Variant&2)?rs4BroadcastFourDword(p.routedDiagonal.data()+offset+j):sfBroadcastFour(p.routedDiagonal.data()+offset+j);
                _mm512_storeu_si512(region+j,_mm512_gf2p8mul_epi8(value,diagonal));
            }
        }
        if constexpr(Timed)phases->route+=std::chrono::duration<double,std::micro>(Clock::now()-mid).count();
    }
    Clock::time_point start;if constexpr(Timed)start=Clock::now();
    if((reinterpret_cast<std::uintptr_t>(output)&63)==0)rs4OuterNext<D,true,Rs4OuterNextVariant::FusedMultiply>(p,planes,output);
    else rs4OuterNext<D,false,Rs4OuterNextVariant::FusedMultiply>(p,planes,output);
    if constexpr(Timed)phases->outer+=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
}
} // namespace spin::research::char2
