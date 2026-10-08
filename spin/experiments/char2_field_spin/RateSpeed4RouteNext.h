#pragma once
#include "RateSpeed4Code.h"
#include "RateSpeed4InnerNext.h"
#include "RateSpeed4OuterNext.h"

namespace spin::research::char2 {
template<unsigned D=40>struct RateSpeed4RoutedPlan:RateSpeed4Plan<D>{
    std::vector<std::uint8_t> routedDiagonal;
    explicit RateSpeed4RoutedPlan(std::size_t k,std::uint64_t seed):RateSpeed4Plan<D>(k,seed),routedDiagonal(this->inputN){
        for(std::size_t i=0;i!=this->inputN;++i)routedDiagonal[this->route[i]]=this->inputDiagonal[i];
        this->planeStride=this->groups+4;
    }
};
// Broad comparison harness only. The selected routed plan above does not
// allocate the gauge/pair coefficients needed by rejected alternatives.
template<unsigned D=40>struct RateSpeed4SchedulePlan:RateSpeed4RoutedPlan<D>{
    std::vector<std::uint8_t> gaugeCoefficients,pairCoefficients;
    explicit RateSpeed4SchedulePlan(std::size_t k,std::uint64_t seed):RateSpeed4RoutedPlan<D>(k,seed){
        gaugeCoefficients.resize(2*this->inputN);pairCoefficients.resize(3*this->inputN/2);
        rs4PrepareGauge(this->stateDiagonal.data(),this->inputN,gaugeCoefficients.data());
        rs4PreparePairs(this->stateDiagonal.data(),this->inputN,pairCoefficients.data());
    }
};

inline __m512i rs4GatherHardware(const Block* input,const std::uint16_t* index){
    const auto compact=_mm_loadl_epi64(reinterpret_cast<const __m128i*>(index));
    auto words=_mm256_cvtepu16_epi32(_mm_unpacklo_epi16(compact,compact));
    words=_mm256_add_epi32(_mm256_slli_epi32(words,1),_mm256_setr_epi32(0,1,0,1,0,1,0,1));
    return _mm512_i32gather_epi64(words,input,8);
}

// Exact-map scheduling controls. Moving Dx across the permutation changes
// coefficient order only: routedDiagonal[route[i]] = inputDiagonal[i].
template<unsigned D,unsigned Prefetch=0,bool Hardware=false,bool RouteDiagonal=false,unsigned InnerVariant=0,bool BatchedRoute=false,bool DirectOutput=false,class Plan>
inline void rs4RegionRouteNext(const Plan& p,const Block* input,Block* planes,Block* temp){
    static_assert(D==Plan::dimension);
    static_assert(Prefetch%4==0,"route prefetch reads complete four-record groups");
    static_assert(!BatchedRoute||RouteDiagonal);
    static_assert(InnerVariant!=4||(RouteDiagonal&&DirectOutput));
    auto state=_mm512_setzero_si512();alignas(64) Block packet[32];
    for(unsigned r=64;r--;){
        const auto offset=r*p.groups;
        for(auto step=(offset+p.groups)/4;step!=offset/4;){
            step-=InnerVariant==4?16:8;const auto begin=4*step;
            auto* target=RouteDiagonal&&DirectOutput?temp+begin-offset:packet;
            if constexpr(InnerVariant==4)rs4InnerTransposeSixteenPipeline<true,true>(input+12*(step/4),target,state,p.stateDiagonal.data()+begin);
            else if constexpr(InnerVariant==1)rs4InnerTransposeEightGauge<true,true>(input+12*(step/4),target,state,p.gaugeCoefficients.data()+2*begin);
            else if constexpr(InnerVariant==2)rs4InnerTransposeEightPairs<true,true,false>(input+12*(step/4),target,state,p.stateDiagonal.data()+begin);
            else if constexpr(InnerVariant==3)rs4InnerTransposeEightPairs<true,true,true>(input+12*(step/4),target,state,p.pairCoefficients.data()+48*(begin/32));
            else rs4InnerTransposeEight<true,true,true>(input+12*(step/4),target,state,p.stateDiagonal.data()+begin);
            __m512i packed,high;
            if constexpr(!RouteDiagonal){
                packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p.packedInputDiagonal.data()+begin)));
                high=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p.packedInputDiagonal.data()+begin+16)));
            }
            if constexpr(!(RouteDiagonal&&DirectOutput))c8Unroll(std::make_index_sequence<8>{},[&](auto chunk){
                auto value=_mm512_load_si512(packet+4*chunk);
                if constexpr(!RouteDiagonal)value=_mm512_gf2p8mul_epi8(value,sf16ExpandRow<unsigned(chunk)%4>(unsigned(chunk)<4?packed:high));
                _mm512_storeu_si512(temp+begin-offset+4*chunk,value);
            });
        }
        auto* region=planes+r*p.planeStride;
        if(!p.foldedInverseLocalRoute.empty()){
            const auto* indices=p.foldedInverseLocalRoute.data()+offset;
            if constexpr(BatchedRoute){
                for(std::size_t j=0;j!=p.groups;j+=16){
                    const auto diagonal=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p.routedDiagonal.data()+offset+j)));
                    const auto v0=c8GatherLocal4(temp,indices+j),v1=c8GatherLocal4(temp,indices+j+4);
                    const auto v2=c8GatherLocal4(temp,indices+j+8),v3=c8GatherLocal4(temp,indices+j+12);
                    _mm512_storeu_si512(region+j,_mm512_gf2p8mul_epi8(v0,sf16ExpandRow<0>(diagonal)));
                    _mm512_storeu_si512(region+j+4,_mm512_gf2p8mul_epi8(v1,sf16ExpandRow<1>(diagonal)));
                    _mm512_storeu_si512(region+j+8,_mm512_gf2p8mul_epi8(v2,sf16ExpandRow<2>(diagonal)));
                    _mm512_storeu_si512(region+j+12,_mm512_gf2p8mul_epi8(v3,sf16ExpandRow<3>(diagonal)));
                }
            }else for(std::size_t j=0;j!=p.groups;j+=4){
                if constexpr(Prefetch){if(j+Prefetch<p.groups){
                    _mm_prefetch(reinterpret_cast<const char*>(temp+indices[j+Prefetch]),_MM_HINT_T0);
                    _mm_prefetch(reinterpret_cast<const char*>(temp+indices[j+Prefetch+1]),_MM_HINT_T0);
                    _mm_prefetch(reinterpret_cast<const char*>(temp+indices[j+Prefetch+2]),_MM_HINT_T0);
                    _mm_prefetch(reinterpret_cast<const char*>(temp+indices[j+Prefetch+3]),_MM_HINT_T0);
                }}
                auto value=Hardware?rs4GatherHardware(temp,indices+j):c8GatherLocal4(temp,indices+j);
                if constexpr(RouteDiagonal)value=_mm512_gf2p8mul_epi8(value,sfBroadcastFour(p.routedDiagonal.data()+offset+j));
                _mm512_storeu_si512(region+j,value);
            }
        }else{
            for(std::size_t j=0;j!=p.groups;j+=4){
                auto value=sf16GatherRegion4(temp,p.foldedInverseRoute.data()+offset+j,static_cast<std::uint32_t>(offset));
                if constexpr(RouteDiagonal)value=_mm512_gf2p8mul_epi8(value,sfBroadcastFour(p.routedDiagonal.data()+offset+j));
                _mm512_storeu_si512(region+j,value);
            }
        }
    }
}
// Selected K20 research schedule. No gauge or pair tables are required.
// The plan adds one cache line to the natural column stride, making the
// stride an odd number of cache lines. Callers supply disjoint planes/temp
// scratch; input/output may
// alias exactly. The whole input is consumed before the outer writes output.
template<unsigned D>
inline void rs4TransposeRouted16(const RateSpeed4RoutedPlan<D>& p,const Block* input,Block* output,Block* planes,Block* temp){
    rs4RegionRouteNext<D,0,false,true,4,false,true>(p,input,planes,temp);
    if((reinterpret_cast<std::uintptr_t>(output)&63)==0)rs4OuterNext<D,true,Rs4OuterNextVariant::FusedMultiply>(p,planes,output);
    else rs4OuterNext<D,false,Rs4OuterNextVariant::FusedMultiply>(p,planes,output);
}
} // namespace spin::research::char2
