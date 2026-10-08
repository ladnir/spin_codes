#pragma once
#include "HighDistanceWide.h"
#include "HdFastOuter.h"

namespace spin::research::char2::hdfast {

#ifndef HDF_OUTER_TILE
#define HDF_OUTER_TILE 64
#endif
#ifndef HDF_SPLIT_ROUTE
#define HDF_SPLIT_ROUTE 1
#endif
#ifndef HDF_PREFETCH_REGION
#define HDF_PREFETCH_REGION 0
#endif
#ifndef HDF_HARDWARE_GATHER
#define HDF_HARDWARE_GATHER 0
#endif
#ifndef HDF_ROUTE_COPY
#define HDF_ROUTE_COPY 1
#endif

// Research-only forward counterpart of the unchanged [256,64]/keep2 map.
// The full GF(2^16) setup distribution and output coordinates are unchanged.
struct ForwardPlan : hdwide::Plan<16,2,256> {
    std::vector<wide16::Coefficient> routedInput;
    std::vector<std::uint32_t> foldedRoute;
    std::vector<std::uint16_t> localRoute;
    explicit ForwardPlan(std::size_t k,std::uint64_t seed):hdwide::Plan<16,2,256>(k,seed) {
        routedInput.resize(inputN);foldedRoute.resize(inputN);
        for(std::size_t i=0;i!=inputN;++i){
            const auto chronological=16*(i/16)+4*(i%4)+(i%16)/4;
            routedInput[i]=wide16::prepareCoefficient(dx[chronological]);foldedRoute[i]=route[chronological];
        }
        if(groups<=65536){
            localRoute.resize(inputN);
            for(std::size_t i=0;i!=inputN;++i)localRoute[i]=static_cast<std::uint16_t>(foldedRoute[i]-(i/groups)*groups);
        }
    }
};

#if defined(_MSC_VER)
#define HDF_INLINE __forceinline
#define HDF_NOINLINE __declspec(noinline)
#else
#define HDF_INLINE inline __attribute__((always_inline))
#define HDF_NOINLINE __attribute__((noinline))
#endif

template<bool Local>HDF_INLINE __m512i routedFour(const ForwardPlan& p,const Block* plane,
                                                std::size_t begin,std::uint32_t regionOffset){
    const auto* index=Local?p.localRoute.data()+begin:nullptr;
    const auto* full=p.foldedRoute.data()+begin;
    if constexpr(Local&&HDF_HARDWARE_GATHER){
        const auto packed=_mm_loadl_epi64(reinterpret_cast<const __m128i*>(index));
        const auto duplicate=_mm_unpacklo_epi16(packed,packed);
        const auto word=_mm256_slli_epi32(_mm256_cvtepu16_epi32(duplicate),1);
        const auto slots=_mm256_or_si256(word,_mm256_setr_epi32(0,1,0,1,0,1,0,1));
        return wide16::multiply512(_mm512_i32gather_epi64(slots,plane,8),p.routedInput.data()+begin);
    }
    const auto i0=Local?index[0]:full[0]-regionOffset,i1=Local?index[1]:full[1]-regionOffset;
    const auto i2=Local?index[2]:full[2]-regionOffset,i3=Local?index[3]:full[3]-regionOffset;
    // Explicit record loads avoid scalar lo/hi reconstruction of each Block.
    auto value=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(plane+i0)));
    value=_mm512_inserti32x4(value,_mm_loadu_si128(reinterpret_cast<const __m128i*>(plane+i1)),1);
    value=_mm512_inserti32x4(value,_mm_loadu_si128(reinterpret_cast<const __m128i*>(plane+i2)),2);
    value=_mm512_inserti32x4(value,_mm_loadu_si128(reinterpret_cast<const __m128i*>(plane+i3)),3);
    return wide16::multiply512(value,p.routedInput.data()+begin);
}

template<bool Local,bool Stored=false>HDF_INLINE void prepareFour(const ForwardPlan& p,const Block* plane,
        std::size_t begin,std::uint32_t regionOffset,__m512i& z0,__m512i& z1,__m512i& z2,__m512i& z3){
    if constexpr(Stored){
        const auto j=begin-regionOffset;
        z0=_mm512_loadu_si512(plane+j);z1=_mm512_loadu_si512(plane+j+4);
        z2=_mm512_loadu_si512(plane+j+8);z3=_mm512_loadu_si512(plane+j+12);
    }else{
        z0=routedFour<Local>(p,plane,begin,regionOffset);
        z1=routedFour<Local>(p,plane,begin+4,regionOffset);
        z2=routedFour<Local>(p,plane,begin+8,regionOffset);
        z3=routedFour<Local>(p,plane,begin+12,regionOffset);
    }
    rs4Circulant<false,true>(z0,z1,z2,z3);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
}

HDF_INLINE void forwardStep(__m512i& moment,__m512i& state,const wide16::Coefficient* coefficient){
    const auto refreshed=wide16::multiply512(state,coefficient);
    state=_mm512_xor_si512(refreshed,moment);
    moment=_mm512_xor_si512(refreshed,rs8ConstantMultiply<true,rs4Scale>(moment));
}

HDF_INLINE void finishFour(__m512i z0,__m512i z1,__m512i z2,__m512i z3,
        __m512i& state,const wide16::Coefficient* ds,Block* output){
    forwardStep(z0,state,ds);forwardStep(z1,state,ds+4);
    forwardStep(z2,state,ds+8);forwardStep(z3,state,ds+12);
    c8TransposeRecords(z0,z1,z2,z3,z0,z1,z2,z3);
    rs4Circulant<true,true>(z0,z1,z2,z3);
    _mm512_storeu_si512(output,z0);_mm512_storeu_si512(output+4,z1);
}

template<bool Local,bool Timed=false>HDF_NOINLINE void routedInner(const ForwardPlan& p,const Block* planes,Block* output,Block* temp,hdwide::PhaseTimes* times){
    auto state=_mm512_setzero_si512();
    for(unsigned region=0;region!=256;++region){
        const auto offset=region*p.groups;
        const auto* plane=planes+region*p.planeStride;
        if constexpr(HDF_PREFETCH_REGION)for(std::size_t j=0;j!=p.groups;j+=4)
            _mm_prefetch(reinterpret_cast<const char*>(plane+j),_MM_HINT_T1);
        hdwide::Clock::time_point start;
        if constexpr(HDF_SPLIT_ROUTE){
            if constexpr(Timed)start=hdwide::Clock::now();
            if constexpr(HDF_ROUTE_COPY){std::memcpy(temp+p.groups,plane,p.groups*sizeof(Block));plane=temp+p.groups;}
            for(std::size_t j=0;j!=p.groups;j+=4)_mm512_storeu_si512(temp+j,routedFour<Local>(p,plane,offset+j,static_cast<std::uint32_t>(offset)));
            plane=temp;
            if constexpr(Timed){times->route+=std::chrono::duration<double,std::micro>(hdwide::Clock::now()-start).count();start=hdwide::Clock::now();}
        }
        for(std::size_t j=0;j!=p.groups;j+=64){
            const auto i=offset+j;auto* destination=output+i/2;
            __m512i z0,z1,z2,z3,w0,w1,w2,w3;
            prepareFour<Local,HDF_SPLIT_ROUTE>(p,plane,i,static_cast<std::uint32_t>(offset),z0,z1,z2,z3);
            prepareFour<Local,HDF_SPLIT_ROUTE>(p,plane,i+16,static_cast<std::uint32_t>(offset),w0,w1,w2,w3);
            finishFour(z0,z1,z2,z3,state,p.stateCoefficient.data()+i,destination);
            prepareFour<Local,HDF_SPLIT_ROUTE>(p,plane,i+32,static_cast<std::uint32_t>(offset),z0,z1,z2,z3);
            finishFour(w0,w1,w2,w3,state,p.stateCoefficient.data()+i+16,destination+8);
            prepareFour<Local,HDF_SPLIT_ROUTE>(p,plane,i+48,static_cast<std::uint32_t>(offset),w0,w1,w2,w3);
            finishFour(z0,z1,z2,z3,state,p.stateCoefficient.data()+i+32,destination+16);
            finishFour(w0,w1,w2,w3,state,p.stateCoefficient.data()+i+48,destination+24);
        }
        if constexpr(Timed&&HDF_SPLIT_ROUTE)times->inner+=std::chrono::duration<double,std::micro>(hdwide::Clock::now()-start).count();
    }
}

// Scratch must be disjoint: planes[planeRecords()] and temp[2*groups]. All
// source is consumed by the outer before any output is written, so exact
// input/output aliasing is supported. Split/copy scheduling is selected by
// measured cost; retained switches are research A/B controls, not code choices.
template<bool Timed=false>HDF_NOINLINE void forward(const ForwardPlan& p,const Block* input,
        Block* output,Block* planes,hdwide::PhaseTimes* times=nullptr,Block* temp=nullptr){
    if constexpr(HDF_SPLIT_ROUTE)if(!temp)throw std::invalid_argument("split route requires regional scratch");
    hdwide::Clock::time_point start;if constexpr(Timed){*times={};start=hdwide::Clock::now();}
    if constexpr(HDF_OUTER_TILE)outerForwardTiled<HDF_OUTER_TILE>(p,input,planes);else outerForward(p,input,planes);
    if constexpr(Timed){times->outer=std::chrono::duration<double,std::micro>(hdwide::Clock::now()-start).count();start=hdwide::Clock::now();}
    if(p.localRoute.empty())routedInner<false,Timed>(p,planes,output,temp,times);else routedInner<true,Timed>(p,planes,output,temp,times);
    if constexpr(Timed&&!HDF_SPLIT_ROUTE)times->inner=std::chrono::duration<double,std::micro>(hdwide::Clock::now()-start).count();
}
#undef HDF_INLINE
#undef HDF_NOINLINE
}
