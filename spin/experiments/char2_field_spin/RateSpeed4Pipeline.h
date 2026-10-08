#pragma once
#include "RateSpeed4RouteNext.h"

namespace spin::research::char2 {

#if defined(_MSC_VER)
#define RS4P_INLINE __forceinline
#define RS4P_NOINLINE __declspec(noinline)
#else
#define RS4P_INLINE inline __attribute__((always_inline))
#define RS4P_NOINLINE __attribute__((noinline))
#endif

namespace rs4_pipeline_detail {

// Route one completed64-record destination interval. This is the selected
// scalar-index/four-record gather and diagonal expansion, without changed
// metadata, hardware gathers, software prefetching, or temporary layout.
template<bool Compact,unsigned Group>
RS4P_INLINE void route4(const Block* source,Block* destination,
                       const std::uint16_t* compact,const std::uint32_t* wide,
                       const std::uint8_t* diagonal,std::uint32_t offset){
    static_assert(Group<16);
    constexpr unsigned j=4*Group;
    __m512i value;
    if constexpr(Compact)value=c8GatherLocal4(source,compact+j);
    else value=sf16GatherRegion4(source,wide+j,offset);
    value=_mm512_gf2p8mul_epi8(value,sfBroadcastFour(diagonal+j));
    _mm512_storeu_si512(destination+j,value);
}

template<bool Compact>
RS4P_INLINE void route64(const Block* source,Block* destination,
                        const std::uint16_t* compact,const std::uint32_t* wide,
                        const std::uint8_t* diagonal,std::uint32_t offset){
    // Explicit calls keep the complete unroll forced-inline. GCC outlined
    // the generic c8Unroll wrapper and repeatedly reloaded captured pointer
    // references from its closure, obscuring this scheduling comparison.
#define RS4P_ROUTE(I) route4<Compact,I>(source,destination,compact,wide,diagonal,offset)
    RS4P_ROUTE(0);RS4P_ROUTE(1);RS4P_ROUTE(2);RS4P_ROUTE(3);
    RS4P_ROUTE(4);RS4P_ROUTE(5);RS4P_ROUTE(6);RS4P_ROUTE(7);
    RS4P_ROUTE(8);RS4P_ROUTE(9);RS4P_ROUTE(10);RS4P_ROUTE(11);
    RS4P_ROUTE(12);RS4P_ROUTE(13);RS4P_ROUTE(14);RS4P_ROUTE(15);
#undef RS4P_ROUTE
}

// Ordinary non-coroutine helper. Two caller-owned regional temporaries let
// the producer's carried-state arithmetic overlap independent random loads
// from the previous completed region. The state itself remains ordinary h,
// and the producer still visits every packet in descending chronological
// order. Exactly one region is produced and one routed at any instant.
template<unsigned D,bool Compact>
RS4P_NOINLINE void inner(const RateSpeed4RoutedPlan<D>& p,const Block* input,
                        Block* planes,Block* temp2){
    auto state=_mm512_setzero_si512();
    auto* previous=temp2;
    auto* next=temp2+p.groups;

    // Warm up the latest chronological region,63. Each call consumes48
    // retained records and emits64 raw records into the folded temporary.
    const auto lastOffset=63*p.groups;
    for(auto local=p.groups;local;){
        local-=64;
        const auto begin=lastOffset+local;
        rs4InnerTransposeSixteenPipeline<true,true>(
            input+12*(begin/16),previous+local,state,p.stateDiagonal.data()+begin);
    }

    for(unsigned region=63;region!=0;){
        --region;
        const auto offset=region*p.groups;
        const auto previousOffset=offset+p.groups;
        auto* destination=planes+(region+1)*p.planeStride;
        const auto* diagonal=p.routedDiagonal.data()+previousOffset;
        const std::uint16_t* compact=nullptr;
        if constexpr(Compact)compact=p.foldedInverseLocalRoute.data()+previousOffset;
        const auto* wide=p.foldedInverseRoute.data()+previousOffset;

        for(std::size_t j=0;j!=p.groups;j+=64){
            const auto local=p.groups-j-64;
            const auto begin=offset+local;
            rs4InnerTransposeSixteenPipeline<true,true>(
                input+12*(begin/16),next+local,state,p.stateDiagonal.data()+begin);
            const std::uint16_t* compactBlock=nullptr;
            if constexpr(Compact)compactBlock=compact+j;
            route64<Compact>(previous,destination+j,compactBlock,wide+j,
                             diagonal+j,static_cast<std::uint32_t>(previousOffset));
        }
        std::swap(previous,next);
    }

    // Drain region0 after its final producer step. Both buffers have always
    // remained disjoint, so no producer write can overwrite a pending read.
    for(std::size_t j=0;j!=p.groups;j+=64){
        const std::uint16_t* compact=nullptr;
        if constexpr(Compact)compact=p.foldedInverseLocalRoute.data()+j;
        route64<Compact>(previous,planes+j,compact,p.foldedInverseRoute.data()+j,
                         p.routedDiagonal.data()+j,0);
    }
}
} // namespace rs4_pipeline_detail

// temp2 contains2*groups records. It must be disjoint from input and planes;
// planeStride has the same >=groups contract as the selected routed plan.
// Compact and wide routes share the same two-buffer schedule.
template<unsigned D>
inline void rs4RegionPipeline(const RateSpeed4RoutedPlan<D>& p,const Block* input,
                              Block* planes,Block* temp2){
    if(!p.foldedInverseLocalRoute.empty())
        rs4_pipeline_detail::inner<D,true>(p,input,planes,temp2);
    else rs4_pipeline_detail::inner<D,false>(p,input,planes,temp2);
}

// Same selected arithmetic, sampler, map, and shortened source convention.
// Caller scratch: planeRecords() plane records plus2*groups temp records.
// input/output may alias exactly: the entire producer/route stage completes
// before the first outer output store. Output alignment selects its stores.
template<unsigned D>
inline void rs4TransposePipeline(const RateSpeed4RoutedPlan<D>& p,const Block* input,
                                 Block* output,Block* planes,Block* temp2){
    rs4RegionPipeline(p,input,planes,temp2);
    if((reinterpret_cast<std::uintptr_t>(output)&63)==0)
        rs4OuterNext<D,true,Rs4OuterNextVariant::FusedMultiply>(p,planes,output);
    else rs4OuterNext<D,false,Rs4OuterNextVariant::FusedMultiply>(p,planes,output);
}

#undef RS4P_INLINE
#undef RS4P_NOINLINE
} // namespace spin::research::char2
