#pragma once
#include "RateSpeed4Code.h"

namespace spin::research::char2 {

// Exact-map research schedules for the GRS[64,39/40] outer. The three
// alternatives isolate ancestor fusion, constant multiplication, and source
// prefetching. No route, diagonal, coefficient basis, or shortening changes.
enum class Rs4OuterNextVariant {
    FusedMultiply,
    FusedAffine,
    FusedAffinePrefetch
};

#if defined(_MSC_VER)
#define RS4N_INLINE __forceinline
#define RS4N_NOINLINE __declspec(noinline)
#else
#define RS4N_INLINE inline __attribute__((always_inline))
#define RS4N_NOINLINE __attribute__((noinline))
#endif

namespace rs4_outer_next_detail {

constexpr std::uint8_t inverse(std::uint8_t value)noexcept{
    std::uint8_t result=1;
    for(unsigned exponent=254;exponent;exponent>>=1){
        if(exponent&1)result=sfMultiply(result,value);
        value=sfMultiply(value,value);
    }
    return result;
}
constexpr std::uint8_t normalized(unsigned level,unsigned offset)noexcept{
    std::uint8_t value=1,scale=1;
    for(unsigned j=0;j!=(1U<<level);++j){
        value=sfMultiply(value,static_cast<std::uint8_t>(offset^j));
        scale=sfMultiply(scale,static_cast<std::uint8_t>((1U<<level)^j));
    }
    return sfMultiply(value,inverse(scale));
}
static_assert(normalized(4,0)==0&&normalized(4,32)==83&&normalized(5,0)==0);

template<bool Affine,unsigned Coefficient>
RS4N_INLINE void pair(__m512i& low,__m512i& high)noexcept{
    low=_mm512_xor_si512(low,high);
    if constexpr(Coefficient==1)high=_mm512_xor_si512(high,low);
    else if constexpr(Coefficient!=0)
        high=_mm512_xor_si512(high,rs8ConstantMultiply<Affine,Coefficient>(low));
}

template<unsigned Offset,unsigned Index,bool Affine,bool Stream,unsigned D>
RS4N_INLINE void emit(__m512i value,__m512i* scratch,Block* output)noexcept{
    if constexpr(Offset==0||Offset==32){
        _mm512_store_si512(scratch+Offset+Index,value);
    }else if constexpr(Offset==16){
        // stage16_0_1: the high result stays unchanged.
        const auto low=_mm512_load_si512(scratch+Index);
        _mm512_store_si512(scratch+Index,_mm512_xor_si512(low,value));
        _mm512_store_si512(scratch+16+Index,value);
    }else{
        static_assert(Offset==48&&D>=32&&D<=48);
        auto low=_mm512_load_si512(scratch+32+Index);
        pair<Affine,83>(low,value);
        // Fuse stage16_83_1 and stage32_0_1 into the retained stores.
        // The last sixteen coefficients are punctured for these profiles.
        sf16Store<Stream>(output+4*Index,
            _mm512_xor_si512(_mm512_load_si512(scratch+Index),low));
        sf16Store<Stream>(output+4*(16+Index),
            _mm512_xor_si512(_mm512_load_si512(scratch+16+Index),value));
        if constexpr(32+Index<D)sf16Store<Stream>(output+4*(32+Index),low);
    }
}

// Explicit sixteen-register leaf in the frozen circuit's butterfly order.
// Its right-leaf outputs feed the ancestors without a scratch round trip.
// PrefetchBytes is either four columns ahead or zero at the plane boundary.
template<unsigned Offset,bool Affine,bool Prefetch,bool Stream,unsigned D>
RS4N_NOINLINE void leaf(const unsigned char* source,std::size_t stride,
                       __m512i* scratch,Block* output,std::size_t prefetchBytes){
    static_assert(Offset==0||Offset==16||Offset==32||Offset==48);
#define RS4N_LOAD(I) \
    if constexpr(Prefetch)_mm_prefetch(reinterpret_cast<const char*>(source+(Offset+I)*stride+prefetchBytes),_MM_HINT_T0); \
    auto x##I=_mm512_loadu_si512(source+(Offset+I)*stride)
    RS4N_LOAD(0);RS4N_LOAD(1);RS4N_LOAD(2);RS4N_LOAD(3);
    RS4N_LOAD(4);RS4N_LOAD(5);RS4N_LOAD(6);RS4N_LOAD(7);
    RS4N_LOAD(8);RS4N_LOAD(9);RS4N_LOAD(10);RS4N_LOAD(11);
    RS4N_LOAD(12);RS4N_LOAD(13);RS4N_LOAD(14);RS4N_LOAD(15);
#undef RS4N_LOAD
#define RS4N_PAIR(A,B,LEVEL,DELTA) pair<Affine,normalized(LEVEL,Offset^DELTA)>(x##A,x##B)
    RS4N_PAIR(14,15,0,14);RS4N_PAIR(12,13,0,12);
    RS4N_PAIR(13,15,1,12);RS4N_PAIR(12,14,1,12);
    RS4N_PAIR(10,11,0,10);RS4N_PAIR(8,9,0,8);
    RS4N_PAIR(9,11,1,8);RS4N_PAIR(8,10,1,8);
    RS4N_PAIR(11,15,2,8);RS4N_PAIR(10,14,2,8);
    RS4N_PAIR(9,13,2,8);RS4N_PAIR(8,12,2,8);
    RS4N_PAIR(6,7,0,6);RS4N_PAIR(4,5,0,4);
    RS4N_PAIR(5,7,1,4);RS4N_PAIR(4,6,1,4);
    RS4N_PAIR(2,3,0,2);RS4N_PAIR(0,1,0,0);
    RS4N_PAIR(1,3,1,0);RS4N_PAIR(0,2,1,0);
    RS4N_PAIR(3,7,2,0);RS4N_PAIR(2,6,2,0);
    RS4N_PAIR(1,5,2,0);RS4N_PAIR(0,4,2,0);
    RS4N_PAIR(7,15,3,0);RS4N_PAIR(6,14,3,0);
    RS4N_PAIR(5,13,3,0);RS4N_PAIR(4,12,3,0);
    RS4N_PAIR(3,11,3,0);RS4N_PAIR(2,10,3,0);
    RS4N_PAIR(1,9,3,0);RS4N_PAIR(0,8,3,0);
#undef RS4N_PAIR
#define RS4N_EMIT(I) emit<Offset,I,Affine,Stream,D>(x##I,scratch,output)
    RS4N_EMIT(0);RS4N_EMIT(1);RS4N_EMIT(2);RS4N_EMIT(3);
    RS4N_EMIT(4);RS4N_EMIT(5);RS4N_EMIT(6);RS4N_EMIT(7);
    RS4N_EMIT(8);RS4N_EMIT(9);RS4N_EMIT(10);RS4N_EMIT(11);
    RS4N_EMIT(12);RS4N_EMIT(13);RS4N_EMIT(14);RS4N_EMIT(15);
#undef RS4N_EMIT
}

// One cold shortened tile uses the frozen prefix, keeping its output-store
// order and exact 1..3-record final memcpy. Its 4KiB scratch is separate
// from the hot caller's 3KiB scratch and never enters a coroutine frame.
template<unsigned D,bool Stream>
RS4N_NOINLINE void tail(const Block* planes,std::size_t stride,Block* output,
                       std::size_t remaining){
    alignas(64) __m512i values[64];
    rs4OuterPrefix(planes,stride,values);
    for(unsigned j=0;j!=D&&4*j<remaining;++j){
        auto value=_mm512_load_si512(values+j);
        if(j<32)value=_mm512_xor_si512(value,_mm512_load_si512(values+32+j));
        const auto left=remaining-4*j;
        if(left>=4)sf16Store<Stream>(output+4*j,value);
        else{
            alignas(64) Block records[4];_mm512_store_si512(records,value);
            std::memcpy(output+4*j,records,left*sizeof(Block));
        }
    }
}
} // namespace rs4_outer_next_detail

// Same arguments, source indexing, and alignment contract as rs4Outer.
// planes/output must be disjoint; Stream=true requires aligned output.
// This is an ordinary non-coroutine helper with caller-independent scratch.
template<unsigned D,bool Stream=true,
         Rs4OuterNextVariant Variant=Rs4OuterNextVariant::FusedMultiply>
RS4N_NOINLINE void rs4OuterNext(const RateSpeed4Plan<D>& p,const Block* planes,
                              Block* output){
    static_assert(D==39||D==40);
    constexpr bool affine=Variant!=Rs4OuterNextVariant::FusedMultiply;
    constexpr bool prefetch=Variant==Rs4OuterNextVariant::FusedAffinePrefetch;
    alignas(64) __m512i scratch[48];
    const auto stride=p.planeStride*sizeof(Block);
    const auto fullRows=(p.requested/(4*D))*4;
    for(std::size_t row=0;row!=fullRows;row+=4){
        const auto* source=reinterpret_cast<const unsigned char*>(planes+row);
        auto* destination=output+row*D;
        const auto ahead=prefetch&&row+16<p.groups?16*sizeof(Block):0;
        rs4_outer_next_detail::leaf<0,affine,prefetch,Stream,D>(source,stride,scratch,destination,ahead);
        rs4_outer_next_detail::leaf<16,affine,prefetch,Stream,D>(source,stride,scratch,destination,ahead);
        rs4_outer_next_detail::leaf<32,affine,prefetch,Stream,D>(source,stride,scratch,destination,ahead);
        rs4_outer_next_detail::leaf<48,affine,prefetch,Stream,D>(source,stride,scratch,destination,ahead);
    }
    const auto base=fullRows*D;
    if(base<p.requested)
        rs4_outer_next_detail::tail<D,Stream>(planes+fullRows,stride,output+base,p.requested-base);
    if constexpr(Stream)_mm_sfence();
}

// Selected cached eight-packet inner followed by one candidate outer. The
// inner completes before output writes, preserving input/output aliasing.
template<unsigned D,Rs4OuterNextVariant Variant=Rs4OuterNextVariant::FusedMultiply>
inline void rs4TransposeCachedNext(const RateSpeed4Plan<D>& p,const Block* input,
                                  Block* output,Block* planes,Block* temp){
    rs4RegionInner<true,D,false,true,false,true>(p,input,planes,temp);
    if((reinterpret_cast<std::uintptr_t>(output)&63)==0)
        rs4OuterNext<D,true,Variant>(p,planes,output);
    else rs4OuterNext<D,false,Variant>(p,planes,output);
}

#undef RS4N_INLINE
#undef RS4N_NOINLINE
} // namespace spin::research::char2
