#pragma once
#include "RateClimb8Code.h"

namespace spin::research::char2 {

// Research A/B: finish the ancestors of each rightmost16-point leaf while
// its outputs are still in registers. The scalar map and source basis are
// unchanged. Frozen left leaves are reused; no generated header is edited.
constexpr std::uint8_t r8sInverse(std::uint8_t value)noexcept{
    std::uint8_t result=1;
    for(unsigned exponent=254;exponent;exponent>>=1){
        if(exponent&1)result=sfMultiply(result,value);
        value=sfMultiply(value,value);
    }
    return result;
}
constexpr std::uint8_t r8sNormalized(unsigned level,unsigned offset)noexcept{
    std::uint8_t value=1,scale=1;
    for(unsigned j=0;j!=(1U<<level);++j){
        value=sfMultiply(value,static_cast<std::uint8_t>(offset^j));
        scale=sfMultiply(scale,static_cast<std::uint8_t>((1U<<level)^j));
    }
    return sfMultiply(value,r8sInverse(scale));
}
static_assert(r8sNormalized(4,32)==83&&r8sNormalized(4,64)==215&&
              r8sNormalized(4,96)==132&&r8sNormalized(5,64)==183);

#if defined(_MSC_VER)
#define R8S_INLINE __forceinline
#define R8S_NOINLINE __declspec(noinline)
#else
#define R8S_INLINE inline __attribute__((always_inline))
#define R8S_NOINLINE __attribute__((noinline))
#endif

template<unsigned Coefficient>
R8S_INLINE void r8sPair(__m512i& low,__m512i& high)noexcept{
    low=_mm512_xor_si512(low,high);
    if constexpr(Coefficient==1)high=_mm512_xor_si512(high,low);
    else if constexpr(Coefficient!=0)
        high=_mm512_xor_si512(high,_mm512_gf2p8mul_epi8(low,
             _mm512_set1_epi8(static_cast<char>(Coefficient))));
}

template<unsigned Index,bool Stream,bool Full,unsigned Dimension>
R8S_INLINE void r8sStore(Block* output,std::size_t remaining,__m512i value)noexcept{
    if constexpr(Index<Dimension){
        if constexpr(Full)sf16Store<Stream>(output+4*Index,value);
        else if(4*Index<remaining){
            const auto count=remaining-4*Index;
            if(count>=4)sf16Store<Stream>(output+4*Index,value);
            else{
                alignas(64) Block tail[4];_mm512_store_si512(tail,value);
                std::memcpy(output+4*Index,tail,count*sizeof(Block));
            }
        }
    }
}

template<unsigned Offset,unsigned Index,bool Stream,bool Full,unsigned Dimension>
R8S_INLINE void r8sEmit(__m512i high,__m512i* scratch,Block* output,
                       std::size_t remaining)noexcept{
    static_assert(Index<16&&(Offset==16||Offset==48||Offset==80||Offset==112));
    if constexpr(Offset==16||Offset==80){
        auto low=_mm512_load_si512(scratch+Offset-16+Index);
        r8sPair<r8sNormalized(4,Offset-16)>(low,high);
        _mm512_store_si512(scratch+Offset-16+Index,low);
        _mm512_store_si512(scratch+Offset+Index,high);
    }else if constexpr(Offset==48){
        auto rightLow=_mm512_load_si512(scratch+32+Index);
        r8sPair<83>(rightLow,high);
        const auto leftLow=_mm512_load_si512(scratch+Index);
        const auto leftHigh=_mm512_load_si512(scratch+16+Index);
        _mm512_store_si512(scratch+Index,_mm512_xor_si512(leftLow,rightLow));
        _mm512_store_si512(scratch+16+Index,_mm512_xor_si512(leftHigh,high));
        _mm512_store_si512(scratch+32+Index,rightLow);
        _mm512_store_si512(scratch+48+Index,high);
    }else{
        auto rightLow=_mm512_load_si512(scratch+96+Index);
        r8sPair<132>(rightLow,high);
        auto low0=_mm512_load_si512(scratch+64+Index);
        r8sPair<183>(low0,rightLow);
        auto low1=_mm512_load_si512(scratch+80+Index);
        r8sPair<183>(low1,high);
        r8sStore<Index,Stream,Full,Dimension>(output,remaining,
            _mm512_xor_si512(_mm512_load_si512(scratch+Index),low0));
        r8sStore<16+Index,Stream,Full,Dimension>(output,remaining,
            _mm512_xor_si512(_mm512_load_si512(scratch+16+Index),low1));
        r8sStore<32+Index,Stream,Full,Dimension>(output,remaining,
            _mm512_xor_si512(_mm512_load_si512(scratch+32+Index),rightLow));
        r8sStore<48+Index,Stream,Full,Dimension>(output,remaining,
            _mm512_xor_si512(_mm512_load_si512(scratch+48+Index),high));
        r8sStore<64+Index,Stream,Full,Dimension>(output,remaining,low0);
        r8sStore<80+Index,Stream,Full,Dimension>(output,remaining,low1);
        r8sStore<96+Index,Stream,Full,Dimension>(output,remaining,rightLow);
        r8sStore<112+Index,Stream,Full,Dimension>(output,remaining,high);
    }
}

// Explicit16-register leaf. Every butterfly and coefficient is compile-time;
// no register array, runtime butterfly loop, callback, or heap allocation.
template<unsigned Offset,bool Fused,bool Stream=false,bool Full=true,unsigned Dimension=0>
R8S_NOINLINE void r8sRightLeaf(const Block* planes,std::size_t strideBytes,
                              __m512i* scratch,Block* output=nullptr,
                              std::size_t remaining=0){
    static_assert(Offset==16||Offset==48||Offset==80||Offset==112);
    static_assert(Offset!=112||(Dimension>=64&&Dimension<=128));
    const auto* source=reinterpret_cast<const unsigned char*>(planes);
#define R8S_LOAD(I) auto x##I=Fused ? _mm512_loadu_si512(source+(Offset+I)*strideBytes) : _mm512_load_si512(scratch+Offset+I)
    R8S_LOAD(0);R8S_LOAD(1);R8S_LOAD(2);R8S_LOAD(3);
    R8S_LOAD(4);R8S_LOAD(5);R8S_LOAD(6);R8S_LOAD(7);
    R8S_LOAD(8);R8S_LOAD(9);R8S_LOAD(10);R8S_LOAD(11);
    R8S_LOAD(12);R8S_LOAD(13);R8S_LOAD(14);R8S_LOAD(15);
#undef R8S_LOAD
#define R8S_PAIR(A,B,LEVEL,DELTA) r8sPair<r8sNormalized(LEVEL,Offset^DELTA)>(x##A,x##B)
    R8S_PAIR(14,15,0,14);R8S_PAIR(12,13,0,12);
    R8S_PAIR(13,15,1,12);R8S_PAIR(12,14,1,12);
    R8S_PAIR(10,11,0,10);R8S_PAIR(8,9,0,8);
    R8S_PAIR(9,11,1,8);R8S_PAIR(8,10,1,8);
    R8S_PAIR(11,15,2,8);R8S_PAIR(10,14,2,8);
    R8S_PAIR(9,13,2,8);R8S_PAIR(8,12,2,8);
    R8S_PAIR(6,7,0,6);R8S_PAIR(4,5,0,4);
    R8S_PAIR(5,7,1,4);R8S_PAIR(4,6,1,4);
    R8S_PAIR(2,3,0,2);R8S_PAIR(0,1,0,0);
    R8S_PAIR(1,3,1,0);R8S_PAIR(0,2,1,0);
    R8S_PAIR(3,7,2,0);R8S_PAIR(2,6,2,0);
    R8S_PAIR(1,5,2,0);R8S_PAIR(0,4,2,0);
    R8S_PAIR(7,15,3,0);R8S_PAIR(6,14,3,0);
    R8S_PAIR(5,13,3,0);R8S_PAIR(4,12,3,0);
    R8S_PAIR(3,11,3,0);R8S_PAIR(2,10,3,0);
    R8S_PAIR(1,9,3,0);R8S_PAIR(0,8,3,0);
#undef R8S_PAIR
#define R8S_EMIT(I) r8sEmit<Offset,I,Stream,Full,Dimension>(x##I,scratch,output,remaining)
    R8S_EMIT(0);R8S_EMIT(1);R8S_EMIT(2);R8S_EMIT(3);
    R8S_EMIT(4);R8S_EMIT(5);R8S_EMIT(6);R8S_EMIT(7);
    R8S_EMIT(8);R8S_EMIT(9);R8S_EMIT(10);R8S_EMIT(11);
    R8S_EMIT(12);R8S_EMIT(13);R8S_EMIT(14);R8S_EMIT(15);
#undef R8S_EMIT
}

template<bool Stream,bool Fused,class Profile>
R8S_NOINLINE void r8SpeedOuterFinish(const RateClimb8Plan<Profile>& p,
                                    const Block* planes,Block* out){
    constexpr unsigned dimension=Profile::dimension;
    static_assert(dimension>=64&&dimension<=128);
    using namespace c128_outer_coefficient_detail;
    alignas(64) __m512i scratch[128];
    const auto stride=p.planeStride*sizeof(Block);
    for(std::size_t row=0;row*dimension<p.requested;row+=4){
        const auto* input=planes+row;
        const auto* source=reinterpret_cast<const unsigned char*>(input);
        if constexpr(!Fused)c8Unroll(std::make_index_sequence<128>{},[&](auto j){
            scratch[j]=_mm512_loadu_si512(source+std::size_t(j)*stride);
        });
        if constexpr(Fused)leaf16_0_strided_transpose(source,stride,scratch);
        else leaf16_0_0_1(scratch);
        r8sRightLeaf<16,Fused>(input,stride,scratch);
        if constexpr(Fused)leaf16_32_strided_transpose(source+32*stride,stride,scratch+32);
        else leaf16_32_0_1(scratch+32);
        r8sRightLeaf<48,Fused>(input,stride,scratch);
        if constexpr(Fused)leaf16_64_strided_transpose(source+64*stride,stride,scratch+64);
        else leaf16_64_0_1(scratch+64);
        r8sRightLeaf<80,Fused>(input,stride,scratch);
        if constexpr(Fused)leaf16_96_strided_transpose(source+96*stride,stride,scratch+96);
        else leaf16_96_0_1(scratch+96);
        const auto base=row*dimension,remaining=p.requested-base;
        if(remaining>=4*dimension)
            r8sRightLeaf<112,Fused,Stream,true,dimension>(input,stride,scratch,out+base,remaining);
        else r8sRightLeaf<112,Fused,Stream,false,dimension>(input,stride,scratch,out+base,remaining);
    }
    if constexpr(Stream)_mm_sfence();
}

template<class Profile>
inline void r8SpeedOuterDispatch(const RateClimb8Plan<Profile>& p,const Block* planes,
                                 Block* out,Subfield16Kernel kernel,Cubic8OuterVariant variant){
    const bool stream=kernel==Subfield16Kernel::RegionAllNonTemporal&&
                      (reinterpret_cast<std::uintptr_t>(out)&63)==0;
    if(variant==Cubic8OuterVariant::FusedLoad){
        if(stream)r8SpeedOuterFinish<true,true>(p,planes,out);
        else r8SpeedOuterFinish<false,true>(p,planes,out);
    }else{
        if(stream)r8SpeedOuterFinish<true,false>(p,planes,out);
        else r8SpeedOuterFinish<false,false>(p,planes,out);
    }
}

#undef R8S_INLINE
#undef R8S_NOINLINE
} // namespace spin::research::char2
