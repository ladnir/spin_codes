#pragma once
#include "RateClimb8Code.h"

namespace spin::research::char2 {

// Exact-map A/B prototype. Reuse every frozen FFT leaf/stage except the
// final stage64_0_1. That stage sends (a,b) to (a^b,b); the retained results
// can therefore be stored directly without materializing the final scratch.
// No generated arithmetic is copied, pruned, or changed here.
template<bool Fused>
inline void r8OuterPrefixFinalStore(const Block* planes,std::size_t strideBytes,
                                   __m512i* values){
    using namespace c128_outer_coefficient_detail;
    const auto* source=reinterpret_cast<const unsigned char*>(planes);
    if constexpr(!Fused){
        c8Unroll(std::make_index_sequence<128>{},[&](auto j){
            values[j]=_mm512_loadu_si512(source+std::size_t(j)*strideBytes);
        });
    }
    if constexpr(Fused){
        leaf16_0_strided_transpose(source,strideBytes,values);
        leaf16_16_strided_transpose(source+16*strideBytes,strideBytes,values+16);
    }else{
        leaf16_0_0_1(values);leaf16_16_0_1(values+16);
    }
    stage16_0_1(values);
    if constexpr(Fused){
        leaf16_32_strided_transpose(source+32*strideBytes,strideBytes,values+32);
        leaf16_48_strided_transpose(source+48*strideBytes,strideBytes,values+48);
    }else{
        leaf16_32_0_1(values+32);leaf16_48_0_1(values+48);
    }
    stage16_83_1(values+32);
    stage32_0_1(values);
    if constexpr(Fused){
        leaf16_64_strided_transpose(source+64*strideBytes,strideBytes,values+64);
        leaf16_80_strided_transpose(source+80*strideBytes,strideBytes,values+80);
    }else{
        leaf16_64_0_1(values+64);leaf16_80_0_1(values+80);
    }
    stage16_215_1(values+64);
    if constexpr(Fused){
        leaf16_96_strided_transpose(source+96*strideBytes,strideBytes,values+96);
        leaf16_112_strided_transpose(source+112*strideBytes,strideBytes,values+112);
    }else{
        leaf16_96_0_1(values+96);leaf16_112_0_1(values+112);
    }
    stage16_132_1(values+96);
    stage32_183_1(values+64);
}

template<unsigned Index>
inline __m512i r8OuterFinalCoefficient(const __m512i* values){
    static_assert(Index<128);
    const auto value=_mm512_load_si512(values+Index);
    if constexpr(Index<64)return _mm512_xor_si512(value,_mm512_load_si512(values+64+Index));
    else return value;
}

#if defined(_MSC_VER)
#define R8_FINAL_STORE_NOINLINE __declspec(noinline)
#else
#define R8_FINAL_STORE_NOINLINE __attribute__((noinline))
#endif

// Ordinary non-coroutine helper: 8KiB SIMD scratch remains on its stack.
// The source planes and output must be disjoint. Stream=true additionally
// requires 64-byte-aligned output, as in the existing sf16Store contract.
// Exactly requested records are written, including a shortened final tile.
template<bool Stream,bool Fused,class Profile>
R8_FINAL_STORE_NOINLINE void r8OuterFinishFinalStore(const RateClimb8Plan<Profile>& p,
                                                    const Block* planes,Block* out){
    constexpr unsigned dimension=Profile::dimension;
    static_assert(dimension>=64&&dimension<=128);
    alignas(64) __m512i values[128];
    for(std::size_t row=0;row*dimension<p.requested;row+=4){
        r8OuterPrefixFinalStore<Fused>(planes+row,p.planeStride*sizeof(Block),values);
        const auto base=row*dimension;
        if(base+4*dimension<=p.requested){
            c8Unroll(std::make_index_sequence<dimension>{},[&](auto j){
                sf16Store<Stream>(out+base+4*j,
                    r8OuterFinalCoefficient<static_cast<unsigned>(decltype(j)::value)>(values));
            });
        }else{
            // Only the final shortened tile takes this bounded cold path.
            for(unsigned j=0;j!=dimension&&base+4*j<p.requested;++j){
                auto value=_mm512_load_si512(values+j);
                if(j<64)value=_mm512_xor_si512(value,_mm512_load_si512(values+64+j));
                const auto left=p.requested-base-4*j;
                if(left>=4)sf16Store<Stream>(out+base+4*j,value);
                else{
                    alignas(64) Block tail[4];_mm512_store_si512(tail,value);
                    std::memcpy(out+base+4*j,tail,left*sizeof(Block));
                }
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}

#undef R8_FINAL_STORE_NOINLINE

template<class Profile>
inline void r8OuterDispatchFinalStore(const RateClimb8Plan<Profile>& p,
                                      const Block* planes,Block* out,
                                      Subfield16Kernel kernel,Cubic8OuterVariant variant){
    const bool stream=kernel==Subfield16Kernel::RegionAllNonTemporal&&
                      (reinterpret_cast<std::uintptr_t>(out)&63)==0;
    if(variant==Cubic8OuterVariant::FusedLoad){
        if(stream)r8OuterFinishFinalStore<true,true>(p,planes,out);
        else r8OuterFinishFinalStore<false,true>(p,planes,out);
    }else{
        if(stream)r8OuterFinishFinalStore<true,false>(p,planes,out);
        else r8OuterFinishFinalStore<false,false>(p,planes,out);
    }
}

} // namespace spin::research::char2
