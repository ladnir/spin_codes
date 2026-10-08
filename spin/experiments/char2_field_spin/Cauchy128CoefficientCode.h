#pragma once
#include "Cauchy128.h"
#include "Cauchy128OuterCoefficients.h"
#if defined(_MSC_VER)
#define C128_COEFFICIENT_NOINLINE __declspec(noinline)
#else
#define C128_COEFFICIENT_NOINLINE __attribute__((noinline))
#endif

namespace spin::research::char2 {

enum class Cauchy128CoefficientVariant : std::uint8_t { Gather, FusedLoad };

// Same full code image and setup as Cauchy128Plan. The source coordinates are
// normalized additive-polynomial coefficients, ordered in four-row tiles.
// This source-basis change deliberately changes the transpose's output bytes.
template<bool Stream,bool FusedLoad=false>
C128_COEFFICIENT_NOINLINE void c128CoefficientOuterFinish(
    const Cauchy128Plan& p,const Block* planes,Block* out){
    alignas(64) __m512i values[128];
    for(std::size_t row=0;row*87<p.requested;row+=4){
        if constexpr(FusedLoad){
            c128OuterCoefficientTransposeStrided(planes+row,p.groups*sizeof(Block),values);
        }else{
            c128Unroll(std::make_index_sequence<128>{},[&](auto j){
                values[j]=_mm512_loadu_si512(planes+std::size_t(j)*p.groups+row);
            });
            c128OuterCoefficientTransposeInplace(values);
        }
        const auto base=row*87;
        if(base+348<=p.requested){
            c128Unroll(std::make_index_sequence<87>{},[&](auto j){
                sf16Store<Stream>(out+base+4*j,values[j]);
            });
        }else{
            for(unsigned j=0;j!=87&&base+4*j<p.requested;++j){
                const auto left=p.requested-base-4*j;
                if(left>=4)sf16Store<Stream>(out+base+4*j,values[j]);
                else{
                    alignas(64) Block tail[4];_mm512_store_si512(tail,values[j]);
                    std::memcpy(out+base+4*j,tail,left*sizeof(Block));
                }
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}
#undef C128_COEFFICIENT_NOINLINE

inline void c128CoefficientOuterDispatch(const Cauchy128Plan& p,const Block* planes,
                                         Block* out,Subfield16Kernel kernel,
                                         Cauchy128CoefficientVariant variant=Cauchy128CoefficientVariant::Gather){
    const bool stream=kernel==Subfield16Kernel::RegionAllNonTemporal&&
                      (reinterpret_cast<std::uintptr_t>(out)&63)==0;
    if(variant==Cauchy128CoefficientVariant::FusedLoad){
        if(stream)c128CoefficientOuterFinish<true,true>(p,planes,out);
        else c128CoefficientOuterFinish<false,true>(p,planes,out);
    }else if(stream)
        c128CoefficientOuterFinish<true>(p,planes,out);
    else c128CoefficientOuterFinish<false>(p,planes,out);
}

// Caller-owned scratch. Exact input/output alias is supported; all scratch
// is disjoint. Only requested outputs are written, preserving an N-buffer tail.
inline void c128CoefficientTranspose(const Cauchy128Plan& p,const Block* in,
                                     Block* out,Block* planes,Block* temp,
                                     Subfield16Kernel kernel,
                                     Cauchy128CoefficientVariant variant=Cauchy128CoefficientVariant::Gather){
    c128InnerDispatch(p,in,planes,temp,kernel,Cauchy128Variant::CompactBatch);
    c128CoefficientOuterDispatch(p,planes,out,kernel,variant);
}
}
