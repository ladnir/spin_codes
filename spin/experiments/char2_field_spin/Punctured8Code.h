#pragma once
#include "Punctured8Inner.h"
#include "Cubic8Code.h"

namespace spin::research::char2 {

// Research profile: GRS[128,84], full-Cauchy width8/state8, keep7 outputs.
// The sampled scalar GF256 code is shared across sixteen independent lanes.
// Setup fixtures are outside encoding timings; no production API is changed.
struct Punctured8Plan {
    static constexpr unsigned dimension=84,outer=128,width=8,state=8;
    std::size_t requested,k,inputN,outputN,groups,planeStride;
    std::vector<std::uint32_t> route,inverseRoute;
    std::vector<std::uint16_t> inverseLocalRoute;
    std::vector<std::uint8_t> inputDiagonal,stateDiagonal;
    explicit Punctured8Plan(std::size_t message,std::uint64_t seed):requested(message){
        constexpr std::size_t unit=dimension*64;
        constexpr std::size_t maxRows=(std::numeric_limits<std::uint32_t>::max()/outer/64)*64;
        if(!message||message>maxRows*dimension)
            throw std::invalid_argument("Punctured8 message size outside uint32 route geometry");
        groups=((message+unit-1)/unit)*64;k=groups*dimension;
        inputN=groups*outer;outputN=groups*112;
        // Optional scratch-stride padding is an exact-map benchmark control.
        // It lost the measured K20 A/B, so compact columns remain the default.
        planeStride=groups;
        std::mt19937_64 rr(seed^0xdf4a1b3285ec7096ULL),ir(seed^0x2e7d18c4f0b3a965ULL),
                        sr(seed^0x918b745cf62d03eaULL);
        route.resize(inputN);inverseRoute.resize(inputN);
        inputDiagonal.resize(inputN);stateDiagonal.resize(inputN);
        for(unsigned r=0;r!=outer;++r){
            auto first=route.begin()+r*groups;
            std::iota(first,first+groups,static_cast<std::uint32_t>(r*groups));
            for(std::size_t j=groups-1;j;--j){
                std::uniform_int_distribution<std::size_t> draw(0,j);
                std::swap(first[j],first[draw(rr)]);
            }
        }
        for(std::size_t i=0;i!=inputN;++i)inverseRoute[route[i]]=static_cast<std::uint32_t>(i);
        if(groups<=65536){
            inverseLocalRoute.resize(inputN);
            for(unsigned r=0;r!=outer;++r){
                const auto offset=r*groups;
                for(std::size_t j=0;j!=groups;++j)
                    inverseLocalRoute[offset+j]=static_cast<std::uint16_t>(inverseRoute[offset+j]-offset);
            }
        }
        for(auto& d:inputDiagonal)d=sfSampleNonzero(ir);
        for(auto& d:stateDiagonal)d=sfSampleNonzero(sr);
    }
    static std::size_t sourceIndex(std::size_t row,unsigned coordinate)noexcept{
        return (row/4)*336+4*coordinate+row%4;
    }
    static std::size_t outputIndex(std::size_t step,unsigned coordinate)noexcept{
        return 28*(step/4)+4*coordinate+step%4;
    }
    std::size_t planeRecords()const noexcept{return 128*planeStride;}
};

template<bool Scatter>
inline void p8InnerRange(const Punctured8Plan& p,const Block* input,Block* dest,
                          std::size_t first,std::size_t past,std::size_t base,
                          __m512i& h0,__m512i& h1){
    alignas(64) Block packet[32];
    for(auto step=past;step!=first;){
        step-=4;const auto begin=8*step;
        p8InnerTransposeFour(input+7*step,packet,h0,h1,p.stateDiagonal.data()+begin);
        // Explicitly unrolled vector stores, no dynamic callback in this kernel.
        c8Unroll(std::make_index_sequence<8>{},[&](auto chunk){
            constexpr unsigned shift=4*unsigned(chunk);
            const auto value=_mm512_gf2p8mul_epi8(_mm512_load_si512(packet+shift),
                sfBroadcastFour(p.inputDiagonal.data()+begin+shift));
            if constexpr(Scatter){
                const auto* index=p.route.data()+begin+shift;
                dest[index[0]-base]=block(_mm512_castsi512_si128(value));
                dest[index[1]-base]=block(_mm512_extracti32x4_epi32(value,1));
                dest[index[2]-base]=block(_mm512_extracti32x4_epi32(value,2));
                dest[index[3]-base]=block(_mm512_extracti32x4_epi32(value,3));
            }else _mm512_storeu_si512(dest+begin-base+shift,value);
        });
    }
}
template<bool Stream,bool Compact>
inline void p8RegionInner(const Punctured8Plan& p,const Block* input,Block* planes,Block* temp){
    auto h0=_mm512_setzero_si512(),h1=_mm512_setzero_si512();
    for(unsigned r=128;r--;){
        const auto offset=r*p.groups;
        auto* region=planes+r*p.planeStride;
        p8InnerRange<false>(p,input,temp,offset/8,(offset+p.groups)/8,offset,h0,h1);
        if constexpr(Compact){
            for(std::size_t j=0;j!=p.groups;j+=4)
                sf16Store<Stream>(region+j,c8GatherLocal4(temp,p.inverseLocalRoute.data()+offset+j));
        }else{
            for(std::size_t j=0;j!=p.groups;j+=4)
                sf16Store<Stream>(region+j,sf16GatherRegion4(temp,p.inverseRoute.data()+offset+j,
                    static_cast<std::uint32_t>(offset)));
        }
    }
    if constexpr(Stream)_mm_sfence();
}
inline void p8InnerDispatch(const Punctured8Plan& p,const Block* input,Block* planes,
                            Block* temp,Subfield16Kernel kernel){
    if(kernel==Subfield16Kernel::Scatter){
        auto h0=_mm512_setzero_si512(),h1=_mm512_setzero_si512();
        for(unsigned r=128;r--;){
            const auto offset=r*p.groups;
            p8InnerRange<true>(p,input,planes+r*p.planeStride,offset/8,(offset+p.groups)/8,offset,h0,h1);
        }
    }else{
        if(!temp)throw std::invalid_argument("Punctured8 regional schedule requires temporary storage");
        const bool aligned=(reinterpret_cast<std::uintptr_t>(planes)&63)==0;
        if(p.inverseLocalRoute.empty()){
            if(aligned)p8RegionInner<true,false>(p,input,planes,temp);
            else p8RegionInner<false,false>(p,input,planes,temp);
        }else{
            if(aligned)p8RegionInner<true,true>(p,input,planes,temp);
            else p8RegionInner<false,true>(p,input,planes,temp);
        }
    }
}
#if defined(_MSC_VER)
#define P8_NOINLINE __declspec(noinline)
#else
#define P8_NOINLINE __attribute__((noinline))
#endif
template<bool Stream,bool Fused>
P8_NOINLINE void p8OuterFinish(const Punctured8Plan& p,const Block* planes,Block* out){
    alignas(64) __m512i values[128];
    for(std::size_t row=0;row*84<p.requested;row+=4){
        if constexpr(Fused)c128OuterCoefficientTransposeStrided(planes+row,p.planeStride*sizeof(Block),values);
        else{
            c8Unroll(std::make_index_sequence<128>{},[&](auto j){
                values[j]=_mm512_loadu_si512(planes+std::size_t(j)*p.planeStride+row);
            });
            c128OuterCoefficientTransposeInplace(values);
        }
        const auto base=row*84;
        if(base+336<=p.requested){
            c8Unroll(std::make_index_sequence<84>{},[&](auto j){sf16Store<Stream>(out+base+4*j,values[j]);});
        }else{
            for(unsigned j=0;j!=84&&base+4*j<p.requested;++j){
                const auto left=p.requested-base-4*j;
                if(left>=4)sf16Store<Stream>(out+base+4*j,values[j]);
                else{alignas(64) Block tail[4];_mm512_store_si512(tail,values[j]);
                     std::memcpy(out+base+4*j,tail,left*sizeof(Block));}
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}
#undef P8_NOINLINE
inline void p8OuterDispatch(const Punctured8Plan& p,const Block* planes,Block* out,
                            Subfield16Kernel kernel,Cubic8OuterVariant variant){
    const bool stream=kernel==Subfield16Kernel::RegionAllNonTemporal&&(reinterpret_cast<std::uintptr_t>(out)&63)==0;
    if(variant==Cubic8OuterVariant::FusedLoad){
        if(stream)p8OuterFinish<true,true>(p,planes,out);else p8OuterFinish<false,true>(p,planes,out);
    }else{
        if(stream)p8OuterFinish<true,false>(p,planes,out);else p8OuterFinish<false,false>(p,planes,out);
    }
}
// Full input is consumed into disjoint planes before any output write.
// Exact in/out alias is supported, with the input suffix preserved.
inline void p8Transpose(const Punctured8Plan& p,const Block* input,Block* out,Block* planes,
                         Block* temp,Subfield16Kernel kernel,
                         Cubic8OuterVariant variant=Cubic8OuterVariant::Gather){
    p8InnerDispatch(p,input,planes,temp,kernel);p8OuterDispatch(p,planes,out,kernel,variant);
}
} // namespace spin::research::char2
