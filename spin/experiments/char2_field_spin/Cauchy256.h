#pragma once
#include "Subfield16.h"
#include "Cauchy256Inner.h"
#include "Cauchy256InnerDense.h"
#include "Cauchy256InnerFft.h"
#include "Cauchy256Outer.h"
#include "Cauchy256OuterFft.h"
#include "Cauchy256OuterParityFft.h"
#include <utility>

namespace spin::research::char2 {

enum class Cauchy256Variant : std::uint8_t { Basis, Dense, Fft, DenseFft, DenseFftFused, FftInnerFused };

// Isolated field research profile. The certificate concerns ideal independent
// setup; these reproducible seeded fixtures are used for implementation tests.
struct Cauchy256Plan {
    static constexpr unsigned dimension=176, outer=256, width=16, state=8;
    std::size_t requested, k, n, groups;
    std::vector<std::uint32_t> route,inverseRoute;
    std::vector<std::uint8_t> inputDiagonal,stateDiagonal;

    explicit Cauchy256Plan(std::size_t message,std::uint64_t seed):requested(message) {
        constexpr std::size_t unit=dimension*64;
        constexpr std::size_t maxRows=(std::numeric_limits<std::uint32_t>::max()/256/64)*64;
        if(!message || message>maxRows*dimension)
            throw std::invalid_argument("Cauchy256 message size outside uint32 route geometry");
        groups=((message+unit-1)/unit)*64;
        k=groups*dimension;n=groups*outer;
        std::mt19937_64 routeRng(seed^0xdf4a1b3285ec7096ULL);
        std::mt19937_64 inputRng(seed^0x2e7d18c4f0b3a965ULL);
        std::mt19937_64 stateRng(seed^0x918b745cf62d03eaULL);
        route.resize(n);inverseRoute.resize(n);inputDiagonal.resize(n);stateDiagonal.resize(n/2);
        for(unsigned r=0;r!=outer;++r) {
            auto first=route.begin()+r*groups;
            std::iota(first,first+groups,static_cast<std::uint32_t>(r*groups));
            for(std::size_t j=groups-1;j;--j) {
                std::uniform_int_distribution<std::size_t> draw(0,j);
                std::swap(first[j],first[draw(routeRng)]);
            }
        }
        for(std::size_t i=0;i!=n;++i) inverseRoute[route[i]]=static_cast<std::uint32_t>(i);
        for(auto& d:inputDiagonal)d=sfSampleNonzero(inputRng);
        for(auto& d:stateDiagonal)d=sfSampleNonzero(stateRng);
    }
    std::size_t setupBytes()const noexcept {
        return (route.size()+inverseRoute.size())*4+inputDiagonal.size()+stateDiagonal.size();
    }
};

template<unsigned Batch,bool Scatter>
inline void c256StoreInner(const Cauchy256Plan& p,Block* dest,std::size_t begin,
                          std::size_t base,__m512i value) {
    constexpr unsigned shift=4*Batch;
    const auto i=begin+shift;
    const auto v=_mm512_gf2p8mul_epi8(value,sfBroadcastFour(p.inputDiagonal.data()+i));
    if constexpr(Scatter) {
        dest[p.route[i]]=block(_mm512_castsi512_si128(v));
        dest[p.route[i+1]]=block(_mm512_extracti32x4_epi32(v,1));
        dest[p.route[i+2]]=block(_mm512_extracti32x4_epi32(v,2));
        dest[p.route[i+3]]=block(_mm512_extracti32x4_epi32(v,3));
    } else _mm512_storeu_si512(dest+begin-base+shift,v);
}

template<unsigned Phase,bool Scatter,bool DenseInner,bool FftInner=false>
inline void c256Packet(const Cauchy256Plan& p,const Block* input,Block* dest,
                       std::size_t step,std::size_t base,__m512i* state) {
    static_assert(!DenseInner || !FftInner);
    const auto begin=16*step;
    const __m512i values[4]={_mm512_loadu_si512(input+begin),_mm512_loadu_si512(input+begin+4),
        _mm512_loadu_si512(input+begin+8),_mm512_loadu_si512(input+begin+12)};
    const __m512i diagonal[2]={sfBroadcastFour(p.stateDiagonal.data()+8*step),
        sfBroadcastFour(p.stateDiagonal.data()+8*step+4)};
    __m512i output[4];
    if constexpr(FftInner)c256InnerTransposeFftPacked<Phase>(values,output,state,diagonal);
    else if constexpr(DenseInner)c256InnerTransposeDensePacked<Phase>(values,output,state,diagonal);
    else c256InnerTransposePacked<Phase>(values,output,state,diagonal);
    c256StoreInner<0,Scatter>(p,dest,begin,base,output[0]);
    c256StoreInner<1,Scatter>(p,dest,begin,base,output[1]);
    c256StoreInner<2,Scatter>(p,dest,begin,base,output[2]);
    c256StoreInner<3,Scatter>(p,dest,begin,base,output[3]);
}

template<bool Scatter,bool DenseInner,bool FftInner=false>
inline void c256InnerRange(const Cauchy256Plan& p,const Block* input,Block* dest,
                          std::size_t first,std::size_t past,std::size_t base,__m512i* state) {
    // Region geometry is a multiple of 64 records: explicit four-phase dispatch.
    for(std::size_t step=past;step!=first;) {
        c256Packet<3,Scatter,DenseInner,FftInner>(p,input,dest,--step,base,state);
        c256Packet<2,Scatter,DenseInner,FftInner>(p,input,dest,--step,base,state);
        c256Packet<1,Scatter,DenseInner,FftInner>(p,input,dest,--step,base,state);
        c256Packet<0,Scatter,DenseInner,FftInner>(p,input,dest,--step,base,state);
    }
}

template<bool Stream,bool DenseInner,bool FftInner=false>
inline void c256RegionInner(const Cauchy256Plan& p,const Block* input,Block* planes,Block* temp) {
    __m512i state[2]={_mm512_setzero_si512(),_mm512_setzero_si512()};
    for(unsigned r=256;r--;) {
        const auto offset=r*p.groups;
        c256InnerRange<false,DenseInner,FftInner>(p,input,temp,offset/16,(offset+p.groups)/16,offset,state);
        for(std::size_t j=0;j!=p.groups;j+=4)
            sf16Store<Stream>(planes+offset+j,sf16GatherRegion4(temp,p.inverseRoute.data()+offset+j,
                static_cast<std::uint32_t>(offset)));
    }
    if constexpr(Stream)_mm_sfence();
}

template<bool DenseInner,bool FftInner=false>
inline void c256InnerDispatchVariant(const Cauchy256Plan& p,const Block* input,Block* planes,
                              Block* temp,Subfield16Kernel kernel) {
    if(kernel==Subfield16Kernel::Scatter) {
        __m512i state[2]={_mm512_setzero_si512(),_mm512_setzero_si512()};
        c256InnerRange<true,DenseInner,FftInner>(p,input,planes,0,p.n/16,0,state);
    } else {
        if(!temp)throw std::invalid_argument("Cauchy256 region schedule requires L temporary records");
        if((reinterpret_cast<std::uintptr_t>(planes)&63)==0)c256RegionInner<true,DenseInner,FftInner>(p,input,planes,temp);
        else c256RegionInner<false,DenseInner,FftInner>(p,input,planes,temp);
    }
}

inline void c256InnerDispatch(const Cauchy256Plan& p,const Block* input,Block* planes,
                              Block* temp,Subfield16Kernel kernel,
                              Cauchy256Variant variant=Cauchy256Variant::Basis) {
    switch(variant) {
    case Cauchy256Variant::Basis:
    case Cauchy256Variant::Fft:c256InnerDispatchVariant<false>(p,input,planes,temp,kernel);break;
    case Cauchy256Variant::Dense:
    case Cauchy256Variant::DenseFft:
    case Cauchy256Variant::DenseFftFused:c256InnerDispatchVariant<true>(p,input,planes,temp,kernel);break;
    case Cauchy256Variant::FftInnerFused:c256InnerDispatchVariant<false,true>(p,input,planes,temp,kernel);break;
    default:throw std::invalid_argument("Cauchy256 inner variant");
    }
}

template<class F,std::size_t... I>
inline void c256Unroll(std::index_sequence<I...>,F&& f){(f(std::integral_constant<std::size_t,I>{}),...);}

template<bool Stream,bool FftOuter,bool ParityOnly=false>
inline void c256OuterFinish(const Cauchy256Plan& p,const Block* planes,Block* output) {
    static_assert(!ParityOnly || FftOuter);
    // Fused parity mode stages 5KiB, rather than the full 16KiB input tile.
    alignas(64) __m512i input[ParityOnly?80:256],result[176];
    // Outputs corresponding only to padded message rows need not be computed.
    const std::size_t needed=(p.requested+175)/176;
    for(std::size_t group=0;group<needed;group+=4) {
        c256Unroll(std::make_index_sequence<ParityOnly?80:256>{},[&](auto i){
            constexpr std::size_t coordinate=std::size_t(i)+(ParityOnly?176:0);
            input[i]=_mm512_loadu_si512(planes+coordinate*p.groups+group);
        });
        if constexpr(ParityOnly)c256OuterParityTransposeFft(input,result);
        else if constexpr(FftOuter)c256OuterTransposeFft(input,result);
        else c256OuterTranspose(input,result);
        if((group+4)*176<=p.requested) {
            c256Unroll(std::make_index_sequence<44>{},[&](auto batch){
                constexpr unsigned j=4*batch;
                auto a=result[j],b=result[j+1],c=result[j+2],d=result[j+3];
                if constexpr(ParityOnly){
                    a=_mm512_xor_si512(a,_mm512_loadu_si512(planes+std::size_t(j)*p.groups+group));
                    b=_mm512_xor_si512(b,_mm512_loadu_si512(planes+std::size_t(j+1)*p.groups+group));
                    c=_mm512_xor_si512(c,_mm512_loadu_si512(planes+std::size_t(j+2)*p.groups+group));
                    d=_mm512_xor_si512(d,_mm512_loadu_si512(planes+std::size_t(j+3)*p.groups+group));
                }
                sf16TransposeRecords4(a,b,c,d);
                sf16Store<Stream>(output+176*group+j,a);
                sf16Store<Stream>(output+176*(group+1)+j,b);
                sf16Store<Stream>(output+176*(group+2)+j,c);
                sf16Store<Stream>(output+176*(group+3)+j,d);
            });
        } else {
            // One terminal tile; never overwrite [requested,N) in-place.
            alignas(64) Block tile[176*4];
            for(unsigned j=0;j!=176;j+=4) {
                auto a=result[j],b=result[j+1],c=result[j+2],d=result[j+3];
                if constexpr(ParityOnly){
                    a=_mm512_xor_si512(a,_mm512_loadu_si512(planes+std::size_t(j)*p.groups+group));
                    b=_mm512_xor_si512(b,_mm512_loadu_si512(planes+std::size_t(j+1)*p.groups+group));
                    c=_mm512_xor_si512(c,_mm512_loadu_si512(planes+std::size_t(j+2)*p.groups+group));
                    d=_mm512_xor_si512(d,_mm512_loadu_si512(planes+std::size_t(j+3)*p.groups+group));
                }
                sf16TransposeRecords4(a,b,c,d);
                _mm512_store_si512(tile+j,a);_mm512_store_si512(tile+176+j,b);
                _mm512_store_si512(tile+352+j,c);_mm512_store_si512(tile+528+j,d);
            }
            std::memcpy(output+176*group,tile,(p.requested-176*group)*sizeof(Block));
        }
    }
    if constexpr(Stream)_mm_sfence();
}

template<bool FftOuter,bool ParityOnly=false>
inline void c256OuterDispatchVariant(const Cauchy256Plan& p,const Block* planes,Block* output,
                              Subfield16Kernel kernel) {
    if(kernel==Subfield16Kernel::RegionAllNonTemporal && (reinterpret_cast<std::uintptr_t>(output)&63)==0)
        c256OuterFinish<true,FftOuter,ParityOnly>(p,planes,output);
    else c256OuterFinish<false,FftOuter,ParityOnly>(p,planes,output);
}

inline void c256OuterDispatch(const Cauchy256Plan& p,const Block* planes,Block* output,
                              Subfield16Kernel kernel,
                              Cauchy256Variant variant=Cauchy256Variant::Basis) {
    switch(variant) {
    case Cauchy256Variant::Basis:
    case Cauchy256Variant::Dense:c256OuterDispatchVariant<false>(p,planes,output,kernel);break;
    case Cauchy256Variant::Fft:
    case Cauchy256Variant::DenseFft:c256OuterDispatchVariant<true>(p,planes,output,kernel);break;
    case Cauchy256Variant::DenseFftFused:
    case Cauchy256Variant::FftInnerFused:c256OuterDispatchVariant<true,true>(p,planes,output,kernel);break;
    default:throw std::invalid_argument("Cauchy256 outer variant");
    }
}

// Caller-owned disjoint scratch, no hot allocations. Input/output may be the
// same N-record buffer, otherwise output contains exactly requested records.
template<bool DenseInner,bool FftOuter,bool ParityOnly=false,bool FftInner=false>
inline void c256TransposeVariant(const Cauchy256Plan& p,const Block* input,Block* output,
                                 Block* planes,Block* temp,Subfield16Kernel kernel) {
    c256InnerDispatchVariant<DenseInner,FftInner>(p,input,planes,temp,kernel);
    c256OuterDispatchVariant<FftOuter,ParityOnly>(p,planes,output,kernel);
}

inline void c256Transpose(const Cauchy256Plan& p,const Block* input,Block* output,
                          Block* planes,Block* temp,Subfield16Kernel kernel,
                          Cauchy256Variant variant=Cauchy256Variant::Basis) {
    // One variant dispatch per encode, never per packet or outer row.
    switch(variant) {
    case Cauchy256Variant::Basis:c256TransposeVariant<false,false>(p,input,output,planes,temp,kernel);break;
    case Cauchy256Variant::Dense:c256TransposeVariant<true,false>(p,input,output,planes,temp,kernel);break;
    case Cauchy256Variant::Fft:c256TransposeVariant<false,true>(p,input,output,planes,temp,kernel);break;
    case Cauchy256Variant::DenseFft:c256TransposeVariant<true,true>(p,input,output,planes,temp,kernel);break;
    case Cauchy256Variant::DenseFftFused:c256TransposeVariant<true,true,true>(p,input,output,planes,temp,kernel);break;
    case Cauchy256Variant::FftInnerFused:c256TransposeVariant<false,true,true,true>(p,input,output,planes,temp,kernel);break;
    default:throw std::invalid_argument("Cauchy256 encoder variant");
    }
}
}
