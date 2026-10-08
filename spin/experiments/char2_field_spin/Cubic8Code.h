#pragma once
#include "Cubic8Inner.h"
#include "Cauchy128OuterCoefficients.h"
#include <utility>

namespace spin::research::char2 {

enum class Cubic8OuterVariant : std::uint8_t { Gather,FusedLoad };

// Separate GRS[128,83] / width8,state4 profile. Source order is explicitly
// (row/4)*332+4*coefficient+row%4 in normalized additive-polynomial basis.
// The first83 degrees span exactly GRS[128,83]. Shortening fixes the padded
// suffix of these source coordinates to zero. mt19937_64 is a test fixture,
// not a production implementation of the ideal setup sampler.
struct Cubic8Plan {
    static constexpr unsigned dimension=83,outer=128,width=8,state=4;
    std::size_t requested,k,n,groups;
    std::vector<std::uint32_t> route,inverseRoute;
    std::vector<std::uint16_t> inverseLocalRoute;
    std::vector<std::uint8_t> inputDiagonal,stateDiagonal;
    explicit Cubic8Plan(std::size_t message,std::uint64_t seed):requested(message){
        constexpr std::size_t unit=dimension*64;
        constexpr std::size_t maxRows=(std::numeric_limits<std::uint32_t>::max()/outer/64)*64;
        if(!message||message>maxRows*dimension)
            throw std::invalid_argument("Cubic8 message size outside uint32 route geometry");
        groups=((message+unit-1)/unit)*64;k=groups*dimension;n=groups*outer;
        std::mt19937_64 rr(seed^0xdf4a1b3285ec7096ULL),ir(seed^0x2e7d18c4f0b3a965ULL),
                        sr(seed^0x918b745cf62d03eaULL);
        route.resize(n);inverseRoute.resize(n);inputDiagonal.resize(n);stateDiagonal.resize(n/2);
        for(unsigned r=0;r!=outer;++r){
            auto first=route.begin()+r*groups;
            std::iota(first,first+groups,static_cast<std::uint32_t>(r*groups));
            for(std::size_t j=groups-1;j;--j){
                std::uniform_int_distribution<std::size_t> draw(0,j);
                std::swap(first[j],first[draw(rr)]);
            }
        }
        for(std::size_t i=0;i!=n;++i)inverseRoute[route[i]]=static_cast<std::uint32_t>(i);
        if(groups<=65536){
            inverseLocalRoute.resize(n);
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
        return (row/4)*332+4*coordinate+row%4;
    }
    std::size_t setupBytes()const noexcept{
        return 4*(route.size()+inverseRoute.size())+2*inverseLocalRoute.size()+
               inputDiagonal.size()+stateDiagonal.size();
    }
};

template<unsigned Batch,bool Scatter>
inline void c8StoreInner(const Cubic8Plan& p,Block* dest,std::size_t begin,
                          std::size_t base,__m512i value,__m512i packedDiagonal){
    constexpr unsigned shift=4*Batch;const auto i=begin+shift;
    const auto v=_mm512_gf2p8mul_epi8(value,sf16ExpandRow<Batch>(packedDiagonal));
    if constexpr(Scatter){
        dest[p.route[i]]=block(_mm512_castsi512_si128(v));
        dest[p.route[i+1]]=block(_mm512_extracti32x4_epi32(v,1));
        dest[p.route[i+2]]=block(_mm512_extracti32x4_epi32(v,2));
        dest[p.route[i+3]]=block(_mm512_extracti32x4_epi32(v,3));
    }else _mm512_storeu_si512(dest+begin-base+shift,v);
}
template<unsigned Phase,bool Scatter>
inline void c8Packet(const Cubic8Plan& p,const Block* input,Block* dest,
                      std::size_t step,std::size_t base,__m512i& state){
    const auto begin=8*step;
    const __m512i in[2]={_mm512_loadu_si512(input+begin),_mm512_loadu_si512(input+begin+4)};
    __m512i out[2];
    c8InnerTransposePacked<Phase>(in,out,state,sfBroadcastFour(p.stateDiagonal.data()+4*step));
    // Exactly eight bytes, including at the last packet; no setup overread.
    const auto packed=_mm512_broadcast_i32x4(_mm_loadl_epi64(
        reinterpret_cast<const __m128i*>(p.inputDiagonal.data()+begin)));
    c8StoreInner<0,Scatter>(p,dest,begin,base,out[0],packed);
    c8StoreInner<1,Scatter>(p,dest,begin,base,out[1],packed);
}
template<unsigned FirstPhase,unsigned SecondPhase,bool Scatter>
inline void c8PacketPair(const Cubic8Plan& p,const Block* input,Block* dest,
                          std::size_t lowStep,std::size_t base,__m512i& state){
    const auto begin=8*lowStep;
    // Higher-index packet is processed first by the adjoint.
    const __m512i in[4]={_mm512_loadu_si512(input+begin+8),_mm512_loadu_si512(input+begin+12),
                         _mm512_loadu_si512(input+begin),_mm512_loadu_si512(input+begin+4)};
    const __m512i diagonal[2]={sfBroadcastFour(p.stateDiagonal.data()+4*lowStep+4),
                               sfBroadcastFour(p.stateDiagonal.data()+4*lowStep)};
    __m512i out[4];c8InnerTransposePairedPacked<FirstPhase,SecondPhase>(in,out,state,diagonal);
    const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(
        reinterpret_cast<const __m128i*>(p.inputDiagonal.data()+begin)));
    c8StoreInner<2,Scatter>(p,dest,begin,base,out[0],packed);
    c8StoreInner<3,Scatter>(p,dest,begin,base,out[1],packed);
    c8StoreInner<0,Scatter>(p,dest,begin,base,out[2],packed);
    c8StoreInner<1,Scatter>(p,dest,begin,base,out[3],packed);
}
template<bool Scatter,bool Paired=false>
inline void c8InnerRange(const Cubic8Plan& p,const Block* in,Block* dest,
                          std::size_t first,std::size_t past,std::size_t base,__m512i& state){
    for(std::size_t step=past;step!=first;){
      if constexpr(Paired){
        step-=2;c8PacketPair<7,6,Scatter>(p,in,dest,step,base,state);
        step-=2;c8PacketPair<5,4,Scatter>(p,in,dest,step,base,state);
        step-=2;c8PacketPair<3,2,Scatter>(p,in,dest,step,base,state);
        step-=2;c8PacketPair<1,0,Scatter>(p,in,dest,step,base,state);
      }else{
        c8Packet<7,Scatter>(p,in,dest,--step,base,state);
        c8Packet<6,Scatter>(p,in,dest,--step,base,state);
        c8Packet<5,Scatter>(p,in,dest,--step,base,state);
        c8Packet<4,Scatter>(p,in,dest,--step,base,state);
        c8Packet<3,Scatter>(p,in,dest,--step,base,state);
        c8Packet<2,Scatter>(p,in,dest,--step,base,state);
        c8Packet<1,Scatter>(p,in,dest,--step,base,state);
        c8Packet<0,Scatter>(p,in,dest,--step,base,state);
      }
    }
}
inline __m512i c8GatherLocal4(const Block* input,const std::uint16_t* index){
    auto value=_mm512_castsi128_si512(reg(input[index[0]]));
    value=_mm512_inserti32x4(value,reg(input[index[1]]),1);
    value=_mm512_inserti32x4(value,reg(input[index[2]]),2);
    return _mm512_inserti32x4(value,reg(input[index[3]]),3);
}
template<bool Stream,bool Compact,bool Paired=false>
inline void c8RegionInner(const Cubic8Plan& p,const Block* in,Block* planes,Block* temp){
    auto state=_mm512_setzero_si512();
    for(unsigned r=128;r--;){
        const auto offset=r*p.groups;
        c8InnerRange<false,Paired>(p,in,temp,offset/8,(offset+p.groups)/8,offset,state);
        if constexpr(Compact){
            for(std::size_t j=0;j!=p.groups;j+=4)
                sf16Store<Stream>(planes+offset+j,c8GatherLocal4(temp,p.inverseLocalRoute.data()+offset+j));
        }else{
            for(std::size_t j=0;j!=p.groups;j+=4)
                sf16Store<Stream>(planes+offset+j,sf16GatherRegion4(temp,p.inverseRoute.data()+offset+j,
                    static_cast<std::uint32_t>(offset)));
        }
    }
    if constexpr(Stream)_mm_sfence();
}
template<bool Paired>
inline void c8InnerDispatchVariant(const Cubic8Plan& p,const Block* in,Block* planes,
                             Block* temp,Subfield16Kernel kernel){
    if(kernel==Subfield16Kernel::Scatter){
        auto state=_mm512_setzero_si512();c8InnerRange<true,Paired>(p,in,planes,0,p.n/8,0,state);
    }else{
        if(!temp)throw std::invalid_argument("Cubic8 regional schedule requires temporary storage");
        const bool aligned=(reinterpret_cast<std::uintptr_t>(planes)&63)==0;
        if(p.inverseLocalRoute.empty()){
            if(aligned)c8RegionInner<true,false,Paired>(p,in,planes,temp);
            else c8RegionInner<false,false,Paired>(p,in,planes,temp);
        }else{
            if(aligned)c8RegionInner<true,true,Paired>(p,in,planes,temp);
            else c8RegionInner<false,true,Paired>(p,in,planes,temp);
        }
    }
}
inline void c8InnerDispatch(const Cubic8Plan& p,const Block* in,Block* planes,
                             Block* temp,Subfield16Kernel kernel){
    c8InnerDispatchVariant<false>(p,in,planes,temp,kernel);
}
inline void c8InnerDispatchPaired(const Cubic8Plan& p,const Block* in,Block* planes,
                                   Block* temp,Subfield16Kernel kernel){
    c8InnerDispatchVariant<true>(p,in,planes,temp,kernel);
}

template<class F,std::size_t... I>
inline void c8Unroll(std::index_sequence<I...>,F&& f){(f(std::integral_constant<std::size_t,I>{}),...);}
#if defined(_MSC_VER)
#define C8_NOINLINE __declspec(noinline)
#else
#define C8_NOINLINE __attribute__((noinline))
#endif
template<bool Stream,bool Fused>
C8_NOINLINE void c8OuterFinish(const Cubic8Plan& p,const Block* planes,Block* out){
    alignas(64) __m512i values[128];
    for(std::size_t row=0;row*83<p.requested;row+=4){
        if constexpr(Fused)c128OuterCoefficientTransposeStrided(planes+row,p.groups*sizeof(Block),values);
        else{
            c8Unroll(std::make_index_sequence<128>{},[&](auto j){
                values[j]=_mm512_loadu_si512(planes+std::size_t(j)*p.groups+row);
            });
            c128OuterCoefficientTransposeInplace(values);
        }
        const auto base=row*83;
        if(base+332<=p.requested){
            c8Unroll(std::make_index_sequence<83>{},[&](auto j){sf16Store<Stream>(out+base+4*j,values[j]);});
        }else{
            for(unsigned j=0;j!=83&&base+4*j<p.requested;++j){
                const auto left=p.requested-base-4*j;
                if(left>=4)sf16Store<Stream>(out+base+4*j,values[j]);
                else{alignas(64) Block tail[4];_mm512_store_si512(tail,values[j]);
                     std::memcpy(out+base+4*j,tail,left*sizeof(Block));}
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}
#undef C8_NOINLINE
inline void c8OuterDispatch(const Cubic8Plan& p,const Block* planes,Block* out,
                             Subfield16Kernel kernel,Cubic8OuterVariant variant=Cubic8OuterVariant::FusedLoad){
    const bool stream=kernel==Subfield16Kernel::RegionAllNonTemporal&&(reinterpret_cast<std::uintptr_t>(out)&63)==0;
    if(variant==Cubic8OuterVariant::FusedLoad){
        if(stream)c8OuterFinish<true,true>(p,planes,out);else c8OuterFinish<false,true>(p,planes,out);
    }else{
        if(stream)c8OuterFinish<true,false>(p,planes,out);else c8OuterFinish<false,false>(p,planes,out);
    }
}
// Exact in/out alias is supported; caller-owned planes/temp must be disjoint.
// Writes exactly requested records, preserving the full input buffer suffix.
inline void c8Transpose(const Cubic8Plan& p,const Block* in,Block* out,Block* planes,
                         Block* temp,Subfield16Kernel kernel,
                         Cubic8OuterVariant variant=Cubic8OuterVariant::FusedLoad){
    c8InnerDispatch(p,in,planes,temp,kernel);c8OuterDispatch(p,planes,out,kernel,variant);
}
inline void c8TransposePaired(const Cubic8Plan& p,const Block* in,Block* out,Block* planes,
                               Block* temp,Subfield16Kernel kernel,
                               Cubic8OuterVariant variant=Cubic8OuterVariant::FusedLoad){
    c8InnerDispatchPaired(p,in,planes,temp,kernel);c8OuterDispatch(p,planes,out,kernel,variant);
}
} // namespace spin::research::char2
