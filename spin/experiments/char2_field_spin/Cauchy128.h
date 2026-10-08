#pragma once
#include "Subfield16.h"
#include "Cauchy256InnerDense.h"
#include "Cauchy128InnerLaneFft.h"
#include "Cauchy128OuterFft.h"
#include <utility>

namespace spin::research::char2 {

enum class Cauchy128Variant : std::uint8_t {
    Baseline, Compact, BatchDiagonal, CompactBatch, LaneFft, CompactLaneFft
};

// GRS[128,87] with the existing width16/state8 inner. Source order is
// EXPLICITLY four-row interleaved: (row/4)*348 + 4*coordinate + row%4.
// This fixed source permutation removes the final record transpose. It
// preserves the full code and every shortened-code distance lower bound.
struct Cauchy128Plan {
    static constexpr unsigned dimension=87,outer=128,width=16,state=8;
    std::size_t requested,k,n,groups;
    std::vector<std::uint32_t> route,inverseRoute;
    // Research controls retain the global uint32 table in every variant.
    // At most 65536 rows means every regional-local index fits uint16.
    std::vector<std::uint16_t> inverseLocalRoute;
    std::vector<std::uint8_t> inputDiagonal,stateDiagonal;
    explicit Cauchy128Plan(std::size_t message,std::uint64_t seed):requested(message){
        constexpr std::size_t unit=dimension*64;
        constexpr std::size_t maxRows=(std::numeric_limits<std::uint32_t>::max()/outer/64)*64;
        if(!message||message>maxRows*dimension)
            throw std::invalid_argument("Cauchy128 message size outside uint32 route geometry");
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
        return (row/4)*348+4*coordinate+row%4;
    }
    std::size_t setupBytes()const noexcept{
        // Both route representations are included, even for Baseline.
        return 4*(route.size()+inverseRoute.size())+2*inverseLocalRoute.size()+
               inputDiagonal.size()+stateDiagonal.size();
    }
};

template<unsigned Batch,bool Scatter>
inline void c128StoreInner(const Cauchy128Plan& p,Block* dest,std::size_t begin,
                           std::size_t base,__m512i value){
    constexpr unsigned shift=4*Batch;
    const auto i=begin+shift;
    const auto v=_mm512_gf2p8mul_epi8(value,sfBroadcastFour(p.inputDiagonal.data()+i));
    if constexpr(Scatter){
        dest[p.route[i]]=block(_mm512_castsi512_si128(v));
        dest[p.route[i+1]]=block(_mm512_extracti32x4_epi32(v,1));
        dest[p.route[i+2]]=block(_mm512_extracti32x4_epi32(v,2));
        dest[p.route[i+3]]=block(_mm512_extracti32x4_epi32(v,3));
    }else _mm512_storeu_si512(dest+begin-base+shift,v);
}
template<unsigned Batch,bool Scatter>
inline void c128StoreInnerBatch(const Cauchy128Plan& p,Block* dest,std::size_t begin,
                                std::size_t base,__m512i value,__m512i packedDiagonal){
    constexpr unsigned shift=4*Batch;
    const auto i=begin+shift;
    const auto v=_mm512_gf2p8mul_epi8(value,sf16ExpandRow<Batch>(packedDiagonal));
    if constexpr(Scatter){
        dest[p.route[i]]=block(_mm512_castsi512_si128(v));
        dest[p.route[i+1]]=block(_mm512_extracti32x4_epi32(v,1));
        dest[p.route[i+2]]=block(_mm512_extracti32x4_epi32(v,2));
        dest[p.route[i+3]]=block(_mm512_extracti32x4_epi32(v,3));
    }else _mm512_storeu_si512(dest+begin-base+shift,v);
}
template<unsigned Phase,bool Scatter,bool BatchDiagonal=false,bool LaneFft=false>
inline void c128Packet(const Cauchy128Plan& p,const Block* input,Block* dest,
                        std::size_t step,std::size_t base,__m512i* state){
    const auto begin=16*step;
    const __m512i in[4]={_mm512_loadu_si512(input+begin),_mm512_loadu_si512(input+begin+4),
                        _mm512_loadu_si512(input+begin+8),_mm512_loadu_si512(input+begin+12)};
    __m512i out[4];
    if constexpr(LaneFft){
        // State stays dword-transposed across every packet and region.
        // Form its diagonal in that layout directly; in/out remain record-major.
        const __m512i diagonal[2]={c128LaneDiagonalFour(p.stateDiagonal.data()+8*step),
                                  c128LaneDiagonalFour(p.stateDiagonal.data()+8*step+4)};
        c128InnerTransposeLaneFftPacked<Phase>(in,out,state,diagonal);
    }else{
        const __m512i diagonal[2]={sfBroadcastFour(p.stateDiagonal.data()+8*step),
                                  sfBroadcastFour(p.stateDiagonal.data()+8*step+4)};
        c256InnerTransposeDensePacked<Phase>(in,out,state,diagonal);
    }
    if constexpr(BatchDiagonal){
        // One compact load/broadcast; selectors expand its four coefficient rows.
        // Load after the dense maps to avoid extending its register lifetime.
        const auto packed=_mm512_broadcast_i32x4(_mm_loadu_si128(
            reinterpret_cast<const __m128i*>(p.inputDiagonal.data()+begin)));
        c128StoreInnerBatch<0,Scatter>(p,dest,begin,base,out[0],packed);
        c128StoreInnerBatch<1,Scatter>(p,dest,begin,base,out[1],packed);
        c128StoreInnerBatch<2,Scatter>(p,dest,begin,base,out[2],packed);
        c128StoreInnerBatch<3,Scatter>(p,dest,begin,base,out[3],packed);
    }else{
        c128StoreInner<0,Scatter>(p,dest,begin,base,out[0]);
        c128StoreInner<1,Scatter>(p,dest,begin,base,out[1]);
        c128StoreInner<2,Scatter>(p,dest,begin,base,out[2]);
        c128StoreInner<3,Scatter>(p,dest,begin,base,out[3]);
    }
}
template<bool Scatter,bool BatchDiagonal=false,bool LaneFft=false>
inline void c128InnerRange(const Cauchy128Plan& p,const Block* in,Block* dest,
                           std::size_t first,std::size_t past,std::size_t base,__m512i* state){
    for(std::size_t step=past;step!=first;){
        c128Packet<3,Scatter,BatchDiagonal,LaneFft>(p,in,dest,--step,base,state);
        c128Packet<2,Scatter,BatchDiagonal,LaneFft>(p,in,dest,--step,base,state);
        c128Packet<1,Scatter,BatchDiagonal,LaneFft>(p,in,dest,--step,base,state);
        c128Packet<0,Scatter,BatchDiagonal,LaneFft>(p,in,dest,--step,base,state);
    }
}
inline __m512i c128GatherLocal4(const Block* input,const std::uint16_t* index){
    auto value=_mm512_castsi128_si512(reg(input[index[0]]));
    value=_mm512_inserti32x4(value,reg(input[index[1]]),1);
    value=_mm512_inserti32x4(value,reg(input[index[2]]),2);
    return _mm512_inserti32x4(value,reg(input[index[3]]),3);
}
template<bool Stream,bool Compact=false,bool BatchDiagonal=false,bool LaneFft=false>
inline void c128RegionInner(const Cauchy128Plan& p,const Block* in,Block* planes,Block* temp){
    __m512i state[2]={_mm512_setzero_si512(),_mm512_setzero_si512()};
    for(unsigned r=128;r--;){
        const auto offset=r*p.groups;
        c128InnerRange<false,BatchDiagonal,LaneFft>(p,in,temp,offset/16,(offset+p.groups)/16,offset,state);
        if constexpr(Compact){
            for(std::size_t j=0;j!=p.groups;j+=4)
                sf16Store<Stream>(planes+offset+j,c128GatherLocal4(temp,p.inverseLocalRoute.data()+offset+j));
        }else{
            for(std::size_t j=0;j!=p.groups;j+=4)
                sf16Store<Stream>(planes+offset+j,sf16GatherRegion4(temp,p.inverseRoute.data()+offset+j,
                    static_cast<std::uint32_t>(offset)));
        }
    }
    if constexpr(Stream)_mm_sfence();
}
template<bool Compact,bool BatchDiagonal,bool LaneFft=false>
inline void c128InnerDispatchVariant(const Cauchy128Plan& p,const Block* in,Block* planes,
                                     Block* temp,Subfield16Kernel kernel){
    if(kernel==Subfield16Kernel::Scatter){
        __m512i state[2]={_mm512_setzero_si512(),_mm512_setzero_si512()};
        // Compact deliberately has no effect on the scatter schedule.
        c128InnerRange<true,BatchDiagonal,LaneFft>(p,in,planes,0,p.n/16,0,state);
    }else{
        if(!temp)throw std::invalid_argument("regional schedule requires temporary storage");
        if constexpr(Compact){
            // Checked once per encode, not once per packet or gather.
            if(p.inverseLocalRoute.empty()){
                c128InnerDispatchVariant<false,BatchDiagonal,LaneFft>(p,in,planes,temp,kernel);
                return;
            }
        }
        if((reinterpret_cast<std::uintptr_t>(planes)&63)==0)
            c128RegionInner<true,Compact,BatchDiagonal,LaneFft>(p,in,planes,temp);
        else c128RegionInner<false,Compact,BatchDiagonal,LaneFft>(p,in,planes,temp);
    }
}
inline void c128InnerDispatch(const Cauchy128Plan& p,const Block* in,Block* planes,
                              Block* temp,Subfield16Kernel kernel,
                              Cauchy128Variant variant=Cauchy128Variant::Baseline){
    // One runtime variant dispatch per encode; every inner loop is specialized.
    switch(variant){
    case Cauchy128Variant::Baseline:
        c128InnerDispatchVariant<false,false>(p,in,planes,temp,kernel);break;
    case Cauchy128Variant::Compact:
        c128InnerDispatchVariant<true,false>(p,in,planes,temp,kernel);break;
    case Cauchy128Variant::BatchDiagonal:
        c128InnerDispatchVariant<false,true>(p,in,planes,temp,kernel);break;
    case Cauchy128Variant::CompactBatch:
        c128InnerDispatchVariant<true,true>(p,in,planes,temp,kernel);break;
    case Cauchy128Variant::LaneFft:
        c128InnerDispatchVariant<false,false,true>(p,in,planes,temp,kernel);break;
    case Cauchy128Variant::CompactLaneFft:
        c128InnerDispatchVariant<true,false,true>(p,in,planes,temp,kernel);break;
    default:throw std::invalid_argument("Cauchy128 encoder variant");
    }
}
template<class F,std::size_t... I>
inline void c128Unroll(std::index_sequence<I...>,F&& f){(f(std::integral_constant<std::size_t,I>{}),...);}
template<bool Stream>
inline void c128OuterFinish(const Cauchy128Plan& p,const Block* planes,Block* out){
    alignas(64) __m512i parity[41],result[87];
    for(std::size_t row=0;row*87<p.requested;row+=4){
        c128Unroll(std::make_index_sequence<41>{},[&](auto j){
            parity[j]=_mm512_loadu_si512(planes+(87+std::size_t(j))*p.groups+row);
        });
        c128OuterParityTransposeFft(parity,result);
        const auto base=row*87;
        if(base+348<=p.requested){
            c128Unroll(std::make_index_sequence<87>{},[&](auto j){
                sf16Store<Stream>(out+base+4*j,_mm512_xor_si512(result[j],
                    _mm512_loadu_si512(planes+std::size_t(j)*p.groups+row)));
            });
        }else{
            // Only the final tile is partial. Preserve an in-place N-buffer tail.
            for(unsigned j=0;j!=87&&base+4*j<p.requested;++j){
                const auto v=_mm512_xor_si512(result[j],_mm512_loadu_si512(planes+std::size_t(j)*p.groups+row));
                const auto left=p.requested-base-4*j;
                if(left>=4)sf16Store<Stream>(out+base+4*j,v);
                else{alignas(64) Block tail[4];_mm512_store_si512(tail,v);
                     std::memcpy(out+base+4*j,tail,left*sizeof(Block));}
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}
inline void c128OuterDispatch(const Cauchy128Plan& p,const Block* planes,Block* out,
                              Subfield16Kernel kernel){
    if(kernel==Subfield16Kernel::RegionAllNonTemporal&&(reinterpret_cast<std::uintptr_t>(out)&63)==0)
        c128OuterFinish<true>(p,planes,out);
    else c128OuterFinish<false>(p,planes,out);
}
// Caller-owned scratch; input/output may alias exactly. Scratch must be disjoint.
inline void c128Transpose(const Cauchy128Plan& p,const Block* in,Block* out,
                          Block* planes,Block* temp,Subfield16Kernel kernel,
                          Cauchy128Variant variant=Cauchy128Variant::Baseline){
    c128InnerDispatch(p,in,planes,temp,kernel,variant);c128OuterDispatch(p,planes,out,kernel);
}
}
