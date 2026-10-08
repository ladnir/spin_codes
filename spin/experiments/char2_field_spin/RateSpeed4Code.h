#pragma once
#include "RateSpeed4Inner.h"
#include "RateClimb8Outer.h"

namespace spin::research::char2 {
// New construction, not an execution schedule for the width8 checkpoint.
// GRS[64,D], four-input/four-state full-Cauchy inner, first3 outputs retained.
// Four chronological packets occupy12 coordinate-major retained records.
template<unsigned D=40>
struct RateSpeed4Plan {
    static_assert(D==39||D==40);
    static constexpr unsigned dimension=D,outer=64,width=4,state=4;
    std::size_t requested,k,groups,inputN,outputN,planeStride;
    std::vector<std::uint32_t> route,inverseRoute,foldedInverseRoute;
    std::vector<std::uint16_t> foldedInverseLocalRoute;
    std::vector<std::uint8_t> inputDiagonal,packedInputDiagonal,stateDiagonal;
    explicit RateSpeed4Plan(std::size_t message,std::uint64_t seed):requested(message){
        constexpr std::size_t unit=D*64,maxRows=(std::numeric_limits<std::uint32_t>::max()/64/64)*64;
        if(!message||message>maxRows*D)throw std::invalid_argument("RateSpeed4 uint32 route geometry");
        groups=((message+unit-1)/unit)*64;k=D*groups;inputN=64*groups;outputN=48*groups;
        planeStride=((groups+255)/256)*256;
        std::mt19937_64 rr(seed^0xdf4a1b3285ec7096ULL),ir(seed^0x2e7d18c4f0b3a965ULL),sr(seed^0x918b745cf62d03eaULL);
        route.resize(inputN);inverseRoute.resize(inputN);foldedInverseRoute.resize(inputN);
        inputDiagonal.resize(inputN);packedInputDiagonal.resize(inputN);stateDiagonal.resize(inputN);
        for(unsigned region=0;region!=64;++region){
            auto first=route.begin()+region*groups;
            std::iota(first,first+groups,static_cast<std::uint32_t>(region*groups));
            for(std::size_t j=groups-1;j;--j){std::uniform_int_distribution<std::size_t> draw(0,j);std::swap(first[j],first[draw(rr)]);}
        }
        for(auto& d:inputDiagonal)d=sfSampleNonzero(ir);
        for(auto& d:stateDiagonal)d=sfSampleNonzero(sr);
        for(std::size_t i=0;i!=inputN;++i){
            inverseRoute[route[i]]=static_cast<std::uint32_t>(i);
            const auto packed=16*(i/16)+4*(i%4)+(i%16)/4;
            foldedInverseRoute[route[i]]=static_cast<std::uint32_t>(packed);packedInputDiagonal[packed]=inputDiagonal[i];
        }
        if(groups<=65536){
            foldedInverseLocalRoute.resize(inputN);
            for(std::size_t i=0;i!=inputN;++i)foldedInverseLocalRoute[i]=static_cast<std::uint16_t>(foldedInverseRoute[i]-(i/groups)*groups);
        }
    }
    std::size_t planeRecords()const noexcept{return 64*planeStride;}
    static constexpr std::size_t sourceIndex(std::size_t row,unsigned j)noexcept{return (row/4)*(4*D)+4*j+row%4;}
    static constexpr std::size_t outputIndex(std::size_t step,unsigned j)noexcept{return 12*(step/4)+4*j+step%4;}
};

template<bool Affine,unsigned D,bool Stream=true,bool Packed=false,bool Scatter=false,bool Eight=false>
inline void rs4RegionInner(const RateSpeed4Plan<D>& p,const Block* input,Block* planes,Block* temp){
    auto state=_mm512_setzero_si512();alignas(64) Block packet[32];
    for(unsigned r=64;r--;){
        const auto offset=r*p.groups;
        for(auto step=(offset+p.groups)/4;step!=offset/4;){
            step-=Eight?8:4;const auto begin=4*step;
            if constexpr(Eight)rs4InnerTransposeEight<Affine,Affine,Packed>(input+12*(step/4),packet,state,p.stateDiagonal.data()+begin);
            else rs4InnerTransposeFour<Affine,Affine,Packed>(input+12*(step/4),packet,state,p.stateDiagonal.data()+begin);
            __m512i packed;
            if constexpr(Packed)packed=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p.packedInputDiagonal.data()+begin)));
            __m512i packedHigh;
            if constexpr(Packed&&Eight)packedHigh=_mm512_broadcast_i32x4(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p.packedInputDiagonal.data()+begin+16)));
            c8Unroll(std::make_index_sequence<Eight?8:4>{},[&](auto chunk){
                constexpr unsigned shift=4*unsigned(chunk);
                const auto diagonal=Packed?sf16ExpandRow<unsigned(chunk)%4>(unsigned(chunk)<4?packed:packedHigh):sfBroadcastFour(p.packedInputDiagonal.data()+begin+shift);
                const auto value=_mm512_gf2p8mul_epi8(_mm512_load_si512(packet+shift),diagonal);
                if constexpr(Scatter){
                    const auto first=begin+16*(unsigned(chunk)/4)+unsigned(chunk)%4;
                    auto* destination=planes+r*p.planeStride;
                    destination[p.route[first]-offset]=block(_mm512_castsi512_si128(value));
                    destination[p.route[first+4]-offset]=block(_mm512_extracti32x4_epi32(value,1));
                    destination[p.route[first+8]-offset]=block(_mm512_extracti32x4_epi32(value,2));
                    destination[p.route[first+12]-offset]=block(_mm512_extracti32x4_epi32(value,3));
                }else _mm512_storeu_si512(temp+begin-offset+shift,value);
            });
        }
        auto* region=planes+r*p.planeStride;
        if constexpr(!Scatter){if(!p.foldedInverseLocalRoute.empty()){
            for(std::size_t j=0;j!=p.groups;j+=4)sf16Store<Stream>(region+j,c8GatherLocal4(temp,p.foldedInverseLocalRoute.data()+offset+j));
        }else{
            for(std::size_t j=0;j!=p.groups;j+=4)sf16Store<Stream>(region+j,sf16GatherRegion4(temp,p.foldedInverseRoute.data()+offset+j,static_cast<std::uint32_t>(offset)));
        }}
    }
    if constexpr(Stream)_mm_sfence();
}

// The first64-point subtree of the frozen additive FFT computes the new
// outer exactly. The last stage32_0_1 is fused with the requested stores.
inline void rs4OuterPrefix(const Block* source,std::size_t strideBytes,__m512i* values){
    using namespace c128_outer_coefficient_detail;
    const auto* bytes=reinterpret_cast<const unsigned char*>(source);
    leaf16_0_strided_transpose(bytes,strideBytes,values);
    leaf16_16_strided_transpose(bytes+16*strideBytes,strideBytes,values+16);
    stage16_0_1(values);
    leaf16_32_strided_transpose(bytes+32*strideBytes,strideBytes,values+32);
    leaf16_48_strided_transpose(bytes+48*strideBytes,strideBytes,values+48);
    stage16_83_1(values+32);
}
template<unsigned Index>
inline __m512i rs4OuterCoefficient(const __m512i* values){
    const auto value=_mm512_load_si512(values+Index);
    if constexpr(Index<32)return _mm512_xor_si512(value,_mm512_load_si512(values+32+Index));
    else return value;
}
#if defined(_MSC_VER)
#define RS4_NOINLINE __declspec(noinline)
#else
#define RS4_NOINLINE __attribute__((noinline))
#endif
template<unsigned D,bool Stream=true>
RS4_NOINLINE void rs4Outer(const RateSpeed4Plan<D>& p,const Block* planes,Block* output){
    alignas(64) __m512i values[64];
    for(std::size_t row=0;row*D<p.requested;row+=4){
        rs4OuterPrefix(planes+row,p.planeStride*sizeof(Block),values);
        const auto base=row*D;
        if(base+4*D<=p.requested){
            c8Unroll(std::make_index_sequence<D>{},[&](auto j){sf16Store<Stream>(output+base+4*j,rs4OuterCoefficient<static_cast<unsigned>(decltype(j)::value)>(values));});
        }else{
            for(unsigned j=0;j!=D&&base+4*j<p.requested;++j){
                auto value=_mm512_load_si512(values+j);if(j<32)value=_mm512_xor_si512(value,_mm512_load_si512(values+32+j));
                const auto left=p.requested-base-4*j;
                if(left>=4)sf16Store<Stream>(output+base+4*j,value);
                else{alignas(64) Block tail[4];_mm512_store_si512(tail,value);std::memcpy(output+base+4*j,tail,left*sizeof(Block));}
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}
#undef RS4_NOINLINE
template<bool Affine=false,unsigned D,bool Packed=false,bool Scatter=false,bool Eight=false>
inline void rs4Transpose(const RateSpeed4Plan<D>& p,const Block* input,Block* output,Block* planes,Block* temp){
    if((reinterpret_cast<std::uintptr_t>(planes)&63)==0)rs4RegionInner<Affine,D,true,Packed,Scatter,Eight>(p,input,planes,temp);
    else rs4RegionInner<Affine,D,false,Packed,Scatter,Eight>(p,input,planes,temp);
    if((reinterpret_cast<std::uintptr_t>(output)&63)==0)rs4Outer<D,true>(p,planes,output);
    else rs4Outer<D,false>(p,planes,output);
}
// Selected K20 research schedule. Cached regional stores trade a slightly
// dearer route for a faster following outer pass. Eight-packet scheduling
// and batched compact diagonals do not change the sampled linear map.
template<unsigned D>
inline void rs4TransposeCached(const RateSpeed4Plan<D>& p,const Block* input,Block* output,Block* planes,Block* temp){
    rs4RegionInner<true,D,false,true,false,true>(p,input,planes,temp);
    if((reinterpret_cast<std::uintptr_t>(output)&63)==0)rs4Outer<D,true>(p,planes,output);
    else rs4Outer<D,false>(p,planes,output);
}
} // namespace spin::research::char2
