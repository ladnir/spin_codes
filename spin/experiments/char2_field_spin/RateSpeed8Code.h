#pragma once
#include "RateClimb8Code.h"
#include "RateClimb8Outer.h"
#include "RateSpeed8Inner.h"
#include "RateSpeed8Outer.h"

namespace spin::research::char2 {

template<bool Maps,bool Scale>
inline void s8AffineRegionInner(const RateClimb8Plan<R8Periodic84>& p,const Block* input,Block* planes,Block* temp){
    auto h0=_mm512_setzero_si512(),h1=_mm512_setzero_si512();alignas(64) Block packet[32];
    for(unsigned r=128;r--;){
        const auto offset=r*p.groups;
        for(auto step=(offset+p.groups)/8;step!=offset/8;){
            step-=4;const auto begin=8*step;
            rs8InnerTransposeFour<R8Periodic84,Maps,Scale,true>(input+25*(step/4),packet,h0,h1,p.stateDiagonal.data()+begin);
            c8Unroll(std::make_index_sequence<8>{},[&](auto chunk){
                constexpr unsigned shift=4*unsigned(chunk);
                _mm512_storeu_si512(temp+begin-offset+shift,_mm512_gf2p8mul_epi8(_mm512_load_si512(packet+shift),sfBroadcastFour(p.packedInputDiagonal.data()+begin+shift)));
            });
        }
        for(std::size_t j=0;j!=p.groups;j+=4)sf16Store<true>(planes+r*p.planeStride+j,c8GatherLocal4(temp,p.foldedInverseLocalRoute.data()+offset+j));
    }
    _mm_sfence();
}

// Exact-map research controls. Tile=0 is the frozen column layout. Otherwise
// outer rows are tiled, with all 128 outer coordinates inside each row tile.
// All strides count 128-bit records; only live entries are read.
template<unsigned Tile>
inline std::size_t s8PlaneRecords(const RateClimb8Plan<R8Periodic84>& p){
    if constexpr(Tile==0)return p.planeRecords();
    else return ((p.groups+Tile-1)/Tile)*Tile*128;
}

template<unsigned Tile,unsigned Prefetch=0,bool Stream=true>
inline void s8RegionInner(const RateClimb8Plan<R8Periodic84>& p,const Block* input,Block* planes,Block* temp){
    static_assert(Tile==0||Tile%4==0);
    auto h0=_mm512_setzero_si512(),h1=_mm512_setzero_si512();
    for(unsigned r=128;r--;){
        const auto offset=r*p.groups;
        r8InnerRange<false,true>(p,input,temp,offset/8,(offset+p.groups)/8,offset,h0,h1);
        const auto* indices=p.foldedInverseLocalRoute.data()+offset;
        constexpr unsigned tileSize=Tile==0?256:Tile;
        for(std::size_t base=0;base<p.groups;base+=tileSize){
            auto* destination=Tile==0?planes+r*p.planeStride+base:planes+(base/Tile)*128*Tile+r*Tile;
            const auto count=std::min<std::size_t>(tileSize,p.groups-base);
            for(std::size_t j=0;j<count;j+=4){
                if constexpr(Prefetch!=0){
                    if(base+j+Prefetch+3<p.groups){
                        const auto* ahead=indices+base+j+Prefetch;
                        _mm_prefetch(reinterpret_cast<const char*>(temp+ahead[0]),_MM_HINT_T0);
                        _mm_prefetch(reinterpret_cast<const char*>(temp+ahead[1]),_MM_HINT_T0);
                        _mm_prefetch(reinterpret_cast<const char*>(temp+ahead[2]),_MM_HINT_T0);
                        _mm_prefetch(reinterpret_cast<const char*>(temp+ahead[3]),_MM_HINT_T0);
                    }
                }
                sf16Store<Stream>(destination+j,c8GatherLocal4(temp,indices+base+j));
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}

#if defined(_MSC_VER)
#define S8_NOINLINE __declspec(noinline)
#else
#define S8_NOINLINE __attribute__((noinline))
#endif
template<unsigned Tile,bool Stream=true>
S8_NOINLINE void s8Outer(const RateClimb8Plan<R8Periodic84>& p,const Block* planes,Block* out){
    constexpr unsigned dimension=R8Periodic84::dimension;
    alignas(64) __m512i values[128];
    for(std::size_t row=0;row*dimension<p.requested;row+=4){
        const auto* source=Tile==0?planes+row:planes+(row/Tile)*128*Tile+row%Tile;
        const auto stride=Tile==0?p.planeStride:Tile;
        r8OuterPrefixFinalStore<true>(source,stride*sizeof(Block),values);
        const auto base=row*dimension;
        if(base+4*dimension<=p.requested){
            c8Unroll(std::make_index_sequence<dimension>{},[&](auto j){
                sf16Store<Stream>(out+base+4*j,
                    r8OuterFinalCoefficient<static_cast<unsigned>(decltype(j)::value)>(values));
            });
        }else{
            for(unsigned j=0;j!=dimension&&base+4*j<p.requested;++j){
                auto value=_mm512_load_si512(values+j);
                if(j<64)value=_mm512_xor_si512(value,_mm512_load_si512(values+64+j));
                const auto left=p.requested-base-4*j;
                if(left>=4)sf16Store<Stream>(out+base+4*j,value);
                else{alignas(64) Block tail[4];_mm512_store_si512(tail,value);std::memcpy(out+base+4*j,tail,left*sizeof(Block));}
            }
        }
    }
    if constexpr(Stream)_mm_sfence();
}
#undef S8_NOINLINE

template<unsigned Tile,unsigned Prefetch=0>
inline void s8Transpose(const RateClimb8Plan<R8Periodic84>& p,const Block* input,Block* output,Block* planes,Block* temp){
    if(p.foldedInverseLocalRoute.empty())throw std::invalid_argument("speed control requires compact row indices");
    const bool planeAligned=(reinterpret_cast<std::uintptr_t>(planes)&63)==0;
    if(planeAligned)s8RegionInner<Tile,Prefetch,true>(p,input,planes,temp);
    else s8RegionInner<Tile,Prefetch,false>(p,input,planes,temp);
    if((reinterpret_cast<std::uintptr_t>(output)&63)==0)s8Outer<Tile,true>(p,planes,output);
    else s8Outer<Tile,false>(p,planes,output);
}
} // namespace spin::research::char2
