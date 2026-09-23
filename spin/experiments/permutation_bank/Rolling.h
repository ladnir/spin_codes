#pragma once
#include "Flows.h"
#include "Masks.h"

namespace spin::experimental::bank {
// Same family, different execution schedule. Route addresses are buffered for
// only one region. Split=true additionally separates inner work from scatter.
template<unsigned Count,bool Split=false,bool Vector=false,unsigned Prefetch=0,bool EpochRoute=false,class Route=Routing<Count,Vector>> class RollingCode {
    Route route_;
    std::vector<std::uint32_t> masks_;
    template<bool Four> void run(kernel::block* input,kernel::block* scratch)const {
        alignas(64) std::array<unsigned,2048> addresses;
        if constexpr(EpochRoute) {
            // Same recurrence as innerReverse<Map128S19>, with region routing
            // prepared outside the 128-point unrolled emitter.
            using Map=kernel::Map128S19;
            alignas(32) __m128i state[Map::S]{},syndrome[Map::S],values[Map::T];
            constexpr unsigned epochs=CodeSize/Map::T;
#if SPIN_COMPOSED_BYTE_ROUTE
            constexpr bool byteRoute=requires{route_.outerRegionBytes(0,addresses.data(),Four);};
#else
            constexpr bool byteRoute=false;
#endif
            for(unsigned epoch=epochs;epoch-->0;) {
                const unsigned base=epoch*Map::T;
#if SPIN_COMPOSED_ROUTE_CHUNK
                constexpr unsigned requestedChunk=SPIN_COMPOSED_ROUTE_CHUNK;
#elif SPIN_COMPOSED_EPOCH_ROUTING
                constexpr unsigned requestedChunk=128;
#else
                constexpr unsigned requestedChunk=2048;
#endif
                static_assert(requestedChunk>=128 && requestedChunk<=2048 && std::has_single_bit(requestedChunk));
                constexpr bool shortRoute=requestedChunk<2048 && requires{route_.template outerChunk<requestedChunk>(0,0,addresses.data(),Four);};
                constexpr unsigned routeChunk=shortRoute?requestedChunk:2048;
                static_assert(!byteRoute || !shortRoute);
                if((epoch&15)==15) {
                    if constexpr(byteRoute)route_.outerRegionBytes(base>>11,addresses.data(),Four);
                    else if constexpr(!shortRoute)route_.outerRegion(base>>11,addresses.data(),Four);
                    if constexpr(requires{route_.prefetchNext(0);})if(base>=2048)route_.prefetchNext((base>>11)-1);
                }
                if constexpr(shortRoute)if((base&(routeChunk-1))==routeChunk-128)
                    route_.template outerChunk<routeChunk>(base>>11,(base&2047)&~(routeChunk-1),addresses.data(),Four);
                if constexpr(requires{route_.prefetchRows();})if((epoch&15)==0)route_.prefetchRows();
                const auto* epochAddresses=addresses.data()+(base&(routeChunk-1));
                auto emit=[&](std::size_t i,kernel::block v){
#if SPIN_COMPOSED_OUTPUT_PREFETCH
                    if((i&127)>=SPIN_COMPOSED_OUTPUT_PREFETCH)
                        __builtin_prefetch(reinterpret_cast<char*>(scratch)+(std::size_t(epochAddresses[(i&127)-SPIN_COMPOSED_OUTPUT_PREFETCH])<<(byteRoute?0:4)),1,3);
#endif
                    if constexpr(byteRoute)*reinterpret_cast<kernel::block*>(reinterpret_cast<char*>(scratch)+epochAddresses[i&127])=v;
                    else scratch[epochAddresses[i&127]]=v;
                };
#if SPIN_COMPOSED_OUTPUT_PREFETCH
                for(unsigned p=128;p-->128-SPIN_COMPOSED_OUTPUT_PREFETCH;)
                    __builtin_prefetch(reinterpret_cast<char*>(scratch)+(std::size_t(epochAddresses[p])<<(byteRoute?0:4)),1,3);
#endif
                if(epoch+1==epochs) {
                    for(unsigned p=Map::T;p-->0;){values[p]=input[base+p].mData;emit(base+p,input[base+p]);}
                } else kernel::imtReversePoints<Map>(input+base,values,state,base,emit,std::make_index_sequence<Map::T>{});
                if(epoch==0)break;
                kernel::zeta<Map::T>(values);Map::finish(values,syndrome);
                if(epoch+1==epochs)std::memcpy(state,syndrome,sizeof(state));
                else {
                    auto u=masks_[2*epoch];auto dot=_mm_setzero_si128();
                    while(u){const auto j=std::countr_zero(u);u&=u-1;dot=_mm_xor_si128(dot,state[j]);}
                    auto v=masks_[2*epoch+1];
                    while(v){const auto j=std::countr_zero(v);v&=v-1;state[j]=_mm_xor_si128(state[j],dot);}
                    for(unsigned j=0;j<Map::S;++j)state[j]=_mm_xor_si128(state[j],syndrome[j]);
                }
            }
        }else if constexpr(Split) {
            // raw values are saved by innerReverse before emitting each output,
            // so overwriting the current coordinate is safe.
            kernel::innerReverse<kernel::Map128S19>(input,CodeSize,masks_.data(),[&](std::size_t i,kernel::block v){input[i]=v;});
            for(unsigned region=0;region<256;++region) {
                route_.outerRegion(region,addresses.data(),Four);
                for(unsigned i=0;i<2048;++i) {
#if defined(__GNUC__) || defined(__clang__)
                    if constexpr(Prefetch)if(i+Prefetch<2048)__builtin_prefetch(scratch+addresses[i+Prefetch],1,3);
#endif
                    scratch[addresses[i]]=input[region*2048+i];
                }
            }
        }else {
            kernel::innerReverse<kernel::Map128S19>(input,CodeSize,masks_.data(),[&](std::size_t i,kernel::block v) {
                if((i&2047)==2047)route_.outerRegion(unsigned(i)>>11,addresses.data(),Four);
#if defined(__GNUC__) || defined(__clang__)
                if constexpr(Prefetch)if((i&2047)>=Prefetch)__builtin_prefetch(scratch+addresses[(i&2047)-Prefetch],1,3);
#endif
                scratch[addresses[i&2047]]=v;
            });
        }
        if constexpr(Four)for(unsigned i=0;i<CodeSize;i+=1024)kernel::bchTranspose4(scratch+i,input+i/2);
        else for(unsigned i=0;i<CodeSize;i+=512)kernel::bchTranspose2(scratch+i,scratch+i+256,input+i/2,input+i/2+128);
    }
public:
    RollingCode(const Tables<Count>& tables,std::uint64_t seed,std::uint64_t maskSeed):route_(tables,seed),masks_(fastMasks(maskSeed)) {}
    RollingCode(Route route,std::uint64_t maskSeed):route_(std::move(route)),masks_(fastMasks(maskSeed)) {}
    void encode(kernel::block* input,std::vector<kernel::block>& scratch)const {
        if(scratch.size()!=CodeSize)throw std::invalid_argument("rolling scratch");
#if SPIN_BCH_AVX512
        if(kernel::bchAvx512Available()){run<true>(input,scratch.data());return;}
#endif
        run<false>(input,scratch.data());
    }
};
}
