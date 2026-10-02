#pragma once
#include "Spin.h"
#include <array>
#include <cstring>
#include <utility>

namespace bare_spin::schedule {
template<bool Packed> OC_FORCEINLINE u32 index(const u8* p24,const u32* p32,std::size_t i) {
    if constexpr(Packed) {u32 v;std::memcpy(&v,p24+3*i,4);return v&0xffffffU;}
    else return p32[i];
}

// Decode addresses and issue independent sequential loads before scatter stores.
// U is compile-time fixed. All arrays are bounded scratch in an ordinary function.
template<bool Packed,unsigned U,unsigned PF,unsigned Write,std::size_t... J>
OC_FORCEINLINE void group(const block* values,block* tile,const u8* p24,const u32* p32,
                          std::size_t pos,bool prefetch,std::index_sequence<J...>) {
    if constexpr(PF) if(prefetch)
        (__builtin_prefetch(tile+index<Packed>(p24,p32,pos+PF+J),Write,3),...);
    const std::array<u32,U> offsets{index<Packed>(p24,p32,pos+J)...};
    const std::array<block,U> inputs{values[pos+J]...};
    ((tile[offsets[J]]=inputs[J]),...);
}
template<bool Packed,unsigned U,unsigned PF,unsigned Write>
OC_FORCEINLINE void scatter(const block* values,block* tile,const u8* p24,const u32* p32,
                           std::size_t base,std::size_t count) {
    std::size_t j=0;
    if constexpr(PF) {
        for(;j+U+PF<=count;j+=U)
            group<Packed,U,PF,Write>(values,tile,p24,p32,base+j,true,std::make_index_sequence<U>{});
    }
    for(;j+U<=count;j+=U)
        group<Packed,U,PF,Write>(values,tile,p24,p32,base+j,false,std::make_index_sequence<U>{});
    for(;j<count;++j) tile[index<Packed>(p24,p32,base+j)]=values[base+j];
}

template<bool Packed,bool Decode,unsigned Lead> struct Emitter {
    block* values;
    const u8* p24;
    const u32* p32;
    std::array<u32,Decode?128:1> decoded;
    OC_FORCEINLINE void prepare(std::size_t base) {
        // Reverse traversal: base-Lead is a future epoch. No padding or
        // out-of-domain pointer calculation is used at the beginning of input.
        if constexpr(Lead) if(base>=Lead)
            for(unsigned p=0;p<128;++p)
                __builtin_prefetch(values+index<Packed>(p24,p32,base-Lead+p),1,3);
        if constexpr(Decode)
            for(unsigned p=0;p<128;++p) decoded[p]=index<Packed>(p24,p32,base+p);
    }
    OC_FORCEINLINE void operator()(std::size_t i,block v) {
        if constexpr(Decode) values[decoded[i&127]]=v;
        else values[index<Packed>(p24,p32,i)]=v;
    }
};
}
