#pragma once
#include "Inner.h"

namespace bare_spin {
template<u32 Mask> OC_FORCEINLINE __m128i sparseFixedSum(const __m128i* state) {
    static_assert(Mask!=0);
    constexpr unsigned j=std::countr_zero(Mask);
    if constexpr((Mask&(Mask-1))==0) return state[j];
    else return _mm_xor_si128(state[j],sparseFixedSum<Mask&(Mask-1)>(state));
}

template<class Map,std::size_t R,class Emit> OC_FORCEINLINE void asymmetricPoint(
    const block* in,__m128i* raw,const __m128i* state,std::size_t base,Emit& emit) {
    constexpr unsigned p=Map::T-1-R;
    const auto input=in[p].mData;
    raw[p]=input; // A^T acts on U, not on the emitted V=U+B^T r.
    emit(base+p,block(_mm_xor_si128(input,sparseFixedSum<Map::columns[p]>(state))));
}

template<class Map,class Emit,std::size_t... R> OC_FORCEINLINE void asymmetricPoints(
    const block* in,__m128i* raw,const __m128i* state,std::size_t base,Emit& emit,std::index_sequence<R...>) {
    (asymmetricPoint<Map,R>(in,raw,state,base,emit),...);
}

template<unsigned S,std::size_t... J> OC_FORCEINLINE void asymmetricMaskedUpdate(
    __m128i* state,const __m128i* syndrome,__m128i dot,u32 mask,std::index_sequence<J...>) {
    ((state[J]=_mm_xor_si128(_mm_xor_si128(state[J],syndrome[J]),
        _mm_and_si128(dot,_mm_set1_epi32(-int((mask>>J)&1))))),...);
}

template<class Map,bool Quarter,bool Masked,class Emit> OC_FORCEINLINE void asymmetricReverse(
    const block* in,std::size_t n,const u32* masks,Emit&& emit) {
    if constexpr(!Quarter) innerReverse<Map>(in,n,masks,emit);
    else {
        static_assert(Map::T==128 && Map::S==19);
        alignas(32) __m128i state[Map::S]{},syndrome[Map::S],raw[Map::T];
        const auto epochs=n/Map::T;
        for(std::size_t epoch=epochs;epoch-->0;) {
            const auto base=epoch*Map::T;
            if(epoch+1==epochs) {
                for(unsigned p=Map::T;p-->0;) {raw[p]=in[base+p].mData;emit(base+p,in[base+p]);}
            } else {
#if W5_SHARED_EMISSION==2
                Map::emitGrouped(in+base,raw,state,base,emit);
#elif W5_SHARED_EMISSION==1
                Map::emitShared(in+base,raw,state,base,emit);
#else
                asymmetricPoints<Map>(in+base,raw,state,base,emit,std::make_index_sequence<Map::T>{});
#endif
            }
            if(epoch==0) break;
            zeta<Map::T>(raw);Map::finish(raw,syndrome);
            if(epoch+1==epochs) std::memcpy(state,syndrome,sizeof(state));
            else {
                // Sparse emission needs no nibble tables. Compute the one
                // remaining random parity directly (setup masks are public).
                auto u=masks[2*epoch];auto dot=_mm_setzero_si128();
                while(u) {const auto j=std::countr_zero(u);u&=u-1;dot=_mm_xor_si128(dot,state[j]);}
                auto v=masks[2*epoch+1];
                if constexpr(Masked) asymmetricMaskedUpdate<Map::S>(state,syndrome,dot,v,std::make_index_sequence<Map::S>{});
                else {
                    while(v) {const auto j=std::countr_zero(v);v&=v-1;state[j]=_mm_xor_si128(state[j],dot);}
                    for(unsigned j=0;j<Map::S;++j) state[j]=_mm_xor_si128(state[j],syndrome[j]);
                }
            }
        }
    }
}
}
