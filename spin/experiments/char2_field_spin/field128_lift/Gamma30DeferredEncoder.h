#pragma once
#include "Gamma30Encoder.h"
#include "SparseByteField128ToeplitzDeferred.h"

// Exact-distribution A/B variant. The Plan IS the frozen gamma30 Plan: no
// resampling, transformed keys, basis conversion, or route changes. Only
// the outer Toeplitz arithmetic defers its linear reductions to each output.
namespace spin::research::char2::field128_lift::gamma30_deferred {
namespace baseline=spin::research::char2::field128_lift::gamma30;
namespace literal=baseline::literal;
namespace field=baseline::field;
namespace old=baseline::old;
using Plan=baseline::Plan;
using Element=baseline::Element;
using Block=baseline::Block;
using State=baseline::State;
using Packet=baseline::Packet;
using Refresh=baseline::Refresh;

// Direct aliases, not wrappers or replacement implementations. These also
// let the tagged dense-oracle driver exercise the unchanged fixed maps.
using baseline::emitAdd;
using baseline::feedback;
using baseline::refreshMap;
using baseline::finish;
using baseline::loadPacket;
using baseline::storePacket;
using baseline::inner;

#if defined(_MSC_VER)
#define GAMMA30_DEFER_INLINE __forceinline
#define GAMMA30_DEFER_NOINLINE __declspec(noinline)
#else
#define GAMMA30_DEFER_INLINE inline __attribute__((always_inline))
#define GAMMA30_DEFER_NOINLINE inline __attribute__((noinline))
#endif

// The traversal below is an explicit copy of LiteralEncoder's outer path.
// The sole arithmetic substitution is this deferred::toeplitz4 call.
template<bool Stream>
GAMMA30_DEFER_INLINE void mixStore(__m512i* column,unsigned j,Block* planes,std::size_t stride,const Element* keys){
    field::deferred::toeplitz4(column,keys+28*j);
    c8Unroll(std::make_index_sequence<4>{},[&](auto r){
        auto* out=planes+(4*j+r)*stride;
        if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(out),column[r]);
        else _mm512_storeu_si512(out,column[r]);
    });
}

template<unsigned Offset,bool Stream>
GAMMA30_DEFER_NOINLINE void coset(const Block* source,Block* planes,std::size_t stride,const Element* keys){
    static_assert(Offset==0||Offset==32);
    alignas(64) __m512i values[128];
    for(unsigned r=0;r!=4;++r){
        c8Unroll(std::make_index_sequence<32>{},[&](auto j){values[32*r+j]=_mm512_loadu_si512(source+4*(32*r+j));});
        old::rsPrefix<32,Offset>(values+32*r);
    }
    for(unsigned j=0;j!=32;j+=2){
        alignas(64) __m512i even[4],odd[4];
        const auto coefficient=_mm512_set1_epi8(static_cast<char>(Offset+j));
        c8Unroll(std::make_index_sequence<4>{},[&](auto r){
            const auto a=values[32*r+j],b=values[32*r+j+1];
            even[r]=_mm512_xor_si512(a,_mm512_gf2p8mul_epi8(b,coefficient));
            odd[r]=_mm512_xor_si512(even[r],b);
        });
        mixStore<Stream>(even,Offset+j,planes,stride,keys);
        mixStore<Stream>(odd,Offset+j+1,planes,stride,keys);
    }
}

template<bool Stream>
GAMMA30_DEFER_NOINLINE void outer(const literal::Plan& p,const Block* input,Block* planes){
    if constexpr(Stream)if(reinterpret_cast<std::uintptr_t>(planes)&63U){gamma30_deferred::outer<false>(p,input,planes);return;}
    constexpr std::size_t tile=64;
    for(std::size_t first=0;first<p.rows;first+=tile){
        const auto last=std::min(first+tile,p.rows);
        for(auto g=first;g<last;g+=4)coset<0,Stream>(input+g*128,planes+g,p.stride,p.outerKeys.data()+g*448);
        for(auto g=first;g<last;g+=4)coset<32,Stream>(input+g*128,planes+g,p.stride,p.outerKeys.data()+g*448);
    }
    if constexpr(Stream)_mm_sfence();
    for(unsigned c=0;c!=256;++c){
        auto* plane=planes+c*p.stride;
        plane[p.rows]=plane[0];
        std::memset(plane+p.rows+1,0,3*sizeof(Block));
    }
}

template<bool Stream=true>
inline void forward(const Plan& plan,const Block* input,Block* output,Block* scratch){
    // Same complete-outer-before-inner alias contract and store policies.
    gamma30_deferred::outer<true>(plan.outerPlan,input,scratch);
    baseline::inner<Stream>(plan,scratch,output);
}

#undef GAMMA30_DEFER_INLINE
#undef GAMMA30_DEFER_NOINLINE
} // namespace spin::research::char2::field128_lift::gamma30_deferred
