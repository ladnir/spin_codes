#pragma once
#include "../feistel/Direct.h"
namespace spin::experimental::bank {
// Same SplitMix word stream and orthogonality correction as makeMasksK18.
// A rare rejected u falls back to the scalar stream, preserving consumption.
inline std::vector<std::uint32_t> fastMasks(std::uint64_t seed,bool* fallback=nullptr) {
    if(fallback)*fallback=false;
#if defined(__AVX512DQ__) && defined(__AVX512VPOPCNTDQ__) && defined(__AVX512VL__)
    std::vector<std::uint32_t> masks(8192);
    constexpr std::uint64_t step=0x9e3779b97f4a7c15ULL;
    const auto offsets=_mm512_setr_epi64(step,3*step,5*step,7*step,9*step,11*step,13*step,15*step);
    const auto mask=_mm256_set1_epi32((1U<<19)-1),one=_mm256_set1_epi32(1),zero=_mm256_setzero_si256();
    auto mix=[](__m512i x){
        x=_mm512_mullo_epi64(_mm512_xor_si512(x,_mm512_srli_epi64(x,30)),_mm512_set1_epi64(0xbf58476d1ce4e5b9ULL));
        x=_mm512_mullo_epi64(_mm512_xor_si512(x,_mm512_srli_epi64(x,27)),_mm512_set1_epi64(0x94d049bb133111ebULL));
        return _mm512_xor_si512(x,_mm512_srli_epi64(x,31));
    };
    auto position=seed;
    for(unsigned i=0;i<masks.size();i+=16,position+=16*step) {
        const auto first=_mm512_add_epi64(_mm512_set1_epi64(position),offsets);
        const auto u=_mm256_and_si256(_mm512_cvtepi64_epi32(mix(first)),mask);
        if(_mm256_movemask_epi8(_mm256_cmpeq_epi32(u,zero))) {
            if(fallback)*fallback=true;
            return feistel::makeMasksK18(seed);
        }
        auto v=_mm256_and_si256(_mm512_cvtepi64_epi32(mix(_mm512_add_epi64(first,_mm512_set1_epi64(step)))),mask);
        const auto parity=_mm256_and_si256(_mm256_popcnt_epi32(_mm256_and_si256(u,v)),one);
        v=_mm256_xor_si256(v,_mm256_and_si256(_mm256_sub_epi32(zero,parity),_mm256_and_si256(u,_mm256_sub_epi32(zero,u))));
        const auto lo=_mm256_unpacklo_epi32(u,v),hi=_mm256_unpackhi_epi32(u,v);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(masks.data()+i),_mm256_permute2x128_si256(lo,hi,0x20));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(masks.data()+i+8),_mm256_permute2x128_si256(lo,hi,0x31));
    }
    return masks;
#else
    return feistel::makeMasksK18(seed);
#endif
}
}
