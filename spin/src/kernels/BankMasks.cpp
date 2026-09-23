#include "BankKernel.h"
#include <immintrin.h>

namespace spin::detail::kernel {
bool bankMasks512(std::uint32_t* masks,std::size_t count,unsigned bits,std::uint64_t seed) noexcept {
    if(count%16)return false;
    constexpr std::uint64_t step=0x9e3779b97f4a7c15ULL;
    const auto offsets=_mm512_setr_epi64(step,3*step,5*step,7*step,9*step,11*step,13*step,15*step);
    const auto mask=_mm256_set1_epi32((1U<<bits)-1),one=_mm256_set1_epi32(1),zero=_mm256_setzero_si256();
    auto mix=[](__m512i x){
        x=_mm512_mullo_epi64(_mm512_xor_si512(x,_mm512_srli_epi64(x,30)),_mm512_set1_epi64(0xbf58476d1ce4e5b9ULL));
        x=_mm512_mullo_epi64(_mm512_xor_si512(x,_mm512_srli_epi64(x,27)),_mm512_set1_epi64(0x94d049bb133111ebULL));
        return _mm512_xor_si512(x,_mm512_srli_epi64(x,31));
    };
    for(std::size_t i=0;i<count;i+=16,seed+=16*step) {
        const auto first=_mm512_add_epi64(_mm512_set1_epi64(seed),offsets);
        const auto u=_mm256_and_si256(_mm512_cvtepi64_epi32(mix(first)),mask);
        // A rejected u changes stream consumption. The caller rewrites every
        // entry with the scalar sampler, including this partial prefix.
        if(_mm256_movemask_epi8(_mm256_cmpeq_epi32(u,zero)))return false;
        auto v=_mm256_and_si256(_mm512_cvtepi64_epi32(mix(_mm512_add_epi64(first,_mm512_set1_epi64(step)))),mask);
        const auto parity=_mm256_and_si256(_mm256_popcnt_epi32(_mm256_and_si256(u,v)),one);
        v=_mm256_xor_si256(v,_mm256_and_si256(_mm256_sub_epi32(zero,parity),_mm256_and_si256(u,_mm256_sub_epi32(zero,u))));
        const auto lo=_mm256_unpacklo_epi32(u,v),hi=_mm256_unpackhi_epi32(u,v);
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(masks+i),_mm256_permute2x128_si256(lo,hi,0x20));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(masks+i+8),_mm256_permute2x128_si256(lo,hi,0x31));
    }
    return true;
}
}
