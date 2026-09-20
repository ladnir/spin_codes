// Literal C++ port of assemble<W> in Hypercat examples/spin_wide.rs.
// No bit transpose: these routines reorder whole 128-bit words into columns.
#include "Block.h"
#include <cstddef>
using osuCrypto::block;
template<unsigned W> static void assemble(const block* in,block* out,std::size_t n) {
    static_assert(W==2 || W==4);
    for(std::size_t j=0;j<n;++j) for(unsigned g=0;g<16;g+=4) {
        const auto* src=in+(g/W)*n*W+j*W;
        __m512i value;
        if constexpr(W==4) value=_mm512_loadu_si512(src);
        else value=_mm512_inserti64x4(_mm512_castsi256_si512(_mm256_loadu_si256((const __m256i*)src)),
                                    _mm256_loadu_si256((const __m256i*)(src+n*W)),1);
        _mm512_stream_si512((__m512i*)(out+j*16+g),value);
    }
    _mm_sfence();
}
extern "C" void spin_columns256(const block* in,block* out,std::size_t n) {assemble<2>(in,out,n);}
extern "C" void spin_columns512(const block* in,block* out,std::size_t n) {assemble<4>(in,out,n);}
