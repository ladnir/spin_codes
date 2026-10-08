#pragma once
#include "HdBlockwiseOuter256.h"

namespace spin::research::char2::hdsumrankouter {
// Same full-row FFT schedule as the frozen GRS256 control. The only
// arithmetic addition is the independently sampled GF65536 output scalar.
// Coefficient planes have 'rows' entries; payload planes have 'stride'.
// This research path requires full four-row groups, not shortened messages.
#if defined(_MSC_VER)
#define HDBWO_NOINLINE __declspec(noinline)
#else
#define HDBWO_NOINLINE __attribute__((noinline))
#endif
using hdfast::outer_detail::normalized;
using hdfast::outer_detail::pair;

template<unsigned Offset,bool Affine,bool Stream>
HDBWO_NOINLINE void leaf16Stream(const __m512i* source,Block* planes,std::size_t stride,
                                    const wide16::Coefficient* diagonal,std::size_t rows) {
#define HDBWO_LOAD(I) auto x##I=_mm512_load_si512(source+I)
    HDBWO_LOAD(0);HDBWO_LOAD(1);HDBWO_LOAD(2);HDBWO_LOAD(3);
    HDBWO_LOAD(4);HDBWO_LOAD(5);HDBWO_LOAD(6);HDBWO_LOAD(7);
    HDBWO_LOAD(8);HDBWO_LOAD(9);HDBWO_LOAD(10);HDBWO_LOAD(11);
    HDBWO_LOAD(12);HDBWO_LOAD(13);HDBWO_LOAD(14);HDBWO_LOAD(15);
#undef HDBWO_LOAD
#define HDBWO_PAIR(A,B,LEVEL,DELTA) pair<Affine,normalized(LEVEL,Offset^DELTA)>(x##A,x##B)
    HDBWO_PAIR(0,8,3,0);HDBWO_PAIR(1,9,3,0);
    HDBWO_PAIR(2,10,3,0);HDBWO_PAIR(3,11,3,0);
    HDBWO_PAIR(4,12,3,0);HDBWO_PAIR(5,13,3,0);
    HDBWO_PAIR(6,14,3,0);HDBWO_PAIR(7,15,3,0);
    HDBWO_PAIR(0,4,2,0);HDBWO_PAIR(1,5,2,0);
    HDBWO_PAIR(2,6,2,0);HDBWO_PAIR(3,7,2,0);
    HDBWO_PAIR(8,12,2,8);HDBWO_PAIR(9,13,2,8);
    HDBWO_PAIR(10,14,2,8);HDBWO_PAIR(11,15,2,8);
    HDBWO_PAIR(0,2,1,0);HDBWO_PAIR(1,3,1,0);
    HDBWO_PAIR(4,6,1,4);HDBWO_PAIR(5,7,1,4);
    HDBWO_PAIR(8,10,1,8);HDBWO_PAIR(9,11,1,8);
    HDBWO_PAIR(12,14,1,12);HDBWO_PAIR(13,15,1,12);
    HDBWO_PAIR(0,1,0,0);HDBWO_PAIR(2,3,0,2);
    HDBWO_PAIR(4,5,0,4);HDBWO_PAIR(6,7,0,6);
    HDBWO_PAIR(8,9,0,8);HDBWO_PAIR(10,11,0,10);
    HDBWO_PAIR(12,13,0,12);HDBWO_PAIR(14,15,0,14);
#undef HDBWO_PAIR
#define HDBWO_STORE(I) do { \
    const auto scaled=wide16::multiply512(x##I,diagonal+I*rows); \
    if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(planes+I*stride),scaled); \
    else _mm512_storeu_si512(planes+I*stride,scaled); \
} while(false)
    HDBWO_STORE(0);HDBWO_STORE(1);HDBWO_STORE(2);HDBWO_STORE(3);
    HDBWO_STORE(4);HDBWO_STORE(5);HDBWO_STORE(6);HDBWO_STORE(7);
    HDBWO_STORE(8);HDBWO_STORE(9);HDBWO_STORE(10);HDBWO_STORE(11);
    HDBWO_STORE(12);HDBWO_STORE(13);HDBWO_STORE(14);HDBWO_STORE(15);
#undef HDBWO_STORE
}

template<unsigned Offset,bool Affine,bool Stream=false>
HDBWO_NOINLINE void coset128(const Block* source,Block* planes,std::size_t stride,
                                const wide16::Coefficient* diagonal,std::size_t rows) {
    static_assert(Offset==0||Offset==128);
    alignas(64) __m512i values[128];
    for(unsigned base=0;base!=16;base+=2){
        c8Unroll(std::make_index_sequence<2>{},[&](auto index){
            const auto j=base+unsigned(index);
            auto a=_mm512_loadu_si512(source+4*j);
            auto b=_mm512_loadu_si512(source+4*(16+j));
            auto c=_mm512_loadu_si512(source+4*(32+j));
            auto d=_mm512_loadu_si512(source+4*(48+j));
            auto e=_mm512_loadu_si512(source+4*(64+j));
            auto f=_mm512_loadu_si512(source+4*(80+j));
            auto g=_mm512_loadu_si512(source+4*(96+j));
            auto h=_mm512_loadu_si512(source+4*(112+j));
            pair<Affine,normalized(6,Offset)>(a,e);
            pair<Affine,normalized(6,Offset)>(b,f);
            pair<Affine,normalized(6,Offset)>(c,g);
            pair<Affine,normalized(6,Offset)>(d,h);
            pair<Affine,normalized(5,Offset)>(a,c);
            pair<Affine,normalized(5,Offset)>(b,d);
            pair<Affine,normalized(5,Offset+64)>(e,g);
            pair<Affine,normalized(5,Offset+64)>(f,h);
            pair<Affine,normalized(4,Offset)>(a,b);
            pair<Affine,normalized(4,Offset+32)>(c,d);
            pair<Affine,normalized(4,Offset+64)>(e,f);
            pair<Affine,normalized(4,Offset+96)>(g,h);
            _mm512_store_si512(values+j,a);
            _mm512_store_si512(values+16+j,b);
            _mm512_store_si512(values+32+j,c);
            _mm512_store_si512(values+48+j,d);
            _mm512_store_si512(values+64+j,e);
            _mm512_store_si512(values+80+j,f);
            _mm512_store_si512(values+96+j,g);
            _mm512_store_si512(values+112+j,h);
        });
    }
    leaf16Stream<Offset,Affine,Stream>(values,planes,stride,diagonal,rows);
    leaf16Stream<Offset+16,Affine,Stream>(values+16,planes+16*stride,stride,diagonal+16*rows,rows);
    leaf16Stream<Offset+32,Affine,Stream>(values+32,planes+32*stride,stride,diagonal+32*rows,rows);
    leaf16Stream<Offset+48,Affine,Stream>(values+48,planes+48*stride,stride,diagonal+48*rows,rows);
    leaf16Stream<Offset+64,Affine,Stream>(values+64,planes+64*stride,stride,diagonal+64*rows,rows);
    leaf16Stream<Offset+80,Affine,Stream>(values+80,planes+80*stride,stride,diagonal+80*rows,rows);
    leaf16Stream<Offset+96,Affine,Stream>(values+96,planes+96*stride,stride,diagonal+96*rows,rows);
    leaf16Stream<Offset+112,Affine,Stream>(values+112,planes+112*stride,stride,diagonal+112*rows,rows);
}


#if defined(HD_SUMRANK_SHARED_OUTER_LOADS)
// Exact-map scheduling control: read each source vector once for the two
// evaluation cosets. Keep both top-stage results in explicitly aligned
// scratch; the established sixteen-register leaf circuits are unchanged.
template<unsigned Offset,bool Affine>
inline void top128(__m512i* values,unsigned j,__m512i a,__m512i b,__m512i c,__m512i d,
                   __m512i e,__m512i f,__m512i g,__m512i h){
    pair<Affine,normalized(6,Offset)>(a,e);pair<Affine,normalized(6,Offset)>(b,f);
    pair<Affine,normalized(6,Offset)>(c,g);pair<Affine,normalized(6,Offset)>(d,h);
    pair<Affine,normalized(5,Offset)>(a,c);pair<Affine,normalized(5,Offset)>(b,d);
    pair<Affine,normalized(5,Offset+64)>(e,g);pair<Affine,normalized(5,Offset+64)>(f,h);
    pair<Affine,normalized(4,Offset)>(a,b);pair<Affine,normalized(4,Offset+32)>(c,d);
    pair<Affine,normalized(4,Offset+64)>(e,f);pair<Affine,normalized(4,Offset+96)>(g,h);
    _mm512_store_si512(values+j,a);_mm512_store_si512(values+16+j,b);
    _mm512_store_si512(values+32+j,c);_mm512_store_si512(values+48+j,d);
    _mm512_store_si512(values+64+j,e);_mm512_store_si512(values+80+j,f);
    _mm512_store_si512(values+96+j,g);_mm512_store_si512(values+112+j,h);
}
template<bool Stream>
HDBWO_NOINLINE void sharedCosetLoads(const Block* source,Block* planes,std::size_t stride,
                                     const wide16::Coefficient* diagonal,std::size_t rows){
    alignas(64) __m512i values[256];
    for(unsigned j=0;j!=16;++j){
        const auto a=_mm512_loadu_si512(source+4*j),b=_mm512_loadu_si512(source+4*(16+j));
        const auto c=_mm512_loadu_si512(source+4*(32+j)),d=_mm512_loadu_si512(source+4*(48+j));
        const auto e=_mm512_loadu_si512(source+4*(64+j)),f=_mm512_loadu_si512(source+4*(80+j));
        const auto g=_mm512_loadu_si512(source+4*(96+j)),h=_mm512_loadu_si512(source+4*(112+j));
        top128<0,true>(values,j,a,b,c,d,e,f,g,h);
        top128<128,true>(values+128,j,a,b,c,d,e,f,g,h);
    }
    c8Unroll(std::make_index_sequence<16>{},[&](auto i){
        constexpr unsigned offset=16*decltype(i)::value;
        leaf16Stream<offset,true,Stream>(values+offset,planes+offset*stride,stride,diagonal+offset*rows,rows);
    });
}
#endif

template<std::size_t TileRows=64,bool Stream=true,unsigned PrefetchRows=0>
HDBWO_NOINLINE void forward(std::size_t k,std::size_t stride,
    const Block* input,Block* planes,const wide16::Coefficient* diagonal,std::size_t diagonalStride=0) {
    static_assert(TileRows>0&&TileRows%4==0);
    const auto rows=k/128;
    if(!diagonalStride)diagonalStride=rows;
    if(diagonalStride<rows)throw std::invalid_argument("short coefficient plane stride");
    if(k%512||stride%4||stride<rows)
        throw std::invalid_argument("scaled GRS256 requires full four-row groups");
    if constexpr(Stream) {
        if(reinterpret_cast<std::uintptr_t>(planes)&63U) {
            forward<TileRows,false,PrefetchRows>(k,stride,input,planes,diagonal,diagonalStride);
            return;
        }
    }
#if defined(HD_SUMRANK_SHARED_OUTER_LOADS)
    for(std::size_t row=0;row!=rows;row+=4)
        sharedCosetLoads<Stream>(input+row*128,planes+row,stride,diagonal+row,diagonalStride);
#else
    for(std::size_t first=0;first<rows;first+=TileRows) {
        const auto last=std::min(first+TileRows,rows);
        for(auto row=first;row!=last;row+=4){
            if constexpr(PrefetchRows){
                static_assert(PrefetchRows%4==0);
                if(row+PrefetchRows<rows){
                    const auto* next=reinterpret_cast<const char*>(input+(row+PrefetchRows)*128);
                    for(unsigned line=0;line!=128;line+=8)
                        c8Unroll(std::make_index_sequence<8>{},[&](auto j){
                            _mm_prefetch(next+64*(line+unsigned(j)),PrefetchRows>4?_MM_HINT_T1:_MM_HINT_T0);
                        });
                }
            }
            coset128<0,true,Stream>(input+row*128,planes+row,stride,diagonal+row,diagonalStride);
        }
        for(auto row=first;row!=last;row+=4)
            coset128<128,true,Stream>(input+row*128,planes+128*stride+row,stride,diagonal+128*diagonalStride+row,diagonalStride);
    }
#endif
    for(unsigned c=0;c!=256;++c)
        std::memset(planes+c*stride+rows,0,(stride-rows)*sizeof(Block));
    if constexpr(Stream)_mm_sfence();
}
#undef HDBWO_NOINLINE
} // namespace spin::research::char2::hdsumrankouter
