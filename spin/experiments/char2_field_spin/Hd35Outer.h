#pragma once
#include "HdSumrankOuter.h"

namespace spin::research::char2::hd35 {
#if defined(_MSC_VER)
#define HD35_NOINLINE __declspec(noinline)
#else
#define HD35_NOINLINE __attribute__((noinline))
#endif

// RS[128,64] in the existing novel-polynomial basis. Preserve the tuned
// sixteen-register leaf and fuse its two ancestors, omitting exactly the
// degree-64 level of the RS256 control. Four rows share every SIMD circuit.
template<unsigned Offset,bool Stream>
HD35_NOINLINE void coset64(const Block* source,Block* planes,std::size_t stride,
                          const wide16::Coefficient* diagonal,std::size_t rows){
    static_assert(Offset==0||Offset==64);
    using hdfast::outer_detail::pair;
    using hdfast::outer_detail::normalized;
    alignas(64) __m512i values[64];
    for(unsigned base=0;base!=16;base+=2){
        c8Unroll(std::make_index_sequence<2>{},[&](auto index){
            const auto j=base+unsigned(index);
            auto a=_mm512_loadu_si512(source+4*j);
            auto b=_mm512_loadu_si512(source+4*(16+j));
            auto c=_mm512_loadu_si512(source+4*(32+j));
            auto d=_mm512_loadu_si512(source+4*(48+j));
            pair<true,normalized(5,Offset)>(a,c);
            pair<true,normalized(5,Offset)>(b,d);
            pair<true,normalized(4,Offset)>(a,b);
            pair<true,normalized(4,Offset+32)>(c,d);
            _mm512_store_si512(values+j,a);_mm512_store_si512(values+16+j,b);
            _mm512_store_si512(values+32+j,c);_mm512_store_si512(values+48+j,d);
        });
    }
    hdsumrankouter::leaf16Stream<Offset,true,Stream>(values,planes,stride,diagonal,rows);
    hdsumrankouter::leaf16Stream<Offset+16,true,Stream>(values+16,planes+16*stride,stride,diagonal+16*rows,rows);
    hdsumrankouter::leaf16Stream<Offset+32,true,Stream>(values+32,planes+32*stride,stride,diagonal+32*rows,rows);
    hdsumrankouter::leaf16Stream<Offset+48,true,Stream>(values+48,planes+48*stride,stride,diagonal+48*rows,rows);
}

template<std::size_t TileRows=64,bool Stream=true>
HD35_NOINLINE void outer128(std::size_t k,std::size_t stride,const Block* input,
                           Block* planes,const wide16::Coefficient* diagonal){
    static_assert(TileRows>0&&TileRows%4==0);
    const auto rows=k/64;
    if(!k||k%256||stride%4||stride<rows)
        throw std::invalid_argument("scaled RS128 requires full four-row groups");
    if constexpr(Stream){
        if(reinterpret_cast<std::uintptr_t>(planes)&63U){
            outer128<TileRows,false>(k,stride,input,planes,diagonal);return;
        }
    }
    for(std::size_t first=0;first<rows;first+=TileRows){
        const auto last=std::min(first+TileRows,rows);
        for(auto row=first;row!=last;row+=4)
            coset64<0,Stream>(input+row*64,planes+row,stride,diagonal+row,rows);
        for(auto row=first;row!=last;row+=4)
            coset64<64,Stream>(input+row*64,planes+64*stride+row,stride,diagonal+64*rows+row,rows);
    }
    for(unsigned c=0;c!=128;++c)
        std::memset(planes+c*stride+rows,0,(stride-rows)*sizeof(Block));
    if constexpr(Stream)_mm_sfence();
}
} // namespace spin::research::char2::hd35
