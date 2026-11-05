#pragma once
#include "Transform256.h"
#include <algorithm>
#include <cstring>
#include <utility>

// Exact RS[256,128] novel-basis evaluation. No changes to evaluation points,
// coefficient order, field representation, or routing distribution.
namespace spin::research::one_route::rs_outer {
namespace d = transform_cantor_detail;
using E = spin::research::one_route::Element;
#define RS_INLINE inline __attribute__((always_inline))

template<unsigned Offset>
RS_INLINE void leaf8Values(__m512i& a,__m512i& b,__m512i& c,__m512i& e,
                          __m512i& f,__m512i& g,__m512i& h,__m512i& i) noexcept {
    d::pair<d::twiddle(2,Offset),false>(a,f);
    d::pair<d::twiddle(2,Offset),false>(b,g);
    d::pair<d::twiddle(2,Offset),false>(c,h);
    d::pair<d::twiddle(2,Offset),false>(e,i);
    d::pair<d::twiddle(1,Offset),false>(a,c);
    d::pair<d::twiddle(1,Offset),false>(b,e);
    d::pair<d::twiddle(1,Offset^4),false>(f,h);
    d::pair<d::twiddle(1,Offset^4),false>(g,i);
    d::pair<d::twiddle(0,Offset),false>(a,b);
    d::pair<d::twiddle(0,Offset^2),false>(c,e);
    d::pair<d::twiddle(0,Offset^4),false>(f,g);
    d::pair<d::twiddle(0,Offset^6),false>(h,i);
}
template<unsigned Offset>
RS_INLINE void leaf8Four(E* p) noexcept {
    auto a=_mm512_load_si512(p+0), b=_mm512_load_si512(p+4);
    auto c=_mm512_load_si512(p+8), e=_mm512_load_si512(p+12);
    auto f=_mm512_load_si512(p+16), g=_mm512_load_si512(p+20);
    auto h=_mm512_load_si512(p+24), i=_mm512_load_si512(p+28);
    leaf8Values<Offset>(a,b,c,e,f,g,h,i);
    _mm512_store_si512(p+0,a); _mm512_store_si512(p+4,b);
    _mm512_store_si512(p+8,c); _mm512_store_si512(p+12,e);
    _mm512_store_si512(p+16,f); _mm512_store_si512(p+20,g);
    _mm512_store_si512(p+24,h); _mm512_store_si512(p+28,i);
}
template<unsigned Level,unsigned Offset,bool Finish=true>
RS_INLINE void fftFour(E* p) noexcept {
    if constexpr(Level==3) {if constexpr(Finish)leaf8Four<Offset>(p);}
    else {
        d::stage<Level,Offset,false,true>(p);
        fftFour<Level-1,Offset,Finish>(p);
        fftFour<Level-1,Offset^(1U<<(Level-1)),Finish>(p+4*(1U<<(Level-1)));
    }
}
RS_INLINE __m512i loadFour(const E* p,unsigned j) noexcept {
    auto a=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(p+j)));
    a=_mm512_inserti32x4(a,_mm_loadu_si128(reinterpret_cast<const __m128i*>(p+128+j)),1);
    a=_mm512_inserti32x4(a,_mm_loadu_si128(reinterpret_cast<const __m128i*>(p+256+j)),2);
    return _mm512_inserti32x4(a,_mm_loadu_si128(reinterpret_cast<const __m128i*>(p+384+j)),3);
}
// Mode 0: frozen transform. 1: register leaves. 2: fuse input loading with
// the first two levels, including the known-zero upper coefficient half.
// Mode 3 leaves the last three levels for the fused output kernel.
template<unsigned Mode>
RS_INLINE void encodeFour(const E* input,E* p) noexcept {
    if constexpr(Mode<2) {
        for(unsigned j=0;j<128;++j) _mm512_store_si512(p+4*j,loadFour(input,j));
        std::memset(p+512,0,512*sizeof(E));
        if constexpr(Mode==0) Transform256Cantor::evaluateNovelFour(p);
        else fftFour<8,0>(p);
    } else {
        for(unsigned j=0;j<64;++j) {
            auto a=loadFour(input,j), b=loadFour(input,j+64);
            auto c=a,e=b;
            d::pair<d::twiddle(6,0),false>(a,b);
            d::pair<d::twiddle(6,128),false>(c,e);
            _mm512_store_si512(p+4*j,a); _mm512_store_si512(p+4*(j+64),b);
            _mm512_store_si512(p+4*(j+128),c); _mm512_store_si512(p+4*(j+192),e);
        }
        fftFour<6,0,Mode!=3>(p); fftFour<6,64,Mode!=3>(p+256);
        fftFour<6,128,Mode!=3>(p+512); fftFour<6,192,Mode!=3>(p+768);
    }
}
template<unsigned Offset>
RS_INLINE void finishTile(std::size_t count,std::size_t blocks,const E* tile,E* output) {
    for(std::size_t r=0;r<count;r+=4){
        const auto* p=tile+r*256+4*Offset;
        auto a=_mm512_load_si512(p+0),b=_mm512_load_si512(p+4);
        auto c=_mm512_load_si512(p+8),e=_mm512_load_si512(p+12);
        auto f=_mm512_load_si512(p+16),g=_mm512_load_si512(p+20);
        auto h=_mm512_load_si512(p+24),i=_mm512_load_si512(p+28);
        leaf8Values<Offset>(a,b,c,e,f,g,h,i);
        _mm512_store_si512(output+(Offset+0)*blocks+r,a);_mm512_store_si512(output+(Offset+1)*blocks+r,b);
        _mm512_store_si512(output+(Offset+2)*blocks+r,c);_mm512_store_si512(output+(Offset+3)*blocks+r,e);
        _mm512_store_si512(output+(Offset+4)*blocks+r,f);_mm512_store_si512(output+(Offset+5)*blocks+r,g);
        _mm512_store_si512(output+(Offset+6)*blocks+r,h);_mm512_store_si512(output+(Offset+7)*blocks+r,i);
    }
}
template<std::size_t... J>
RS_INLINE void finishTiles(std::size_t count,std::size_t blocks,const E* tile,E* output,std::index_sequence<J...>){
    (finishTile<8*J>(count,blocks,tile,output),...);
}
template<unsigned Rows>
__attribute__((noinline)) void outerFused(std::size_t blocks,const E* input,E* output){
    static_assert(Rows>=4&&Rows<=256&&Rows%4==0);
    alignas(64) E tile[Rows*256];
    for(std::size_t first=0;first<blocks;first+=Rows){
        const auto count=std::min<std::size_t>(Rows,blocks-first);
        for(std::size_t r=0;r<count;r+=4)encodeFour<3>(input+(first+r)*128,tile+r*256);
        finishTiles(count,blocks,tile,output+first,std::make_index_sequence<32>{});
    }
}

// Block stores without changing the regional layout expected by the route.
// Workspace lives on the ordinary stack; there are no allocations or callbacks.
template<unsigned Rows,unsigned Mode,bool Tiled=false,bool Stream=false>
__attribute__((noinline)) void outer(std::size_t blocks,const E* input,E* output) {
    static_assert(Rows>=4 && Rows<=256 && Rows%4==0);
    alignas(64) E tile[Rows*256];
    for(std::size_t first=0;first<blocks;first+=Rows) {
        const auto count=std::min<std::size_t>(Rows,blocks-first);
        for(std::size_t r=0;r<count;r+=4) encodeFour<Mode>(input+(first+r)*128,tile+r*256);
        for(unsigned j=0;j<256;++j) {
            E* dest=Tiled?output+first*256+j*count:output+j*blocks+first;
            for(std::size_t r=0;r<count;r+=4) {
                const auto v=_mm512_load_si512(tile+r*256+4*j);
                if constexpr(Stream) _mm512_stream_si512(reinterpret_cast<__m512i*>(dest+r),v);
                else _mm512_store_si512(dest+r,v);
            }
        }
    }
    if constexpr(Stream) _mm_sfence();
}
RS_INLINE void encodeOne(const E* input,E* p) noexcept {
    for(unsigned j=0;j<64;j+=4) {
        auto a=_mm512_loadu_si512(input+j),b=_mm512_loadu_si512(input+64+j);
        auto c=a,e=b;
        d::pair<d::twiddle(6,0),false>(a,b);
        d::pair<d::twiddle(6,128),false>(c,e);
        _mm512_store_si512(p+j,a);_mm512_store_si512(p+64+j,b);
        _mm512_store_si512(p+128+j,c);_mm512_store_si512(p+192+j,e);
    }
    d::fft<6,0,false,false>(p);d::fft<6,64,false,false>(p+64);
    d::fft<6,128,false,false>(p+128);d::fft<6,192,false,false>(p+192);
}
template<bool Stream>
RS_INLINE void store(E* dest,__m512i value) noexcept {
    if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(dest),value);
    else _mm512_store_si512(dest,value);
}
template<unsigned Rows,bool Stream=false>
__attribute__((noinline)) void outerOne(std::size_t blocks,const E* input,E* output) {
    static_assert(Rows>=4&&Rows<=256&&Rows%4==0);
    alignas(64) E tile[Rows*256];
    for(std::size_t first=0;first<blocks;first+=Rows) {
        const auto count=std::min<std::size_t>(Rows,blocks-first);
        for(std::size_t r=0;r<count;++r)encodeOne(input+(first+r)*128,tile+r*256);
        for(unsigned j=0;j<256;j+=4)for(std::size_t r=0;r<count;r+=4){
            const auto a=_mm512_load_si512(tile+(r+0)*256+j),b=_mm512_load_si512(tile+(r+1)*256+j);
            const auto c=_mm512_load_si512(tile+(r+2)*256+j),e=_mm512_load_si512(tile+(r+3)*256+j);
            const auto ab0=_mm512_shuffle_i32x4(a,b,0x44),ab1=_mm512_shuffle_i32x4(a,b,0xee);
            const auto ce0=_mm512_shuffle_i32x4(c,e,0x44),ce1=_mm512_shuffle_i32x4(c,e,0xee);
            store<Stream>(output+(j+0)*blocks+first+r,_mm512_shuffle_i32x4(ab0,ce0,0x88));
            store<Stream>(output+(j+1)*blocks+first+r,_mm512_shuffle_i32x4(ab0,ce0,0xdd));
            store<Stream>(output+(j+2)*blocks+first+r,_mm512_shuffle_i32x4(ab1,ce1,0x88));
            store<Stream>(output+(j+3)*blocks+first+r,_mm512_shuffle_i32x4(ab1,ce1,0xdd));
        }
    }
    if constexpr(Stream)_mm_sfence();
}
template<unsigned Rows,bool Tiled>
__attribute__((noinline)) void route(std::size_t blocks,const std::uint32_t* labels,const E* input,E* output) {
    for(unsigned j=0;j<256;++j) {
        auto* dest=output+j*blocks; const auto* ix=labels+j*blocks;
        for(std::size_t i=0;i<blocks;++i) {
            const auto row=ix[i];
            if constexpr(Tiled) {
                const auto first=(row/Rows)*Rows;
                const auto width=std::min<std::size_t>(Rows,blocks-first);
                dest[i]=input[first*256+j*width+(row-first)];
            } else dest[i]=input[j*blocks+row];
        }
    }
}
RS_INLINE __m512i gatherFour(const E* source,const std::uint32_t* labels) noexcept {
    auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[0])));
    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[1])),1);
    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[2])),2);
    return _mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(source+labels[3])),3);
}
template<unsigned Prefetch>
__attribute__((noinline)) void routeBatched(std::size_t blocks,const std::uint32_t* labels,const E* input,E* output) {
    for(unsigned j=0;j<256;++j) {
        const auto* source=input+j*blocks;const auto* ix=labels+j*blocks;auto* dest=output+j*blocks;
        if constexpr(Prefetch==1)
            for(std::size_t i=0;i<blocks;i+=4)_mm_prefetch(reinterpret_cast<const char*>(source+i),_MM_HINT_T1);
        if constexpr(Prefetch==3) {
            // Demand loads populate L2 reliably; a burst of nonbinding prefetch
            // hints can be dropped. Four independent chains preserve load MLP.
            auto a=_mm512_setzero_si512(),b=a,c=a,e=a;
            std::size_t k=0;
            for(;k+16<=blocks;k+=16){
                a=_mm512_xor_si512(a,_mm512_load_si512(source+k));
                b=_mm512_xor_si512(b,_mm512_load_si512(source+k+4));
                c=_mm512_xor_si512(c,_mm512_load_si512(source+k+8));
                e=_mm512_xor_si512(e,_mm512_load_si512(source+k+12));
            }
            asm volatile("" : : "v"(a),"v"(b),"v"(c),"v"(e) : "memory");
        }
        std::size_t i=0;
        for(;i+16<=blocks;i+=16) {
            if constexpr(Prefetch==2) {
                if(i+48<=blocks)for(unsigned k=0;k<16;++k)_mm_prefetch(reinterpret_cast<const char*>(source+ix[i+32+k]),_MM_HINT_T0);
            }
            if constexpr(Prefetch==4) {
                if(i+144<=blocks)for(unsigned k=0;k<16;++k)_mm_prefetch(reinterpret_cast<const char*>(source+ix[i+128+k]),_MM_HINT_T0);
            }
            const auto a=gatherFour(source,ix+i),b=gatherFour(source,ix+i+4);
            const auto c=gatherFour(source,ix+i+8),e=gatherFour(source,ix+i+12);
            _mm512_store_si512(dest+i,a);_mm512_store_si512(dest+i+4,b);
            _mm512_store_si512(dest+i+8,c);_mm512_store_si512(dest+i+12,e);
        }
        for(;i<blocks;++i)dest[i]=source[ix[i]];
    }
}
// Selected K20 path. Contract: blocks is a positive multiple of four;
// input contains blocks*128 elements, output blocks*256; output is 64-byte
// aligned and does not overlap input. The outer retains coordinate-major
// layout. Routing writes a separate 64-byte-aligned output with exactly the
// original per-coordinate labels. These routines do not allocate or sample.
inline void encodeRegional(std::size_t blocks,const E* input,E* output){
    outerFused<128>(blocks,input,output);
}
inline void routeRegional(std::size_t blocks,const std::uint32_t* labels,const E* input,E* output){
    routeBatched<3>(blocks,labels,input,output);
}
#undef RS_INLINE
} // namespace spin::research::one_route::rs_outer
