#pragma once
#include "HdFastOuter.h"

namespace spin::research::char2::hdbwouter {

// Isolated forward GRS[256,128] cost-control outer, not a distance proof.
// Input rows are coefficients in the same normalized additive-polynomial
// basis as HdFastOuter, now retaining degrees0..127. All256 evaluation
// points0..255 are emitted; there is no Cauchy/systematic basis conversion.
inline constexpr unsigned outerLength=256;
inline constexpr unsigned outerDimension=128;
inline constexpr std::size_t sourceIndex(std::size_t row,unsigned degree) noexcept {
    return (row/4)*512+4*degree+row%4;
}

#if defined(_MSC_VER)
#define HDBWO_INLINE __forceinline
#define HDBWO_NOINLINE __declspec(noinline)
#else
#define HDBWO_INLINE inline __attribute__((always_inline))
#define HDBWO_NOINLINE __attribute__((noinline))
#endif

namespace detail {
using hdfast::outer_detail::normalized;
using hdfast::outer_detail::pair;

// Circuit-identical copy of HdFastOuter's sixteen-register leaf, differing
// only in its stores. The caller guarantees planes aligned to64 bytes and
// stride divisible by4 records. Its single fence belongs to the outer call.
template<unsigned Offset,bool Affine>
HDBWO_NOINLINE void leaf16Stream(const __m512i* source,Block* planes,std::size_t stride) {
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
#define HDBWO_STORE(I) _mm512_stream_si512(reinterpret_cast<__m512i*>(planes+I*stride),x##I)
    HDBWO_STORE(0);HDBWO_STORE(1);HDBWO_STORE(2);HDBWO_STORE(3);
    HDBWO_STORE(4);HDBWO_STORE(5);HDBWO_STORE(6);HDBWO_STORE(7);
    HDBWO_STORE(8);HDBWO_STORE(9);HDBWO_STORE(10);HDBWO_STORE(11);
    HDBWO_STORE(12);HDBWO_STORE(13);HDBWO_STORE(14);HDBWO_STORE(15);
#undef HDBWO_STORE
}

template<unsigned Offset,bool Affine,bool Stream>
HDBWO_INLINE void emitLeaf(const __m512i* source,Block* planes,std::size_t stride) {
    if constexpr(Stream)leaf16Stream<Offset,Affine>(source,planes,stride);
    else hdfast::outer_detail::leaf16<Offset,Affine>(source,planes,stride);
}

// Fuse all three ancestors of the existing sixteen-register FFT leaves.
// Eight source values produce corresponding inputs of all eight leaves.
// Only this8KiB aligned array is materialized; there is no separate level6
// output array and no coset64 scratch nested beneath it. Two independent
// eight-vector groups (two j values) are unrolled, bounding register demand.
template<unsigned Offset,bool Affine,bool Stream=false>
HDBWO_NOINLINE void coset128(const Block* source,Block* planes,std::size_t stride) {
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
    emitLeaf<Offset,Affine,Stream>(values,planes,stride);
    emitLeaf<Offset+16,Affine,Stream>(values+16,planes+16*stride,stride);
    emitLeaf<Offset+32,Affine,Stream>(values+32,planes+32*stride,stride);
    emitLeaf<Offset+48,Affine,Stream>(values+48,planes+48*stride,stride);
    emitLeaf<Offset+64,Affine,Stream>(values+64,planes+64*stride,stride);
    emitLeaf<Offset+80,Affine,Stream>(values+80,planes+80*stride,stride);
    emitLeaf<Offset+96,Affine,Stream>(values+96,planes+96*stride,stride);
    emitLeaf<Offset+112,Affine,Stream>(values+112,planes+112*stride,stride);
}

template<bool Affine>
HDBWO_INLINE void tile(const Block* source,Block* planes,std::size_t stride) {
    // Degree<128 makes omitted FFT level7 pure replication, so the two
    // offset128 subtransforms are the exact full256-point evaluation.
    coset128<0,Affine>(source,planes,stride);
    coset128<128,Affine>(source,planes+128*stride,stride);
}

// The shortened input group alone needs a copy. Keep its8KiB zero-padded
// input outside the normal outer stack frame and outside the hot path.
template<bool Affine>
HDBWO_NOINLINE void tail(const Block* source,std::size_t count,
                        Block* planes,std::size_t stride) {
    alignas(64) Block padded[512]{};
    std::memcpy(padded,source,count*sizeof(Block));
    tile<Affine>(padded,planes,stride);
}
} // namespace detail

// input[k] and planes[256*rows] must be disjoint. rows is the plane stride,
// must be a multiple of4, and must fit ceil(k/512) four-row source groups.
// Payload need not be SIMD aligned. Missing source records and all unused
// plane rows are zero-filled, independent of previous output contents.
// For a natural full-length message, pass rows=k/128 (k divisible by512).
// TileRows=64 keeps a128KiB input row tile; Affine selects instruction form
// only, not the map. Stream=true uses NT stores for aligned full row groups,
// retaining cached stores for shortened input and guard rows. Unaligned
// output falls back to the unchanged cached path. The streaming path fences
// once before return. There is no allocation in the encoder.
template<std::size_t TileRows=64,bool Affine=true,bool Stream=false>
HDBWO_NOINLINE void outerForward256(std::size_t k,std::size_t rows,
                                   const Block* input,Block* planes) {
    static_assert(TileRows>0&&TileRows%4==0);
    const auto fullRows=(k/512)*4;
    const auto remaining=k%512;
    if(rows%4||fullRows>rows||(remaining&&rows-fullRows<4))
        throw std::invalid_argument("GRS256/128 outer requires enough four-row plane groups");
    if constexpr(Stream){
        if(reinterpret_cast<std::uintptr_t>(planes)&63U){
            outerForward256<TileRows,Affine,false>(k,rows,input,planes);
            return;
        }
    }
    for(std::size_t first=0;first<fullRows;first+=TileRows){
        const auto last=std::min(first+TileRows,fullRows);
        for(auto row=first;row!=last;row+=4)
            detail::coset128<0,Affine,Stream>(input+row*128,planes+row,rows);
        for(auto row=first;row!=last;row+=4)
            detail::coset128<128,Affine,Stream>(input+row*128,planes+128*rows+row,rows);
    }
    auto assignedRows=fullRows;
    if(remaining){
        detail::tail<Affine>(input+fullRows*128,remaining,planes+fullRows,rows);
        assignedRows+=4;
    }
    if(assignedRows!=rows){
        for(unsigned region=0;region!=256;++region)
            std::memset(planes+region*rows+assignedRows,0,(rows-assignedRows)*sizeof(Block));
    }
    if constexpr(Stream)_mm_sfence();
}

#undef HDBWO_INLINE
#undef HDBWO_NOINLINE
} // namespace spin::research::char2::hdbwouter
