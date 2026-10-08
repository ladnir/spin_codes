#pragma once
#include "HighDistanceWide.h"

namespace spin::research::char2::hdfast {

// Forward evaluation in the unchanged normalized additive-polynomial basis.
// Four rows share a ZMM; the GF256 coefficients act independently on every
// byte of a GF16/GF32 record. This is not the systematic Cauchy256 transform.
#if defined(_MSC_VER)
#define HDFO_INLINE __forceinline
#define HDFO_NOINLINE __declspec(noinline)
#else
#define HDFO_INLINE inline __attribute__((always_inline))
#define HDFO_NOINLINE __attribute__((noinline))
#endif

namespace outer_detail {
using rs4_outer_next_detail::normalized;

template<bool Affine,unsigned Coefficient>
HDFO_INLINE void pair(__m512i& low,__m512i& high) noexcept {
    if constexpr(Coefficient==1)low=_mm512_xor_si512(low,high);
    else if constexpr(Coefficient!=0)
        low=_mm512_xor_si512(low,rs8ConstantMultiply<Affine,Coefficient>(high));
    high=_mm512_xor_si512(high,low);
}

// A leaf keeps exactly sixteen values live and emits final plane columns.
// Its order is the reversed circuit order of hdwide::coset64, with forward
// butterflies rather than transposed butterflies. No leaf-output scratch.
template<unsigned Offset,bool Affine>
HDFO_NOINLINE void leaf16(const __m512i* source,Block* planes,std::size_t stride) {
#define HDFO_LOAD(I) auto x##I=_mm512_load_si512(source+I)
    HDFO_LOAD(0);HDFO_LOAD(1);HDFO_LOAD(2);HDFO_LOAD(3);
    HDFO_LOAD(4);HDFO_LOAD(5);HDFO_LOAD(6);HDFO_LOAD(7);
    HDFO_LOAD(8);HDFO_LOAD(9);HDFO_LOAD(10);HDFO_LOAD(11);
    HDFO_LOAD(12);HDFO_LOAD(13);HDFO_LOAD(14);HDFO_LOAD(15);
#undef HDFO_LOAD
#define HDFO_PAIR(A,B,LEVEL,DELTA) pair<Affine,normalized(LEVEL,Offset^DELTA)>(x##A,x##B)
    HDFO_PAIR(0,8,3,0);HDFO_PAIR(1,9,3,0);
    HDFO_PAIR(2,10,3,0);HDFO_PAIR(3,11,3,0);
    HDFO_PAIR(4,12,3,0);HDFO_PAIR(5,13,3,0);
    HDFO_PAIR(6,14,3,0);HDFO_PAIR(7,15,3,0);
    HDFO_PAIR(0,4,2,0);HDFO_PAIR(1,5,2,0);
    HDFO_PAIR(2,6,2,0);HDFO_PAIR(3,7,2,0);
    HDFO_PAIR(8,12,2,8);HDFO_PAIR(9,13,2,8);
    HDFO_PAIR(10,14,2,8);HDFO_PAIR(11,15,2,8);
    HDFO_PAIR(0,2,1,0);HDFO_PAIR(1,3,1,0);
    HDFO_PAIR(4,6,1,4);HDFO_PAIR(5,7,1,4);
    HDFO_PAIR(8,10,1,8);HDFO_PAIR(9,11,1,8);
    HDFO_PAIR(12,14,1,12);HDFO_PAIR(13,15,1,12);
    HDFO_PAIR(0,1,0,0);HDFO_PAIR(2,3,0,2);
    HDFO_PAIR(4,5,0,4);HDFO_PAIR(6,7,0,6);
    HDFO_PAIR(8,9,0,8);HDFO_PAIR(10,11,0,10);
    HDFO_PAIR(12,13,0,12);HDFO_PAIR(14,15,0,14);
#undef HDFO_PAIR
#define HDFO_STORE(I) _mm512_storeu_si512(planes+I*stride,x##I)
    HDFO_STORE(0);HDFO_STORE(1);HDFO_STORE(2);HDFO_STORE(3);
    HDFO_STORE(4);HDFO_STORE(5);HDFO_STORE(6);HDFO_STORE(7);
    HDFO_STORE(8);HDFO_STORE(9);HDFO_STORE(10);HDFO_STORE(11);
    HDFO_STORE(12);HDFO_STORE(13);HDFO_STORE(14);HDFO_STORE(15);
#undef HDFO_STORE
}

// Fuse both ancestor levels before their only materialization. Four original
// coefficient vectors give the corresponding inputs of all four leaves.
// The bounded quartet schedule avoids trying to hold 64 values in registers.
template<unsigned Offset,bool Affine>
HDFO_NOINLINE void coset64(const Block* source,Block* planes,std::size_t stride) {
    static_assert(Offset==0||Offset==64||Offset==128||Offset==192);
    alignas(64) __m512i values[64];
    for(unsigned base=0;base!=16;base+=4){
        c8Unroll(std::make_index_sequence<4>{},[&](auto index){
            const auto j=base+unsigned(index);
            auto a=_mm512_loadu_si512(source+4*j);
            auto b=_mm512_loadu_si512(source+4*(16+j));
            auto c=_mm512_loadu_si512(source+4*(32+j));
            auto d=_mm512_loadu_si512(source+4*(48+j));
            pair<Affine,normalized(5,Offset)>(a,c);
            pair<Affine,normalized(5,Offset)>(b,d);
            pair<Affine,normalized(4,Offset)>(a,b);
            pair<Affine,normalized(4,Offset+32)>(c,d);
            _mm512_store_si512(values+j,a);
            _mm512_store_si512(values+16+j,b);
            _mm512_store_si512(values+32+j,c);
            _mm512_store_si512(values+48+j,d);
        });
    }
    leaf16<Offset,Affine>(values,planes,stride);
    leaf16<Offset+16,Affine>(values+16,planes+16*stride,stride);
    leaf16<Offset+32,Affine>(values+32,planes+32*stride,stride);
    leaf16<Offset+48,Affine>(values+48,planes+48*stride,stride);
}

template<bool Affine>
HDFO_INLINE void tile(const Block* source,Block* planes,std::size_t stride) {
    // Degree <64 makes the two omitted high FFT levels pure replication.
    // Evaluating four 64-point cosets therefore retains the exact 256-point
    // evaluation map without ever materializing 256 coefficient vectors.
    coset64<0,Affine>(source,planes,stride);
    coset64<64,Affine>(source,planes+64*stride,stride);
    coset64<128,Affine>(source,planes+128*stride,stride);
    coset64<192,Affine>(source,planes+192*stride,stride);
}

// Only the single shortened source tile needs a copy and zero fill. Keep its
// aligned 4KiB scratch outside the hot outer caller's stack frame.
template<bool Affine>
HDFO_NOINLINE void tail(const Block* source,std::size_t count,
                       Block* planes,std::size_t stride) {
    alignas(64) Block padded[256]{};
    std::memcpy(padded,source,count*sizeof(Block));
    tile<Affine>(padded,planes,stride);
}
} // namespace outer_detail

// input[requested] and planes[256*planeStride] must be disjoint. Neither needs
// SIMD alignment. Input indexing is Plan::sourceIndex; each padded plane is
// contiguous in row order. Both shortened records and the four guard rows
// are assigned zero, so prior contents of planes are immaterial.
// Affine=false is an exact-map instruction-selection control, not a new map.
template<bool Affine=true,unsigned Bits>
HDFO_NOINLINE void outerForward(const hdwide::Plan<Bits,2,256>& p,
                                const Block* input,Block* planes) {
    const auto fullRows=(p.requested/256)*4;
    for(std::size_t row=0;row!=fullRows;row+=4)
        outer_detail::tile<Affine>(input+row*64,planes+row,p.planeStride);
    auto assignedRows=fullRows;
    const auto remaining=p.requested-fullRows*64;
    if(remaining){
        outer_detail::tail<Affine>(input+fullRows*64,remaining,
                                  planes+fullRows,p.planeStride);
        assignedRows+=4;
    }
    for(unsigned region=0;region!=256;++region)
        std::memset(planes+region*p.planeStride+assignedRows,0,
                    (p.planeStride-assignedRows)*sizeof(Block));
}

// Exact-map loop-interchange control: keep a bounded source-row tile hot,
// while completing its first 64 planes before moving to the next coset.
// TileRows=64 uses a 64KiB input tile; TileRows=256 uses 256KiB. The per-coset
// kernels and their register/scratch schedule are unchanged from outerForward.
template<std::size_t TileRows=64,bool Affine=true,unsigned Bits>
HDFO_NOINLINE void outerForwardTiled(const hdwide::Plan<Bits,2,256>& p,
                                     const Block* input,Block* planes) {
    static_assert(TileRows>0&&TileRows%4==0);
    const auto fullRows=(p.requested/256)*4;
    for(std::size_t first=0;first<fullRows;first+=TileRows){
        const auto last=std::min(first+TileRows,fullRows);
        for(auto row=first;row!=last;row+=4)
            outer_detail::coset64<0,Affine>(input+row*64,planes+row,p.planeStride);
        for(auto row=first;row!=last;row+=4)
            outer_detail::coset64<64,Affine>(input+row*64,
                planes+64*p.planeStride+row,p.planeStride);
        for(auto row=first;row!=last;row+=4)
            outer_detail::coset64<128,Affine>(input+row*64,
                planes+128*p.planeStride+row,p.planeStride);
        for(auto row=first;row!=last;row+=4)
            outer_detail::coset64<192,Affine>(input+row*64,
                planes+192*p.planeStride+row,p.planeStride);
    }
    auto assignedRows=fullRows;
    const auto remaining=p.requested-fullRows*64;
    if(remaining){
        outer_detail::tail<Affine>(input+fullRows*64,remaining,
                                  planes+fullRows,p.planeStride);
        assignedRows+=4;
    }
    for(unsigned region=0;region!=256;++region)
        std::memset(planes+region*p.planeStride+assignedRows,0,
                    (p.planeStride-assignedRows)*sizeof(Block));
}

#undef HDFO_INLINE
#undef HDFO_NOINLINE
} // namespace spin::research::char2::hdfast
