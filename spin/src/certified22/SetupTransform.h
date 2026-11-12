// Shared pieces of the fast exact inner setups (SetupAvx2.cpp,
// SetupPortable.cpp): the SplitMix word, and the per-epoch basis transforms
// of a sampled GL22 matrix L into the packed reverse table
// pack(basis * L * inverseBasis) and forward table pack(basis * L^T *
// inverseBasis), exactly as Setup.cpp's literal prepareInner computes them.
// Constant left factors are XOR patterns, constant right factors byte tables.
#pragma once
#include "productionMap.h"
#include "../kernels/Block.h"
#include <array>
#include <cstdint>
#include <utility>

namespace spin::detail::certified22::transform {
using Matrix=std::array<std::uint32_t,22>;
inline constexpr std::uint64_t step=0x9e3779b97f4a7c15ULL;
// Word t of kernel::setup::Words(seed) is mix(seed+step*(t+1)).
constexpr std::uint64_t mix(std::uint64_t z) {
    z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
    z=(z^(z>>27))*0x94d049bb133111ebULL;
    return z^(z>>31);
}
constexpr Matrix transposeOf(const Matrix& m) {
    Matrix result{};
    for(unsigned r=0;r<22;++r)
        for(unsigned c=0;c<22;++c) result[c]|=((m[r]>>c)&1U)<<r;
    return result;
}
inline constexpr Matrix inverseTranspose=transposeOf(map::inverseBasis);
inline constexpr Matrix basisTranspose=transposeOf(map::basis);

// (M*X)[r] = XOR of X[b] over the bits b of the constant row M[r].
template<std::uint32_t Mask,std::size_t... B>
SPIN_FORCEINLINE std::uint32_t combine(const std::uint32_t* x,std::index_sequence<B...>) {
    return (((Mask>>B)&1U?x[B]:0U)^...);
}
template<const Matrix& M,std::size_t... R>
SPIN_FORCEINLINE void leftConstant(const std::uint32_t* x,std::uint32_t* out,std::index_sequence<R...>) {
    ((out[R]=combine<M[R]>(x,std::make_index_sequence<22>{})),...);
}
// (X*C)[r] = XOR of C[c] over the bits c of X[r], three byte lookups.
struct RightTables { std::uint32_t t[3][256]; };
constexpr RightTables rightTables(const Matrix& c) {
    RightTables result{};
    for(unsigned i=0;i<3;++i)
        for(unsigned v=0;v<256;++v)
            for(unsigned b=0;b<8;++b)
                if(((v>>b)&1U) && 8*i+b<22) result.t[i][v]^=c[8*i+b];
    return result;
}
alignas(64) inline constexpr RightTables byInverse=rightTables(map::inverseBasis);
alignas(64) inline constexpr RightTables byBasisTranspose=rightTables(basisTranspose);
SPIN_FORCEINLINE std::uint32_t right(const RightTables& t,std::uint32_t x) {
    return t.t[0][x&255U]^t.t[1][(x>>8)&255U]^t.t[2][x>>16];
}
// u = basis*L*inverseBasis and g = (basis*L^T*inverseBasis)^T =
// inverseBasis^T*L*basis^T, as rows (24 entries; rows 22, 23 zero).
SPIN_FORCEINLINE void transforms(const std::uint32_t* literal,std::uint32_t* u,std::uint32_t* g) {
    std::uint32_t a[22];
    leftConstant<map::basis>(literal,a,std::make_index_sequence<22>{});
    for(unsigned r=0;r<22;++r) u[r]=right(byInverse,a[r]);
    leftConstant<inverseTranspose>(literal,a,std::make_index_sequence<22>{});
    for(unsigned r=0;r<22;++r) g[r]=right(byBasisTranspose,a[r]);
    u[22]=u[23]=g[22]=g[23]=0;
}
// Setup.cpp's packBlocks: qword 3o+i, byte j is byte i of row 8o+j.
SPIN_FORCEINLINE std::uint64_t packBlock(const std::uint32_t* rows,unsigned o,unsigned i) {
    std::uint64_t q=0;
    for(unsigned j=0;j<8;++j) q|=std::uint64_t((rows[8*o+j]>>(8*i))&255U)<<(8*j);
    return q;
}
// 8x8 bit transpose of a qword: byte j bit k <-> byte k bit j.
constexpr std::uint64_t bitTranspose(std::uint64_t x) {
    std::uint64_t t=(x^(x>>7))&0x00aa00aa00aa00aaULL; x^=t^(t<<7);
    t=(x^(x>>14))&0x0000cccc0000ccccULL; x^=t^(t<<14);
    t=(x^(x>>28))&0x00000000f0f0f0f0ULL; x^=t^(t<<28);
    return x;
}
// Portable packing of both tables. Forward block (o,i) is the bit transpose
// of block (i,o) of g.
SPIN_FORCEINLINE void pack(const std::uint32_t* u,const std::uint32_t* g,std::uint64_t* reverse,std::uint64_t* forward) {
    for(unsigned o=0;o<3;++o)
        for(unsigned i=0;i<3;++i) {
            reverse[3*o+i]=packBlock(u,o,i);
            forward[3*o+i]=bitTranspose(packBlock(g,i,o));
        }
}
}
