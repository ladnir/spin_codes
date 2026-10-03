#pragma once
#include "../../src/kernels/Block.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::research::packet8wide24 {
using Block=detail::storage::block;
inline constexpr std::size_t groupStride=516;
struct Update {
    // Forward update is the binary adjoint of multiplication by this scalar
    // in GF256[z]/(z^3+z+1). Reverse evaluation uses ordinary multiplication.
    std::uint32_t scalar=0;
    std::array<std::uint8_t,6> coefficients{}; // r0,r1,r2,r0^r1,r0^r2,r1^r2
};
struct Plan {
    Plan(std::size_t messageBits,std::uint64_t seed);
    std::size_t k=0,n=0,groups=0;
    std::vector<std::uint32_t> route;
    std::vector<Update> updates;
    std::vector<std::array<std::uint8_t,9>> outerCoefficients;
    // Reverse symbol rows: logical RS coordinate c receives the indicated
    // ordinary, contiguous physical symbol input bits. No nibble repacking.
    std::vector<std::array<std::uint32_t,32>> outerRows;
    std::size_t scratchBlocks() const noexcept { return groups*groupStride; }
};
constexpr std::uint8_t multiply8(std::uint8_t a,std::uint8_t b) {
    unsigned x=a,y=b,result=0;
    while(y){if(y&1)result^=x;y>>=1;x=(x<<1)^((x&128)?0x11b:0);}
    return std::uint8_t(result);
}
constexpr std::uint64_t adjointMatrix(std::uint8_t scalar) {
    std::uint64_t result=0;
    for(unsigned i=0;i<8;++i)
        result|=std::uint64_t(multiply8(scalar,std::uint8_t(1U<<i)))<<(8*(7-i));
    return result;
}
constexpr std::uint32_t multiply24(std::uint32_t x,std::uint32_t y) {
    std::uint8_t c[5]{};
    for(unsigned i=0;i<3;++i)for(unsigned j=0;j<3;++j)
        c[i+j]^=multiply8(std::uint8_t(x>>(8*i)),std::uint8_t(y>>(8*j)));
    for(unsigned degree=4;degree>=3;--degree){c[degree-3]^=c[degree];c[degree-2]^=c[degree];}
    return std::uint32_t(c[0])|(std::uint32_t(c[1])<<8)|(std::uint32_t(c[2])<<16);
}

// Input/output are ordinary 128-bit XOR elements, with 16-byte alignment.
// Scratch is a disjoint 64-byte-aligned array of scratchBlocks() records.
// Each contiguous eight-coordinate packet occupies two packed ZMM vectors.
// Hot functions allocate nothing. Full transpose allows output==input.
void reverseRouteFast(const Block*,Block*,const Plan&);
void outerFast(const Block*,Block*,const Plan&);
void transposeFast(const Block*,Block*,Block*,const Plan&);
void routeOnlyFast(const Block*,Block*,const Plan&);

// Independent literal binary-coordinate validation oracles, not timed.
void reverseRouteScalar(const Block*,Block*,const Plan&);
void outerScalar(const Block*,Block*,const Plan&);
void transposeScalar(const Block*,Block*,Block*,const Plan&);
void forwardScalar(const Block*,Block*,const Plan&);
void packScratchScalar(const Block*,Block*,const Plan&);
void unpackScratchScalar(const Block*,Block*,const Plan&);
}
