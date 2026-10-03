#pragma once
#include "../k16_codesign_100us/kernel/K16NativeOuter.h"
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace spin::research::packet8 {
using Block=rs::Block;
struct Update {
    // Literal forward GF256 matrix in row-major order.
    std::array<std::uint8_t,4> forward{};
    // Binary adjoint of each scalar product, encoded for GFNI affine.
    std::array<std::uint64_t,4> adjoint{};
};
struct Plan {
    explicit Plan(std::size_t k,std::uint64_t seed);
    rs::Plan outer;
    k16codesign::NativeOuterTables native;
    // Sequential inner byte packet -> first raw coordinate in padded outer group.
    // This is a fresh 32-region route, not paired entries of a 64-region route.
    std::vector<std::uint32_t> route;
    std::vector<Update> updates;
    std::size_t k() const noexcept { return outer.k; }
    std::size_t n() const noexcept { return outer.n; }
    std::size_t scratchBlocks() const noexcept { return outer.scratchBlocks(); }
};

constexpr std::uint8_t multiply(std::uint8_t a,std::uint8_t b) {
    unsigned x=a,y=b,result=0;
    while(y) { if(y&1)result^=x; y>>=1; x=(x<<1)^((x&128)?0x11b:0); }
    return static_cast<std::uint8_t>(result);
}
constexpr std::uint64_t adjointMatrix(std::uint8_t scalar) {
    std::uint64_t result=0;
    for(unsigned i=0;i<8;++i)
        result|=std::uint64_t(multiply(scalar,std::uint8_t(1U<<i)))<<(8*(7-i));
    return result;
}

// Input/output: ordinary arrays of 128-bit XOR elements, only 16-byte aligned.
// Scratch: 64-byte aligned, padded groupStride; each eight-coordinate packet
// is two adjacent 64-byte bitplane vectors. It is disjoint from input/output.
// Hot functions allocate nothing. Output may overwrite the input prefix.
void reverseRouteFast(const Block*,Block* packedScratch,const Plan&);
void outerFast(const Block* packedScratch,Block*,const Plan&);
void outerFastShared(const Block* packedScratch,Block*,const Plan&);
void outerFastHalf(const Block* packedScratch,Block*,const Plan&);
void transposeFast(const Block*,Block*,Block* packedScratch,const Plan&);
void routeOnlyFast(const Block*,Block* packedScratch,const Plan&);
void reverseRoutePrefetch8(const Block*,Block*,const Plan&);
void reverseRoutePrefetch16(const Block*,Block*,const Plan&);
void reverseRoutePrefetch32(const Block*,Block*,const Plan&);
void routeOnlyPrefetch8(const Block*,Block*,const Plan&);
void routeOnlyPrefetch16(const Block*,Block*,const Plan&);
void routeOnlyPrefetch32(const Block*,Block*,const Plan&);

// Independent literal binary-coordinate oracles. Their scratch is RAW.
void reverseRouteScalar(const Block*,Block* rawScratch,const Plan&);
void transposeScalar(const Block*,Block*,Block* rawScratch,const Plan&);
void forwardScalar(const Block*,Block*,const Plan&);
void packScratchScalar(const Block* rawScratch,Block* packedScratch,const Plan&);
void unpackScratchScalar(const Block* packedScratch,Block* rawScratch,const Plan&);
}
