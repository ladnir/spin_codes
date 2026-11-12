// Constant byte maps of the RS8 outer as nibble lookup tables (two 16-entry
// tables per GF(2)-linear byte map f: f(x) = lo[x&15] ^ hi[x>>4]), shared by
// the byte-plane kernels (Avx2.cpp, Neon.cpp): multiplication by each GF(256)
// constant, its adjoint, and the RS parity network's fixed nibble maps.
#pragma once
#include <bit>
#include <cstdint>

namespace spin::detail::rs8::outer_maps {
inline constexpr std::uint8_t mul8(std::uint8_t a,std::uint8_t b) {
    std::uint8_t result=0;
    for(unsigned i=0;i<8;++i) {
        if(b&1U) result^=a;
        a=std::uint8_t((unsigned(a)<<1)^((a&128U)?0x11bU:0U));
        b>>=1;
    }
    return result;
}
// GFNI affine semantics: bit i of the result is parity(byte 7-i of m & x).
inline constexpr std::uint8_t affine(std::uint64_t m,std::uint8_t x) {
    std::uint8_t result=0;
    for(unsigned i=0;i<8;++i)
        result|=std::uint8_t((std::popcount(unsigned((m>>(8*(7-i)))&x)&255U)&1U)<<i);
    return result;
}
inline constexpr std::uint64_t transposeAffine(std::uint64_t m) {
    std::uint64_t result=0;
    for(unsigned i=0;i<8;++i) for(unsigned j=0;j<8;++j)
        result|=((m>>(8*(7-i)+j))&1ULL)<<(8*(7-j)+i);
    return result;
}
// A GF(2)-linear byte map f as two nibble tables: f(x)=lo[x&15]^hi[x>>4].
struct alignas(32) NibbleMap { std::uint8_t lo[16]{},hi[16]{}; };
template<class F> constexpr NibbleMap nibbleMap(F f) {
    NibbleMap result{};
    for(unsigned v=0;v<16;++v) {
        result.lo[v]=f(std::uint8_t(v));
        result.hi[v]=f(std::uint8_t(v<<4));
    }
    return result;
}
struct Maps { NibbleMap multiply[256],adjoint[256]; };
// Built from the eight columns c*e_b of each multiplication (cheap enough for
// default constexpr step limits).
inline constexpr Maps makeMaps() {
    Maps maps{};
    for(unsigned c=0;c<256;++c) {
        std::uint8_t column[8]{};
        for(unsigned b=0;b<8;++b) column[b]=mul8(std::uint8_t(c),std::uint8_t(1U<<b));
        const auto times=[&](unsigned x) {
            std::uint8_t r=0;
            for(unsigned b=0;b<8;++b) if((x>>b)&1U) r^=column[b];
            return r;
        };
        // Adjoint of multiplication by c: bit b is parity(c*e_b & y).
        const auto adjoint=[&](unsigned y) {
            std::uint8_t r=0;
            for(unsigned b=0;b<8;++b) r|=std::uint8_t((std::popcount(unsigned(column[b]&y))&1U)<<b);
            return r;
        };
        for(unsigned v=0;v<16;++v) {
            maps.multiply[c].lo[v]=times(v);maps.multiply[c].hi[v]=times(v<<4);
            maps.adjoint[c].lo[v]=adjoint(v);maps.adjoint[c].hi[v]=adjoint(v<<4);
        }
    }
    return maps;
}
alignas(64) inline constexpr Maps maps=makeMaps();
// The parity network's adjoint nibble multiplications (rs8/Fast.cpp) and
// their transposes (rs8/ForwardFast.cpp).
inline constexpr std::uint64_t m6=0x060c0b0560c0b050ULL,m8=0x0803060c803060c0ULL;
struct Affine {
    std::uint64_t m;
    constexpr std::uint8_t operator()(std::uint8_t x) const { return affine(m,x); }
};
alignas(32) inline constexpr NibbleMap a6=nibbleMap(Affine{m6});
alignas(32) inline constexpr NibbleMap a8=nibbleMap(Affine{m8});
alignas(32) inline constexpr NibbleMap a6t=nibbleMap(Affine{transposeAffine(m6)});
alignas(32) inline constexpr NibbleMap a8t=nibbleMap(Affine{transposeAffine(m8)});
// The same maps as GFNI affine matrices (byte 7-i of the matrix is the mask
// of input bits that make output bit i): multiplication by c and its adjoint.
struct AffineMaps { std::uint64_t multiply[256],adjoint[256]; };
inline constexpr AffineMaps makeAffineMaps() {
    AffineMaps result{};
    for(unsigned c=0;c<256;++c)
        for(unsigned b=0;b<8;++b) {
            const auto column=mul8(std::uint8_t(c),std::uint8_t(1U<<b));
            // Adjoint: output bit b is parity(c*e_b & y).
            result.adjoint[c]|=std::uint64_t(column)<<(8*(7-b));
            // Multiply: output bit i takes input bit b when bit i of c*e_b is set.
            for(unsigned i=0;i<8;++i)
                if((column>>i)&1U) result.multiply[c]|=std::uint64_t(1U<<b)<<(8*(7-i));
        }
    return result;
}
alignas(64) inline constexpr AffineMaps affineMaps=makeAffineMaps();
}
