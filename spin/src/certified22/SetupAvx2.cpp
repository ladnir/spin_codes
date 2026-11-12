// AVX2 replica of the exact certified-s22 inner setup (Setup.cpp prepareInner),
// for CPUs without the AVX-512 packet ISA. It consumes the same SplitMix word
// stream and accepts the same candidates, so every table is bit-identical
// (checked by spin_turbo_setup_test). Candidate c is words 22c..22c+21 of the
// stream; epoch e takes the e-th full-rank candidate. Eight candidates are
// rank-tested per vector, one per 32-bit lane. The basis transforms are those
// of SetupTransform.h, with vector packing.
#include "Plan.h"
#include "SetupTransform.h"
#include "../kernels/Block.h"
#include <immintrin.h>
#include <array>
#include <utility>

namespace spin::detail::certified22 {
namespace {
using namespace transform;

// Low 64 bits of a*c for a constant c, from three 32x32 multiplies.
static SPIN_FORCEINLINE __m256i multiply(__m256i a,__m256i c,__m256i cHigh) {
    const auto low=_mm256_mul_epu32(a,c);
    const auto cross=_mm256_add_epi64(_mm256_mul_epu32(_mm256_srli_epi64(a,32),c),_mm256_mul_epu32(a,cHigh));
    return _mm256_add_epi64(low,_mm256_slli_epi64(cross,32));
}
static SPIN_FORCEINLINE __m256i mix(__m256i v) {
    const auto c1=_mm256_set1_epi64x(std::int64_t(0xbf58476d1ce4e5b9ULL)),h1=_mm256_srli_epi64(c1,32);
    const auto c2=_mm256_set1_epi64x(std::int64_t(0x94d049bb133111ebULL)),h2=_mm256_srli_epi64(c2,32);
    v=multiply(_mm256_xor_si256(v,_mm256_srli_epi64(v,30)),c1,h1);
    v=multiply(_mm256_xor_si256(v,_mm256_srli_epi64(v,27)),c2,h2);
    return _mm256_xor_si256(v,_mm256_srli_epi64(v,31));
}

struct Batch {
    alignas(32) std::uint32_t rows[22][8];
    unsigned alive;
};

// Rank-tests candidates first..first+8*NB-1, exactly as SetupFast.cpp: each
// new row is reduced by the accepted rows in insertion order (pivot = lowest
// set bit); a lane is singular exactly when some row reduces to zero.
template<unsigned NB>
static void rankBatches(std::uint64_t seed,std::uint64_t first,Batch* out) {
    __m256i x[NB][2],basis[NB][22],pivot[NB][22],alive[NB];
    const auto spread=_mm256_setr_epi64x(0,std::int64_t(22*step),std::int64_t(2*22*step),std::int64_t(3*22*step));
    for(unsigned b=0;b<NB;++b) {
        const auto start=seed+step*(22*(first+8*b)+1);
        x[b][0]=_mm256_add_epi64(_mm256_set1_epi64x(std::int64_t(start)),spread);
        x[b][1]=_mm256_add_epi64(x[b][0],_mm256_set1_epi64x(std::int64_t(4*22*step)));
        alive[b]=_mm256_set1_epi32(-1);
    }
    const auto low22=_mm256_set1_epi32(0x3fffff);
    const auto advance=_mm256_set1_epi64x(std::int64_t(step));
    const auto even=_mm256_setr_epi32(0,2,4,6,1,3,5,7);
    const auto zero=_mm256_setzero_si256();
    for(unsigned r=0;r<22;++r) {
        __m256i v[NB];
        for(unsigned b=0;b<NB;++b) {
            const auto lo=_mm256_permutevar8x32_epi32(mix(x[b][0]),even);
            const auto hi=_mm256_permutevar8x32_epi32(mix(x[b][1]),even);
            x[b][0]=_mm256_add_epi64(x[b][0],advance);
            x[b][1]=_mm256_add_epi64(x[b][1],advance);
            v[b]=_mm256_and_si256(_mm256_permute2x128_si256(lo,hi,0x20),low22);
            _mm256_store_si256(reinterpret_cast<__m256i*>(out[b].rows[r]),v[b]);
        }
        for(unsigned i=0;i<r;++i)
            for(unsigned b=0;b<NB;++b) {
                const auto clear=_mm256_cmpeq_epi32(_mm256_and_si256(v[b],pivot[b][i]),zero);
                v[b]=_mm256_xor_si256(v[b],_mm256_andnot_si256(clear,basis[b][i]));
            }
        for(unsigned b=0;b<NB;++b) {
            alive[b]=_mm256_andnot_si256(_mm256_cmpeq_epi32(v[b],zero),alive[b]);
            basis[b][r]=v[b];
            pivot[b][r]=_mm256_and_si256(v[b],_mm256_sub_epi32(zero,v[b]));
        }
    }
    for(unsigned b=0;b<NB;++b) out[b].alive=unsigned(_mm256_movemask_ps(_mm256_castsi256_ps(alive[b])));
}

// Setup.cpp's packBlocks: qword 3o+i, byte j is byte i of row 8o+j. Rows of
// block o give their bytes 0..3 in qwords 0..3 of one vector (rows 22, 23
// are zero padding).
static SPIN_FORCEINLINE __m256i packRows(const std::uint32_t* rows) {
    const auto bytes=_mm256_setr_epi8(0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15,
                                      0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15);
    const auto x=_mm256_shuffle_epi8(_mm256_load_si256(reinterpret_cast<const __m256i*>(rows)),bytes);
    return _mm256_permutevar8x32_epi32(x,_mm256_setr_epi32(0,4,1,5,2,6,3,7));
}
// 8x8 bit transpose in every qword: byte j bit k <-> byte k bit j.
static SPIN_FORCEINLINE __m256i bitTranspose(__m256i x) {
    auto t=_mm256_and_si256(_mm256_xor_si256(x,_mm256_srli_epi64(x,7)),_mm256_set1_epi64x(0x00aa00aa00aa00aaLL));
    x=_mm256_xor_si256(x,_mm256_xor_si256(t,_mm256_slli_epi64(t,7)));
    t=_mm256_and_si256(_mm256_xor_si256(x,_mm256_srli_epi64(x,14)),_mm256_set1_epi64x(0x0000cccc0000ccccLL));
    x=_mm256_xor_si256(x,_mm256_xor_si256(t,_mm256_slli_epi64(t,14)));
    t=_mm256_and_si256(_mm256_xor_si256(x,_mm256_srli_epi64(x,28)),_mm256_set1_epi64x(0x00000000f0f0f0f0LL));
    return _mm256_xor_si256(x,_mm256_xor_si256(t,_mm256_slli_epi64(t,28)));
}

// Reverse: pack(basis * L * inverseBasis). Forward: F = basis * L^T *
// inverseBasis, whose transpose is G = inverseBasis^T * L * basis^T; block
// (o,i) of F is the bit transpose of block (i,o) of G.
static void transformEpoch(const std::uint32_t* literal,std::uint64_t* reverse,std::uint64_t* forward) {
    alignas(32) std::uint32_t u[24],g[24];
    transforms(literal,u,g);
    const auto three=_mm256_setr_epi64x(-1,-1,-1,0);
    alignas(32) std::uint64_t blocks[12];
    for(unsigned o=0;o<3;++o) {
        _mm256_maskstore_epi64(reinterpret_cast<long long*>(reverse+3*o),three,packRows(u+8*o));
        _mm256_store_si256(reinterpret_cast<__m256i*>(blocks+4*o),bitTranspose(packRows(g+8*o)));
    }
    for(unsigned o=0;o<3;++o)
        for(unsigned i=0;i<3;++i) forward[3*o+i]=blocks[4*i+o];
}
}

void prepareInnerAvx2(Plan& p,std::uint64_t seed) {
    const auto epochs=p.n/64;
    p.reverseMatrices.resize(epochs);
    p.updates.resize(9*epochs);
    p.forwardUpdates.resize(9*epochs);
    // Setup.cpp draws from Words(seed^0x3f625a92); word t is mix(s+step*(t+1)).
    const auto streamSeed=seed^0x3f625a92ULL;
    constexpr unsigned NB=4;
    Batch batches[NB];
    std::size_t epoch=0;
    for(std::uint64_t first=0;epoch<epochs;first+=8*NB) {
        rankBatches<NB>(streamSeed,first,batches);
        for(unsigned b=0;b<NB && epoch<epochs;++b)
            for(unsigned lane=0;lane<8 && epoch<epochs;++lane) {
                if(!((batches[b].alive>>lane)&1U)) continue;
                auto& literal=p.reverseMatrices[epoch];
                for(unsigned r=0;r<22;++r) literal[r]=batches[b].rows[r][lane];
                transformEpoch(literal.data(),p.updates.data()+9*epoch,p.forwardUpdates.data()+9*epoch);
                ++epoch;
            }
    }
}
}
