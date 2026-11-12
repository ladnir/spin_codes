// AVX-512 replica of the exact certified-s22 inner setup. It consumes the same
// SplitMix word stream as Setup.cpp's scalar rejection sampler and accepts the
// same candidates, so every table is bit-identical (checked by spin_setup_test).
// Candidate c is words 22c..22c+21 of the stream; epoch e takes the e-th
// full-rank candidate. Sixteen candidates are rank-tested per vector, one per
// lane, by incremental reduction against the rows already accepted.
#include "Plan.h"
#include "productionMap.h"
#include <immintrin.h>

namespace spin::detail::certified22 {
namespace {
constexpr std::uint64_t step=0x9e3779b97f4a7c15ULL;

inline __m512i mix(__m512i v) {
    v=_mm512_mullo_epi64(_mm512_xor_si512(v,_mm512_srli_epi64(v,30)),_mm512_set1_epi64(0xbf58476d1ce4e5b9LL));
    v=_mm512_mullo_epi64(_mm512_xor_si512(v,_mm512_srli_epi64(v,27)),_mm512_set1_epi64(0x94d049bb133111ebLL));
    return _mm512_xor_si512(v,_mm512_srli_epi64(v,31));
}

struct Batch {
    alignas(64) std::uint32_t rows[22][16];
    __mmask16 alive;
};

// Rank-tests candidates first..first+16*NB-1. Each new row is reduced by the
// accepted rows in insertion order; row i was reduced by rows 0..i-1 and so is
// zero at their pivots, which keeps every earlier pivot cleared. A lane is
// singular exactly when some row reduces to zero.
template<unsigned NB>
void rankBatches(std::uint64_t seed,std::uint64_t first,Batch* out) {
    __m512i x[NB][2],basis[NB][22],pivot[NB][22];
    __mmask16 alive[NB];
    const auto spread=_mm512_mullo_epi64(_mm512_setr_epi64(0,1,2,3,4,5,6,7),_mm512_set1_epi64(std::int64_t(22*step)));
    for(unsigned b=0; b<NB; ++b) {
        const auto start=seed+step*(22*(first+16*b)+1);
        x[b][0]=_mm512_add_epi64(_mm512_set1_epi64(std::int64_t(start)),spread);
        x[b][1]=_mm512_add_epi64(x[b][0],_mm512_set1_epi64(std::int64_t(8*22*step)));
        alive[b]=0xffff;
    }
    const auto low22=_mm512_set1_epi32(0x3fffff);
    const auto advance=_mm512_set1_epi64(std::int64_t(step));
    for(unsigned r=0; r<22; ++r) {
        __m512i v[NB];
        for(unsigned b=0; b<NB; ++b) {
            const auto lo=_mm512_cvtepi64_epi32(mix(x[b][0])),hi=_mm512_cvtepi64_epi32(mix(x[b][1]));
            x[b][0]=_mm512_add_epi64(x[b][0],advance);
            x[b][1]=_mm512_add_epi64(x[b][1],advance);
            v[b]=_mm512_and_si512(_mm512_inserti64x4(_mm512_castsi256_si512(lo),hi,1),low22);
            _mm512_store_si512(out[b].rows[r],v[b]);
        }
        for(unsigned i=0; i<r; ++i)
            for(unsigned b=0; b<NB; ++b)
                v[b]=_mm512_mask_xor_epi32(v[b],_mm512_test_epi32_mask(v[b],pivot[b][i]),v[b],basis[b][i]);
        for(unsigned b=0; b<NB; ++b) {
            alive[b]&=_mm512_test_epi32_mask(v[b],v[b]);
            basis[b][r]=v[b];
            pivot[b][r]=_mm512_and_si512(v[b],_mm512_sub_epi32(_mm512_setzero_si512(),v[b]));
        }
    }
    for(unsigned b=0; b<NB; ++b) out[b].alive=alive[b];
}

// Row masks of a constant left factor: bit r of lo[b]/hi[b] says row r (or
// 16+r) of the matrix has column b.
struct LeftMasks { __mmask16 lo[22],hi[22]; };
constexpr LeftMasks leftMasks(const std::array<std::uint32_t,22>& m) {
    LeftMasks result{};
    for(unsigned b=0; b<22; ++b)
        for(unsigned r=0; r<22; ++r)
            if((m[r]>>b)&1U) {
                if(r<16) result.lo[b]=__mmask16(result.lo[b]|(1U<<r));
                else result.hi[b]=__mmask16(result.hi[b]|(1U<<(r-16)));
            }
    return result;
}
constexpr auto transposeOf(const std::array<std::uint32_t,22>& m) {
    std::array<std::uint32_t,22> result{};
    for(unsigned r=0; r<22; ++r)
        for(unsigned c=0; c<22; ++c) result[c]|=((m[r]>>c)&1U)<<r;
    return result;
}
constexpr auto basisMasks=leftMasks(map::basis);
constexpr auto inverseTransposeMasks=leftMasks(transposeOf(map::inverseBasis));

// Rows 0..15 in lo, 16..21 in hi lanes 0..5 (other lanes zero).
struct Rows { __m512i lo,hi; };
inline Rows leftConstant(const LeftMasks& m,const std::uint32_t* right) {
    auto lo=_mm512_setzero_si512(),hi=_mm512_setzero_si512();
    for(unsigned b=0; b<22; ++b) {
        const auto row=_mm512_set1_epi32(int(right[b]));
        lo=_mm512_mask_xor_epi32(lo,m.lo[b],lo,row);
        hi=_mm512_mask_xor_epi32(hi,m.hi[b],hi,row);
    }
    return {lo,hi};
}
inline Rows rightConstant(Rows x,const std::array<std::uint32_t,22>& right) {
    auto lo=_mm512_setzero_si512(),hi=_mm512_setzero_si512();
    for(unsigned b=0; b<22; ++b) {
        const auto bit=_mm512_set1_epi32(int(1U<<b)),row=_mm512_set1_epi32(int(right[b]));
        lo=_mm512_mask_xor_epi32(lo,_mm512_test_epi32_mask(x.lo,bit),lo,row);
        hi=_mm512_mask_xor_epi32(hi,_mm512_test_epi32_mask(x.hi,bit),hi,row);
    }
    return {lo,hi};
}
inline void transposeRows(Rows x,std::uint32_t* out) {
    for(unsigned c=0; c<22; ++c) {
        const auto bit=_mm512_set1_epi32(int(1U<<c));
        out[c]=std::uint32_t(_mm512_test_epi32_mask(x.lo,bit))|(std::uint32_t(_mm512_test_epi32_mask(x.hi,bit))<<16);
    }
}
// Setup.cpp's packBlocks: qword 3o+i, byte j is byte i of row 8o+j.
constexpr auto packIndex=[] {
    std::array<std::uint8_t,64> result{};
    for(unsigned q=0; q<8; ++q)
        for(unsigned j=0; j<8; ++j) result[8*q+j]=std::uint8_t(4*(8*(q/3)+j)+q%3);
    return result;
}();
constexpr auto packLastIndex=[] {
    std::array<std::uint8_t,64> result{};
    for(unsigned j=0; j<8; ++j) result[j]=std::uint8_t(4*j+2);
    return result;
}();
inline void pack(Rows x,std::uint64_t* out) {
    const auto first=_mm512_loadu_si512(packIndex.data()),last=_mm512_loadu_si512(packLastIndex.data());
    _mm512_storeu_si512(out,_mm512_permutex2var_epi8(x.lo,first,x.hi));
    out[8]=std::uint64_t(_mm_cvtsi128_si64(_mm512_castsi512_si128(_mm512_permutexvar_epi8(last,x.hi))));
}

void transformEpoch(const std::uint32_t* literal,std::uint64_t* reverse,std::uint64_t* forward) {
    // Reverse: basis * R * inverseBasis.
    pack(rightConstant(leftConstant(basisMasks,literal),map::inverseBasis),reverse);
    // Forward: basis * R^T * inverseBasis = basis * (inverseBasis^T * R)^T.
    alignas(64) std::uint32_t transposed[22];
    transposeRows(leftConstant(inverseTransposeMasks,literal),transposed);
    pack(leftConstant(basisMasks,transposed),forward);
}
}

void fillWords(std::uint64_t seed,std::uint64_t first,std::uint64_t* out,std::size_t count) noexcept {
    auto x=_mm512_add_epi64(_mm512_set1_epi64(std::int64_t(seed+step*(first+1))),
        _mm512_mullo_epi64(_mm512_setr_epi64(0,1,2,3,4,5,6,7),_mm512_set1_epi64(std::int64_t(step))));
    const auto advance=_mm512_set1_epi64(std::int64_t(8*step));
    for(std::size_t i=0; i<count; i+=8) {
        _mm512_storeu_si512(out+i,mix(x));
        x=_mm512_add_epi64(x,advance);
    }
}

void prepareInnerFast(Plan& p,std::uint64_t seed) {
    const auto epochs=p.n/64;
    p.reverseMatrices.resize(epochs);
    p.updates.resize(9*epochs);
    p.forwardUpdates.resize(9*epochs);
    // Setup.cpp draws from Words(seed^0x3f625a92); word t is mix(s+step*(t+1)).
    const auto streamSeed=seed^0x3f625a92ULL;
#ifndef SPIN_SETUP_RANK_BATCHES
#define SPIN_SETUP_RANK_BATCHES 4
#endif
    constexpr unsigned NB=SPIN_SETUP_RANK_BATCHES;
    Batch batches[NB];
    std::size_t epoch=0;
    for(std::uint64_t first=0; epoch<epochs; first+=16*NB) {
        rankBatches<NB>(streamSeed,first,batches);
        for(unsigned b=0; b<NB && epoch<epochs; ++b)
            for(unsigned lane=0; lane<16 && epoch<epochs; ++lane) {
                if(!((batches[b].alive>>lane)&1U)) continue;
                auto& literal=p.reverseMatrices[epoch];
                for(unsigned r=0; r<22; ++r) literal[r]=batches[b].rows[r][lane];
                transformEpoch(literal.data(),p.updates.data()+9*epoch,p.forwardUpdates.data()+9*epoch);
                ++epoch;
            }
    }
}
}
