#pragma once
#include "Gfni.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstring>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace spin::research::char2 {
struct alignas(16) Block { std::uint64_t lo{}, hi{}; bool operator==(const Block&) const = default; };
inline Block bxor(Block a, Block b) { return {a.lo^b.lo,a.hi^b.hi}; }
inline __m128i reg(Block a) { return std::bit_cast<__m128i>(a); }
inline Block block(__m128i a) { return std::bit_cast<Block>(a); }
inline Block xtime(Block a) {
    return {(a.lo<<1)^((0-(a.hi>>63))&0x87), (a.hi<<1)|(a.lo>>63)};
}
inline Block mulSmallOracle(Block a, std::uint32_t c) {
    Block r{};
    while(c) { if(c&1) r=bxor(r,a); a=xtime(a); c>>=1; }
    return r;
}
inline Block mulOracle(Block a, Block b) {
    Block r{};
    for(unsigned i=0;i<128;++i) {
        if((i<64?b.lo>>(i):b.hi>>(i-64))&1) r=bxor(r,a);
        a=xtime(a);
    }
    return r;
}
constexpr std::uint32_t pmul(std::uint32_t a,std::uint32_t b) {
    std::uint32_t r=0;
    while(b) { if(b&1) r^=a; a<<=1; b>>=1; }
    return r;
}
// Unnormalized subspace polynomial S_j(X)=product_(v<2^j)(X+v).
// Butterfly constants fit in 16 bits; complete basis evaluations fit in
// 32 bits. These are polynomial coefficients in GF(2^128), not integers.
constexpr std::uint32_t subspace(unsigned j,unsigned x) {
    std::uint32_t r=1;
    for(unsigned v=0;v<(1U<<j);++v) r=pmul(r,x^v);
    return r;
}
constexpr std::uint32_t basis(unsigned j,unsigned x) {
    std::uint32_t r=1;
    for(unsigned d=0;j;++d,j>>=1) if(j&1)r=pmul(r,subspace(d,x));
    return r;
}

inline __m128i multiplyX8(__m128i x) {
    const auto overflow=_mm_srli_si128(x,15);
    const auto lo=_mm_gf2p8affine_epi64_epi8(overflow,
        _mm_set1_epi64x(static_cast<long long>(affineMatrix(0x87,false))),0);
    const auto hi=_mm_gf2p8affine_epi64_epi8(overflow,
        _mm_set1_epi64x(static_cast<long long>(affineMatrix(0x87,true))),0);
    return _mm_xor_si128(_mm_slli_si128(x,1),_mm_xor_si128(lo,_mm_slli_si128(hi,1)));
}
template<std::uint32_t C> inline __m128i constantMul(__m128i x) {
    static_assert(C<65536);
    if constexpr(C==0)return _mm_setzero_si128();
    else if constexpr(C==1)return x;
    else if constexpr(C<256)return smallmul128(x,Coeff128{affineMatrix(C,false),affineMatrix(C,true)});
    else return _mm_xor_si128(constantMul<(C&255)>(x),multiplyX8(constantMul<(C>>8)>(x)));
}
template<unsigned Count,unsigned Offset> inline void evalTranspose(__m128i* v) {
    if constexpr(Count>1) {
        constexpr unsigned half=Count/2, degree=std::countr_zero(half);
        evalTranspose<half,Offset>(v);
        evalTranspose<half,Offset^half>(v+half);
        for(unsigned j=0;j<half;++j) {
            const auto sum=_mm_xor_si128(v[j],v[j+half]);
            v[j+half]=_mm_xor_si128(constantMul<subspace(degree,Offset)>(sum),
                                   constantMul<subspace(degree,half)>(v[j+half]));
            v[j]=sum;
        }
    }
}
template<unsigned B> inline void outerTranspose(const Block* source,Block* destination,std::size_t groups) {
    constexpr unsigned k=B/2;
    for(std::size_t group=0;group<groups;++group) {
        __m128i v[B];
        for(unsigned j=0;j<B;++j)v[j]=reg(source[group*B+j]);
        evalTranspose<k,0>(v);
        evalTranspose<k,k>(v+k);
        for(unsigned j=0;j<k;++j)destination[group*k+j]=block(_mm_xor_si128(v[j],v[j+k]));
    }
}

struct FullPlan {
    std::size_t k,n;
    unsigned outer;
    std::vector<std::uint32_t> permutation[2];
    std::vector<std::uint8_t> coefficient[2];
    std::vector<Coeff4> prepared[2];
    std::vector<std::uint32_t> inverse[2];
    std::vector<std::uint8_t> gatheredCoefficient[2];
    FullPlan(std::size_t message,unsigned b,std::uint64_t seed,bool prepare=false):k(message),n(2*message),outer(b) {
        if((b!=8&&b!=16)||!k||k%(b/2)||k>std::numeric_limits<std::uint32_t>::max()/2)
            throw std::invalid_argument("natural message size or outer size");
        std::mt19937_64 rng(seed);
        for(unsigned stage=0;stage<2;++stage) {
            auto& p=permutation[stage]; p.resize(n);std::iota(p.begin(),p.end(),std::uint32_t{0});
            for(std::size_t i=n-1;i;--i) {
                std::uniform_int_distribution<std::size_t> draw(0,i);
                std::swap(p[i],p[draw(rng)]);
            }
            auto& c=coefficient[stage];c.resize(n);
            for(auto& x:c)x=static_cast<std::uint8_t>(1+(rng()&127));
            if(prepare) {
                prepared[stage].resize(n/4);
                for(std::size_t i=0;i<n;i+=4)prepared[stage][i/4]=prepare4(c.data()+i);
            }
        }
    }
    std::size_t setupBytes()const {
        return 2*n*(sizeof(std::uint32_t)+1)+
            (prepared[0].size()+prepared[1].size())*sizeof(Coeff4)+
            (inverse[0].size()+inverse[1].size())*sizeof(std::uint32_t)+
            gatheredCoefficient[0].size()+gatheredCoefficient[1].size();
    }
    void prepareGather() {
        for(unsigned stage=0;stage<2;++stage) {
            inverse[stage].resize(n);gatheredCoefficient[stage].resize(n);
            for(std::size_t i=0;i<n;++i) {
                inverse[stage][permutation[stage][i]]=static_cast<std::uint32_t>(i);
                gatheredCoefficient[stage][permutation[stage][i]]=coefficient[stage][i];
            }
        }
    }
};

inline std::vector<Block> forwardOracle(const FullPlan& plan,const std::vector<Block>& input) {
    if(input.size()!=plan.k)throw std::invalid_argument("oracle input size");
    std::vector<Block> output(plan.n),temp(plan.n);
    const unsigned dim=plan.outer/2;
    for(std::size_t group=0;group<plan.k/dim;++group)
        for(unsigned x=0;x<plan.outer;++x)
            for(unsigned j=0;j<dim;++j)
                output[group*plan.outer+x]=bxor(output[group*plan.outer+x],
                    mulSmallOracle(input[group*dim+j],basis(j,x)));
    for(unsigned stage=0;stage<2;++stage) {
        Block state{};
        for(std::size_t i=0;i<plan.n;++i) {
            state=bxor(state,mulSmallOracle(output[plan.permutation[stage][i]],plan.coefficient[stage][i]));
            temp[i]=state;
        }
        output.swap(temp);
    }
    return output;
}
inline std::vector<Block> transposeOracle(const FullPlan& plan,const std::vector<Block>& input) {
    if(input.size()!=plan.n)throw std::invalid_argument("oracle input size");
    auto x=input;std::vector<Block> tmp(plan.n),out(plan.k);
    for(int stage=1;stage>=0;--stage) {
        Block state{};
        for(std::size_t i=plan.n;i--;) {
            state=bxor(state,x[i]);
            tmp[plan.permutation[stage][i]]=mulSmallOracle(state,plan.coefficient[stage][i]);
        }
        x.swap(tmp);
    }
    const unsigned dim=plan.outer/2;
    for(std::size_t group=0;group<plan.k/dim;++group)
        for(unsigned j=0;j<dim;++j)
            for(unsigned e=0;e<plan.outer;++e)
                out[group*dim+j]=bxor(out[group*dim+j],mulSmallOracle(x[group*plan.outer+e],basis(j,e)));
    return out;
}

inline __m512i suffix4(__m512i x,__m128i carry) {
    x=_mm512_xor_si512(x,_mm512_maskz_permutexvar_epi64(0x3f,_mm512_setr_epi64(2,3,4,5,6,7,0,1),x));
    x=_mm512_xor_si512(x,_mm512_maskz_permutexvar_epi64(0x0f,_mm512_setr_epi64(4,5,6,7,0,1,2,3),x));
    return _mm512_xor_si512(x,_mm512_broadcast_i32x4(carry));
}
template<bool Prepared,unsigned Prefetch=128>
inline void scanTranspose(const Block* source,Block* destination,const std::uint32_t* permutation,
                          const std::uint8_t* coefficients,const Coeff4* prepared,std::size_t n) {
    auto carry=_mm_setzero_si128();
    for(std::size_t end=n;end;end-=4) {
        const auto i=end-4;
        if constexpr(Prefetch)if(i>=Prefetch) {
            _mm_prefetch(reinterpret_cast<const char*>(destination+permutation[i-Prefetch]),_MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(destination+permutation[i-Prefetch+1]),_MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(destination+permutation[i-Prefetch+2]),_MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(destination+permutation[i-Prefetch+3]),_MM_HINT_T0);
        }
        auto x=suffix4(_mm512_loadu_si512(source+i),carry);
        carry=_mm512_castsi512_si128(x);
        if constexpr(Prepared)x=smallmul4(x,prepared[i/4]);
        else x=smallmul4(x,coefficients+i);
        destination[permutation[i]]=block(_mm512_castsi512_si128(x));
        destination[permutation[i+1]]=block(_mm512_extracti32x4_epi32(x,1));
        destination[permutation[i+2]]=block(_mm512_extracti32x4_epi32(x,2));
        destination[permutation[i+3]]=block(_mm512_extracti32x4_epi32(x,3));
    }
}
template<bool Prepared> inline void transposeFast(const FullPlan& plan,const Block* input,Block* output,
                                                 Block* scratch1,Block* scratch2) {
    scanTranspose<Prepared>(input,scratch1,plan.permutation[1].data(),plan.coefficient[1].data(),
                            plan.prepared[1].data(),plan.n);
    scanTranspose<Prepared>(scratch1,scratch2,plan.permutation[0].data(),plan.coefficient[0].data(),
                            plan.prepared[0].data(),plan.n);
    if(plan.outer==8)outerTranspose<8>(scratch2,output,plan.n/8);
    else outerTranspose<16>(scratch2,output,plan.n/16);
}
inline void suffixWrite(const Block* input,Block* output,std::size_t n) {
    auto carry=_mm_setzero_si128();
    for(std::size_t end=n;end;end-=4) {
        auto x=suffix4(_mm512_loadu_si512(input+end-4),carry);
        carry=_mm512_castsi512_si128(x);
        _mm512_storeu_si512(output+end-4,x);
    }
}
template<unsigned Prefetch=64>
inline void gatherMultiply(const Block* input,Block* output,const std::uint32_t* inverse,
                           const std::uint8_t* coefficients,std::size_t n) {
    for(std::size_t i=0;i<n;i+=4) {
        if constexpr(Prefetch)if(i+Prefetch+3<n) {
            _mm_prefetch(reinterpret_cast<const char*>(input+inverse[i+Prefetch]),_MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(input+inverse[i+Prefetch+1]),_MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(input+inverse[i+Prefetch+2]),_MM_HINT_T0);
            _mm_prefetch(reinterpret_cast<const char*>(input+inverse[i+Prefetch+3]),_MM_HINT_T0);
        }
        auto x=_mm512_castsi128_si512(reg(input[inverse[i]]));
        x=_mm512_inserti32x4(x,reg(input[inverse[i+1]]),1);
        x=_mm512_inserti32x4(x,reg(input[inverse[i+2]]),2);
        x=_mm512_inserti32x4(x,reg(input[inverse[i+3]]),3);
        _mm512_storeu_si512(output+i,smallmul4(x,coefficients+i));
    }
}
inline void transposeGather(const FullPlan& plan,const Block* input,Block* output,
                            Block* scratch1,Block* scratch2) {
    suffixWrite(input,scratch1,plan.n);
    gatherMultiply(scratch1,scratch2,plan.inverse[1].data(),plan.gatheredCoefficient[1].data(),plan.n);
    suffixWrite(scratch2,scratch2,plan.n);
    gatherMultiply(scratch2,scratch1,plan.inverse[0].data(),plan.gatheredCoefficient[0].data(),plan.n);
    if(plan.outer==8)outerTranspose<8>(scratch1,output,plan.n/8);
    else outerTranspose<16>(scratch1,output,plan.n/16);
}
} // namespace spin::research::char2
