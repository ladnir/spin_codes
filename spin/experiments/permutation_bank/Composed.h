#pragma once
#include "Flows.h"

namespace spin::experimental::bank {
inline unsigned packDestination(unsigned x){return (x&~1023U)|((x&255U)<<2)|((x>>8)&3U);}
inline unsigned unpackDestination(unsigned x){return (x&~1023U)|((x&3U)<<8)|((x>>2)&255U);}

// A different heuristic family: fixed independent row permutations, independent
// banks for each region, and fresh input wrappers on the chosen region routes.
// Each table directly contains four-row-packed outer destinations.
template<unsigned Count> struct ComposedBank {
    static_assert(std::has_single_bit(Count));
    std::vector<unsigned> routes;
#if SPIN_COMPOSED_ROW_INDEX
    // Experimental redundant layout; the original array remains a scalar
    // reference. Both arrays are reusable bank material, never per-code setup.
    std::vector<unsigned> indexed;
#endif
    explicit ComposedBank(std::uint64_t seed):routes(std::size_t(Count)*CodeSize) {
        Bank rows(8,2048,seed^0x726f7773ULL);
        Words seeds(seed^0x72656773ULL);
        for(unsigned g=0;g<256;++g) {
            Bank regions(11,Count,seeds());
            for(unsigned j=0;j<Count;++j)for(unsigned p=0;p<2048;++p) {
                const unsigned row=regions.p[j*2048+p];
                routes[(g*Count+j)*2048+p]=packDestination(row*256+rows.inverse[row*256+g]);
            }
        }
#if SPIN_COMPOSED_ROW_INDEX
        indexed.resize(routes.size());
        for(std::size_t i=0;i<routes.size();++i) {
            const unsigned x=routes[i],row=((x>>8)&~3U)|(x&3U);
            indexed[i]=x|(row<<19);
        }
#endif
    }
    std::size_t bytes()const{return routes.size()*sizeof(unsigned)
#if SPIN_COMPOSED_ROW_INDEX
        +indexed.size()*sizeof(unsigned)
#endif
        ;}
};

template<unsigned Count,bool ShiftOnly=false,bool RowShift=false,bool RowAffine=false,bool RowRotate=false> class ComposedRouting {
    static_assert(!RowAffine || (RowShift && ShiftOnly));
    static_assert(!RowRotate || RowAffine);
    const ComposedBank<Count>& bank_;
    struct Parameters {unsigned entry,a,b,r;};
    std::array<Parameters,256> params_;
    std::array<unsigned,256> region_;
    std::array<unsigned,RowShift?2048:0> rowDelta_;
    const unsigned* shiftedTable(unsigned g)const {
        const auto offset=(region_[g]*Count+params_[g].entry)*2048;
#if SPIN_COMPOSED_ROW_INDEX
        if constexpr(RowShift)return bank_.indexed.data()+offset;
#endif
        return bank_.routes.data()+offset;
    }
    static unsigned applyColumn(unsigned x,unsigned d) {
        if constexpr(RowAffine) {
            unsigned c=(x>>2)&255;
#if SPIN_COMPOSED_ROTATE_KEYS
            if constexpr(RowRotate) {
                c=std::rotr(std::uint8_t(c),int((d>>26)-2));
                return (x&((CodeSize-1)&~1020U))|((((c*((d>>2)&255)+(d>>10))^(d>>18))&255)<<2);
            }
#endif
            if constexpr(RowRotate)c=std::rotr(std::uint8_t(c),int(d>>24));
            return (x&((CodeSize-1)&~1020U))|((((c*(d&255)+(d>>8))^(d>>16))&255)<<2);
        }
        else return (x&((CodeSize-1)&~1020U))|((x+d)&1020U);
    }
#if defined(__AVX2__)
    static __m256i applyColumn(__m256i x,__m256i d) {
        const auto mask=_mm256_set1_epi32(1020);
        if constexpr(RowAffine) {
            const auto byte=_mm256_set1_epi32(255);
            auto c=_mm256_and_si256(x,mask);
#if SPIN_COMPOSED_ROTATE_KEYS
            if constexpr(RowRotate) {
                c=_mm256_srlv_epi32(_mm256_or_si256(c,_mm256_slli_epi32(c,8)),_mm256_srli_epi32(d,26));
                auto y=_mm256_mullo_epi16(c,_mm256_and_si256(d,mask));
                y=_mm256_and_si256(_mm256_xor_si256(_mm256_add_epi32(y,_mm256_srli_epi32(d,8)),_mm256_srli_epi32(d,16)),mask);
                return _mm256_or_si256(_mm256_and_si256(_mm256_set1_epi32((CodeSize-1)&~1020U),x),y);
            }
#endif
            if constexpr(RowRotate)c=_mm256_and_si256(_mm256_srlv_epi32(_mm256_or_si256(c,_mm256_slli_epi32(c,8)),_mm256_srli_epi32(d,24)),mask);
#if SPIN_BANK_AFFINE16
            auto y=_mm256_mullo_epi16(c,_mm256_and_si256(d,byte));
#else
            auto y=_mm256_mullo_epi32(c,_mm256_and_si256(d,byte));
#endif
            y=_mm256_and_si256(_mm256_xor_si256(_mm256_add_epi32(y,_mm256_srli_epi32(d,6)),_mm256_srli_epi32(d,14)),mask);
            return _mm256_or_si256(_mm256_and_si256(_mm256_set1_epi32((CodeSize-1)&~1020U),x),y);
        } else return _mm256_or_si256(_mm256_and_si256(_mm256_set1_epi32((CodeSize-1)&~1020U),x),_mm256_and_si256(_mm256_add_epi32(x,d),mask));
    }
#endif
#if defined(__AVX512F__)
    static __m512i applyColumn(__m512i x,__m512i d) {
        const auto mask=_mm512_set1_epi32(1020);
        if constexpr(RowAffine) {
            const auto byte=_mm512_set1_epi32(255);
            auto c=_mm512_and_si512(x,mask);
#if SPIN_COMPOSED_ROTATE_KEYS
            if constexpr(RowRotate) {
                c=_mm512_srlv_epi32(_mm512_or_si512(c,_mm512_slli_epi32(c,8)),_mm512_srli_epi32(d,26));
#if defined(__AVX512BW__)
                auto y=_mm512_mullo_epi16(c,_mm512_and_si512(d,mask));
#else
                auto y=_mm512_mullo_epi32(c,_mm512_and_si512(d,mask));
#endif
                y=_mm512_and_si512(_mm512_xor_si512(_mm512_add_epi32(y,_mm512_srli_epi32(d,8)),_mm512_srli_epi32(d,16)),mask);
                return _mm512_or_si512(_mm512_and_si512(_mm512_set1_epi32((CodeSize-1)&~1020U),x),y);
            }
#endif
            if constexpr(RowRotate)c=_mm512_and_si512(_mm512_srlv_epi32(_mm512_or_si512(c,_mm512_slli_epi32(c,8)),_mm512_srli_epi32(d,24)),mask);
#if defined(__AVX512BW__) && SPIN_BANK_AFFINE16
            auto y=_mm512_mullo_epi16(c,_mm512_and_si512(d,byte));
#else
            auto y=_mm512_mullo_epi32(c,_mm512_and_si512(d,byte));
#endif
            y=_mm512_and_si512(_mm512_xor_si512(_mm512_add_epi32(y,_mm512_srli_epi32(d,6)),_mm512_srli_epi32(d,14)),mask);
            return _mm512_or_si512(_mm512_and_si512(_mm512_set1_epi32((CodeSize-1)&~1020U),x),y);
        } else return _mm512_or_si512(_mm512_and_si512(_mm512_set1_epi32((CodeSize-1)&~1020U),x),_mm512_and_si512(_mm512_add_epi32(x,d),mask));
    }
#endif
    template<bool Bytes=false>
#if SPIN_COMPOSED_INLINE_ROUTE
    SPIN_FORCEINLINE
#endif
    void copyShifted(const unsigned* src,unsigned* dst,unsigned n,bool packed)const {
        unsigned i=0;
#if defined(__AVX512F__) && SPIN_COMPOSED_SHIFT512
#if SPIN_COMPOSED_SHIFT_UNROLL
#ifndef SPIN_COMPOSED_UNROLL_VECTORS
#define SPIN_COMPOSED_UNROLL_VECTORS 4
#endif
        constexpr unsigned batch=SPIN_COMPOSED_UNROLL_VECTORS;
        static_assert(batch>0 && batch<=16);
        for(;i+16*batch<=n;i+=16*batch) {
            __m512i x[batch],d[batch];
            for(unsigned j=0;j<batch;++j) {
                x[j]=_mm512_loadu_si512(src+i+16*j);
                const auto row=
#if SPIN_COMPOSED_ROW_INDEX
                    _mm512_srli_epi32(x[j],19);
#else
                    _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi32(x[j],8),_mm512_set1_epi32(~3U)),_mm512_and_si512(x[j],_mm512_set1_epi32(3)));
#endif
                d[j]=_mm512_i32gather_epi32(row,rowDelta_.data(),4);
            }
            for(unsigned j=0;j<batch;++j) {
                auto y=applyColumn(x[j],d[j]);
                if(!packed)y=_mm512_or_si512(_mm512_and_si512(y,_mm512_set1_epi32(~1023U)),_mm512_or_si512(_mm512_slli_epi32(_mm512_and_si512(y,_mm512_set1_epi32(3)),8),_mm512_and_si512(_mm512_srli_epi32(y,2),_mm512_set1_epi32(255))));
                if constexpr(Bytes)y=_mm512_slli_epi32(y,4);
                _mm512_storeu_si512(dst+i+16*j,y);
            }
        }
#endif
        for(;i+16<=n;i+=16) {
            auto x=_mm512_loadu_si512(src+i);
            const auto row=
#if SPIN_COMPOSED_ROW_INDEX
                _mm512_srli_epi32(x,19);
#else
                _mm512_or_si512(_mm512_and_si512(_mm512_srli_epi32(x,8),_mm512_set1_epi32(~3U)),_mm512_and_si512(x,_mm512_set1_epi32(3)));
#endif
            const auto d=_mm512_i32gather_epi32(row,rowDelta_.data(),4);
            x=applyColumn(x,d);
            if(!packed)x=_mm512_or_si512(_mm512_and_si512(x,_mm512_set1_epi32(~1023U)),_mm512_or_si512(_mm512_slli_epi32(_mm512_and_si512(x,_mm512_set1_epi32(3)),8),_mm512_and_si512(_mm512_srli_epi32(x,2),_mm512_set1_epi32(255))));
            if constexpr(Bytes)x=_mm512_slli_epi32(x,4);
            _mm512_storeu_si512(dst+i,x);
        }
#elif defined(__AVX2__)
#if SPIN_COMPOSED_SHIFT256_UNROLL
        constexpr unsigned batch256=SPIN_COMPOSED_SHIFT256_UNROLL;
        static_assert(batch256>0 && batch256<=8);
        for(;i+8*batch256<=n;i+=8*batch256) {
            __m256i x[batch256],d[batch256];
            for(unsigned j=0;j<batch256;++j) {
                x[j]=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src+i+8*j));
                const auto row=
#if SPIN_COMPOSED_ROW_INDEX
                    _mm256_srli_epi32(x[j],19);
#else
                    _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi32(x[j],8),_mm256_set1_epi32(~3U)),_mm256_and_si256(x[j],_mm256_set1_epi32(3)));
#endif
                d[j]=_mm256_i32gather_epi32(reinterpret_cast<const int*>(rowDelta_.data()),row,4);
            }
            for(unsigned j=0;j<batch256;++j) {
                auto y=applyColumn(x[j],d[j]);
                if(!packed)y=_mm256_or_si256(_mm256_and_si256(y,_mm256_set1_epi32(~1023U)),_mm256_or_si256(_mm256_slli_epi32(_mm256_and_si256(y,_mm256_set1_epi32(3)),8),_mm256_and_si256(_mm256_srli_epi32(y,2),_mm256_set1_epi32(255))));
                if constexpr(Bytes)y=_mm256_slli_epi32(y,4);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst+i+8*j),y);
            }
        }
#endif
        for(;i+8<=n;i+=8) {
            auto x=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src+i));
            const auto row=
#if SPIN_COMPOSED_ROW_INDEX
                _mm256_srli_epi32(x,19);
#else
                _mm256_or_si256(_mm256_and_si256(_mm256_srli_epi32(x,8),_mm256_set1_epi32(~3U)),_mm256_and_si256(x,_mm256_set1_epi32(3)));
#endif
            const auto d=_mm256_i32gather_epi32(reinterpret_cast<const int*>(rowDelta_.data()),row,4);
            x=applyColumn(x,d);
            if(!packed)x=_mm256_or_si256(_mm256_and_si256(x,_mm256_set1_epi32(~1023U)),_mm256_or_si256(_mm256_slli_epi32(_mm256_and_si256(x,_mm256_set1_epi32(3)),8),_mm256_and_si256(_mm256_srli_epi32(x,2),_mm256_set1_epi32(255))));
            if constexpr(Bytes)x=_mm256_slli_epi32(x,4);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst+i),x);
        }
#endif
        for(;i<n;++i){unsigned x=src[i];const unsigned row=((x>>8)&2044U)|(x&3U);x=applyColumn(x,rowDelta_[row]);dst[i]=(packed?x:unpackDestination(x))<<(Bytes?4:0);}
    }
    void shiftRows(unsigned* out,bool packed)const {
        if constexpr(RowShift) {
#if defined(__AVX2__)
            for(unsigned i=0;i<2048;i+=8) {
                auto x=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out+i));
                auto row=packed?_mm256_or_si256(_mm256_and_si256(_mm256_srli_epi32(x,8),_mm256_set1_epi32(~3U)),_mm256_and_si256(x,_mm256_set1_epi32(3))):_mm256_srli_epi32(x,8);
                auto d=_mm256_i32gather_epi32(reinterpret_cast<const int*>(rowDelta_.data()),row,4);
                if(!packed)d=_mm256_srli_epi32(d,2);
                const auto mask=_mm256_set1_epi32(packed?1020:255);
                x=_mm256_or_si256(_mm256_andnot_si256(mask,x),_mm256_and_si256(_mm256_add_epi32(x,d),mask));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+i),x);
            }
#else
            for(unsigned i=0;i<2048;++i){const unsigned x=out[i],row=packed?((x>>8)&~3U)|(x&3U):x>>8,mask=packed?1020:255;out[i]=(x&~mask)|((x+(packed?rowDelta_[row]:rowDelta_[row]>>2))&mask);}
#endif
        }
    }
public:
    void prefetchNext(unsigned g)const {
#if SPIN_COMPOSED_LOOKAHEAD
        const auto* table=shiftedTable(g);
        for(unsigned i=0;i<2048;i+=16)_mm_prefetch(reinterpret_cast<const char*>(table+i),_MM_HINT_T1);
#else
        (void)g;
#endif
    }
    void prefetchRows()const {
#if SPIN_COMPOSED_LOOKAHEAD
        if constexpr(RowShift)for(unsigned i=0;i<2048;i+=16)_mm_prefetch(reinterpret_cast<const char*>(rowDelta_.data()+i),_MM_HINT_T0);
#endif
    }
    ComposedRouting(const ComposedBank<Count>& bank,std::uint64_t seed):bank_(bank) {
        Words w(seed);detail::kernel::setup::Divisor rotations(11);
        for(auto& p:params_)p={unsigned(w())&(Count-1),unsigned(w())&2047,unsigned(w())&2047,unsigned(rotations.sample(w,11))};
        // Fresh shared region relabeling. Uniform Fisher--Yates, not a full route.
        std::iota(region_.begin(),region_.end(),0);
        for(unsigned n=256;n>1;--n){detail::kernel::setup::Divisor d(n);std::swap(region_[n-1],region_[d.sample(w,n)]);}
        if constexpr(RowShift)for(auto& d:rowDelta_)d=RowAffine?((unsigned(w())&(RowRotate?0x7ffffff:0xffffff))|1):((unsigned(w())&255)<<2);
#if SPIN_COMPOSED_ROTATE_KEYS
        if constexpr(RowRotate)for(auto& d:rowDelta_)d=((d&0xffffff)<<2)|(((d>>24)+2)<<26);
#endif
    }
    unsigned packed(unsigned inner)const {
        const unsigned g=inner>>11;const auto& p=params_[g];
        const unsigned q=ShiftOnly?((inner+p.b)&2047):((rotate((inner&2047)^p.a,p.r,11)+p.b)&2047);
        unsigned x=bank_.routes[(region_[g]*Count+p.entry)*2048+q];
        if constexpr(RowShift){const unsigned row=((x>>8)&~3U)|(x&3U);x=applyColumn(x,rowDelta_[row]);}
        return x;
    }
    unsigned outer(unsigned inner)const{return unpackDestination(packed(inner));}
    template<unsigned Size>
#if SPIN_COMPOSED_INLINE_ROUTE
    SPIN_FORCEINLINE
#endif
    void outerChunk(unsigned g,unsigned offset,unsigned* out,bool packed)const requires(ShiftOnly) {
        static_assert(Size>=128 && Size<=2048 && std::has_single_bit(Size));
        const auto& p=params_[g];const auto* table=shiftedTable(g);
        const unsigned q=(offset+p.b)&2047,first=std::min(Size,2048-q);
        if constexpr(RowShift) {
            copyShifted(table+q,out,first,packed);
            if(first!=Size)copyShifted(table,out+first,Size-first,packed);
        }else {
            std::memcpy(out,table+q,first*sizeof(unsigned));
            if(first!=Size)std::memcpy(out+first,table,(Size-first)*sizeof(unsigned));
            if(!packed)for(unsigned i=0;i<Size;++i)out[i]=unpackDestination(out[i]);
        }
    }
    void outerEpoch(unsigned g,unsigned offset,unsigned* out,bool packed)const requires(ShiftOnly) {
        outerChunk<128>(g,offset,out,packed);
    }
    SPIN_NOINLINE void outerRegionBytes(unsigned g,unsigned* out,bool packed)const requires(ShiftOnly && RowShift) {
        const auto& p=params_[g];const auto* table=shiftedTable(g);
        copyShifted<true>(table+p.b,out,2048-p.b,packed);
        copyShifted<true>(table,out+2048-p.b,p.b,packed);
    }
#if SPIN_COMPOSED_INLINE_REGION
    SPIN_FORCEINLINE
#elif defined(_MSC_VER)
    __declspec(noinline)
#else
    __attribute__((noinline))
#endif
    void outerRegion(unsigned g,unsigned* out,bool packed=false)const {
        const auto& p=params_[g];const auto* table=bank_.routes.data()+(region_[g]*Count+p.entry)*2048;
        if constexpr(ShiftOnly) {
            table=shiftedTable(g);
            if constexpr(RowShift) {
                copyShifted(table+p.b,out,2048-p.b,packed);
                copyShifted(table,out+2048-p.b,p.b,packed);
                return;
            }
            std::memcpy(out,table+p.b,(2048-p.b)*sizeof(unsigned));
            std::memcpy(out+2048-p.b,table,p.b*sizeof(unsigned));
            if(!packed)for(unsigned i=0;i<2048;++i)out[i]=unpackDestination(out[i]);
            shiftRows(out,packed);
            return;
        }
#if SPIN_COMPOSED_PREFETCH
        for(unsigned i=0;i<2048;i+=16)_mm_prefetch(reinterpret_cast<const char*>(table+i),_MM_HINT_T0);
#endif
#if defined(__AVX512F__) && SPIN_COMPOSED_SCATTER512
        const auto a=_mm512_set1_epi32(p.a),b=_mm512_set1_epi32(p.b),mask=_mm512_set1_epi32(2047);
        for(unsigned i=0;i<2048;i+=16) {
            auto q=_mm512_and_si512(_mm512_sub_epi32(_mm512_add_epi32(_mm512_set1_epi32(i),_mm512_setr_epi32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15)),b),mask);
            q=_mm512_xor_si512(_mm512_and_si512(_mm512_or_si512(_mm512_srl_epi32(q,_mm_cvtsi32_si128(p.r)),_mm512_sll_epi32(q,_mm_cvtsi32_si128(11-p.r))),mask),a);
            auto x=_mm512_loadu_si512(table+i);
            if(!packed)x=_mm512_or_si512(_mm512_and_si512(x,_mm512_set1_epi32(~1023U)),_mm512_or_si512(_mm512_slli_epi32(_mm512_and_si512(x,_mm512_set1_epi32(3)),8),_mm512_and_si512(_mm512_srli_epi32(x,2),_mm512_set1_epi32(255))));
            _mm512_i32scatter_epi32(out,q,x,4);
        }
#elif defined(__AVX2__)
        const auto a=_mm256_set1_epi32(p.a),b=_mm256_set1_epi32(p.b),mask=_mm256_set1_epi32(2047);
        for(unsigned i=0;i<2048;i+=8) {
#if defined(__AVX512VL__) && SPIN_COMPOSED_SCATTER
            auto q=_mm256_and_si256(_mm256_sub_epi32(_mm256_add_epi32(_mm256_set1_epi32(i),_mm256_setr_epi32(0,1,2,3,4,5,6,7)),b),mask);
            q=_mm256_xor_si256(_mm256_and_si256(_mm256_or_si256(_mm256_srl_epi32(q,_mm_cvtsi32_si128(p.r)),_mm256_sll_epi32(q,_mm_cvtsi32_si128(11-p.r))),mask),a);
            auto x=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(table+i));
#else
            auto x=_mm256_xor_si256(_mm256_add_epi32(_mm256_set1_epi32(i),_mm256_setr_epi32(0,1,2,3,4,5,6,7)),a);
            x=_mm256_or_si256(_mm256_sll_epi32(x,_mm_cvtsi32_si128(p.r)),_mm256_srl_epi32(x,_mm_cvtsi32_si128(11-p.r)));
            x=_mm256_and_si256(_mm256_add_epi32(x,b),mask);
            x=_mm256_i32gather_epi32(reinterpret_cast<const int*>(table),x,4);
#endif
            if(!packed)x=_mm256_or_si256(_mm256_and_si256(x,_mm256_set1_epi32(~1023U)),_mm256_or_si256(_mm256_slli_epi32(_mm256_and_si256(x,_mm256_set1_epi32(3)),8),_mm256_and_si256(_mm256_srli_epi32(x,2),_mm256_set1_epi32(255))));
#if defined(__AVX512VL__) && SPIN_COMPOSED_SCATTER
            _mm256_i32scatter_epi32(out,q,x,4);
#else
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+i),x);
#endif
        }
#else
        for(unsigned i=0;i<2048;++i){const unsigned x=table[(rotate(i^p.a,p.r,11)+p.b)&2047];out[i]=packed?x:unpackDestination(x);}
#endif
        shiftRows(out,packed);
    }
};

template<unsigned Count> class DirectComposed {
    ComposedRouting<Count> route_;
    std::vector<std::uint32_t> masks_;
    template<bool Four> void run(kernel::block* input,kernel::block* scratch)const {
        kernel::innerReverse<kernel::Map128S19>(input,CodeSize,masks_.data(),[&](std::size_t i,kernel::block v){
            if constexpr(Four)scratch[route_.packed(unsigned(i))]=v;
            else scratch[route_.outer(unsigned(i))]=v;
        });
        if constexpr(Four)for(unsigned i=0;i<CodeSize;i+=1024)kernel::bchTranspose4(scratch+i,input+i/2);
        else for(unsigned i=0;i<CodeSize;i+=512)kernel::bchTranspose2(scratch+i,scratch+i+256,input+i/2,input+i/2+128);
    }
public:
    DirectComposed(const ComposedBank<Count>& bank,std::uint64_t seed,std::uint64_t maskSeed):route_(bank,seed),masks_(feistel::makeMasksK18(maskSeed)){}
    void encode(kernel::block* input,std::vector<kernel::block>& scratch)const {
        if(scratch.size()!=CodeSize)throw std::invalid_argument("composed scratch");
#if SPIN_BCH_AVX512
        if(kernel::bchAvx512Available()){run<true>(input,scratch.data());return;}
#endif
        run<false>(input,scratch.data());
    }
};
}
