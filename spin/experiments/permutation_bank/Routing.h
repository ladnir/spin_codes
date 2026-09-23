#pragma once
#include "Bank.h"
#include <array>
#if defined(__AVX2__)
#include <immintrin.h>
#endif
namespace spin::experimental::bank {
// Experimental K=2^18 inverse routing. Preprocessed banks persist across codes.
template<unsigned Count> struct Tables {
    static_assert(Count==16 || Count==64);
    std::vector<std::uint8_t> row;
    std::vector<std::uint16_t> region;
    explicit Tables(std::uint64_t seed):row(256*Count+3),region(2048*Count+1) {
        Bank r(8,Count,seed^0x726f7773ULL),g(11,Count,seed^0x72656773ULL);
        for(unsigned i=0;i<256*Count;++i)row[i]=std::uint8_t(r.inverse[i]);
        for(unsigned i=0;i<2048*Count;++i)region[i]=std::uint16_t(g.inverse[i]);
    }
    std::size_t bytes()const{return row.size()+2*region.size();}
};
struct PackedParameters {
    std::uint16_t entry,a,b,c,e;
    std::uint8_t ri,ro;
};
template<unsigned Bits> inline unsigned undoRotation(unsigned x,unsigned r) {
    if constexpr(Bits==8)return std::rotr(std::uint8_t(x),int(r));
    else return ((x>>r)|(x<<(Bits-r)))&((1U<<Bits)-1);
}
template<unsigned Bits> inline unsigned undoOutput(unsigned x,const PackedParameters& k) {
    return undoRotation<Bits>((x-k.e)&((1U<<Bits)-1),k.ro)^k.c;
}
template<unsigned Bits> inline unsigned undoInput(unsigned x,const PackedParameters& k) {
    return undoRotation<Bits>((x-k.b)&((1U<<Bits)-1),k.ri)^k.a;
}
template<unsigned Count,bool Vector=false> class Routing {
    const Tables<Count>& tables_;
    std::vector<PackedParameters> rows_,regions_;
    std::vector<std::uint32_t> rowAB_,rowCE_;
    template<unsigned Bits> static PackedParameters sample(Words& w,const detail::kernel::setup::Divisor& rd) {
        constexpr unsigned mask=(1U<<Bits)-1;
        return {std::uint16_t(w()&(Count-1)),std::uint16_t(w()&mask),std::uint16_t(w()&mask),
            std::uint16_t(w()&mask),std::uint16_t(w()&mask),std::uint8_t(rd.sample(w,Bits)),std::uint8_t(rd.sample(w,Bits))};
    }
public:
    static constexpr unsigned Size=1U<<19;
    Routing(const Tables<Count>& tables,std::uint64_t seed):tables_(tables),rows_(2048),regions_(256),rowAB_(Vector?2048:0),rowCE_(Vector?2048:0) {
        Words rw(seed^0x726f7773ULL),gw(seed^0x72656773ULL);
        detail::kernel::setup::Divisor dr(8),dg(11);
        for(auto& k:rows_)k=sample<8>(rw,dr);
        for(auto& k:regions_)k=sample<11>(gw,dg);
        if constexpr(Vector)for(unsigned row=0;row<2048;++row) {
            const auto& k=rows_[row];rowAB_[row]=k.a|(k.b<<8)|(k.ri<<16);
            rowCE_[row]=k.c|(k.e<<8)|(k.ro<<16)|(k.entry<<20);
        }
    }
    unsigned outer(unsigned inner)const {
        const unsigned region=inner>>11;const auto& g=regions_[region];
        const unsigned row=undoInput<11>(tables_.region[g.entry*2048+undoOutput<11>(inner&2047,g)],g);
        const auto& r=rows_[row];
        return row*256+undoInput<8>(tables_.row[r.entry*256+undoOutput<8>(region,r)],r);
    }
    void outerBatch(unsigned first,unsigned* out)const {
        const unsigned region=first>>11;const auto& g=regions_[region];
        const auto* gp=tables_.region.data()+g.entry*2048;
        std::array<unsigned,16> row,at;
        // Stage independent lookups to expose memory-level parallelism.
        for(unsigned i=0;i<16;++i)row[i]=undoInput<11>(gp[undoOutput<11>((first+i)&2047,g)],g);
        for(unsigned i=0;i<16;++i){const auto& r=rows_[row[i]];at[i]=r.entry*256+undoOutput<8>(region,r);}
        for(unsigned i=0;i<16;++i)at[i]=tables_.row[at[i]];
        for(unsigned i=0;i<16;++i)out[i]=row[i]*256+undoInput<8>(at[i],rows_[row[i]]);
    }
#if defined(_MSC_VER)
    __declspec(noinline)
#else
    __attribute__((noinline))
#endif
    void outerRegion(unsigned region,unsigned* out,bool packed=false)const {
        // One 8 KiB region, not a full route. Stream row parameters in order
        // before applying the region permutation, avoiding random parameter
        // loads on the destination-address dependency chain.
        alignas(64) std::array<unsigned,2048> coordinate;
#if defined(__AVX512F__) && SPIN_BANK_ROUTE512
        if constexpr(Vector) {
            const auto mask=_mm512_set1_epi32(255),eight=_mm512_set1_epi32(8),seven=_mm512_set1_epi32(7);
            const auto rg=_mm512_set1_epi32(region);
            for(unsigned row=0;row<2048;row+=16) {
                const auto ab=_mm512_loadu_si512(reinterpret_cast<const __m512i*>(rowAB_.data()+row));
                const auto ce=_mm512_loadu_si512(reinterpret_cast<const __m512i*>(rowCE_.data()+row));
                const auto ro=_mm512_and_si512(_mm512_srli_epi32(ce,16),seven);
                const auto ri=_mm512_and_si512(_mm512_srli_epi32(ab,16),seven);
                auto x=_mm512_and_si512(_mm512_sub_epi32(rg,_mm512_and_si512(_mm512_srli_epi32(ce,8),mask)),mask);
                x=_mm512_and_si512(_mm512_or_si512(_mm512_srlv_epi32(x,ro),_mm512_sllv_epi32(x,_mm512_sub_epi32(eight,ro))),mask);
                x=_mm512_xor_si512(x,_mm512_and_si512(ce,mask));
                const auto at=_mm512_add_epi32(x,_mm512_slli_epi32(_mm512_srli_epi32(ce,20),8));
                // Three zero padding bytes permit a 32-bit gather at the last byte.
                x=_mm512_and_si512(_mm512_i32gather_epi32(at,tables_.row.data(),1),mask);
                x=_mm512_and_si512(_mm512_sub_epi32(x,_mm512_and_si512(_mm512_srli_epi32(ab,8),mask)),mask);
                x=_mm512_and_si512(_mm512_or_si512(_mm512_srlv_epi32(x,ri),_mm512_sllv_epi32(x,_mm512_sub_epi32(eight,ri))),mask);
                x=_mm512_xor_si512(x,_mm512_and_si512(ab,mask));
                const auto ids=_mm512_add_epi32(_mm512_set1_epi32(row),_mm512_setr_epi32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15));
                if(packed)x=_mm512_or_si512(_mm512_slli_epi32(x,2),_mm512_or_si512(_mm512_slli_epi32(_mm512_and_si512(ids,_mm512_set1_epi32(~3U)),8),_mm512_and_si512(ids,_mm512_set1_epi32(3))));
                else x=_mm512_or_si512(x,_mm512_slli_epi32(ids,8));
                _mm512_store_si512(reinterpret_cast<__m512i*>(coordinate.data()+row),x);
            }
            const auto& g=regions_[region];const auto* gp=tables_.region.data()+g.entry*2048;
            const auto mask11=_mm512_set1_epi32(2047),a=_mm512_set1_epi32(g.a),b=_mm512_set1_epi32(g.b),c=_mm512_set1_epi32(g.c),e=_mm512_set1_epi32(g.e);
            for(unsigned i=0;i<2048;i+=16) {
                auto x=_mm512_add_epi32(_mm512_set1_epi32(i),_mm512_setr_epi32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15));
                x=_mm512_and_si512(_mm512_sub_epi32(x,e),mask11);
                x=_mm512_and_si512(_mm512_or_si512(_mm512_srl_epi32(x,_mm_cvtsi32_si128(g.ro)),_mm512_sll_epi32(x,_mm_cvtsi32_si128(11-g.ro))),mask11);
                x=_mm512_xor_si512(x,c);
                x=_mm512_and_si512(_mm512_i32gather_epi32(x,gp,2),mask11);
                x=_mm512_and_si512(_mm512_sub_epi32(x,b),mask11);
                x=_mm512_and_si512(_mm512_or_si512(_mm512_srl_epi32(x,_mm_cvtsi32_si128(g.ri)),_mm512_sll_epi32(x,_mm_cvtsi32_si128(11-g.ri))),mask11);
                x=_mm512_xor_si512(x,a);
                x=_mm512_i32gather_epi32(x,coordinate.data(),4);
                _mm512_storeu_si512(reinterpret_cast<__m512i*>(out+i),x);
            }
            return;
        }
#endif
#if defined(__AVX2__)
        if constexpr(Vector) {
            const auto mask=_mm256_set1_epi32(255),eight=_mm256_set1_epi32(8),seven=_mm256_set1_epi32(7);
            const auto rg=_mm256_set1_epi32(region);
            for(unsigned row=0;row<2048;row+=8) {
                const auto ab=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(rowAB_.data()+row));
                const auto ce=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(rowCE_.data()+row));
                const auto ro=_mm256_and_si256(_mm256_srli_epi32(ce,16),seven);
                const auto ri=_mm256_and_si256(_mm256_srli_epi32(ab,16),seven);
                auto x=_mm256_and_si256(_mm256_sub_epi32(rg,_mm256_and_si256(_mm256_srli_epi32(ce,8),mask)),mask);
                x=_mm256_and_si256(_mm256_or_si256(_mm256_srlv_epi32(x,ro),_mm256_sllv_epi32(x,_mm256_sub_epi32(eight,ro))),mask);
                x=_mm256_xor_si256(x,_mm256_and_si256(ce,mask));
                const auto at=_mm256_add_epi32(x,_mm256_slli_epi32(_mm256_srli_epi32(ce,20),8));
                // Three zero padding bytes permit a 32-bit gather at the last byte.
                x=_mm256_and_si256(_mm256_i32gather_epi32(reinterpret_cast<const int*>(tables_.row.data()),at,1),mask);
                x=_mm256_and_si256(_mm256_sub_epi32(x,_mm256_and_si256(_mm256_srli_epi32(ab,8),mask)),mask);
                x=_mm256_and_si256(_mm256_or_si256(_mm256_srlv_epi32(x,ri),_mm256_sllv_epi32(x,_mm256_sub_epi32(eight,ri))),mask);
                x=_mm256_xor_si256(x,_mm256_and_si256(ab,mask));
                const auto ids=_mm256_add_epi32(_mm256_set1_epi32(row),_mm256_setr_epi32(0,1,2,3,4,5,6,7));
                if(packed)x=_mm256_or_si256(_mm256_slli_epi32(x,2),_mm256_or_si256(_mm256_slli_epi32(_mm256_and_si256(ids,_mm256_set1_epi32(~3U)),8),_mm256_and_si256(ids,_mm256_set1_epi32(3))));
                else x=_mm256_or_si256(x,_mm256_slli_epi32(ids,8));
                _mm256_store_si256(reinterpret_cast<__m256i*>(coordinate.data()+row),x);
            }
            const auto& g=regions_[region];const auto* gp=tables_.region.data()+g.entry*2048;
            const auto mask11=_mm256_set1_epi32(2047),a=_mm256_set1_epi32(g.a),b=_mm256_set1_epi32(g.b),c=_mm256_set1_epi32(g.c),e=_mm256_set1_epi32(g.e);
            for(unsigned i=0;i<2048;i+=8) {
                auto x=_mm256_add_epi32(_mm256_set1_epi32(i),_mm256_setr_epi32(0,1,2,3,4,5,6,7));
                x=_mm256_and_si256(_mm256_sub_epi32(x,e),mask11);
                x=_mm256_and_si256(_mm256_or_si256(_mm256_srl_epi32(x,_mm_cvtsi32_si128(g.ro)),_mm256_sll_epi32(x,_mm_cvtsi32_si128(11-g.ro))),mask11);
                x=_mm256_xor_si256(x,c);
                x=_mm256_and_si256(_mm256_i32gather_epi32(reinterpret_cast<const int*>(gp),x,2),mask11);
                x=_mm256_and_si256(_mm256_sub_epi32(x,b),mask11);
                x=_mm256_and_si256(_mm256_or_si256(_mm256_srl_epi32(x,_mm_cvtsi32_si128(g.ri)),_mm256_sll_epi32(x,_mm_cvtsi32_si128(11-g.ri))),mask11);
                x=_mm256_xor_si256(x,a);
                x=_mm256_i32gather_epi32(reinterpret_cast<const int*>(coordinate.data()),x,4);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+i),x);
            }
            return;
        }
#endif
        for(unsigned row=0;row<2048;++row) {
            const auto& r=rows_[row];
            const unsigned c=undoInput<8>(tables_.row[r.entry*256+undoOutput<8>(region,r)],r);
            coordinate[row]=packed?((row&~3U)*256+(c<<2)+(row&3)):row*256+c;
        }
        const auto& g=regions_[region];const auto* gp=tables_.region.data()+g.entry*2048;
        for(unsigned i=0;i<2048;++i) {
            const unsigned row=undoInput<11>(gp[undoOutput<11>(i,g)],g);
            out[i]=coordinate[row];
        }
    }
    std::size_t bytes()const{return sizeof(PackedParameters)*(rows_.size()+regions_.size())+4*(rowAB_.size()+rowCE_.size());}
};
}
