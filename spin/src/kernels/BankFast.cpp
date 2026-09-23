#include "BankKernel.h"
#include "Inner.h"

namespace spin::detail::kernel {
void bchTranspose4(const block*,block*);

// Preserve the prototype's fixed-width batching and byte-address scatter.
SPIN_NOINLINE static void bankRegionK18(const BankState& b,unsigned g,u32* dst) {
    const auto* table=b.indexed.data()+b.region[g]*2048;
    auto copy=[&](const u32* src,unsigned count) {
        unsigned i=0;
        auto batch=[&]<unsigned Batch>() {
            __m256i x[Batch],d[Batch];
            for(unsigned j=0;j<Batch;++j) {
                x[j]=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src+i+8*j));
                d[j]=_mm256_i32gather_epi32(reinterpret_cast<const int*>(b.packed_keys.data()),_mm256_srli_epi32(x[j],19),4);
            }
            for(unsigned j=0;j<Batch;++j) {
                const auto mask=_mm256_set1_epi32(1020);
                auto c=_mm256_and_si256(x[j],mask);
                c=_mm256_srlv_epi32(_mm256_or_si256(c,_mm256_slli_epi32(c,8)),_mm256_srli_epi32(d[j],26));
                auto y=_mm256_mullo_epi16(c,_mm256_and_si256(d[j],mask));
                y=_mm256_and_si256(_mm256_xor_si256(_mm256_add_epi32(y,_mm256_srli_epi32(d[j],8)),_mm256_srli_epi32(d[j],16)),mask);
                y=_mm256_or_si256(y,_mm256_and_si256(x[j],_mm256_set1_epi32(((1U<<19)-1)&~1020U)));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst+i+8*j),_mm256_slli_epi32(y,4));
            }
        };
        for(;i+32<=count;i+=32)batch.template operator()<4>();
        for(;i+8<=count;i+=8)batch.template operator()<1>();
        for(;i<count;++i) {
            const auto x=src[i],d=b.packed_keys[x>>19];
            const auto c=std::rotr(std::uint8_t((x>>2)&255),int((d>>26)-2));
            const auto y=(x&(((1U<<19)-1)&~1020U))|((((c*((d>>2)&255)+(d>>10))^(d>>18))&255)<<2);
            dst[i]=y<<4;
        }
        dst+=count;
    };
    copy(table+b.shift[g],2048-b.shift[g]);copy(table,b.shift[g]);
}
void bankRunK18(const BankState& b,const block* input,block* output,block* scratch,u32*) {
    // A private fixed-size route buffer also tells the compiler these addresses
    // cannot alias input, output, or the scatter destination. This is an ordinary
    // non-coroutine helper; the aligned scratch never crosses a suspension.
    alignas(64) u32 addresses[2048];
    using Map=Map128S19;
    alignas(32) __m128i state[Map::S]{},syndrome[Map::S],values[Map::T];
    constexpr unsigned n=1U<<19,epochs=n/Map::T;
    for(unsigned epoch=epochs;epoch-->0;) {
        const unsigned base=epoch*Map::T;
        if((epoch&15)==15) {
            bankRegionK18(b,base>>11,addresses);
            if(base>=2048) {
                const auto* next=b.indexed.data()+b.region[(base>>11)-1]*2048;
                for(unsigned i=0;i<2048;i+=16)_mm_prefetch(reinterpret_cast<const char*>(next+i),_MM_HINT_T1);
            }
        }
        if((epoch&15)==0)for(unsigned i=0;i<2048;i+=16)_mm_prefetch(reinterpret_cast<const char*>(b.packed_keys.data()+i),_MM_HINT_T0);
        const auto* epochAddresses=addresses+(base&2047);
        auto emit=[&](std::size_t i,block v){*reinterpret_cast<block*>(reinterpret_cast<char*>(scratch)+epochAddresses[i&127])=v;};
        if(epoch+1==epochs) {
            for(unsigned p=Map::T;p-->0;){values[p]=input[base+p].mData;emit(base+p,input[base+p]);}
        } else imtReversePoints<Map>(input+base,values,state,base,emit,std::make_index_sequence<Map::T>{});
        if(epoch==0)break;
        zeta<Map::T>(values);Map::finish(values,syndrome);
        if(epoch+1==epochs)std::memcpy(state,syndrome,sizeof(state));
        else {
            imtStep(state,b.masks[2*epoch],b.masks[2*epoch+1]);
            for(unsigned j=0;j<Map::S;++j)state[j]=_mm_xor_si128(state[j],syndrome[j]);
        }
    }
    for(unsigned i=0;i<n;i+=1024)bchTranspose4(scratch+i,output+i/2);
}
}
