#include "BankKernel.h"
#include "Inner.h"
#include "generated/BchCircuit.h"

namespace spin::detail::kernel {
#if SPIN_BCH_AVX512
void bchTranspose4(const block*,block*);
void bankRunK18(const BankState&,const block*,block*,block*,u32*);
#endif
template<bool Four,bool Bytes> SPIN_NOINLINE void bankRegion(const BankState& b,unsigned g,u32* dst) {
    const auto* table=b.bank.data()+b.region[g]*b.rows;
    const auto copy=[&](const u32* src,std::size_t count) {
        std::size_t i=0;
        // Four independent gathers, matching the tuned prototype's batching.
        const auto batch=[&]<unsigned Batch>() {
            __m256i x[Batch],d[Batch];
            for(unsigned j=0;j<Batch;++j) {
                x[j]=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(src+i+8*j));
                d[j]=_mm256_i32gather_epi32(reinterpret_cast<const int*>(b.row_keys.data()),_mm256_srli_epi32(x[j],8),4);
            }
            for(unsigned j=0;j<Batch;++j) {
                const auto byte=_mm256_set1_epi32(255);
                auto c=_mm256_and_si256(x[j],byte);
                c=_mm256_and_si256(_mm256_srlv_epi32(_mm256_or_si256(c,_mm256_slli_epi32(c,8)),_mm256_srli_epi32(d[j],24)),byte);
                auto y=_mm256_mullo_epi16(c,_mm256_and_si256(d[j],byte));
                y=_mm256_and_si256(_mm256_xor_si256(_mm256_add_epi32(y,_mm256_srli_epi32(d[j],8)),_mm256_srli_epi32(d[j],16)),byte);
                if constexpr(Four) {
                    y=_mm256_or_si256(_mm256_slli_epi32(y,2),_mm256_and_si256(_mm256_srli_epi32(x[j],8),_mm256_set1_epi32(3)));
                    y=_mm256_or_si256(y,_mm256_and_si256(x[j],_mm256_set1_epi32(~1023U)));
                } else y=_mm256_or_si256(y,_mm256_andnot_si256(byte,x[j]));
                if constexpr(Bytes)y=_mm256_slli_epi32(y,4);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(dst+i+8*j),y);
            }
        };
        for(;i+32<=count;i+=32)batch.template operator()<4>();
        for(;i+8<=count;i+=8)batch.template operator()<1>();
        for(;i<count;++i) {
            const auto x=src[i];
            auto y=(x&~255U)|BankState::column(x&255,b.row_keys[x>>8]);
            if constexpr(Four)y=(y&~1023U)|((y&255U)<<2)|((y>>8)&3U);
            dst[i]=Bytes?y<<4:y;
        }
        dst+=count;
    };
    copy(table+b.shift[g],b.rows-b.shift[g]);copy(table,b.shift[g]);
}
template<class Map,bool Four,bool Bytes> void bankRun(const BankState& b,const block* input,block* output,block* scratch,u32* addresses) {
    alignas(32) __m128i state[Map::S]{},syndrome[Map::S],values[Map::T];
    const auto n=2*b.spec.message_size,epochs=n/Map::T;
    for(unsigned g=256;g-->0;) {
        bankRegion<Four,Bytes>(b,g,addresses);
        for(std::size_t offset=b.rows;offset;offset-=Map::T) {
            const auto base=g*b.rows+offset-Map::T,epoch=base/Map::T;
            const auto* epochAddresses=addresses+offset-Map::T;
            auto emit=[&](std::size_t i,block v) {
                if constexpr(Bytes)*reinterpret_cast<block*>(reinterpret_cast<char*>(scratch)+epochAddresses[i%Map::T])=v;
                else scratch[epochAddresses[i%Map::T]]=v;
            };
            if(epoch+1==epochs) {
                for(unsigned p=Map::T;p-->0;){values[p]=input[base+p].mData;emit(base+p,input[base+p]);}
            } else if constexpr(std::is_same_v<Map,Map128S19>)
                imtReversePoints<Map>(input+base,values,state,base,emit,std::make_index_sequence<Map::T>{});
            else Map::emitShared(input+base,values,state,base,emit);
            if(epoch==0)break;
            zeta<Map::T>(values);Map::finish(values,syndrome);
            if(epoch+1==epochs)std::memcpy(state,syndrome,sizeof(state));
            else {
                if constexpr(isTwoRoundMap<Map>) {
                    imtStep(state,b.masks[4*epoch+2],b.masks[4*epoch+3]);
                    imtStep(state,b.masks[4*epoch],b.masks[4*epoch+1]);
                } else imtStep(state,b.masks[2*epoch],b.masks[2*epoch+1]);
                for(unsigned j=0;j<Map::S;++j)state[j]=_mm_xor_si128(state[j],syndrome[j]);
            }
        }
    }
#if SPIN_BCH_AVX512
    if constexpr(Four)for(std::size_t i=0;i<n;i+=1024)bchTranspose4(scratch+i,output+i/2);
    else
#endif
    for(std::size_t i=0;i<n;i+=512)bchTranspose2(scratch+i,scratch+i+256,output+i/2,output+i/2+128);
}
template<class Map> void bankDispatch(const BankState& b,bool four,const block* in,block* out,block* scratch,u32* addresses) {
    const bool bytes=2*b.spec.message_size<=std::size_t(1)<<28;
#if SPIN_BCH_AVX512
    if(four) {
        if(bytes)bankRun<Map,true,true>(b,in,out,scratch,addresses);
        else bankRun<Map,true,false>(b,in,out,scratch,addresses);
        return;
    }
#else
    (void)four;
#endif
    if(bytes)bankRun<Map,false,true>(b,in,out,scratch,addresses);
    else bankRun<Map,false,false>(b,in,out,scratch,addresses);
}
void bankTranspose(const BankState& b,bool four,const void* in,void* out,block* scratch,u32* addresses) {
    const auto* input=static_cast<const block*>(in);auto* output=static_cast<block*>(out);
#if SPIN_BCH_AVX512
    if(four && !b.indexed.empty()){bankRunK18(b,input,output,scratch,addresses);return;}
#endif
    switch(b.spec.parameters) {
    case Parameters::T128S19:bankDispatch<Map128S19>(b,four,input,output,scratch,addresses);break;
    case Parameters::T64S12:bankDispatch<Map64S12>(b,four,input,output,scratch,addresses);break;
    case Parameters::T64S12R2:bankDispatch<Map64S12R2>(b,four,input,output,scratch,addresses);break;
    }
}
}
