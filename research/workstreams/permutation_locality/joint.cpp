// Research-only joint packet/inner screen. No production defaults change.
#include "Spin.h"
#include "Inner.h"
#include "SetupRandom.h"
#include "WorkspaceRouting.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>

namespace k=spin::detail::kernel;
namespace spin::detail::kernel {
struct JointMapR2 : Map128S19 {
    template<class Emit> static SPIN_FORCEINLINE void emitShared(
        const block* in,__m128i* raw,const __m128i* state,std::size_t base,Emit& emit) {
        imtReversePoints<Map128S19>(in,raw,state,base,emit,std::make_index_sequence<128>{});
    }
};
template<> inline constexpr bool isTwoRoundMap<JointMapR2> = true;
void bchTranspose4(const block*,block*);
void bchTranspose4GfniBlend(const block*,block*);
}
using Clock=std::chrono::steady_clock;
using Map=k::Map128S19;
constexpr unsigned tileStride=1028;

// Fixed-width update sequence: preserve the original unrolled emission and
// zeta circuit, changing only the number of independent state updates.
template<unsigned R> SPIN_FORCEINLINE void reverseUpdates(__m128i* state,const unsigned* masks) {
    k::imtStep(state,masks[2*R-2],masks[2*R-1]);
    if constexpr(R>1)reverseUpdates<R-1>(state,masks);
}
template<unsigned R,class Emit> SPIN_FORCEINLINE void reverseInner(
    const k::block* in,std::size_t n,const unsigned* masks,Emit&& emit) {
    static_assert(R==2 || R==3);
    alignas(32) __m128i state[19]{},values[128],syndrome[19];
    for(std::size_t epoch=n/128;epoch-->0;) {
        const auto base=128*epoch;
        if(epoch+1==n/128) {
            for(unsigned p=128;p-->0;){values[p]=in[base+p].mData;emit(base+p,in[base+p]);}
        } else k::imtReversePoints<Map>(in+base,values,state,base,emit,std::make_index_sequence<128>{});
        if(!epoch)break;
        k::zeta<128>(values);Map::finish(values,syndrome);
        if(epoch+1==n/128)std::memcpy(state,syndrome,sizeof(state));
        else {
            reverseUpdates<R>(state,masks+2*R*epoch);
            for(unsigned j=0;j<19;++j)state[j]=_mm_xor_si128(state[j],syndrome[j]);
        }
    }
}

struct Sampler {
    k::setup::Words words;
    std::vector<k::setup::Divisor> divisors;
    Sampler(unsigned rows,std::uint64_t seed):words(seed),divisors(std::max(rows,256U)+1) {
        for(unsigned i=2;i<divisors.size();++i)divisors[i]=k::setup::Divisor(i);
    }
    template<class V> void shuffle(V& v) {
        std::iota(v.begin(),v.end(),0U);
        for(unsigned i=unsigned(v.size());i>1;--i)std::swap(v[i-1],v[divisors[i].sample(words,i)]);
    }
};

template<unsigned G> struct Packets {
    static_assert(G==2 || G==4);
    std::vector<unsigned> inverse,bases;
    std::vector<std::uint64_t> controls;
    std::vector<std::uint16_t> offsets;
    Packets(unsigned rows,std::uint64_t seed):inverse(rows*256),bases(inverse.size()/G),
        controls(bases.size()),offsets(inverse.size(),65535) {
        Sampler rng(rows,seed);
        std::vector<unsigned> coordinates(inverse.size()),groups(rows/G);
        std::array<unsigned,256> permutation;
        std::array<unsigned,G> lanes;
        for(unsigned row=0;row<rows;++row) {
            rng.shuffle(permutation);
            for(unsigned c=0;c<256;++c)coordinates[row*256+c]=permutation[c];
        }
        // Every row shuffle, regional packet permutation, and lane shuffle
        // uses a separate draw. There is no shared row permutation.
        for(unsigned region=0;region<256;++region) {
            rng.shuffle(groups);
            for(unsigned group=0;group<rows/G;++group) {
                rng.shuffle(lanes);
                for(unsigned lane=0;lane<G;++lane) {
                    const auto row=group*G+lane;
                    inverse[region*rows+groups[group]*G+lanes[lane]]=row*256+coordinates[row*256+region];
                }
            }
        }
        std::vector<bool> seen(inverse.size());
        for(unsigned region=0;region<256;++region) {
            std::vector<bool> rowSeen(rows);
            for(unsigned pos=0;pos<rows;pos+=G) {
                const auto i=region*rows+pos;
                const auto packetGroup=inverse[i]/(256*G),tile=packetGroup/(4/G);
                const auto localBase=(packetGroup%(4/G))*256*G;
                unsigned anchor=256,laneSeen=0;
                for(unsigned l=0;l<G;++l)if((inverse[i+l]/256)%G==0)anchor=inverse[i+l]%256;
                if(anchor==256)throw std::runtime_error("packet lacks anchor row");
                bases[i/G]=tile*tileStride+localBase+G*anchor;
                for(unsigned l=0;l<G;++l) {
                    const auto x=inverse[i+l],row=x/256,column=x%256,lane=row%G;
                    if(x>=inverse.size() || seen[x] || rowSeen[row] || row/G!=packetGroup || (laneSeen&(1U<<lane)))
                        throw std::runtime_error("route/packet bijection");
                    seen[x]=true;rowSeen[row]=true;laneSeen|=1U<<lane;
                    auto& offset=offsets[tile*1024+4*column+row%4];
                    if(offset!=65535)throw std::runtime_error("duplicate row offset");
                    offset=std::uint16_t(localBase+G*anchor+lane);
                    for(unsigned half=0;half<2;++half)
                        controls[i/G]|=std::uint64_t(2*l+half)<<(8*(2*lane+half));
                }
            }
        }
        for(auto offset:offsets)if(offset>=1024)throw std::runtime_error("missing row offset");
    }
};

template<unsigned G> static SPIN_NOINLINE void repackBch(
    const k::block* scratch,k::block* out,std::size_t n,const std::uint16_t* offsets) {
    alignas(64) k::block local[1024],canonical[1024];
    for(std::size_t tile=0;tile<n/1024;++tile) {
        const auto* source=scratch+tile*tileStride;
        // Copy cold routed data sequentially before the local permutation.
        for(unsigned i=0;i<1024;i+=16) {
            const auto a=_mm512_load_si512(source+i),b=_mm512_load_si512(source+i+4);
            const auto c=_mm512_load_si512(source+i+8),d=_mm512_load_si512(source+i+12);
            _mm512_store_si512(local+i,a);_mm512_store_si512(local+i+4,b);
            _mm512_store_si512(local+i+8,c);_mm512_store_si512(local+i+12,d);
        }
        const auto* p=offsets+tile*1024;
#if defined(__GNUC__)
#pragma GCC unroll 4
#endif
        for(unsigned c=0;c<256;++c) {
            auto x=_mm512_castsi128_si512(_mm_load_si128(reinterpret_cast<const __m128i*>(local+G*c)));
            x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(local+p[4*c+1])),1);
            if constexpr(G==2)x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(local+512+2*c)),2);
            else x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(local+p[4*c+2])),2);
            x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(local+p[4*c+3])),3);
            _mm512_store_si512(canonical+4*c,x);
        }
        k::bchTranspose4GfniBlend(canonical,out+tile*512);
    }
}

template<unsigned G,unsigned R,bool Streaming> static SPIN_NOINLINE void encode(
    k::block* in,k::block* scratch,std::size_t n,const unsigned* masks,const Packets<G>& packets,double* phases) {
    __m128i v1{},v2{},v3{};
    auto start=Clock::time_point{};if(phases)start=Clock::now();
    reverseInner<R>(in,n,masks,[&](std::size_t i,k::block v) {
        if constexpr(G==2) {
            if(i&1)v1=v.mData;
            else {
                auto x=_mm256_inserti128_si256(_mm256_castsi128_si256(v.mData),v1,1);
                const auto control=_mm256_cvtepu8_epi64(_mm_cvtsi32_si128(int(packets.controls[i/2])));
                x=_mm256_permutexvar_epi64(control,x);
                auto* dst=reinterpret_cast<__m256i*>(scratch+packets.bases[i/2]);
                if constexpr(Streaming)_mm256_stream_si256(dst,x);else _mm256_store_si256(dst,x);
            }
        } else {
            switch(i&3) {
                case 3:v3=v.mData;break;
                case 2:v2=v.mData;break;
                case 1:v1=v.mData;break;
                case 0:{
                    auto x=_mm512_castsi128_si512(v.mData);
                    x=_mm512_inserti32x4(x,v1,1);x=_mm512_inserti32x4(x,v2,2);x=_mm512_inserti32x4(x,v3,3);
                    const auto control=_mm512_cvtepu8_epi64(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(packets.controls.data()+i/4)));
                    x=_mm512_permutexvar_epi64(control,x);
                    auto* dst=reinterpret_cast<__m512i*>(scratch+packets.bases[i/4]);
                    if constexpr(Streaming)_mm512_stream_si512(dst,x);else _mm512_store_si512(dst,x);
                    break;
                }
            }
        }
    });
    if constexpr(Streaming)_mm_sfence();
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    repackBch<G>(scratch,in,n,packets.offsets.data());
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}

// Dense matrix-row reference, independent of sparse SIMD mixing and zeta.
template<unsigned R> static void reverseReference(const k::block* in,k::block* out,std::size_t n,const unsigned* masks) {
    std::array<k::block,19> state{},next{};
    for(std::size_t epoch=n/128;epoch-->0;) {
        for(unsigned p=0;p<128;++p) {
            auto v=in[128*epoch+p].mData;
            for(unsigned j=0;j<19;++j)if((Map::feedbackColumns[p]>>j)&1)v=_mm_xor_si128(v,state[j].mData);
            out[128*epoch+p]=k::block(v);
        }
        for(unsigned j=0;j<19;++j) {
            auto row=1U<<j;
            for(unsigned r=0;r<R;++r)if(std::popcount(row&masks[2*(R*epoch+r)+1])&1)row^=masks[2*(R*epoch+r)];
            auto v=_mm_setzero_si128();
            for(unsigned i=0;i<19;++i)if((row>>i)&1)v=_mm_xor_si128(v,state[i].mData);
            for(unsigned p=0;p<128;++p)if((Map::columns[p]>>j)&1)v=_mm_xor_si128(v,in[128*epoch+p].mData);
            next[j]=k::block(v);
        }
        state=next;
    }
}
template<unsigned R> static void forwardReference(const k::block* in,k::block* out,std::size_t n,const unsigned* masks) {
    std::array<k::block,19> state{},syndrome{};
    for(std::size_t epoch=0;epoch<n/128;++epoch) {
        syndrome={};
        for(unsigned p=0;p<128;++p) {
            auto v=in[128*epoch+p].mData;
            for(unsigned j=0;j<19;++j) {
                if((Map::columns[p]>>j)&1)v=_mm_xor_si128(v,state[j].mData);
                if((Map::feedbackColumns[p]>>j)&1)syndrome[j].mData=_mm_xor_si128(syndrome[j].mData,in[128*epoch+p].mData);
            }
            out[128*epoch+p]=k::block(v);
        }
        for(unsigned r=0;r<R;++r) {
            const auto u=masks[2*(R*epoch+r)],v=masks[2*(R*epoch+r)+1];
            auto dot=_mm_setzero_si128();
            for(unsigned j=0;j<19;++j)if((v>>j)&1)dot=_mm_xor_si128(dot,state[j].mData);
            for(unsigned j=0;j<19;++j)if((u>>j)&1)state[j].mData=_mm_xor_si128(state[j].mData,dot);
        }
        for(unsigned j=0;j<19;++j)state[j].mData=_mm_xor_si128(state[j].mData,syndrome[j].mData);
    }
}
static void fill(std::vector<k::block>& x,std::uint64_t seed) {
    k::setup::Words rng(seed);for(auto& v:x){auto a=rng(),b=rng();v=k::block(a,b);}
}
static std::uint64_t hash(const std::vector<k::block>& x) {
    std::uint64_t h=0;for(const auto& v:x){std::uint64_t halves[2];std::memcpy(halves,&v,16);for(auto a:halves)h=(h^a)*0x100000001b3ULL;}return h;
}
template<unsigned G,unsigned R,bool Streaming> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    const unsigned n=2U<<exponent;
    const Packets<G> packets(n/256,seed);
    std::vector<unsigned> masks(2*R*(n/128));
    k::setup::Words rng(seed^0x3f625a92ULL);
    for(std::size_t i=0;i<masks.size();i+=2) {
        unsigned u;do{u=unsigned(rng())&((1U<<19)-1);}while(!u);
        unsigned v=unsigned(rng())&((1U<<19)-1);
        if(std::popcount(u&v)&1)v^=u&-u;
        masks[i]=u;masks[i+1]=v;
    }
    std::vector<k::block> input(n),inner(n),actual(n),other(n),forward(n),expected(n),materialized(n);
    fill(input,913);fill(other,971);expected=input;
    reverseReference<R>(input.data(),inner.data(),n,masks.data());
    reverseInner<R>(input.data(),n,masks.data(),[&](std::size_t i,k::block v){actual[i]=v;});
    if(std::memcmp(actual.data(),inner.data(),n*16))throw std::runtime_error("dense inner mismatch");
    if constexpr(R==2) {
        k::innerReverse<k::JointMapR2>(input.data(),n,masks.data(),[&](std::size_t i,k::block v){actual[i]=v;});
        if(std::memcmp(actual.data(),inner.data(),n*16))throw std::runtime_error("original two-update kernel mismatch");
    }
    forwardReference<R>(other.data(),forward.data(),n,masks.data());
    auto lhs=_mm_setzero_si128(),rhs=lhs;
    for(std::size_t i=0;i<n;++i){lhs=_mm_xor_si128(lhs,_mm_and_si128(input[i].mData,forward[i].mData));rhs=_mm_xor_si128(rhs,_mm_and_si128(other[i].mData,inner[i].mData));}
    if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535)throw std::runtime_error("inner adjoint mismatch");
    for(std::size_t i=0;i<n;++i){auto x=packets.inverse[i];materialized[(x/1024)*1024+4*(x%256)+(x/256)%4]=inner[i];}
    for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(materialized.data()+i,expected.data()+i/2);
    const auto scratchBlocks=(n/1024)*tileStride;
    std::vector<k::block> storage(scratchBlocks+4);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63));
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    encode<G,R,Streaming>(input.data(),scratch,n,masks.data(),packets,nullptr);
    if(std::memcmp(input.data(),expected.data(),n*16))throw std::runtime_error("full encoder/suffix mismatch");
    std::cout<<"checks: bijection, regional occupancy, packet membership, dense inner, adjoint, full encoder, suffix\n";
    if(!calls)return;
    fill(input,913);double phases[2]{};
    auto run=[&]{encode<G,R,Streaming>(input.data(),scratch,n,masks.data(),packets,profile?phases:nullptr);};
    for(unsigned i=0;i<3;++i)run();
    std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i){const auto start=Clock::now();run();times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());}
    std::sort(times.begin(),times.end());
    std::cout<<"exponent,packet,t,s,updates,stream,seed,calls,median_ms,checksum\n"<<exponent<<','<<G<<",128,19,"<<R<<','<<Streaming<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)<<times[calls/2]<<','<<std::hex<<hash(input)<<std::dec<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)<<",repack_bch,"<<phases[1]/(calls+3)<<'\n';
}
template<unsigned G,unsigned R> static void dispatch(bool stream,unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    if(stream)experiment<G,R,true>(exponent,seed,calls,profile);else experiment<G,R,false>(exponent,seed,calls,profile);
}
int main(int argc,char** argv) {try {
    if(argc<7 || argc>8)throw std::invalid_argument("usage: joint exponent packet updates stream seed calls [profile]; calls=0 checks only");
    const unsigned exponent=std::stoul(argv[1]),packet=std::stoul(argv[2]),updates=std::stoul(argv[3]),stream=std::stoul(argv[4]);
    const auto seed=std::stoull(argv[5]);const unsigned calls=std::stoul(argv[6]);const bool profile=argc==8 && std::stoul(argv[7]);
    if(exponent<14 || exponent>22 || (packet!=2 && packet!=4) || (updates!=2 && updates!=3) || stream>1)throw std::invalid_argument("unsupported screen parameters");
#if defined(__GNUC__)
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") || !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))throw std::runtime_error("requires AVX512F/VL/BW and GFNI");
#endif
    if(packet==2){if(updates==2)dispatch<2,2>(stream,exponent,seed,calls,profile);else dispatch<2,3>(stream,exponent,seed,calls,profile);}
    else {if(updates==2)dispatch<4,2>(stream,exponent,seed,calls,profile);else dispatch<4,3>(stream,exponent,seed,calls,profile);}
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
