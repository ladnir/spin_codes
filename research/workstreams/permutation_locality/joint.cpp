// Research-only joint packet/inner screen. No production defaults change.
#include "Spin.h"
#include "Inner.h"
#include "SetupRandom.h"
#include "WorkspaceRouting.h"
#include "TwoBitTiled.h"
#include "Gf16Packet.h"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
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
void bchTranspose4GfniIndependent(const block*,block*,const std::uint16_t*);
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
    static_assert(R>=1 && R<=4);
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
    std::vector<std::uint8_t> multipliers;
    Packets(unsigned rows,std::uint64_t seed,bool field=false,bool shared=false,bool pairwise=false):inverse(rows*256),bases(inverse.size()/G),
        controls(bases.size()),offsets(inverse.size(),65535) {
        if(field && G!=4)throw std::invalid_argument("GF16 needs four-bit packets");
        if(shared && G!=4)throw std::invalid_argument("shared route needs four-row groups");
        if(pairwise && (G!=4 || !field || shared))throw std::invalid_argument("pairwise route needs GF16 and distinct pair shuffles");
        if(field)multipliers.resize(bases.size());
        Sampler rng(rows,seed);
        std::vector<unsigned> coordinates(inverse.size()),groups(rows/G);
        std::array<unsigned,256> permutation;
        std::array<unsigned,G> lanes;
        for(unsigned row=0;row<rows;++row) {
            if(pairwise ? row%2==0 : (!shared || row%G==0))rng.shuffle(permutation);
            for(unsigned c=0;c<256;++c)coordinates[row*256+c]=permutation[c];
        }
        // Default: independent row shuffles. The explicit shared experiment
        // draws one coordinate permutation per group; regions remain independent.
        for(unsigned region=0;region<256;++region) {
            rng.shuffle(groups);
            for(unsigned group=0;group<rows/G;++group) {
                if(field)std::iota(lanes.begin(),lanes.end(),0U);else rng.shuffle(lanes);
                if(field)multipliers[(region*rows)/G+groups[group]]=std::uint8_t(1+rng.divisors[15].sample(rng.words,15));
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
                    if(shared && column!=anchor)throw std::runtime_error("shared packet column mismatch");
                    if(pairwise && column!=inverse[i+(lane&~1U)]%256)
                        throw std::runtime_error("paired packet column mismatch");
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
        if(field)for(std::size_t i=0;i<controls.size();++i)controls[i]=spin::research::gf16::transposeControl(multipliers[i]);
    }
};

template<unsigned G,bool Pairwise=false,bool Direct=false,bool Fused=false> static SPIN_NOINLINE void repackBch(
    const k::block* scratch,k::block* out,std::size_t n,const std::uint16_t* offsets) {
    static_assert(!Direct || Pairwise || Fused);
    alignas(64) k::block local[1024],canonical[1024];
    for(std::size_t tile=0;tile<n/1024;++tile) {
        const auto* source=scratch+tile*tileStride;
        // The independent route retains its sequential cold-data copy.
        // Pairwise Direct reads the 16 KiB routed tile in place instead.
        if constexpr(!Direct)for(unsigned i=0;i<1024;i+=16) {
            const auto a=_mm512_load_si512(source+i),b=_mm512_load_si512(source+i+4);
            const auto c=_mm512_load_si512(source+i+8),d=_mm512_load_si512(source+i+12);
            _mm512_store_si512(local+i,a);_mm512_store_si512(local+i+4,b);
            _mm512_store_si512(local+i+8,c);_mm512_store_si512(local+i+12,d);
        }
        const auto* pairedSource=Direct?source:local;
        const auto* p=offsets+tile*1024;
        if constexpr(Fused) {
            static_assert(G==4 && !Pairwise);
            // Load each independently permuted dense coordinate directly into
            // GFNI's bit-transpose preparation. No canonical dense tile is
            // written and reloaded; the original BCH map is unchanged.
            k::bchTranspose4GfniIndependent(pairedSource,out+tile*512,p);
            continue;
        }
#if defined(__GNUC__)
#pragma GCC unroll 4
#endif
        for(unsigned c=0;c<256;++c) {
            if constexpr(Pairwise) {
                static_assert(G==4);
                // Each pair already has its two BCH lanes together. Keep
                // full-line routing, but restore columns with two 256-bit
                // loads instead of four independent 128-bit lane loads.
                const auto low=_mm256_load_si256(reinterpret_cast<const __m256i*>(pairedSource+4*c));
                const auto high=_mm256_load_si256(reinterpret_cast<const __m256i*>(pairedSource+p[4*c+2]));
                const auto x=_mm512_inserti64x4(_mm512_castsi256_si512(low),high,1);
                _mm512_store_si512(canonical+4*c,x);
            } else {
            auto x=_mm512_castsi128_si512(_mm_load_si128(reinterpret_cast<const __m128i*>(local+G*c)));
            x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(local+p[4*c+1])),1);
            if constexpr(G==2)x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(local+512+2*c)),2);
            else x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(local+p[4*c+2])),2);
            x=_mm512_inserti32x4(x,_mm_load_si128(reinterpret_cast<const __m128i*>(local+p[4*c+3])),3);
            _mm512_store_si512(canonical+4*c,x);
            }
        }
        k::bchTranspose4GfniBlend(canonical,out+tile*512);
    }
}

template<unsigned G,unsigned R,bool Streaming,bool Field=false,bool Shared=false,bool Pairwise=false,bool Direct=false,bool Fused=false> static SPIN_NOINLINE void encode(
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
                    if constexpr(Field)x=spin::research::gf16::transpose(x,std::uint32_t(packets.controls[i/4]));
                    else {
                        const auto control=_mm512_cvtepu8_epi64(_mm_loadl_epi64(reinterpret_cast<const __m128i*>(packets.controls.data()+i/4)));
                        x=_mm512_permutexvar_epi64(control,x);
                    }
                    auto* dst=reinterpret_cast<__m512i*>(scratch+packets.bases[i/4]);
                    if constexpr(Streaming)_mm512_stream_si512(dst,x);else _mm512_store_si512(dst,x);
                    break;
                }
            }
        }
    });
    if constexpr(Streaming)_mm_sfence();
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    if constexpr(Shared) {
        static_assert(G==4);
        // Shared coordinates put all four lanes directly in canonical BCH
        // order. No tile copy, coordinate gather, or intermediate array.
        for(std::size_t tile=0;tile<n/1024;++tile)
            k::bchTranspose4GfniBlend(scratch+tile*tileStride,in+tile*512);
    } else repackBch<G,Pairwise,Direct,Fused>(scratch,in,n,packets.offsets.data());
    if(phases){phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();}
}

// Stable bucket routing, as in the production fast path. Consecutive inner
// pairs go to consecutive slots within their destination bucket. This turns
// the long-range random stores into a few sequential streams. All randomness
// stays in the sampled inverse route; only the order of memory accesses changes.
template<unsigned R,bool PairStores,bool Gfni,bool Gather=false> static SPIN_NOINLINE void encodeTiled(
    k::block* in,k::block* buckets,k::block* tile,std::size_t n,const unsigned* masks,
    const spin::research::TwoBitTiled& layout,double* phases) {
    auto start=Clock::time_point{};if(phases)start=Clock::now();
    __m128i high{};
    const auto* slots=layout.pairSlots.data();
    reverseInner<R>(in,n,masks,[&](std::size_t i,k::block v) {
        if constexpr(PairStores) {
            if(i&1)high=v.mData;
            else {
                const auto pair=_mm256_inserti128_si256(_mm256_castsi128_si256(v.mData),high,1);
                _mm256_store_si256(reinterpret_cast<__m256i*>(buckets+slots[i/2]),pair);
            }
        } else buckets[slots[i/2]+(i&1)]=v;
    });
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    const auto tileSize=layout.tileBlocks;
    for(std::size_t base=0;base<n;base+=tileSize) {
        const auto* source=buckets+(base/tileSize)*layout.bucketStride;
        if constexpr(Gather) {
            const auto* indices=layout.gatherSources.data()+base;
            // Read four independent 128-bit values, then write one full line.
            // Only a four-row canonical tile is materialized at a time.
            for(unsigned group=0;group<tileSize;group+=1024) {
                for(unsigned j=0;j<1024;j+=4) {
                    const auto* p=indices+group+j;
                    auto v=_mm512_castsi128_si512(source[p[0]].mData);
                    v=_mm512_inserti32x4(v,source[p[1]].mData,1);
                    v=_mm512_inserti32x4(v,source[p[2]].mData,2);
                    v=_mm512_inserti32x4(v,source[p[3]].mData,3);
                    _mm512_store_si512(tile+j,v);
                }
                k::bchTranspose4GfniBlend(tile,in+(base+group)/2);
            }
        } else {
        const auto* offsets=layout.offsets.data()+base;
        // Same lookahead and scalar cache-local scatter as Fast.cpp. Four-row
        // BCH remains fixed-width SIMD; no generalized loop replaces its circuit.
        for(unsigned j=0;j<tileSize;++j) {
            if(j+32<tileSize)
                _mm_prefetch(reinterpret_cast<const char*>(tile+offsets[j+32]),_MM_HINT_T0);
            tile[offsets[j]]=source[j];
        }
        for(unsigned j=0;j<tileSize;j+=1024) {
            if constexpr(Gfni)k::bchTranspose4GfniBlend(tile+j,in+(base+j)/2);
            else k::bchTranspose4(tile+j,in+(base+j)/2);
        }
        }
    }
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

static void checkTiledLayout() {
    std::vector<unsigned> route(2048);
    for(unsigned pair=0;pair<4;++pair)for(unsigned c=0;c<256;++c) {
        route[512*pair+2*c]=512*pair+c;
        route[512*pair+2*c+1]=512*pair+256+c;
    }
    for(bool padded:{false,true})for(unsigned rows:{4U,8U,16U}) {
        const spin::research::TwoBitTiled layout(route,rows,padded);
        const spin::research::TwoBitTiled gather(route,rows,padded,true);
        std::vector<bool> seen((route.size()/layout.tileBlocks)*layout.bucketStride);
        for(unsigned i=0;i<route.size();++i) {
            const auto slot=layout.pairSlots[i/2]+(i&1),bucket=slot/layout.bucketStride;
            const auto local=slot%layout.bucketStride;
            if(slot>=seen.size() || seen[slot] || local>=layout.tileBlocks ||
               bucket*layout.tileBlocks+layout.offsets[bucket*layout.tileBlocks+local]!=spin::research::TwoBitTiled::packed(route[i]))
                throw std::runtime_error("tiled layout coordinate replay");
            seen[slot]=true;
            if(gather.gatherSources[bucket*layout.tileBlocks+layout.offsets[bucket*layout.tileBlocks+local]]!=local)
                throw std::runtime_error("tiled gather coordinate replay");
        }
    }
    auto rejected=[](const std::vector<unsigned>& inverse,unsigned rows) {
        try {const spin::research::TwoBitTiled layout(inverse,rows,true);}
        catch(const std::invalid_argument&) {return;}
        throw std::runtime_error("invalid tiled layout accepted");
    };
    rejected({},4);rejected(route,2);rejected(route,6);
    auto bad=route;bad[1]=bad[0];rejected(bad,4);
    bad=route;bad[1]=bad[2];rejected(bad,4);
    bad=route;bad[1]=unsigned(bad.size());rejected(bad,4);
    bad=route;std::swap(bad[1],bad[513]);rejected(bad,4);
    bad=route;bad.pop_back();rejected(bad,4);
}
// Modes 0/1 preserve the earlier direct-store baselines. Modes 2/3/4 select
// tiled scalar stores, tiled pair stores, and pair stores with original BCH.
// Mode 5 pads each bucket by one cache line to avoid power-of-two set aliasing.
static void checkGf16() {
    for(unsigned a=1;a<16;++a) {
        unsigned outputs=0;
        for(unsigned x=0;x<16;++x) {
            const auto y=spin::research::gf16::multiply(a,x);outputs|=1U<<y;
            for(unsigned z=0;z<16;++z)if(spin::research::gf16::multiply(a,x^z)!=(y^spin::research::gf16::multiply(a,z)))
                throw std::runtime_error("GF16 linearity");
            alignas(64) std::uint64_t lanes[8];
            for(unsigned j=0;j<4;++j)lanes[2*j]=lanes[2*j+1]=0-std::uint64_t((x>>j)&1);
            const auto actual=spin::research::gf16::transpose(_mm512_load_si512(lanes),spin::research::gf16::transposeControl(a));
            _mm512_store_si512(lanes,actual);
            for(unsigned j=0;j<4;++j) {
                const auto expected=0-std::uint64_t(std::popcount(x&spin::research::gf16::multiply(a,1U<<j))&1);
                if(lanes[2*j]!=expected || lanes[2*j+1]!=expected)throw std::runtime_error("GF16 SIMD adjoint");
            }
        }
        if(outputs!=65535)throw std::runtime_error("GF16 not bijective");
    }
    for(unsigned x=1;x<16;++x) {
        unsigned outputs=0;
        for(unsigned a=1;a<16;++a)outputs|=1U<<spin::research::gf16::multiply(a,x);
        if(outputs!=65534)throw std::runtime_error("GF16 nonzero action");
    }
}
template<unsigned G,unsigned R,unsigned Mode> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile,unsigned tileRows) {
    const unsigned n=2U<<exponent;
    constexpr bool pairwise=Mode>=12 && Mode<=15;
    constexpr bool field=Mode==7 || Mode==8 || Mode==9 || Mode==10 || pairwise || Mode==16 || Mode==17;
    constexpr bool shared=Mode>=9 && Mode<=11;
    const Packets<G> packets(n/256,seed,field,shared,pairwise);
    // Allocate the same additional workspace in all two-bit comparisons.
    // Setup/validation and allocation are outside every timed call.
    const auto tiled=G==2?std::make_unique<spin::research::TwoBitTiled>(packets.inverse,tileRows,Mode>=5,Mode==6):nullptr;
    std::vector<k::block> tileStorage(tiled?tiled->tileBlocks+4:0);
    auto* tile=tiled?reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(tileStorage.data())+63)&~std::uintptr_t(63)):nullptr;
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
    if constexpr(R==1) {
        k::innerReverse<Map>(input.data(),n,masks.data(),[&](std::size_t i,k::block v){actual[i]=v;});
        if(std::memcmp(actual.data(),inner.data(),n*16))throw std::runtime_error("original one-update kernel mismatch");
    }
    forwardReference<R>(other.data(),forward.data(),n,masks.data());
    auto lhs=_mm_setzero_si128(),rhs=lhs;
    for(std::size_t i=0;i<n;++i){lhs=_mm_xor_si128(lhs,_mm_and_si128(input[i].mData,forward[i].mData));rhs=_mm_xor_si128(rhs,_mm_and_si128(other[i].mData,inner[i].mData));}
    if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535)throw std::runtime_error("inner adjoint mismatch");
    auto materialize=[&] {
        for(std::size_t i=0;i<n;++i) {
            auto value=inner[i].mData;
            if constexpr(field) {
                value=_mm_setzero_si128();
                const auto column=spin::research::gf16::multiply(packets.multipliers[i/4],1U<<(i&3));
                for(unsigned bit=0;bit<4;++bit)if((column>>bit)&1)value=_mm_xor_si128(value,inner[(i&~std::size_t(3))+bit].mData);
            }
            const auto x=packets.inverse[i];materialized[(x/1024)*1024+4*(x%256)+(x/256)%4]=k::block(value);
        }
    };
    materialize();
    for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(materialized.data()+i,expected.data()+i/2);
    const auto scratchBlocks=(n/1024)*tileStride;
    std::vector<k::block> storage(scratchBlocks+4);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63));
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    auto encodeSelected=[&](double* phases) {
        if constexpr(Mode==16 || Mode==17)encode<G,R,true,true,false,false,Mode==17,true>(input.data(),scratch,n,masks.data(),packets,phases);
        else if constexpr(pairwise)encode<G,R,Mode==12 || Mode==14,true,false,true,Mode>=14>(input.data(),scratch,n,masks.data(),packets,phases);
        else if constexpr(shared)encode<G,R,Mode!=10,field,true>(input.data(),scratch,n,masks.data(),packets,phases);
        else if constexpr(field)encode<G,R,Mode==7,true>(input.data(),scratch,n,masks.data(),packets,phases);
        else if constexpr(Mode<2)encode<G,R,Mode==1>(input.data(),scratch,n,masks.data(),packets,phases);
        else encodeTiled<R,Mode!=2,Mode!=4,Mode==6>(input.data(),scratch,tile,n,masks.data(),*tiled,phases);
    };
    encodeSelected(nullptr);
    if(std::memcmp(input.data(),expected.data(),n*16))throw std::runtime_error("full encoder/suffix mismatch");
    std::cout<<"checks: bijection, regional occupancy, packet membership, dense inner, adjoint, full encoder, suffix\n";
    if(!calls) {
        // Replay a sparse input and a different dense input through the complete
        // encoder, not just the inner. Verify in-place suffix preservation too.
        for(unsigned pattern=0;pattern<2;++pattern) {
            if(pattern)fill(input,997);else {std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
            expected=input;
            reverseReference<R>(input.data(),inner.data(),n,masks.data());
            materialize();
            for(std::size_t i=0;i<n;i+=1024)k::bchTranspose4(materialized.data()+i,expected.data()+i/2);
            encodeSelected(nullptr);
            if(std::memcmp(input.data(),expected.data(),n*16))throw std::runtime_error("sparse/dense full replay mismatch");
        }
        return;
    }
    fill(input,913);double phases[2]{};
    auto run=[&]{encodeSelected(profile?phases:nullptr);};
    for(unsigned i=0;i<3;++i)run();
    std::vector<double> times;times.reserve(calls);
    for(unsigned i=0;i<calls;++i){const auto start=Clock::now();run();times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());}
    std::sort(times.begin(),times.end());
    std::cout<<"exponent,packet,t,s,updates,mode,tile_rows,seed,calls,median_ms,checksum\n"<<exponent<<','<<G<<",128,19,"<<R<<','<<Mode<<','<<(tiled?tiled->tileBlocks/256:0)<<','<<seed<<','<<calls<<','<<std::fixed<<std::setprecision(6)<<times[calls/2]<<','<<std::hex<<hash(input)<<std::dec<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)<<",repack_bch,"<<phases[1]/(calls+3)<<'\n';
}
template<unsigned G,unsigned R> static void dispatch(unsigned mode,unsigned exponent,std::uint64_t seed,unsigned calls,bool profile,unsigned tileRows) {
    if(mode==0)experiment<G,R,0>(exponent,seed,calls,profile,tileRows);
    else if(mode==1)experiment<G,R,1>(exponent,seed,calls,profile,tileRows);
    else if constexpr(G==4) {
        if(mode==7)experiment<G,R,7>(exponent,seed,calls,profile,tileRows);
        else if(mode==8)experiment<G,R,8>(exponent,seed,calls,profile,tileRows);
        else if(mode==9)experiment<G,R,9>(exponent,seed,calls,profile,tileRows);
        else if(mode==10)experiment<G,R,10>(exponent,seed,calls,profile,tileRows);
        else if(mode==11)experiment<G,R,11>(exponent,seed,calls,profile,tileRows);
        else if(mode==12)experiment<G,R,12>(exponent,seed,calls,profile,tileRows);
        else if(mode==13)experiment<G,R,13>(exponent,seed,calls,profile,tileRows);
        else if(mode==14)experiment<G,R,14>(exponent,seed,calls,profile,tileRows);
        else if(mode==15)experiment<G,R,15>(exponent,seed,calls,profile,tileRows);
        else if(mode==16)experiment<G,R,16>(exponent,seed,calls,profile,tileRows);
        else if(mode==17)experiment<G,R,17>(exponent,seed,calls,profile,tileRows);
    }
    else if constexpr(G==2) {
        switch(mode) {
            case 2:experiment<G,R,2>(exponent,seed,calls,profile,tileRows);break;
            case 3:experiment<G,R,3>(exponent,seed,calls,profile,tileRows);break;
            case 4:experiment<G,R,4>(exponent,seed,calls,profile,tileRows);break;
            case 5:experiment<G,R,5>(exponent,seed,calls,profile,tileRows);break;
            case 6:if constexpr(R==2)experiment<G,R,6>(exponent,seed,calls,profile,tileRows);break;
        }
    }
}
#ifndef SPIN_JOINT_NO_MAIN
int main(int argc,char** argv) {try {
    if(argc<7 || argc>9)throw std::invalid_argument("usage: joint exponent packet updates mode seed calls [profile [tile_rows]]; mode=0 cached,1 stream,2 tiled-scalar,3 tiled-pair,4 tiled-original-BCH,5 tiled-padded,6 tiled-gather,7 GF16-stream,8 GF16-cached,9 shared-GF16-stream,10 shared-GF16-cached,11 shared-lane-stream,12 pairwise-GF16-stream,13 pairwise-GF16-cached,14 direct-pairwise-stream,15 direct-pairwise-cached,16 copied-fused-BCH,17 direct-fused-BCH; calls=0 checks only");
    const unsigned exponent=std::stoul(argv[1]),packet=std::stoul(argv[2]),updates=std::stoul(argv[3]),stream=std::stoul(argv[4]);
    const auto seed=std::stoull(argv[5]);const unsigned calls=std::stoul(argv[6]);const bool profile=argc>=8 && std::stoul(argv[7]);
    const unsigned tileRows=argc==9?std::stoul(argv[8]):16;
    if(exponent<14 || exponent>22 || (packet!=2 && packet!=4) || updates<1 || updates>4 || stream>17 || (packet==4 && ((stream>1 && stream<7) || updates==1)) || (packet==2 && stream>6) || (stream==6 && updates!=2) || (updates==4 && (packet!=4 || stream<7)))throw std::invalid_argument("unsupported screen parameters");
#if defined(__GNUC__)
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") || !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))throw std::runtime_error("requires AVX512F/VL/BW and GFNI");
#endif
    if(!calls)checkTiledLayout();
    if(stream>=7)checkGf16();
    if(packet==2){if(updates==1)dispatch<2,1>(stream,exponent,seed,calls,profile,tileRows);else if(updates==2)dispatch<2,2>(stream,exponent,seed,calls,profile,tileRows);else dispatch<2,3>(stream,exponent,seed,calls,profile,tileRows);}
    else {if(updates==2)dispatch<4,2>(stream,exponent,seed,calls,profile,tileRows);else if(updates==3)dispatch<4,3>(stream,exponent,seed,calls,profile,tileRows);else dispatch<4,4>(stream,exponent,seed,calls,profile,tileRows);}
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
#endif
