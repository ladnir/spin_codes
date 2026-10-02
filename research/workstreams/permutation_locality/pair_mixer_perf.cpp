// Research-only matched-packet mixer. Reuse the existing checked inner/BCH
// kernels; all new randomness and buffers are prepared outside timed calls.
#define SPIN_JOINT_NO_MAIN
#include "joint.cpp"
#undef SPIN_JOINT_NO_MAIN

namespace pm {
namespace gf=spin::research::gf16;
struct Pair {
    std::uint16_t x,y;
    std::uint8_t a,b,e;
    std::uint32_t ca,cb,ce;
};
struct Setup {
    Packets<4> route;
    std::vector<Pair> pairs;
    Setup(unsigned rows,std::uint64_t seed):route(rows,seed,true,true) {
        Sampler rng(rows,seed^0x58fc01679a34ULL);
        std::array<unsigned,256> order;
        pairs.reserve(rows/4*128);
        for(unsigned group=0;group<rows/4;++group) {
            rng.shuffle(order);
            const auto begin=pairs.size();
            for(unsigned j=0;j<256;j+=2) {
                auto x=order[j],y=order[j+1];
                if(x>y)std::swap(x,y);
                const auto a=1+rng.divisors[15].sample(rng.words,15);
                const auto b=1+rng.divisors[15].sample(rng.words,15);
                const auto e=1+rng.divisors[15].sample(rng.words,15);
                pairs.push_back(Pair{std::uint16_t(4*x),std::uint16_t(4*y),
                    std::uint8_t(a),std::uint8_t(b),std::uint8_t(e),
                    gf::transposeControl(a),gf::transposeControl(b),gf::transposeControl(e)});
            }
            // The random matching is unchanged. Order disjoint work by its
            // first cache-line address to make the local traversal predictable.
            std::sort(pairs.begin()+begin,pairs.end(),[](const Pair& a,const Pair& b){return a.x<b.x;});
        }
    }
};

static void fieldReference(const k::block* in,k::block* out,unsigned a,bool transpose) {
    for(unsigned j=0;j<4;++j) {
        auto v=_mm_setzero_si128();
        for(unsigned i=0;i<4;++i) {
            const auto bit=transpose ? ((gf::multiply(a,1U<<j)>>i)&1U)
                                     : ((gf::multiply(a,1U<<i)>>j)&1U);
            if(bit)v=_mm_xor_si128(v,in[i].mData);
        }
        out[j]=k::block(v);
    }
}

static void localReference(k::block* tile,const Pair* pairs,bool transpose) {
    for(unsigned j=0;j<128;++j) {
        const auto& p=pairs[j];
        k::block u[4],v[4],w[4],x[4],y[4];
        if(transpose) {
            fieldReference(tile+p.y,v,p.e,true);
            fieldReference(v,w,2,true);
            for(unsigned b=0;b<4;++b) {
                u[b]=k::block(_mm_xor_si128(tile[p.x+b].mData,v[b].mData));
                y[b]=k::block(_mm_xor_si128(tile[p.x+b].mData,w[b].mData));
            }
            fieldReference(u,x,p.a,true);fieldReference(y,v,p.b,true);
            std::memcpy(tile+p.x,x,sizeof(x));std::memcpy(tile+p.y,v,sizeof(v));
        } else {
            fieldReference(tile+p.x,u,p.a,false);fieldReference(tile+p.y,v,p.b,false);
            fieldReference(v,w,2,false);
            for(unsigned b=0;b<4;++b) {
                x[b]=k::block(_mm_xor_si128(u[b].mData,v[b].mData));
                y[b]=k::block(_mm_xor_si128(u[b].mData,w[b].mData));
            }
            fieldReference(y,v,p.e,false);
            std::memcpy(tile+p.x,x,sizeof(x));std::memcpy(tile+p.y,v,sizeof(v));
        }
    }
}

static SPIN_FORCEINLINE void localTranspose(k::block* tile,const Pair* pairs) {
    for(unsigned j=0;j<128;++j) {
        const auto& p=pairs[j];
        const auto x=_mm512_load_si512(tile+p.x);
        const auto v=gf::transpose(_mm512_load_si512(tile+p.y),p.ce);
        const auto a=gf::transpose(_mm512_xor_si512(x,v),p.ca);
        // In the polynomial basis, multiplication by alpha=2 transposes to
        // (v1,v2,v3,v0+v1). Keep its two fixed shuffles explicit.
        const auto rotated=_mm512_permutexvar_epi64(_mm512_setr_epi64(2,3,4,5,6,7,0,1),v);
        const auto extra=_mm512_maskz_permutexvar_epi64(0xc0,_mm512_setr_epi64(0,1,0,1,0,1,2,3),v);
        const auto b=gf::transpose(_mm512_xor_si512(x,_mm512_xor_si512(rotated,extra)),p.cb);
        _mm512_store_si512(tile+p.x,a);_mm512_store_si512(tile+p.y,b);
    }
}

template<unsigned R,bool LocalCopy> static SPIN_NOINLINE void encode(
    k::block* input,k::block* scratch,unsigned n,const unsigned* masks,const Setup& setup,double* phases) {
    auto start=Clock::time_point{};if(phases)start=Clock::now();
    __m128i v1{},v2{},v3{};
    reverseInner<R>(input,n,masks,[&](std::size_t i,k::block v) {
        switch(i&3) {
            case 3:v3=v.mData;break;
            case 2:v2=v.mData;break;
            case 1:v1=v.mData;break;
            case 0:{
                auto x=_mm512_castsi128_si512(v.mData);
                x=_mm512_inserti32x4(x,v1,1);x=_mm512_inserti32x4(x,v2,2);x=_mm512_inserti32x4(x,v3,3);
                _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+setup.route.bases[i/4]),x);
                break;
            }
        }
    });
    _mm_sfence();
    auto middle=Clock::time_point{};if(phases)middle=Clock::now();
    alignas(64) k::block local[1024];
    for(unsigned tile=0;tile<n/1024;++tile) {
        auto* block=scratch+tile*tileStride;
        if constexpr(LocalCopy) {
            for(unsigned i=0;i<1024;i+=16) {
                const auto a=_mm512_load_si512(block+i),b=_mm512_load_si512(block+i+4);
                const auto c=_mm512_load_si512(block+i+8),d=_mm512_load_si512(block+i+12);
                _mm512_store_si512(local+i,a);_mm512_store_si512(local+i+4,b);
                _mm512_store_si512(local+i+8,c);_mm512_store_si512(local+i+12,d);
            }
        }
        auto* source=LocalCopy?local:block;
        localTranspose(source,setup.pairs.data()+128*tile);
        k::bchTranspose4GfniBlend(source,input+512*tile);
    }
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}

template<unsigned R> static void experiment(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile,bool localCopy) {
    const unsigned n=2U<<exponent;
    Setup setup(n/256,seed);
    std::vector<unsigned> masks(2*R*(n/128));
    k::setup::Words rng(seed^0x3f625a92ULL);
    for(std::size_t i=0;i<masks.size();i+=2) {
        unsigned u;do{u=unsigned(rng())&((1U<<19)-1);}while(!u);
        unsigned v=unsigned(rng())&((1U<<19)-1);
        if(std::popcount(u&v)&1)v^=u&-u;
        masks[i]=u;masks[i+1]=v;
    }
    std::vector<k::block> input(n),expected(n),inner(n),materialized(n),adjoint(1024),probe(1024),forward(1024);
    const auto scratchBlocks=(n/1024)*tileStride;
    std::vector<k::block> storage(scratchBlocks+4);
    auto* scratch=reinterpret_cast<k::block*>((reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63));
    k::workspace_routing::adviseOwned(scratch,scratchBlocks*16);
    auto encodeSelected=[&](double* phases) {
        if(localCopy)encode<R,true>(input.data(),scratch,n,masks.data(),setup,phases);
        else encode<R,false>(input.data(),scratch,n,masks.data(),setup,phases);
    };
    fill(probe,115);fill(adjoint,113);forward=probe;
    localReference(forward.data(),setup.pairs.data(),false);
    auto lhs=_mm_setzero_si128(),rhs=lhs;
    for(unsigned i=0;i<1024;++i)lhs=_mm_xor_si128(lhs,_mm_and_si128(adjoint[i].mData,forward[i].mData));
    localReference(adjoint.data(),setup.pairs.data(),true);
    for(unsigned i=0;i<1024;++i)rhs=_mm_xor_si128(rhs,_mm_and_si128(adjoint[i].mData,probe[i].mData));
    if(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs,rhs))!=65535)throw std::runtime_error("pair mixer adjoint mismatch");
    for(unsigned pattern=0;pattern<(calls?1U:3U);++pattern) {
        if(pattern==1) {std::fill(input.begin(),input.end(),k::block{});input[n-3]=k::block(1,3);}
        else fill(input,pattern?997:913);
        expected=input;
        reverseReference<R>(input.data(),inner.data(),n,masks.data());
        for(unsigned i=0;i<n;++i) {
            const auto x=setup.route.inverse[i];
            materialized[(x/1024)*1024+4*(x%256)+(x/256)%4]=inner[i];
        }
        for(unsigned tile=0;tile<n/1024;++tile) {
            localReference(materialized.data()+1024*tile,setup.pairs.data()+128*tile,true);
            k::bchTranspose4(materialized.data()+1024*tile,expected.data()+512*tile);
        }
        encodeSelected(nullptr);
        if(std::memcmp(input.data(),expected.data(),n*16))throw std::runtime_error("pair mixer full encoder/suffix mismatch");
    }
    std::cout<<"checks: matching bijection by shuffle, local adjoint, independent scalar full encoder, suffix\n";
    if(!calls)return;
    fill(input,913);double phases[2]{};
    auto run=[&]{encodeSelected(profile?phases:nullptr);};
    for(unsigned j=0;j<3;++j)run();
    std::vector<double> times;times.reserve(calls);
    for(unsigned j=0;j<calls;++j) {
        const auto start=Clock::now();run();
        times.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
    }
    std::sort(times.begin(),times.end());
    std::cout<<"exponent,updates,seed,calls,local_copy,median_ms,checksum\n"<<exponent<<','<<R<<','<<seed<<','<<calls<<','<<localCopy<<','
        <<std::fixed<<std::setprecision(6)<<times[calls/2]<<','<<std::hex<<hash(input)<<std::dec<<'\n';
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)
        <<",mixer_bch,"<<phases[1]/(calls+3)<<'\n';
}
}

int main(int argc,char** argv) {try {
    if(argc<5 || argc>7)throw std::invalid_argument("usage: pair_mixer exponent updates seed calls [profile [local_copy=1]]");
    const unsigned exponent=std::stoul(argv[1]),updates=std::stoul(argv[2]),calls=std::stoul(argv[4]);
    const auto seed=std::stoull(argv[3]);const bool profile=argc>=6 && std::stoul(argv[5]);
    const bool localCopy=argc<7 || std::stoul(argv[6]);
    if(exponent<14 || exponent>22 || updates<1 || updates>4)throw std::invalid_argument("unsupported parameters");
#if defined(__GNUC__)
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") || !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))
        throw std::runtime_error("requires AVX512F/VL/BW and GFNI");
#endif
    if(updates==1)pm::experiment<1>(exponent,seed,calls,profile,localCopy);
    else if(updates==2)pm::experiment<2>(exponent,seed,calls,profile,localCopy);
    else if(updates==3)pm::experiment<3>(exponent,seed,calls,profile,localCopy);
    else pm::experiment<4>(exponent,seed,calls,profile,localCopy);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
