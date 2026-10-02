// Isolated BCH and code-equivalent full-encoder comparison. Research only.
#define SPIN_JOINT_NO_MAIN
#include "joint.cpp"
#undef SPIN_JOINT_NO_MAIN
#include <x86intrin.h>

namespace spin::detail::kernel {
void bchTranspose4Restrict(const block*,block*);
void bchTranspose4AlgebraicExact(const block*,block*);
void bchTranspose4AlgebraicRaw(const block*,block*);
void bchTranspose4AlgebraicRawReference(const block*,block*);
}
namespace bc {
template<unsigned Backend> SPIN_FORCEINLINE void bch(const k::block* in,k::block* out) {
    static_assert(Backend<6);
    if constexpr(Backend==0)k::bchTranspose4(in,out);
    else if constexpr(Backend==1)k::bchTranspose4Restrict(in,out);
    else if constexpr(Backend<4)k::bchTranspose4GfniBlend(in,out);
    else if constexpr(Backend==4)k::bchTranspose4AlgebraicExact(in,out);
    else k::bchTranspose4AlgebraicRaw(in,out);
}
struct Buffer {
    std::vector<k::block> storage;
    k::block* data;
    explicit Buffer(std::size_t n):storage(n+4),data(reinterpret_cast<k::block*>(
        (reinterpret_cast<std::uintptr_t>(storage.data())+63)&~std::uintptr_t(63))) {
        k::workspace_routing::adviseOwned(data,n*sizeof(k::block));
    }
};
static std::uint64_t checksum(const k::block* x,std::size_t n) {
    std::uint64_t h=0;
    for(std::size_t i=0;i<n;++i) {
        std::uint64_t v[2];std::memcpy(v,x+i,16);
        for(auto a:v)h=(h^a)*0x100000001b3ULL;
    }
    return h;
}
static void fillBlocks(k::block* x,std::size_t n,std::uint64_t seed) {
    k::setup::Words rng(seed);
    for(std::size_t i=0;i<n;++i){const auto a=rng(),b=rng();x[i]=k::block(a,b);}
}
template<class Run> static void measure(const std::string& kind,unsigned size,unsigned backend,
    std::uint64_t seed,unsigned calls,std::size_t tilesPerCall,Run&& run,const k::block* out,std::size_t outSize) {
    for(unsigned i=0;i<3;++i)run();
    std::vector<double> ms,ticks;ms.reserve(calls);ticks.reserve(calls);
    for(unsigned i=0;i<calls;++i) {
        unsigned aux;
        const auto start=Clock::now();
        _mm_lfence();const auto t0=__rdtscp(&aux);_mm_lfence();
        run();
        _mm_lfence();const auto t1=__rdtscp(&aux);_mm_lfence();
        ms.push_back(std::chrono::duration<double,std::milli>(Clock::now()-start).count());
        ticks.push_back(double(t1-t0));
    }
    std::sort(ms.begin(),ms.end());std::sort(ticks.begin(),ticks.end());
    std::cout<<"kind,size,backend,seed,calls,tiles_per_call,median_ms,p10_ms,p90_ms,ns_per_tile,tsc_ticks_per_tile,checksum\n"
        <<kind<<','<<size<<','<<backend<<','<<seed<<','<<calls<<','<<tilesPerCall<<','<<std::fixed<<std::setprecision(6)
        <<ms[calls/2]<<','<<ms[calls/10]<<','<<ms[(9*calls)/10]<<','<<1e6*ms[calls/2]/tilesPerCall<<','
        <<ticks[calls/2]/tilesPerCall<<','<<std::hex<<checksum(out,outSize)<<std::dec<<'\n';
}
template<unsigned Backend> static SPIN_NOINLINE void bulk(const k::block* in,k::block* out,std::size_t tiles) {
    for(std::size_t tile=0;tile<tiles;++tile)bch<Backend>(in+tile*tileStride,out+tile*512);
}
template<unsigned Backend> static void isolated(bool hot,unsigned size,std::uint64_t seed,unsigned calls) {
    const std::size_t tiles=hot?size:(std::size_t(2)<<size)/1024;
    const std::size_t repeats=hot?std::max<std::size_t>(1,2048/tiles):1;
    Buffer in(tiles*tileStride),out(tiles*512),reference(tiles*512);
    fillBlocks(in.data,tiles*tileStride,seed);
    bulk<0>(in.data,reference.data,tiles);bulk<Backend>(in.data,out.data,tiles);
    if(std::memcmp(out.data,reference.data,tiles*512*16))throw std::runtime_error("isolated BCH mismatch");
    std::cout<<"checks: full isolated BCH output equality\n";
    if(!calls)return;
    measure(hot?"hot":"bulk",size,Backend,seed,calls,tiles*repeats,[&]{
        for(std::size_t repeat=0;repeat<repeats;++repeat)bulk<Backend>(in.data,out.data,tiles);
    },out.data,tiles*512);
}
template<unsigned Backend> static SPIN_NOINLINE void fullEncode(k::block* input,k::block* scratch,std::size_t n,
    const unsigned* masks,const Packets<4>& packets,double* phases) {
    const auto start=phases?Clock::now():Clock::time_point{};
    __m128i v1{},v2{},v3{};
    reverseInner<2>(input,n,masks,[&](std::size_t i,k::block v) {
        switch(i&3) {
            case 3:v3=v.mData;break;
            case 2:v2=v.mData;break;
            case 1:v1=v.mData;break;
            case 0: {
                auto x=_mm512_castsi128_si512(v.mData);
                x=_mm512_inserti32x4(x,v1,1);x=_mm512_inserti32x4(x,v2,2);x=_mm512_inserti32x4(x,v3,3);
                x=spin::research::gf16::transpose(x,std::uint32_t(packets.controls[i/4]));
                _mm512_stream_si512(reinterpret_cast<__m512i*>(scratch+packets.bases[i/4]),x);
                break;
            }
        }
    });
    _mm_sfence();
    const auto middle=phases?Clock::now():Clock::time_point{};
    bulk<Backend>(scratch,input,n/1024);
    if(phases) {
        phases[0]+=std::chrono::duration<double,std::milli>(middle-start).count();
        phases[1]+=std::chrono::duration<double,std::milli>(Clock::now()-middle).count();
    }
}
template<unsigned Backend> static void full(unsigned exponent,std::uint64_t seed,unsigned calls,bool profile) {
    const std::size_t n=std::size_t(2)<<exponent;
    const Packets<4> packets(unsigned(n/256),seed,true,true);
    std::vector<unsigned> masks(4*(n/128));
    k::setup::Words rng(seed^0x3f625a92ULL);
    for(std::size_t i=0;i<masks.size();i+=2) {
        unsigned u;do{u=unsigned(rng())&((1U<<19)-1);}while(!u);
        unsigned v=unsigned(rng())&((1U<<19)-1);
        if(std::popcount(u&v)&1)v^=u&-u;
        masks[i]=u;masks[i+1]=v;
    }
    Buffer input(n),reference(n),scratch(n/1024*tileStride);
    auto selected=[&](double* phases) {
        if constexpr(Backend==3)::encode<4,2,true,true,true>(input.data,scratch.data,n,masks.data(),packets,phases);
        else fullEncode<Backend>(input.data,scratch.data,n,masks.data(),packets,phases);
    };
    for(unsigned pattern=0;pattern<(calls?1U:3U);++pattern) {
        fillBlocks(input.data,n,913+pattern*84);
        if(pattern==1){std::memset(input.data,0,n*16);input.data[n-3]=k::block(1,3);}
        std::memcpy(reference.data,input.data,n*16);
        // Original shared-GF16 R2 implementation remains the code-equivalence oracle.
        ::encode<4,2,true,true,true>(reference.data,scratch.data,n,masks.data(),packets,nullptr);
        selected(nullptr);
        if(std::memcmp(input.data,reference.data,n*16))throw std::runtime_error("full encoder/suffix mismatch");
    }
    std::cout<<"checks: original shared-GF16 R2 encoder equality including suffix\n";
    if(!calls)return;
    fillBlocks(input.data,n,913);double phases[2]{};
    measure(profile?"full-profile":"full",exponent,Backend,seed,calls,n/1024,[&]{
        selected(profile?phases:nullptr);
    },input.data,n);
    if(profile)std::cout<<"phase_means_including_warmups,inner_route,"<<phases[0]/(calls+3)
        <<",bch,"<<phases[1]/(calls+3)<<'\n';
}
static void basisCheck() {
    Buffer input(1024),expected(512),actual(512);
    std::memset(input.data,0,1024*16);
    std::size_t checked=0;
    for(unsigned coordinate=0;coordinate<1024;++coordinate) {
        for(unsigned bit=0;bit<128;++bit) {
            input.data[coordinate]=bit<64?k::block(0,1ULL<<bit):k::block(1ULL<<(bit-64),0);
            bch<0>(input.data,expected.data);
            bch<1>(input.data,actual.data);
            if(std::memcmp(actual.data,expected.data,512*16))throw std::runtime_error("Restrict basis mismatch");
            bch<2>(input.data,actual.data);
            if(std::memcmp(actual.data,expected.data,512*16))throw std::runtime_error("GFNI basis mismatch");
            ++checked;
        }
        input.data[coordinate]=k::block{};
    }
    std::cout<<"Passed "<<checked<<" complete four-row, 128-bit-element basis vectors for both alternate BCH kernels.\n";
}
template<unsigned Backend> static void dispatch(const std::string& kind,unsigned size,std::uint64_t seed,unsigned calls) {
    if(kind=="hot")isolated<Backend>(true,size,seed,calls);
    else if(kind=="bulk")isolated<Backend>(false,size,seed,calls);
    else full<Backend>(size,seed,calls,kind=="full-profile");
}
}
int main(int argc,char** argv) {try {
    if(!__builtin_cpu_supports("avx512f") || !__builtin_cpu_supports("avx512vl") ||
       !__builtin_cpu_supports("avx512bw") || !__builtin_cpu_supports("gfni"))throw std::runtime_error("AVX512/GFNI required");
    if(argc==2 && std::string(argv[1])=="check"){bc::basisCheck();return 0;}
    if(argc!=6)throw std::invalid_argument("usage: bch_compare hot|bulk|full|full-profile size backend seed calls; or check");
    const std::string kind=argv[1];const unsigned size=std::stoul(argv[2]),backend=std::stoul(argv[3]),calls=std::stoul(argv[5]);
    const auto seed=std::stoull(argv[4]);
    if((kind!="hot" && kind!="bulk" && kind!="full" && kind!="full-profile") || backend>3 ||
        (kind=="hot" ? (size!=1 && size!=4 && size!=16) : (size<14 || size>22)))throw std::invalid_argument("unsupported comparison");
    if(backend==0)bc::dispatch<0>(kind,size,seed,calls);
    else if(backend==1)bc::dispatch<1>(kind,size,seed,calls);
    else if(backend==2)bc::dispatch<2>(kind,size,seed,calls);
    else bc::dispatch<3>(kind,size,seed,calls);
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
