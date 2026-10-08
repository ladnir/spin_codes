#include "Full.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstdlib>
#if defined(__linux__)
#include <sys/mman.h>
#endif

using namespace spin::research::char2;
using Clock=std::chrono::steady_clock;
static void fill(std::vector<Block>& x,std::mt19937_64& rng) {for(auto& b:x)b={rng(),rng()};}
struct Payload {
    Block* p{};std::size_t count,allocated;
    explicit Payload(std::size_t n):count(n) {
        const bool huge=n*sizeof(Block)>=16*1024*1024;
        const std::size_t alignment=huge?2*1024*1024:64;
        allocated=(n*sizeof(Block)+alignment-1)&~(alignment-1);
#if defined(_MSC_VER)
        p=static_cast<Block*>(_aligned_malloc(allocated,alignment));
#else
        p=static_cast<Block*>(std::aligned_alloc(alignment,allocated));
#endif
        if(!p)throw std::bad_alloc();
#if defined(__linux__)
        if(huge)madvise(p,allocated,MADV_HUGEPAGE);
#endif
        std::memset(p,0,allocated);
    }
    ~Payload(){
#if defined(_MSC_VER)
        _aligned_free(p);
#else
        std::free(p);
#endif
    }
    Payload(const Payload&)=delete;Payload& operator=(const Payload&)=delete;
};
static void check() {
    std::mt19937_64 rng(20261004);
    for(auto b:{8U,16U})for(auto k:{16U,24U,64U,256U})for(auto seed:{1U,91U}) {
        FullPlan plan(k,b,seed,true);plan.prepareGather();
        std::vector<Block> x(k),z(2*k),fast(k),s1(2*k),s2(2*k);
        fill(x,rng);fill(z,rng);
        auto expected=transposeOracle(plan,z);
        transposeFast<false>(plan,z.data(),fast.data(),s1.data(),s2.data());
        if(fast!=expected)throw std::runtime_error("lookup transpose mismatch");
        transposeFast<true>(plan,z.data(),fast.data(),s1.data(),s2.data());
        if(fast!=expected)throw std::runtime_error("prepared transpose mismatch");
        transposeGather(plan,z.data(),fast.data(),s1.data(),s2.data());
        if(fast!=expected)throw std::runtime_error("gather transpose mismatch");
        const auto y=forwardOracle(plan,x);Block lhs{},rhs{};
        for(std::size_t i=0;i<2*k;++i)lhs=bxor(lhs,mulOracle(y[i],z[i]));
        for(std::size_t i=0;i<k;++i)rhs=bxor(rhs,mulOracle(x[i],fast[i]));
        if(lhs!=rhs)throw std::runtime_error("field adjoint mismatch");
        auto aliased=z;
        transposeFast<false>(plan,aliased.data(),aliased.data(),s1.data(),s2.data());
        if(!std::equal(expected.begin(),expected.end(),aliased.begin())||
           !std::equal(z.begin()+k,z.end(),aliased.begin()+k))throw std::runtime_error("in-place suffix");
        // Check F-linearity with a full-width scalar, not only binary scalars.
        const Block scalar{rng(),rng()};auto scaled=x;
        for(auto& a:scaled)a=mulOracle(a,scalar);
        const auto ys=forwardOracle(plan,scaled);
        for(std::size_t i=0;i<2*k;++i)if(ys[i]!=mulOracle(y[i],scalar))throw std::runtime_error("field linearity");
    }
    std::cout<<"full encoder: scalar, SIMD, field adjoint, field linearity, alias checks pass\n";
}
template<bool Prepared,bool Gather=false>static void bench(std::size_t k,unsigned b,std::uint64_t seed,unsigned calls) {
    const auto setupStart=Clock::now();FullPlan plan(k,b,seed,Prepared);
    if constexpr(Gather)plan.prepareGather();
    const auto setupMs=std::chrono::duration<double,std::milli>(Clock::now()-setupStart).count();
    Payload input(plan.n),s1(plan.n),s2(plan.n);std::mt19937_64 rng(913);
    for(std::size_t i=0;i<plan.n;++i)input.p[i]={rng(),rng()};
    auto encode=[&]{
        if constexpr(Gather)transposeGather(plan,input.p,input.p,s1.p,s2.p);
        else transposeFast<Prepared>(plan,input.p,input.p,s1.p,s2.p);
    };
    for(unsigned w=0;w<5;++w)encode();
    std::vector<double> times(calls);
    for(auto& t:times) {const auto begin=Clock::now();encode();t=std::chrono::duration<double,std::micro>(Clock::now()-begin).count();}
    std::sort(times.begin(),times.end());std::uint64_t hash=0;
    for(std::size_t i=0;i<k;++i){hash=(hash^input.p[i].lo)*0x100000001b3ULL;hash=(hash^input.p[i].hi)*0x100000001b3ULL;}
    std::cout<<"outer,K,seed,calls,coefficient_mode,median_us,p10_us,p90_us,setup_ms,setup_bytes,scratch_bytes,checksum\n"
        <<b<<','<<k<<','<<seed<<','<<calls<<','<<(Gather?"gather":Prepared?"prepared":"lookup")<<','<<std::fixed<<std::setprecision(4)
        <<times[calls/2]<<','<<times[calls/10]<<','<<times[9*calls/10]<<','<<setupMs<<','<<plan.setupBytes()<<','
        <<2*plan.n*sizeof(Block)<<','<<std::hex<<hash<<std::dec<<'\n';
}
static void phases(std::size_t k,unsigned b,std::uint64_t seed,unsigned calls) {
    FullPlan plan(k,b,seed);Payload input(plan.n),s1(plan.n),s2(plan.n);
    std::mt19937_64 rng(913);for(std::size_t i=0;i<plan.n;++i)input.p[i]={rng(),rng()};
    for(unsigned i=0;i<5;++i)transposeFast<false>(plan,input.p,input.p,s1.p,s2.p);
    std::vector<double> first(calls),second(calls),outer(calls),total(calls);
    auto us=[](auto a,auto b){return std::chrono::duration<double,std::micro>(b-a).count();};
    for(unsigned i=0;i<calls;++i) {
        auto a=Clock::now();
        scanTranspose<false>(input.p,s1.p,plan.permutation[1].data(),plan.coefficient[1].data(),nullptr,plan.n);
        auto c=Clock::now();
        scanTranspose<false>(s1.p,s2.p,plan.permutation[0].data(),plan.coefficient[0].data(),nullptr,plan.n);
        auto d=Clock::now();
        if(b==8)outerTranspose<8>(s2.p,input.p,plan.n/8);else outerTranspose<16>(s2.p,input.p,plan.n/16);
        auto e=Clock::now();first[i]=us(a,c);second[i]=us(c,d);outer[i]=us(d,e);total[i]=us(a,e);
    }
    for(auto* values:{&first,&second,&outer,&total})std::sort(values->begin(),values->end());
    std::cout<<"phases,"<<b<<','<<k<<','<<seed<<','<<std::fixed<<std::setprecision(4)
        <<first[calls/2]<<','<<second[calls/2]<<','<<outer[calls/2]<<','<<total[calls/2]<<'\n';
}
int main(int argc,char** argv) {try {
    if(!gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc<4||argc>6)throw std::invalid_argument("usage: full OUTER K lookup|prepared [seed=113] [calls=51]");
    const auto b=std::stoul(argv[1]);const auto k=std::stoull(argv[2]);const std::string mode=argv[3];
    const auto seed=argc>4?std::stoull(argv[4]):113;const auto calls=argc>5?std::stoul(argv[5]):51;
    if(calls<1||calls>1000000)throw std::invalid_argument("call count");
    if(mode=="lookup")bench<false>(k,b,seed,calls);
    else if(mode=="prepared")bench<true>(k,b,seed,calls);
    else if(mode=="phases")phases(k,b,seed,calls);
    else if(mode=="gather")bench<false,true>(k,b,seed,calls);
    else throw std::invalid_argument("coefficient mode");
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
