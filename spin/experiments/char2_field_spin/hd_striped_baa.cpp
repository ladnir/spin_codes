#include "HdFastOuter.h"
#include "BenchPayload.h"
#include "RateSpeed4Oracle.h"
#include <atomic>
#include <iomanip>
#include <iostream>

// Rejected research control: grouped-check degree4 limits K20 distance to
// <=832321/2097152, below40%, for every setup. Kept for route cost experiments.
// One true global
// permutation, a column-major GRS[128,64] outer, and two GF16 input-weighted
// prefix scans. Setup, coefficients, and caller-owned workspaces are reused.
namespace c2=spin::research::char2;
namespace out=c2::hdfast::outer_detail;
using Block=c2::Block;
using Clock=std::chrono::steady_clock;
#ifndef SB_PREFETCH
#define SB_PREFETCH 0
#endif
#ifndef SB_SPLIT
#define SB_SPLIT 0
#endif
#ifndef SB_BUCKET_LOG
#define SB_BUCKET_LOG 0
#endif
#if defined(_MSC_VER)
#define SB_NOINLINE __declspec(noinline)
#else
#define SB_NOINLINE __attribute__((noinline))
#endif
struct Plan {
    std::size_t k,n,rows;
    std::vector<std::uint32_t> route;
    std::vector<std::uint32_t> stageIndex,localRoute;
    std::vector<std::uint16_t> d0,d1;
    std::vector<c2::wide16::Coefficient> c0,c1;
    explicit Plan(std::size_t size,std::uint64_t seed):k(size),n(2*size),rows(size/64){
        if(!size||size%256||size>(std::uint64_t(1)<<30))throw std::invalid_argument("natural K multiple256 required");
        std::mt19937_64 rng(seed);route.resize(n);std::iota(route.begin(),route.end(),0);
        for(auto i=n-1;i;--i){std::uniform_int_distribution<std::size_t> draw(0,i);std::swap(route[i],route[draw(rng)]);}
        if constexpr(SB_BUCKET_LOG>0){
            const std::size_t bucketSize=std::size_t(1)<<SB_BUCKET_LOG;
            std::vector<std::uint32_t> inverse(n),next((n+bucketSize-1)/bucketSize);
            for(std::size_t i=0;i<n;++i)inverse[route[i]]=std::uint32_t(i);
            for(std::size_t b=0;b<next.size();++b)next[b]=std::uint32_t(b*bucketSize);
            stageIndex.resize(n);localRoute.resize(n);
            for(std::size_t i=0;i<n;++i){const auto dest=inverse[i];stageIndex[i]=next[dest>>SB_BUCKET_LOG]++;localRoute[dest]=stageIndex[i];}
        }
        for(auto* d:{&d0,&d1}){d->resize(n);for(auto& a:*d)do{a=std::uint16_t(rng());}while(!a);}
        c0.resize(n);c1.resize(n);c2::wide16::prepare(d0.data(),c0.data(),n);c2::wide16::prepare(d1.data(),c1.data(),n);
    }
};
SB_NOINLINE static void outer(const Plan& p,const Block* input,Block* planes){
    for(std::size_t first=0;first<p.rows;first+=64){
        const auto last=std::min(first+64,p.rows);
        for(auto row=first;row!=last;row+=4)out::coset64<0,true>(input+64*row,planes+row,p.rows);
        for(auto row=first;row!=last;row+=4)out::coset64<64,true>(input+64*row,planes+64*p.rows+row,p.rows);
    }
}
// Prefix within four records, then add the state broadcast. The two fixed
// lane permutations have compile-time masks; no serial chain of four XORs.
static inline __m512i prefix(__m512i x,__m512i& carry){
    const auto shift1=_mm512_setr_epi64(0,1,0,1,2,3,4,5);
    const auto shift2=_mm512_setr_epi64(0,1,0,1,0,1,2,3);
    x=_mm512_xor_si512(x,_mm512_maskz_permutexvar_epi64(0xfc,shift1,x));
    x=_mm512_xor_si512(x,_mm512_maskz_permutexvar_epi64(0xf0,shift2,x));
    x=_mm512_xor_si512(x,carry);
    carry=_mm512_shuffle_i64x2(x,x,0xff);
    return x;
}
static inline void stageStore(Block* data,const std::uint32_t* index,__m512i x){
    _mm_storeu_si128(reinterpret_cast<__m128i*>(data+index[0]),_mm512_castsi512_si128(x));
    _mm_storeu_si128(reinterpret_cast<__m128i*>(data+index[1]),_mm512_extracti32x4_epi32(x,1));
    _mm_storeu_si128(reinterpret_cast<__m128i*>(data+index[2]),_mm512_extracti32x4_epi32(x,2));
    _mm_storeu_si128(reinterpret_cast<__m128i*>(data+index[3]),_mm512_extracti32x4_epi32(x,3));
}
SB_NOINLINE static void scan(const Plan& p,Block* data,Block* stage){
    auto carry=_mm512_setzero_si512();
    for(std::size_t i=0;i!=p.n;i+=16){
        // Independent coefficient work is scheduled before recurrence work.
        auto a=c2::wide16::multiply512(_mm512_loadu_si512(data+i),p.c0.data()+i);
        auto b=c2::wide16::multiply512(_mm512_loadu_si512(data+i+4),p.c0.data()+i+4);
        auto c=c2::wide16::multiply512(_mm512_loadu_si512(data+i+8),p.c0.data()+i+8);
        auto d=c2::wide16::multiply512(_mm512_loadu_si512(data+i+12),p.c0.data()+i+12);
        a=prefix(a,carry);b=prefix(b,carry);c=prefix(c,carry);d=prefix(d,carry);
        if constexpr(SB_BUCKET_LOG>0){
            stageStore(stage,p.stageIndex.data()+i,a);stageStore(stage,p.stageIndex.data()+i+4,b);
            stageStore(stage,p.stageIndex.data()+i+8,c);stageStore(stage,p.stageIndex.data()+i+12,d);
        }else{
            _mm512_storeu_si512(data+i,a);_mm512_storeu_si512(data+i+4,b);
            _mm512_storeu_si512(data+i+8,c);_mm512_storeu_si512(data+i+12,d);
        }
    }
}
static inline __m512i gather4(const Block* data,const std::uint32_t* indices){
    auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(data+indices[0])));
    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(data+indices[1])),1);
    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(data+indices[2])),2);
    return _mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(data+indices[3])),3);
}
SB_NOINLINE static void routedScan(const Plan& p,const Block* data,Block* output){
    auto carry=_mm512_setzero_si512();
    const auto* indices=SB_BUCKET_LOG?p.localRoute.data():p.route.data();
    for(std::size_t i=0;i!=p.n;i+=16){
        if constexpr(SB_PREFETCH>0){if(i+SB_PREFETCH+16<=p.n){
            c2::c8Unroll(std::make_index_sequence<16>{},[&](auto j){
                _mm_prefetch(reinterpret_cast<const char*>(data+indices[i+SB_PREFETCH+unsigned(j)]),_MM_HINT_T0);
            });
        }}
        auto a=c2::wide16::multiply512(gather4(data,indices+i),p.c1.data()+i);
        auto b=c2::wide16::multiply512(gather4(data,indices+i+4),p.c1.data()+i+4);
        auto c=c2::wide16::multiply512(gather4(data,indices+i+8),p.c1.data()+i+8);
        auto d=c2::wide16::multiply512(gather4(data,indices+i+12),p.c1.data()+i+12);
        if constexpr(!SB_SPLIT){a=prefix(a,carry);b=prefix(b,carry);c=prefix(c,carry);d=prefix(d,carry);}
        _mm512_storeu_si512(output+i,a);_mm512_storeu_si512(output+i+4,b);
        _mm512_storeu_si512(output+i+8,c);_mm512_storeu_si512(output+i+12,d);
    }
    if constexpr(SB_SPLIT)for(std::size_t i=0;i!=p.n;i+=4)
        _mm512_storeu_si512(output+i,prefix(_mm512_loadu_si512(output+i),carry));
}
struct Times{double outer,pre,post;};
template<bool Timed=false> SB_NOINLINE static void forward(const Plan& p,const Block* input,Block* output,Block* scratch,Times* t=nullptr,Block* stage=nullptr){
    if constexpr(SB_BUCKET_LOG>0)if(!stage)throw std::invalid_argument("bucket workspace missing");
    auto a=Clock::time_point{};if constexpr(Timed)a=Clock::now();outer(p,input,scratch);
    if constexpr(Timed){auto b=Clock::now();t->outer=std::chrono::duration<double,std::micro>(b-a).count();a=b;}
    scan(p,scratch,stage);
    if constexpr(Timed){auto b=Clock::now();t->pre=std::chrono::duration<double,std::micro>(b-a).count();a=b;}
    routedScan(p,SB_BUCKET_LOG?stage:scratch,output);
    if constexpr(Timed)t->post=std::chrono::duration<double,std::micro>(Clock::now()-a).count();
}
static std::uint16_t mul(std::uint16_t a,std::uint16_t b){
    const auto& t=c2::p8Oracle().mul;auto high=t[a>>8][b>>8];
    return std::uint16_t(t[a&255][b&255]^t[32][high])|
        (std::uint16_t(t[a&255][b>>8]^t[a>>8][b&255]^high)<<8);
}
static Block scale(Block x,std::uint16_t c){std::uint16_t a[8];std::memcpy(a,&x,16);for(auto& v:a)v=mul(v,c);std::memcpy(&x,a,16);return x;}
static void add(Block& a,Block b){a.lo^=b.lo;a.hi^=b.hi;}
static const auto& evaluations(){
    static const auto table=[] {
        std::array<std::array<std::uint8_t,64>,128> r{};const auto& o=c2::p8Oracle();
        for(unsigned p=0;p!=128;++p){std::uint8_t v[6]{};
            for(unsigned bit=0;bit!=6;++bit){std::uint8_t n=1,d=1;for(unsigned j=0;j!=(1U<<bit);++j){n=o.mul[n][p^j];d=o.mul[d][(1U<<bit)^j];}v[bit]=o.mul[n][o.inverse[d]];}
            for(unsigned j=0;j!=64;++j){r[p][j]=1;for(unsigned bit=0;bit!=6;++bit)if(j&(1U<<bit))r[p][j]=o.mul[r[p][j]][v[bit]];}
        }return r;
    }();return table;
}
static std::vector<Block> oracle(const Plan& p,const Block* x){
    std::vector<Block> work(p.n),result(p.n);const auto& e=evaluations();
    for(std::size_t i=0;i!=p.k;++i){const auto row=(i/256)*4+i%4,degree=(i%256)/4;
        for(unsigned j=0;j!=128;++j)add(work[j*p.rows+row],scale(x[i],e[j][degree]));
    }
    Block s{};for(std::size_t i=0;i!=p.n;++i){add(s,scale(work[i],p.d0[i]));work[i]=s;}
    s={};for(std::size_t i=0;i!=p.n;++i){add(s,scale(work[p.route[i]],p.d1[i]));result[i]=s;}return result;
}
static void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
static void check(){
    std::mt19937_64 rng(17);const Block guard{113,917};
    if constexpr(SB_BUCKET_LOG>0){
        Plan p(1U<<16,7919);std::vector<bool> seen(p.n);
        for(std::size_t i=0;i<p.n;++i){
            require(p.localRoute[i]<p.n&&!seen[p.localRoute[i]],"bucket bijection");seen[p.localRoute[i]]=true;
            require(p.localRoute[i]==p.stageIndex[p.route[i]],"bucket exact permutation");
            require((p.localRoute[i]>>SB_BUCKET_LOG)==(i>>SB_BUCKET_LOG),"bucket locality");
        }
        std::cout<<"check,K,65536,multi-bucket-exact-map,pass\n";
    }
    for(auto k:{256U,512U,1280U}){
        Plan p(k,113);std::vector<Block> source(k);for(auto& x:source)x={rng(),rng()};const auto want=oracle(p,source.data());
        for(unsigned offset=0;offset!=4;++offset){
            std::vector<Block> input(p.n+8,guard),output(p.n+8,guard),scratch(p.n+8,guard),stage(p.n+8,guard);std::copy(source.begin(),source.end(),input.begin()+offset);
            forward(p,input.data()+offset,output.data()+offset,scratch.data()+offset,nullptr,stage.data()+offset);
            require(std::equal(want.begin(),want.end(),output.begin()+offset),"scalar equality");
            require(std::equal(source.begin(),source.end(),input.begin()+offset),"input preserved");
            Times t{};forward<true>(p,input.data()+offset,input.data()+offset,scratch.data()+offset,&t,stage.data()+offset);
            require(std::equal(want.begin(),want.end(),input.begin()+offset),"alias timed equality");
            for(const auto* v:{&input,&output,&scratch,&stage}){for(unsigned i=0;i<offset;++i)require((*v)[i]==guard,"prefix guard");for(auto i=offset+p.n;i<v->size();++i)require((*v)[i]==guard,"suffix guard");}
        }std::cout<<"check,K,"<<k<<",scalar-alias-guards-offsets,pass\n";
    }
}
static std::uint64_t hash(const Block* x,std::size_t n){std::uint64_t h=0;for(std::size_t i=0;i<n;++i){h=(h^x[i].lo)*0x100000001b3ULL;h=(h^x[i].hi)*0x100000001b3ULL;}return h;}
static void produce(Block* p,std::size_t n){std::uint64_t x=913,y=113;for(std::size_t i=0;i<n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;_mm_stream_si128(reinterpret_cast<__m128i*>(p+i),_mm_set_epi64x(y,x));}_mm_sfence();}
static void stats(const char* label,std::vector<double> v){std::sort(v.begin(),v.end());std::cout<<label<<','<<v[v.size()/2]<<','<<v[v.size()/10]<<','<<v[9*v.size()/10]<<'\n';}
static void bench(unsigned calls,bool phases){
    Plan p(1U<<20,113);c2::BenchPayload input(p.k),output(p.n),scratch(p.n),stage(SB_BUCKET_LOG?p.n:1);produce(input.p,p.k);const auto inHash=hash(input.p,p.k);
    for(unsigned j=0;j<3;++j)forward(p,input.p,output.p,scratch.p,nullptr,stage.p);
    std::vector<double> total(calls),ot(calls),pre(calls),post(calls);std::uint64_t consume=0;
    for(unsigned j=0;j<calls;++j){produce(input.p,p.k);Times t{};auto start=Clock::now();
        if(phases)forward<true>(p,input.p,output.p,scratch.p,&t,stage.p);else forward(p,input.p,output.p,scratch.p,nullptr,stage.p);
        total[j]=std::chrono::duration<double,std::micro>(Clock::now()-start).count();ot[j]=t.outer;pre[j]=t.pre;post[j]=t.post;
        std::atomic_signal_fence(std::memory_order_seq_cst);consume^=output.p[(j*104729U)%p.n].lo;
    }
    require(inHash==hash(input.p,p.k),"benchmark source modified");
    std::cout<<"profile,striped-grs128-64-two-weighted-scans,K,"<<p.k<<",N,"<<p.n<<",calls,"<<calls<<",phase_clocks,"<<phases<<",prefetch,"<<SB_PREFETCH<<",split,"<<SB_SPLIT<<",bucket_log,"<<SB_BUCKET_LOG<<",checksum,"<<hash(output.p,p.n)<<",consume,"<<consume<<"\nname,median_us,p10_us,p90_us\n"<<std::fixed<<std::setprecision(3);
    stats("full",total);if(phases){stats("outer",ot);stats("pre_scan",pre);stats("global_route_post_scan",post);}
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc<2||std::string(argv[1])!="--bench")throw std::invalid_argument("--check | --bench [calls=31] [full|phases]");
    const unsigned calls=argc>2?std::stoul(argv[2]):31;const std::string mode=argc>3?argv[3]:"full";
    if(!calls||calls>100000||(mode!="full"&&mode!="phases"))throw std::invalid_argument("benchmark arguments");
    bench(calls,mode=="phases");return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
