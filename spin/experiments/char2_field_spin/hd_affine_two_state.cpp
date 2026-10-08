#include "HdFastOuter.h"
#include "BenchPayload.h"
#include "RateSpeed4Oracle.h"
#include <atomic>
#include <iomanip>
#include <iostream>

// Research cost control, NOT a distance-certified code. GRS[64,32], fixed
// affine placement, one weighted prefix, one regional route, two additive
// states. Setup, allocations, input restoration and hashing are not timed.
namespace c2=spin::research::char2;
namespace out=c2::hdfast::outer_detail;
using Block=c2::Block;
using Coefficient=c2::wide16::Coefficient;
using Clock=std::chrono::steady_clock;
#ifndef AF_SPLIT_ROUTE
#define AF_SPLIT_ROUTE 0
#endif
#ifndef AF_GROUP
#define AF_GROUP 1
#endif
#ifndef AF_NT
#define AF_NT 0
#endif
#ifndef AF_TILE
#define AF_TILE 64
#endif
static_assert(AF_GROUP==1||AF_GROUP==2||AF_GROUP==4);
static_assert(AF_TILE>=4&&AF_TILE%4==0);
#if defined(_MSC_VER)
#define AF_NOINLINE __declspec(noinline)
#define AF_INLINE __forceinline
#else
#define AF_NOINLINE __attribute__((noinline))
#define AF_INLINE inline __attribute__((always_inline))
#endif
static constexpr std::array<unsigned,64> strides={
  12623,13781,1327,8485,15923,13269,9939,15617,11733,7157,4563,9235,4579,3107,8209,4815,
  10163,3237,2417,10819,15471,3299,11593,14227,10361,6701,15631,14507,8535,2041,461,3057,
  13069,37,16173,10917,7993,10657,2063,6261,7265,7819,4669,14679,2989,2637,10487,16033,
  3573,9877,9539,4089,10903,6659,9425,14581,3003,12613,10389,7933,9513,6025,6205,6119};
static constexpr unsigned inverseOdd(unsigned a){
    unsigned x=a;for(unsigned i=0;i!=5;++i)x*=2-a*x;return x;
}
static std::size_t placed(unsigned column,std::size_t row,std::size_t rows){
    return ((inverseOdd(strides[column])*(row/AF_GROUP))&(rows/AF_GROUP-1))*AF_GROUP+row%AF_GROUP;
}
static std::uint16_t mul(std::uint16_t a,std::uint16_t b){
    return c2::medium32::quadraticMultiply(a,b);
}
struct Plan {
    std::size_t k,n,rows;
    std::vector<std::uint16_t> route;
    std::array<std::vector<std::uint16_t>,5> scalar;
    std::array<std::vector<Coefficient>,5> coefficient;
    explicit Plan(std::size_t size,std::uint64_t seed):k(size),n(2*size),rows(size/32){
        if(size<128||size%128||(rows&(rows-1))||rows>65536)
            throw std::invalid_argument("research geometry: K=power of two, 128<=K<=2^21");
        std::mt19937_64 rng(seed);route.resize(n);
        for(unsigned c=0;c!=64;++c){
            auto* r=route.data()+c*rows;std::iota(r,r+rows,std::uint16_t(0));
            for(auto i=rows-1;i;--i){std::uniform_int_distribution<std::size_t> d(0,i);std::swap(r[i],r[d(rng)]);}
        }
        for(auto& s:scalar)s.resize(n);
        auto nonzero=[&](){std::uint16_t a;do{a=std::uint16_t(rng());}while(!a);return a;};
        for(std::size_t i=0;i!=n;++i){
            scalar[0][i]=nonzero();
            std::uint16_t d1,d2,f,cp;
            do{d1=nonzero();d2=nonzero();f=nonzero();cp=nonzero()^d1^mul(f,d2);}while(!cp);
            scalar[1][i]=d1;scalar[2][i]=d2;scalar[3][i]=cp;scalar[4][i]=f;
        }
        for(unsigned j=0;j!=5;++j){coefficient[j].resize(n);c2::wide16::prepare(scalar[j].data(),coefficient[j].data(),n);}
    }
};

template<unsigned Column,bool Affine,bool Stream>
AF_INLINE void storeColumn(__m512i x,Block* planes,std::size_t rows,std::size_t row){
    if constexpr(!Affine){_mm512_storeu_si512(planes+Column*rows+row,x);}
    else if constexpr(AF_GROUP==4){
        constexpr auto a=inverseOdd(strides[Column]);
        auto* target=planes+Column*rows+((a*(row/4))&(rows/4-1))*4;
        if constexpr(Stream)_mm512_stream_si512(reinterpret_cast<__m512i*>(target),x);
        else _mm512_storeu_si512(target,x);
    }else if constexpr(AF_GROUP==2){
        constexpr auto a=inverseOdd(strides[Column]);auto* plane=planes+Column*rows;
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(plane+((a*(row/2))&(rows/2-1))*2),_mm512_castsi512_si256(x));
        _mm256_storeu_si256(reinterpret_cast<__m256i*>(plane+((a*(row/2+1))&(rows/2-1))*2),_mm512_extracti64x4_epi64(x,1));
    }
    else{
        constexpr auto a=inverseOdd(strides[Column]);const auto mask=rows-1;
        auto* plane=planes+Column*rows;
        _mm_storeu_si128(reinterpret_cast<__m128i*>(plane+((a*(row+0))&mask)),_mm512_castsi512_si128(x));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(plane+((a*(row+1))&mask)),_mm512_extracti32x4_epi32(x,1));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(plane+((a*(row+2))&mask)),_mm512_extracti32x4_epi32(x,2));
        _mm_storeu_si128(reinterpret_cast<__m128i*>(plane+((a*(row+3))&mask)),_mm512_extracti32x4_epi32(x,3));
    }
}
template<unsigned Offset,bool Affine,bool Stream>
AF_NOINLINE void leaf16(const __m512i* source,Block* planes,std::size_t rows,std::size_t row){
#define AF_LOAD(I) auto x##I=_mm512_load_si512(source+I)
    AF_LOAD(0);AF_LOAD(1);AF_LOAD(2);AF_LOAD(3);AF_LOAD(4);AF_LOAD(5);AF_LOAD(6);AF_LOAD(7);
    AF_LOAD(8);AF_LOAD(9);AF_LOAD(10);AF_LOAD(11);AF_LOAD(12);AF_LOAD(13);AF_LOAD(14);AF_LOAD(15);
#undef AF_LOAD
#define AF_PAIR(A,B,L,D) out::pair<true,out::normalized(L,Offset^D)>(x##A,x##B)
    AF_PAIR(0,8,3,0);AF_PAIR(1,9,3,0);AF_PAIR(2,10,3,0);AF_PAIR(3,11,3,0);
    AF_PAIR(4,12,3,0);AF_PAIR(5,13,3,0);AF_PAIR(6,14,3,0);AF_PAIR(7,15,3,0);
    AF_PAIR(0,4,2,0);AF_PAIR(1,5,2,0);AF_PAIR(2,6,2,0);AF_PAIR(3,7,2,0);
    AF_PAIR(8,12,2,8);AF_PAIR(9,13,2,8);AF_PAIR(10,14,2,8);AF_PAIR(11,15,2,8);
    AF_PAIR(0,2,1,0);AF_PAIR(1,3,1,0);AF_PAIR(4,6,1,4);AF_PAIR(5,7,1,4);
    AF_PAIR(8,10,1,8);AF_PAIR(9,11,1,8);AF_PAIR(12,14,1,12);AF_PAIR(13,15,1,12);
    AF_PAIR(0,1,0,0);AF_PAIR(2,3,0,2);AF_PAIR(4,5,0,4);AF_PAIR(6,7,0,6);
    AF_PAIR(8,9,0,8);AF_PAIR(10,11,0,10);AF_PAIR(12,13,0,12);AF_PAIR(14,15,0,14);
#undef AF_PAIR
#define AF_STORE(I) storeColumn<Offset+I,Affine,Stream>(x##I,planes,rows,row)
    AF_STORE(0);AF_STORE(1);AF_STORE(2);AF_STORE(3);AF_STORE(4);AF_STORE(5);AF_STORE(6);AF_STORE(7);
    AF_STORE(8);AF_STORE(9);AF_STORE(10);AF_STORE(11);AF_STORE(12);AF_STORE(13);AF_STORE(14);AF_STORE(15);
#undef AF_STORE
}
template<unsigned Offset,bool Affine,bool Stream>
AF_NOINLINE void coset32(const Block* source,Block* planes,std::size_t rows,std::size_t row){
    alignas(64) __m512i v[32];
    c2::c8Unroll(std::make_index_sequence<16>{},[&](auto j){
        auto a=_mm512_loadu_si512(source+4*j),b=_mm512_loadu_si512(source+4*(16+j));
        out::pair<true,out::normalized(4,Offset)>(a,b);_mm512_store_si512(v+j,a);_mm512_store_si512(v+16+j,b);
    });
    leaf16<Offset,Affine,Stream>(v,planes,rows,row);leaf16<Offset+16,Affine,Stream>(v+16,planes,rows,row);
}
template<bool Affine,bool Stream>
AF_NOINLINE void outerImpl(const Plan& p,const Block* input,Block* scratch){
    for(std::size_t first=0;first<p.rows;first+=AF_TILE){
        const auto last=std::min(first+AF_TILE,p.rows);
        for(auto r=first;r!=last;r+=4)coset32<0,Affine,Stream>(input+32*r,scratch,p.rows,r);
        for(auto r=first;r!=last;r+=4)coset32<32,Affine,Stream>(input+32*r,scratch,p.rows,r);
    }
    if constexpr(Stream)_mm_sfence();
}
template<bool Affine=true>
AF_NOINLINE void outer(const Plan& p,const Block* input,Block* scratch){
    if constexpr(AF_NT&&Affine&&AF_GROUP==4){
        if((reinterpret_cast<std::uintptr_t>(scratch)&63)==0)outerImpl<Affine,true>(p,input,scratch);
        else outerImpl<Affine,false>(p,input,scratch);
    }else outerImpl<Affine,false>(p,input,scratch);
}
AF_INLINE __m512i prefix(__m512i x,__m512i& state){
    const auto a=_mm512_setr_epi64(0,1,0,1,2,3,4,5),b=_mm512_setr_epi64(0,1,0,1,0,1,2,3);
    x=_mm512_xor_si512(x,_mm512_maskz_permutexvar_epi64(0xfc,a,x));
    x=_mm512_xor_si512(x,_mm512_maskz_permutexvar_epi64(0xf0,b,x));
    x=_mm512_xor_si512(x,state);state=_mm512_shuffle_i64x2(x,x,0xff);return x;
}
AF_NOINLINE void prescan(const Block* input,Block* stage,const Coefficient* coefficient,std::size_t n,__m512i& state){
    std::size_t i=0;
    for(;i+16<=n;i+=16){
        auto a=c2::wide16::multiply512(_mm512_loadu_si512(input+i),coefficient+i);
        auto b=c2::wide16::multiply512(_mm512_loadu_si512(input+i+4),coefficient+i+4);
        auto c=c2::wide16::multiply512(_mm512_loadu_si512(input+i+8),coefficient+i+8);
        auto d=c2::wide16::multiply512(_mm512_loadu_si512(input+i+12),coefficient+i+12);
        a=prefix(a,state);b=prefix(b,state);c=prefix(c,state);d=prefix(d,state);
        _mm512_storeu_si512(stage+i,a);_mm512_storeu_si512(stage+i+4,b);
        _mm512_storeu_si512(stage+i+8,c);_mm512_storeu_si512(stage+i+12,d);
    }
    for(;i!=n;i+=4)
        _mm512_storeu_si512(stage+i,prefix(c2::wide16::multiply512(_mm512_loadu_si512(input+i),coefficient+i),state));
}
AF_INLINE __m512i gather4(const Block* input,const std::uint16_t* r){
    auto x=_mm512_castsi128_si512(_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+r[0])));
    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+r[1])),1);
    x=_mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+r[2])),2);
    return _mm512_inserti32x4(x,_mm_loadu_si128(reinterpret_cast<const __m128i*>(input+r[3])),3);
}
template<bool Routed>
AF_INLINE __m512i postInput(const Plan& p,const Block* stage,std::size_t j,std::size_t start){
    if constexpr(Routed)return _mm512_loadu_si512(stage+j);
    else return gather4(stage,p.route.data()+start+j);
}
template<bool Routed=false>
AF_NOINLINE void postscan(const Plan& p,const Block* stage,Block* output,std::size_t start,__m512i& s1,__m512i& s2){
    std::size_t j=0;
    for(;j+16<=p.rows;j+=16){
        const auto i=start+j;
#define AF_GATHER(J) auto x##J=postInput<Routed>(p,stage,j+4*J,start)
        AF_GATHER(0);AF_GATHER(1);AF_GATHER(2);AF_GATHER(3);
#undef AF_GATHER
#define AF_PRODUCT(J) auto a##J=c2::wide16::multiply512(x##J,p.coefficient[1].data()+i+4*J); \
        auto b##J=c2::wide16::multiply512(x##J,p.coefficient[2].data()+i+4*J); \
        auto y##J=c2::wide16::multiply512(x##J,p.coefficient[3].data()+i+4*J)
        AF_PRODUCT(0);AF_PRODUCT(1);AF_PRODUCT(2);AF_PRODUCT(3);
#undef AF_PRODUCT
        a0=prefix(a0,s1);a1=prefix(a1,s1);a2=prefix(a2,s1);a3=prefix(a3,s1);
        b0=prefix(b0,s2);b1=prefix(b1,s2);b2=prefix(b2,s2);b3=prefix(b3,s2);
#define AF_OUTPUT(J) _mm512_storeu_si512(output+i+4*J,_mm512_xor_si512(y##J, \
        _mm512_xor_si512(a##J,c2::wide16::multiply512(b##J,p.coefficient[4].data()+i+4*J))))
        AF_OUTPUT(0);AF_OUTPUT(1);AF_OUTPUT(2);AF_OUTPUT(3);
#undef AF_OUTPUT
    }
    for(;j!=p.rows;j+=4){
        const auto i=start+j;auto x=postInput<Routed>(p,stage,j,start);
        auto a=c2::wide16::multiply512(x,p.coefficient[1].data()+i);
        auto b=c2::wide16::multiply512(x,p.coefficient[2].data()+i);
        auto y=c2::wide16::multiply512(x,p.coefficient[3].data()+i);
        a=prefix(a,s1);b=prefix(b,s2);
        y=_mm512_xor_si512(y,_mm512_xor_si512(a,c2::wide16::multiply512(b,p.coefficient[4].data()+i)));
        _mm512_storeu_si512(output+i,y);
    }
}
struct Times{double outer=0,pre=0,post=0;};
template<bool Timed=false>
AF_NOINLINE void forward(const Plan& p,const Block* input,Block* output,Block* scratch,Block* stage,Times* t=nullptr){
    auto start=Clock::time_point{};if constexpr(Timed)start=Clock::now();outer(p,input,scratch);
    if constexpr(Timed)t->outer=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
    auto pre=_mm512_setzero_si512(),s1=pre,s2=pre;
    for(unsigned c=0;c!=64;++c){
        const auto begin=c*p.rows;if constexpr(Timed)start=Clock::now();
        prescan(scratch+begin,stage,p.coefficient[0].data()+begin,p.rows,pre);
        if constexpr(Timed){auto now=Clock::now();t->pre+=std::chrono::duration<double,std::micro>(now-start).count();start=now;}
        if constexpr(AF_SPLIT_ROUTE){
            for(std::size_t j=0;j!=p.rows;j+=4)
                _mm512_storeu_si512(output+begin+j,gather4(stage,p.route.data()+begin+j));
            postscan<true>(p,output+begin,output,begin,s1,s2);
        }else postscan(p,stage,output,begin,s1,s2);
        if constexpr(Timed)t->post+=std::chrono::duration<double,std::micro>(Clock::now()-start).count();
    }
}
static Block scale(Block x,std::uint16_t c){std::uint16_t a[8];std::memcpy(a,&x,16);for(auto& v:a)v=mul(v,c);std::memcpy(&x,a,16);return x;}
static void add(Block& a,Block b){a.lo^=b.lo;a.hi^=b.hi;}
static void require(bool b,const char* s){if(!b)throw std::runtime_error(s);}
static const auto& evaluations(){
    static const auto result=[](){std::array<std::array<std::uint8_t,32>,64> e{};const auto& o=c2::p8Oracle();
        for(unsigned point=0;point!=64;++point){std::uint8_t v[5]{};
            for(unsigned bit=0;bit!=5;++bit){std::uint8_t n=1,d=1;for(unsigned j=0;j!=(1U<<bit);++j){n=o.mul[n][point^j];d=o.mul[d][(1U<<bit)^j];}v[bit]=o.mul[n][o.inverse[d]];}
            for(unsigned j=0;j!=32;++j){e[point][j]=1;for(unsigned bit=0;bit!=5;++bit)if(j&(1U<<bit))e[point][j]=o.mul[e[point][j]][v[bit]];}
        }return e;}();return result;
}
static std::vector<Block> oracle(const Plan& p,const Block* input){
    std::vector<Block> work(p.n),output(p.n);const auto& e=evaluations();
    for(std::size_t i=0;i!=p.k;++i){const auto row=(i/128)*4+i%4,degree=(i%128)/4;
        for(unsigned c=0;c!=64;++c)add(work[c*p.rows+placed(c,row,p.rows)],scale(input[i],e[c][degree]));
    }
    Block pre{},s1{},s2{};
    for(std::size_t i=0;i!=p.n;++i){add(pre,scale(work[i],p.scalar[0][i]));work[i]=pre;}
    for(std::size_t i=0;i!=p.n;++i){auto x=work[(i/p.rows)*p.rows+p.route[i]];
        const auto c=p.scalar[3][i]^p.scalar[1][i]^mul(p.scalar[4][i],p.scalar[2][i]);require(c!=0,"nonzero original diagonal");
        auto y=scale(x,c);add(y,s1);add(y,scale(s2,p.scalar[4][i]));output[i]=y;
        add(s1,scale(x,p.scalar[1][i]));add(s2,scale(x,p.scalar[2][i]));
    }return output;
}
static void check(){
    std::mt19937_64 rng(17);const Block guard{119,811};
    for(auto k:{128U,256U,1024U,4096U}){
        Plan p(k,113);std::vector<Block> source(k);for(auto& x:source)x={rng(),rng()};auto want=oracle(p,source.data());
        for(unsigned offset=0;offset!=4;++offset){
            std::vector<Block> input(p.n+8,guard),output(p.n+8,guard),scratch(p.n+8,guard),stage(p.rows+8,guard),plain(p.n);
            std::copy(source.begin(),source.end(),input.begin()+offset);
            forward(p,input.data()+offset,output.data()+offset,scratch.data()+offset,stage.data()+offset);
            require(std::equal(want.begin(),want.end(),output.begin()+offset),"scalar encoder equality");
            require(std::equal(source.begin(),source.end(),input.begin()+offset),"source preserved");
            outer<false>(p,input.data()+offset,plain.data());
            for(unsigned c=0;c!=64;++c)for(std::size_t r=0;r!=p.rows;++r)
                require(plain[c*p.rows+r]==scratch[c*p.rows+placed(c,r,p.rows)+offset],"affine outer equality");
            Times times;forward<true>(p,input.data()+offset,input.data()+offset,scratch.data()+offset,stage.data()+offset,&times);
            require(std::equal(want.begin(),want.end(),input.begin()+offset),"exact alias");
            for(auto* v:{&input,&output,&scratch}){for(unsigned i=0;i!=offset;++i)require((*v)[i]==guard,"prefix guard");for(auto i=offset+p.n;i!=v->size();++i)require((*v)[i]==guard,"suffix guard");}
            for(unsigned i=0;i!=offset;++i)require(stage[i]==guard,"stage prefix guard");
            for(auto i=offset+p.rows;i!=stage.size();++i)require(stage[i]==guard,"stage suffix guard");
        }std::cout<<"check,K,"<<k<<",scalar-affine-alias-guards-offsets,pass\n";
    }
}
static std::uint64_t hash(const Block* p,std::size_t n){std::uint64_t h=0;for(std::size_t i=0;i!=n;++i){h=(h^p[i].lo)*0x100000001b3ULL;h=(h^p[i].hi)*0x100000001b3ULL;}return h;}
static void produce(Block* p,std::size_t n){std::uint64_t x=913,y=113;for(std::size_t i=0;i!=n;++i){x^=x<<13;x^=x>>7;x^=x<<17;y+=0x9e3779b97f4a7c15ULL;_mm_stream_si128(reinterpret_cast<__m128i*>(p+i),_mm_set_epi64x(y,x));}_mm_sfence();}
static void stats(const char* label,std::vector<double> v){std::sort(v.begin(),v.end());std::cout<<label<<','<<v[v.size()/2]<<','<<v[v.size()/10]<<','<<v[9*v.size()/10]<<'\n';}
static void bench(unsigned calls,const std::string& mode){
    Plan p(1U<<20,113);c2::BenchPayload input(p.k),output(p.n),scratch(p.n),stage(p.rows);produce(input.p,p.k);
    auto run=[&](Times* t){if(mode=="outer-affine")outer(p,input.p,output.p);else if(mode=="outer-linear")outer<false>(p,input.p,output.p);
        else if(mode=="phases")forward<true>(p,input.p,output.p,scratch.p,stage.p,t);else forward(p,input.p,output.p,scratch.p,stage.p);};
    for(unsigned j=0;j!=3;++j){Times t;run(&t);}
    std::vector<double> total(calls),ot(calls),pt(calls),it(calls);std::uint64_t consume=0;
    for(unsigned j=0;j!=calls;++j){produce(input.p,p.k);Times t;auto start=Clock::now();run(&t);
        total[j]=std::chrono::duration<double,std::micro>(Clock::now()-start).count();ot[j]=t.outer;pt[j]=t.pre;it[j]=t.post;
        std::atomic_signal_fence(std::memory_order_seq_cst);consume^=output.p[(j*104729U)%p.n].lo;
    }
    std::cout<<"profile,affine-grs64-32-pre1-post2,K,"<<p.k<<",N,"<<p.n<<",calls,"<<calls<<",mode,"<<mode<<",placement_group,"<<AF_GROUP<<",split_route,"<<AF_SPLIT_ROUTE<<",outer_nt,"<<AF_NT<<",outer_tile,"<<AF_TILE<<",checksum,"<<hash(output.p,p.n)<<",consume,"<<consume<<"\nname,median_us,p10_us,p90_us\n"<<std::fixed<<std::setprecision(3);
    stats(mode.c_str(),total);if(mode=="phases"){stats("outer",ot);stats("pre_scan",pt);stats("regional_route_post",it);}
}
int main(int argc,char** argv){try{
    if(!c2::gfni512Available())return 77;
    if(argc==2&&std::string(argv[1])=="--check"){check();return 0;}
    if(argc<2||std::string(argv[1])!="--bench")throw std::invalid_argument("--check | --bench [calls=31] [full|phases|outer-affine|outer-linear]");
    const unsigned calls=argc>2?std::stoul(argv[2]):31;const std::string mode=argc>3?argv[3]:"full";
    if(!calls||calls>100000||(mode!="full"&&mode!="phases"&&mode!="outer-affine"&&mode!="outer-linear"))throw std::invalid_argument("benchmark arguments");
    bench(calls,mode);return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
